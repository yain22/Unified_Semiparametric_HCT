# Matched-score inference for the population causal-gap, using the baseline draws.
# From the replication package root:
#   Rscript analysis/simulations/check_gap_inference.R [repetitions] [output.csv]
args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args)) as.integer(args[1]) else 500L
out <- if (length(args) >= 2L) args[2] else file.path("results", "gap_inference_summary.csv")
stopifnot(length(args) <= 2L, !is.na(reps), reps >= 2L, reps <= 5000L)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(here, "lib", "fitted_helpers.R"))
helpers <- fg_load(here)

# Load only the pure generator, without launching the full-grid simulation.
seed_base <- 52461004L
designs <- data.frame(N_trial = 250L, NT = 167L, NC = 83L, N0 = 500L)
baseline <- parse(file.path(here, "run_baseline.R"))
generator <- Filter(function(x) is.call(x) && identical(x[[1L]], as.name("<-")) &&
                      identical(x[[2L]], as.name("generate")), as.list(baseline))
stopifnot(length(generator) == 1L)
eval(generator[[1L]])
paired <- read.csv(file.path(here, "..", "comparators", "inputs", "paired_reference.csv"))
paired <- paired[paired$N_trial == 250L, ]
gaps <- c(0, .15)

# For delta(x) = intercept + x1, the gap is intercept*K0 + K1.
external_weight <- function(x) plogis(.8 * x - log(designs$NC / designs$N0))
K0 <- integrate(function(x) external_weight(x) * dnorm(x), -Inf, Inf, rel.tol = 1e-12)$value
K1 <- integrate(function(x) x * external_weight(x) * dnorm(x), -Inf, Inf, rel.tol = 1e-12)$value
amplitudes <- (gaps - K1) / K0

estimate_gap <- function(data) {
  fit <- fg_fit(data, helpers, mean_basis = "quadratic", variance_basis = "x1",
                methods = "U", keep = TRUE)[["U"]]
  d <- fit$prepared
  state <- fit$state
  trial_raw <- d$S * (state$mT - state$mC) +
    d$T / d$e * (data$Y - state$mT) -
    d$C / (1 - d$e) * (data$Y - state$mC)
  trial <- sum(trial_raw) / d$nR
  estimate <- sum(state$raw - trial_raw) / d$nR
  stopifnot(abs(estimate - (fit$estimate - trial)) < 1e-10)

  # Differentiate the matched trial estimating equation with respect to the
  # SAME nuisance vector as the HCT equation before taking their difference.
  nuisance <- setdiff(seq_along(fit$parameters), d$ix$theta)
  derivative <- numeric(length(fit$parameters))
  derivative[d$ix$bT] <- colMeans(d$z * (d$S - d$T / d$e))
  derivative[d$ix$bC] <- colMeans(d$z * (-d$S + d$C / (1 - d$e)))
  adjustment <- solve(t(state$jacobian[nuisance, nuisance]), derivative[nuisance])
  trial_if <- drop((trial_raw - d$S * trial -
                      state$equations[, nuisance, drop = FALSE] %*% adjustment) / d$q)
  gap_if <- fit$influence - trial_if
  se <- helpers$u$hct_unrestricted_group_variance(gap_if, d$group)$se

  # An independent covariance calculation checks that pairing has been retained.
  joint <- matrix(0, 2, 2)
  for (group in levels(d$group)) {
    selected <- d$group == group
    joint <- joint + sum(selected) * cov(cbind(fit$influence[selected], trial_if[selected])) / d$N^2
  }
  paired_variance <- joint[1, 1] + joint[2, 2] - 2 * joint[1, 2]
  stopifnot(abs(se^2 - paired_variance) < 1e-10, abs(mean(gap_if)) < 1e-7)
  c(estimate = estimate, se = se)
}

estimates <- standard_errors <- matrix(NA_real_, reps, length(gaps))
for (b in seq_len(reps)) {
  original <- generate(250L, b)
  record <- paired[match(b, paired$replicate), ]
  stopifnot(original$seed == record$seed,
            digest::digest(original$data, algo = "sha256") == record$full_sample_sha256)
  for (j in seq_along(gaps)) {
    data <- original$data
    data$Y <- data$Y - (1 - data$S) * amplitudes[j]
    result <- estimate_gap(data)
    estimates[b, j] <- result["estimate"]
    standard_errors[b, j] <- result["se"]
  }
}

results <- do.call(rbind, lapply(seq_along(gaps), function(j) {
  error <- estimates[, j] - gaps[j]
  covered <- abs(error) <= qnorm(.975) * standard_errors[, j]
  data.frame(n_treated = designs$NT, n_internal = designs$NC, n_external = designs$N0,
             repetitions = reps, gap = gaps[j], discrepancy_intercept = amplitudes[j],
             bias = mean(error), bias_mcse = sd(error) / sqrt(reps), empirical_sd = sd(error),
             mean_se = mean(standard_errors[, j]), covered = sum(covered), coverage = mean(covered),
             coverage_mcse = sd(as.integer(covered)) / sqrt(reps))
}))
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
write.csv(results, out, row.names = FALSE)
print(results, digits = 5, row.names = FALSE)
