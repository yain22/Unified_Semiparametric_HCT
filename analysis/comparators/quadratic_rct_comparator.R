# Supplementary, paired quadratic ANCOVA comparator for the baseline simulation.
# Run from any directory; --help lists output and repetition options.
# Uses bundled simulated-data hashes and never overwrites published reference results.
args <- commandArgs(FALSE)
script <- normalizePath(sub('^--file=', '', grep('^--file=', args, value = TRUE)[1]), winslash = '/')
here <- dirname(script)
bio <- normalizePath(file.path(here, '../..'), winslash = '/')
cli <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% cli) {
  cat("Rscript quadratic_rct_comparator.R [--repetitions 5000] [--output-dir PATH]\n")
  quit(status = 0)
}
allowed <- c("--repetitions", "--output-dir")
if (length(cli) %% 2L || (length(cli) && any(!cli[seq(1L, length(cli), 2L)] %in% allowed))) stop("Use --help for valid options")
option <- function(key, default) {
  pos <- match(key, cli)
  if (is.na(pos)) default else cli[[pos + 1L]]
}
B <- suppressWarnings(as.integer(option("--repetitions", "5000")))
if (is.na(B) || B < 1L || B > 5000L) stop("repetitions must be between 1 and 5000")
output <- option("--output-dir", file.path(bio, "results/quadratic"))
dir.create(output, recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(output, winslash = "/", mustWork = TRUE)
if (identical(output, here)) stop("Use a separate output directory; published references must not be overwritten")
library(data.table)
setDTthreads(1L)
source(file.path(here, 'inputs/trial_hc2_helper.R'))
pairing <- fread(file.path(here, 'inputs/paired_reference.csv'))
old <- pairing[, .(N_trial, replicate, estimate = linear_estimate, se = linear_se)]
source_files <- file.path(here, 'inputs', c('paired_reference.csv', 'original-baseline-protocol.json', 'trial_hc2_helper.R'))
hashes <- vapply(source_files, function(f) digest::digest(file = f, algo = 'sha256'), character(1))
linear <- as.formula(paste('~', paste(paste0('X', 1:10), collapse = '+')))
quadratic <- as.formula(paste('~', paste(c(paste0('X', 1:10), paste0('I(X', 1:10, '^2)')), collapse = '+')))
generate <- function(n, b) {
  nT <- if (n == 250L) 167L else 333L
  nC <- n - nT
  nE <- 2L * n
  seed <- 52461004L + 100000L * (n == 500L) + b
  set.seed(seed)
  xe <- matrix(rnorm(nE * 10L), ncol = 10L)
  epsE <- rnorm(nE)
  xt <- matrix(rnorm(nT * 10L), ncol = 10L)
  xc <- matrix(rnorm(nC * 10L), ncol = 10L)
  epsT <- rnorm(nT)
  epsC <- rnorm(nC)
  X <- rbind(xt, xc, xe)
  f <- .5 + (.6 / sqrt(5)) * rowSums(X[, 1:5]) + (.25 / sqrt(5)) * rowSums(X[, 1:5]^2)
  error <- c(epsT, exp(.2 * xc[, 1]) * epsC, exp(-.2 * xe[, 1]) * epsE)
  data <- setNames(as.data.frame(X), paste0('X', 1:10))
  data$S <- c(rep(1L, n), rep(0L, nE))
  data$A <- c(rep(1L, nT), rep(0L, nC + nE))
  data$Y <- f + error - (1 - data$S) * X[, 1]
  list(data = data, seed = seed)
}
rows <- vector('list', 2L * B)
validations <- vector('list', 2L * B)
start <- Sys.time()
index <- 0L
for (n in c(250L, 500L)) for (b in seq_len(B)) {
  index <- index + 1L
  d <- generate(n, b)
  hash <- digest::digest(d$data, algo = 'sha256')
  saved <- pairing[N_trial == n & replicate == b]
  stopifnot(nrow(saved) == 1L, saved$seed == d$seed, hash == saved$full_sample_sha256)
  qfit <- hct_trial_wald_fit(d$data, quadratic)
  # Reproduce the archived linear result for every sample as an independent
  # check of outcome generation, treatment coding, HC2, and interval convention.
  lfit <- hct_trial_wald_fit(d$data, linear)
  prev <- old[N_trial == n & replicate == b]
  esterr <- abs(lfit$estimate - prev$estimate)
  seerr <- abs(lfit$se - prev$se)
  stopifnot(nrow(prev) == 1L, esterr < 1e-10, seerr < 1e-10)
  rows[[index]] <- data.table(
    N_trial = n, replicate = b, seed = d$seed, method = 'Quadratic ANCOVA (RCT)',
    estimate = qfit$estimate, se = qfit$se, coverage = as.integer(abs(qfit$estimate) <= qnorm(.975) * qfit$se),
    squared_error = qfit$estimate^2, width = 2 * qnorm(.975) * qfit$se,
    linear_estimate = lfit$estimate, linear_se = lfit$se,
    paired_squared_error_difference = qfit$estimate^2 - lfit$estimate^2,
    paired_width_difference = 2 * qnorm(.975) * (qfit$se - lfit$se)
  )
  validations[[index]] <- data.table(
    N_trial = n, replicate = b, seed = d$seed, full_sample_sha256 = hash,
    linear_estimate_error = esterr, linear_se_error = seerr
  )
  if (b %% 1000L == 0L) {
    cat('Completed nRCT=', n, ', repetitions=', b, '\n', sep = '')
    flush.console()
  }
}
r <- rbindlist(rows)
v <- rbindlist(validations)
s <- r[,
  {
    a <- list(
      n_reps = .N, bias = mean(estimate), empirical_sd = sd(estimate), mean_se = mean(se),
      coverage = mean(coverage), mse = mean(squared_error), width = mean(width)
    )

    for (z in c('estimate', 'coverage', 'squared_error', 'width', 'paired_squared_error_difference', 'paired_width_difference')) {
      a[[paste0(z, '_mean')]] <- mean(get(z))
      a[[paste0(z, '_mcse')]] <- sd(get(z)) / sqrt(.N)
    }
    a
  },
  by = .(N_trial, method)
]
fwrite(r, file.path(output, 'quadratic_rct_replicates.csv.gz'), compress = 'gzip')
fwrite(s, file.path(output, 'quadratic_rct_summary.csv'))
qa_dir <- output
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(v, file.path(qa_dir, 'quadratic_rct_pairing_validation.csv'))
protocol <- list(
  repetitions_per_size = B, seed = '52461004 + 100000*(nRCT==500) + replicate',
  draw_order = 'Exact original external-X, external-error, treated-X, control-X, treated-error, control-error order',
  method = 'OLS Y on intercept, treatment, X1,...,X10, X1^2,...,X10^2 in trial records only; treatment coefficient, HC2 standard error, normal 95% Wald interval',
  comparison = 'Same outcome basis as the HCT and ECT nuisance regressions; retains original ANOVA and linear ANCOVA. Common quadratic main effects are correct for this DGP; it does not claim complete equality of all estimation/uncertainty procedures.',
  gap = 'RCT outcomes and estimates are invariant across the external-discrepancy grid; true trial causal gap is zero',
  samples_hashed = nrow(v), all_hashes_match = TRUE, max_linear_estimate_error = max(v$linear_estimate_error),
  max_linear_se_error = max(v$linear_se_error), source_files = basename(source_files), source_sha256 = as.list(setNames(hashes, basename(source_files))),
  elapsed_seconds = as.numeric(difftime(Sys.time(), start, units = 'secs'))
)
stopifnot(identical(hashes, vapply(source_files, function(f) digest::digest(file = f, algo = 'sha256'), character(1))))
jsonlite::write_json(protocol, file.path(output, 'quadratic_rct_protocol.json'), pretty = TRUE, auto_unbox = TRUE, digits = 16)
print(s)
