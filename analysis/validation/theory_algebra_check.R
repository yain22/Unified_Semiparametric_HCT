# On a finite observed law, calculate target changes directly from probabilities
# and compare them to the derivatives represented by the proposed gradients.
f_x <- c(.4, .6)
p_trial <- c(.3, .7)
e <- c(.25, .8)
outcomes <- c(-2, 1, 4, 7)
# Rows are treatment, internal control, and external control within each stratum.
outcome_prob <- list(rbind(c(.1,.2,.4,.3), c(.3,.4,.2,.1), c(.4,.1,.2,.3)),
                     rbind(c(.4,.3,.2,.1), c(.15,.25,.45,.15), c(.3,.4,.2,.1)))

share <- function(p, vi, ve) p * (1-e) / (p * (1-e) + (1-p) * vi/ve)

score_numerator <- function(records, p, m, v) {
  x <- records$x
  g <- records$group
  s <- as.integer(g != 3L)
  lambda <- share(p, v[,2], v[,3])[x]
  px <- p[x]
  ex <- e[x]
  resid <- records$y - m[cbind(x, g)]
  contrast <- m[x,1] - lambda * m[x,2] - (1-lambda) * m[x,3]
  core <- s * contrast + (g == 1L) * resid / ex -
    (g == 2L) * lambda * resid / (1-ex) -
    (g == 3L) * px * (1-lambda) * resid / (1-px)
  derivative <- -(s-px)/(1-px) +
    (g == 2L)/(1-ex) * (resid^2/v[x,2] - 1) -
    (g == 3L) * px/(1-px) * (resid^2/v[x,3] - 1)
  core + (m[x,2] - m[x,3]) * lambda * (1-lambda) * derivative
}

observed_law <- function(f, p, conditional) {
  means <- t(vapply(conditional, function(pr) drop(pr %*% outcomes), numeric(3)))
  variances <- t(vapply(seq_along(f), function(x) {
    centered <- sweep(matrix(outcomes, nrow = 3, ncol = 4, byrow = TRUE), 1, means[x,])
    rowSums(conditional[[x]] * centered^2)
  }, numeric(3)))
  q <- sum(f*p)
  lambda <- share(p, variances[,2], variances[,3])
  tce <- sum(f*p * (means[,1] - lambda * means[,2] - (1-lambda) * means[,3])) / q
  trial <- sum(f*p * (means[,1] - means[,2])) / q
  cells <- expand.grid(j = seq_along(outcomes), group = 1:3, x = seq_along(f))
  x <- cells$x
  g <- cells$group
  cells$y <- outcomes[cells$j]
  prob_group <- cbind(p*e, p*(1-e), 1-p)
  prob_outcome <- mapply(function(x,g,j) conditional[[x]][g,j], x, g, cells$j)
  cells$prob <- f[x] * prob_group[cbind(x,g)] * prob_outcome
  s <- as.integer(g != 3L)
  phi <- (score_numerator(cells, p, means, variances) - s*tce) / q
  phi_trial <- (s*(means[x,1] - means[x,2] - trial) +
    (g == 1L)/e[x] * (cells$y - means[x,1]) -
    (g == 2L)/(1-e[x]) * (cells$y - means[x,2])) / q
  list(q = q, tce = tce, trial = trial, gap = tce-trial, means = means, variances = variances,
       shares = lambda, records = cells, gradients = cbind(tce = phi, trial = phi_trial, gap = phi - phi_trial))
}

tilt <- function(probability, direction, amount) {
  w <- probability * exp(amount * direction)
  w / sum(w)
}
perturb <- function(component, amount) {
  f <- f_x
  p <- p_trial
  conditional <- outcome_prob
  if (component == -2L) {
    f <- tilt(f, c(.5, -.2), amount)
  } else if (component == -1L) {
    p <- plogis(qlogis(p) + amount*c(.4, -.3))
  } else {
    x <- component %/% 3L + 1L
    g <- component %% 3L + 1L
    conditional[[x]][g,] <- tilt(conditional[[x]][g,], c(.4,-.6,.2,.8), amount)
  }
  observed_law(f, p, conditional)
}

# The supplement expresses E[H_fit] - q*TCE in terms of nuisance errors.
# Compare that identity to direct summation over all possible records.
exact_drift <- function(base, p_fit, m_fit, v_fit) {
  p_error <- p_fit - p_trial
  m_error <- m_fit - base$means
  v_error <- v_fit - base$variances
  lambda_fit <- share(p_fit, v_fit[,2], v_fit[,3])
  weight_fit <- (1-p_trial)*p_fit/(1-p_fit)
  component <- -p_trial*(base$means[,2]-base$means[,3])*(lambda_fit-base$shares) +
    (1-lambda_fit)*p_error*m_error[,3]/(1-p_fit) +
    (m_fit[,2]-m_fit[,3])*lambda_fit*(1-lambda_fit) * (
      p_error/(1-p_fit) - p_trial*v_error[,2]/v_fit[,2] + weight_fit*v_error[,3]/v_fit[,3] +
      p_trial*m_error[,2]^2/v_fit[,2] - weight_fit*m_error[,3]^2/v_fit[,3])
  sum(f_x * component)
}

theory_check <- function() {
  base <- observed_law(f_x, p_trial, outcome_prob)
  step <- 1e-6
  labels <- c("covariate exponential tilt", "source logistic tilt",
              sprintf("outcome tilt: stratum %d, group %d", rep(0:1, each = 3), rep(0:2, 2)))
  targets <- c("tce", "trial", "gap")
  checks <- lapply(-2:5, function(component) {
    plus <- perturb(component, step)
    minus <- perturb(component, -step)
    derivative <- (unlist(plus[targets]) - unlist(minus[targets]))/(2*step)
    dp <- (plus$records$prob - minus$records$prob)/(2*step)
    represented <- drop(crossprod(dp, base$gradients))
    data.frame(path = labels[component+3L], target = targets, finite_difference = derivative,
               eif_inner_product = represented, absolute_error = abs(derivative-represented),
               row.names = NULL)
  })
  checks <- do.call(rbind, checks)
  set.seed(20261007, kind = "Mersenne-Twister")
  drift_errors <- replicate(500L, {
    p_fit <- runif(2, .1, .9)
    m_fit <- base$means + matrix(runif(6, -1.5, 1.5), nrow = 2)
    v_fit <- base$variances * exp(matrix(runif(6, -.8, .8), nrow = 2))
    direct <- sum(base$records$prob * score_numerator(base$records, p_fit, m_fit, v_fit)) - base$q*base$tce
    abs(direct - exact_drift(base, p_fit, m_fit, v_fit))
  })
  gradient_means <- drop(crossprod(base$records$prob, base$gradients))
  stopifnot(max(checks$absolute_error) < 1e-8, max(drift_errors) < 1e-10,
            max(abs(gradient_means)) < 1e-12)
  list(status = "passed", purpose = "Independent finite-law EIF and exact nuisance-drift identity checks",
       limits = "Algebraic verification on a finite law. General regularity and asymptotics require separate arguments.",
       inputs = list(covariate_probabilities = f_x, source_probabilities = p_trial,
                     known_allocation = e, outcomes = outcomes, outcome_probabilities = outcome_prob),
       baseline = base[c("q", targets, "means", "variances", "shares")],
       finite_difference_step = step, factor_paths = length(labels), target_path_checks = nrow(checks),
       eif_means = as.list(gradient_means), max_derivative_absolute_error = max(checks$absolute_error),
       exact_drift_random_seed = 20261007L, exact_drift_tuples = length(drift_errors),
       max_exact_drift_absolute_error = max(drift_errors), derivative_checks = checks)
}
if (sys.nframe() == 0L) {
  script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  root <- dirname(dirname(dirname(normalizePath(script))))
  args <- commandArgs(TRUE)
  usage <- "Rscript theory_algebra_check.R [--output-dir DIR]"
  if (identical(args, "--help")) { cat(usage, "\n"); quit(status = 0L) }
  if (length(args) && (length(args) != 2 || args[1] != "--output-dir")) stop(usage, call. = FALSE)
  out <- if (length(args)) args[2] else file.path(root, "results", "theory-check")
  report <- theory_check()
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(report, file.path(out, "theory_algebra_results.json"), auto_unbox = TRUE,
                       pretty = TRUE, digits = NA)
  cat(jsonlite::toJSON(report[c("status", "factor_paths", "target_path_checks", "eif_means",
      "max_derivative_absolute_error", "exact_drift_tuples", "max_exact_drift_absolute_error")],
      auto_unbox = TRUE, pretty = TRUE, digits = NA), "\n")
}
