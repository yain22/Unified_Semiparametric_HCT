# Trial-only regression comparators on the observed outcome.
# No known prognostic function is supplied to either estimator.
# Both methods use HC2 standard errors and two-sided normal Wald intervals.

hct_trial_wald_fit <- function(data, adjustment_terms = ~1, alpha = .05) {
  if (!is.data.frame(data) || !all(c('Y', 'R', 'A') %in% names(data)) ||
    !is.numeric(data$Y) || any(!is.finite(data$Y)) ||
    anyNA(data) || !all(data$R %in% c(0, 1)) || !all(data$A %in% c(0, 1)))
    stop('data must contain finite numeric Y and binary R and A, without missing values.')
  if (!inherits(adjustment_terms, 'formula') || length(adjustment_terms) != 2L ||
    any(c('Y', 'R', 'A') %in% all.vars(adjustment_terms)))
    stop('adjustment_terms must be a one-sided formula containing only baseline covariates.')
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1)
    stop('alpha must be strictly between zero and one.')
  trial <- data[data$R == 1, , drop = FALSE]
  if (any(tabulate(trial$A + 1L, nbins = 2L) < 2L))
    stop('At least two observations are required in each randomized trial arm.')
  rhs <- call('+', as.name('A'), adjustment_terms[[2L]])
  form <- as.formula(call('~', rhs), env = environment(adjustment_terms))
  X <- model.matrix(form, data = trial)
  fit <- lm.fit(X, trial$Y)
  if (fit$rank != ncol(X) || nrow(X) <= ncol(X) || anyNA(fit$coefficients))
    stop('Trial regression is rank deficient or has no residual degrees of freedom.')
  bread <- chol2inv(chol(crossprod(X)))
  leverages <- rowSums((X %*% bread) * X)
  if (any(!is.finite(leverages)) || any(leverages >= 1))
    stop('HC2 requires every regression leverage to be less than one.')
  treatment <- match('A', colnames(X))
  if (is.na(treatment)) stop('Treatment coefficient is missing from the trial regression.')
  coefficient_weight <- drop(X %*% bread[, treatment])
  variance <- sum(coefficient_weight^2 * fit$residuals^2 / (1 - leverages))
  estimate <- unname(fit$coefficients[treatment])
  se <- sqrt(variance)
  if (!is.finite(estimate) || !is.finite(se) || se <= 0)
    stop('Trial estimate or HC2 standard error is invalid.')
  list(
    estimate = estimate, se = se, ci = estimate + qnorm(c(alpha / 2, 1 - alpha / 2)) * se,
    p_value = pnorm(estimate / se, lower.tail = FALSE), success = TRUE,
    variance_type = 'HC2', critical_value_type = 'normal', n_trial = nrow(X),
    rank = fit$rank, maximum_leverage = max(leverages),
    residual_df = nrow(X) - fit$rank
  )
}

hct_trial_wald_comparators <- function(data, adjustment_terms, alpha = .05) {
  list(
    'ANOVA (RCT)' = hct_trial_wald_fit(data, ~1, alpha),
    'ANCOVA (RCT)' = hct_trial_wald_fit(data, adjustment_terms, alpha)
  )
}
