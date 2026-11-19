# ECT fits use only the treated and external-control records.

scalar_ect_import_functions <- function(path, names) {
  env <- new.env(parent = .GlobalEnv)
  for (ex in parse(path)) if (is.call(ex) && identical(ex[[1]], as.name('<-')) &&
    is.symbol(ex[[2]]) && as.character(ex[[2]]) %in% names) eval(ex, envir = env)
  if (!all(vapply(names, exists, logical(1), envir = env, inherits = FALSE)))
    stop('Missing requested read-only source functions.')
  env
}

scalar_ect_load <- function(paper, sim_source) {
  sim_names <- c(
    'truncate_prob', 'as_numeric_matrix', 'make_quad_features',
    'fit_or_model', 'predict_or_model', 'fit_ps_model', 'predict_ps_model',
    'predict_qhat', 'att_hat_dr', 'eif_wald_ci_dr'
  )
  list(sim = scalar_ect_import_functions(sim_source, sim_names))
}

scalar_ect_fit <- function(data, sim, alpha = .05, eps = .01) {
  # data is the original complete HCT draw, before any shift. This explicit
  # split guarantees internal-control observations enter neither ECT fit.
  D1 <- data[data$R == 1 & data$A == 1, , drop = FALSE]
  D0 <- data[data$R == 0, , drop = FALSE]
  n1 <- nrow(D1)
  n0 <- nrow(D0)
  N <- n1 + n0
  if (n1 < 22L || n0 < 22L) stop('ECT outcome-basis support is insufficient.')
  warnings_seen <- character()
  fitted <- withCallingHandlers(
    {
      or <- sim$fit_or_model(D0, or_method = 'quad_lm')
      ps <- sim$fit_ps_model(D1, D0, ps_method = 'glm')
      if (or$fit$rank != length(coef(or$fit)) || any(!is.finite(coef(or$fit))))
        stop('ECT external outcome regression is rank deficient.')
      if (!isTRUE(ps$fit$converged) || isTRUE(ps$fit$boundary) ||
        ps$fit$rank != length(coef(ps$fit)) || any(!is.finite(coef(ps$fit))))
        stop('ECT source logistic fit failed; no fallback.')
      estimate <- sim$att_hat_dr(D1, D0, or, ps, eps = eps)
      inference <- sim$eif_wald_ci_dr(D1, D0, or, ps, estimate, alpha = alpha, eps = eps)
      f1 <- sim$predict_or_model(or, D1)
      f0 <- sim$predict_or_model(or, D0)
      raw_p1 <- sim$predict_ps_model(ps, D1)
      raw_p0 <- sim$predict_ps_model(ps, D0)
      q0 <- sim$predict_qhat(ps, D0, eps = eps)
      phi1 <- (D1$Y - f1 - estimate) / (n1 / N)
      phi0 <- -q0 * (D0$Y - f0) / (n1 / N)
      # Diagnostic only. The primary interval reproduces the SIM plug-in EIF
      # variance, not the HCT full nuisance-stacked HC2 sandwich.
      fixed_se <- sqrt((n1 * var(phi1) + n0 * var(phi0)) / N^2)
      list(
        estimate = estimate, se = inference$se, se_fixed_group = fixed_se,
        lower = unname(inference$ci['lower']), upper = unname(inference$ci['upper']),
        influence = c(phi1, phi0), n1 = n1, n0 = n0,
        influence_mean = mean(c(phi1, phi0)),
        external_clipped = sum(raw_p0 < eps | raw_p0 > 1 - eps),
        treated_clipped = sum(raw_p1 < eps | raw_p1 > 1 - eps),
        raw_p_min = min(c(raw_p1, raw_p0)), raw_p_max = max(c(raw_p1, raw_p0)),
        external_ess = sum(q0)^2 / sum(q0^2),
        source_iterations = ps$fit$iter,
        source_nobs = nobs(ps$fit), outcome_nobs = nobs(or$fit),
        formula_error = abs(estimate - (mean(D1$Y - f1) - sum(q0 * (D0$Y - f0)) / n1)),
        variance_error = abs(inference$se - sqrt(mean(c(phi1, phi0)^2) / N))
      )
    },
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart('muffleWarning')
    }
  )
  if (any(grepl('numerically 0 or 1|did not converge', warnings_seen)))
    stop('ECT logistic fit reached a numerical boundary.')
  if (any(!is.finite(c(fitted$estimate, fitted$se, fitted$se_fixed_group))) ||
    fitted$se <= 0 || fitted$se_fixed_group <= 0 || abs(fitted$influence_mean) > 1e-10 ||
    fitted$source_nobs != N || fitted$outcome_nobs != n0)
    stop('ECT fit failed its estimating-equation/sample check.')
  fitted$warnings <- unique(warnings_seen)
  fitted
}
