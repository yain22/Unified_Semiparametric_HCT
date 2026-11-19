# Exact function extracts; see source_manifest.json. Sourcing performs no I/O.

cgm_est_variance <- function(z, y, mean_fit, u) {
  residual <- y - drop(z %*% mean_fit$beta)
  response <- residual^2 / (1 - mean_fit$leverage)
  if (any(!is.finite(response)) || any(response < 0) || sum(response) <= 0)
    stop('Invalid HC2 training squared residuals for conditional variance.')
  center <- colMeans(z)
  scale <- apply(z, 2, sd)
  center[1] <- 0
  scale[1] <- 1
  if (any(!is.finite(scale)) || any(scale <= 0)) stop('Degenerate variance basis.')
  zs <- sweep(sweep(z, 2, center, '-'), 2, scale, '/')
  # Quasi-Gamma estimating equations specify only the conditional mean and
  # working variance proportional to mu^2, not a Gamma outcome distribution.
  # Exact convex quasi-Gamma objective. Damped Newton is used uniformly;
  # ordinary Fisher scoring can oscillate in participant bootstrap samples.
  coefficients <- c(log(mean(response)), rep(0, ncol(z) - 1L))
  state <- function(b) {
    eta <- drop(zs %*% b)
    ratio <- response * exp(-eta)
    if (any(!is.finite(c(eta, ratio)))) return(NULL)
    list(objective = sum(eta + ratio), gradient = drop(crossprod(zs, 1 - ratio)), ratio = ratio)
  }
  for (iter in seq_len(100L)) {
    current <- state(coefficients)
    if (is.null(current)) stop('Invalid conditional variance objective.')
    if (max(abs(current$gradient)) / nrow(z) < 1e-10) break
    hessian <- crossprod(zs, zs * current$ratio)
    step <- tryCatch(solve(hessian, current$gradient), error = function(e)
      stop('Singular conditional variance Hessian; no fallback.'))
    decrement <- sum(step * current$gradient)
    rate <- 1
    repeat {
      next_state <- state(coefficients - rate * step)
      if (!is.null(next_state) && next_state$objective <= current$objective -
        1e-4 * rate * decrement + 1e-12 * max(1, abs(current$objective))) break
      rate <- rate / 2
      if (rate < 2^-35) stop('Conditional variance line search failed; no fallback.')
    }
    coefficients <- coefficients - rate * step
  }
  current <- state(coefficients)
  fitted <- exp(drop(zs %*% coefficients))
  if (is.null(current) || max(abs(current$gradient)) / nrow(z) >= 1e-10 ||
    any(!is.finite(fitted)) || any(fitted <= 0))
    stop('Conditional variance quasi-GLM did not converge; no fallback or clipping.')
  beta <- coefficients / scale
  beta[1] <- beta[1] - sum(center[-1] * beta[-1])
  list(beta = beta, diagnostics = c(
    n = nrow(z), rank = mean_fit$rank,
    mean_df = mean_fit$df, max_mean_leverage = max(mean_fit$leverage),
    hc2_squared_residual_min = min(response), hc2_squared_residual_mean = mean(response),
    fitted_variance_min = min(fitted), fitted_variance_max = max(fitted),
    variance_iterations = iter, variance_equation_error = max(abs(current$gradient)) / nrow(z)
  ), warnings = character())
}

cgm_est_nuisance <- function(data, u, source = 'logit', fold_id = NULL) {
  z <- model.matrix(cgm6m_basis, data)
  R <- data$R
  group <- ifelse(R == 0, 'E', ifelse(data$A == 1, 'T', 'C'))
  n <- nrow(data)
  if (is.null(fold_id)) fold_id <- rep(0L, n)
  if (length(fold_id) != n || anyNA(fold_id)) stop('Invalid original-fold vector.')
  predictions <- matrix(NA_real_, n, 7L, dimnames = list(NULL, c('mT', 'mC', 'mE', 'p', 'odds', 'vC', 'vE')))
  fit_diagnostics <- list()
  models <- list()
  warnings <- character()
  for (fold in sort(unique(fold_id))) {
    validation <- which(fold_id == fold)
    training <- if (fold == 0L) seq_len(n) else which(fold_id != fold)
    means <- lapply(c('T', 'C', 'E'), function(g) {
      ii <- training[group[training] == g]
      u$hct_unrestricted_ols(z[ii, , drop = FALSE], data$Y[ii])
    })
    names(means) <- c('T', 'C', 'E')
    vars <- lapply(c('C', 'E'), function(g) {
      ii <- training[group[training] == g]
      cgm_est_variance(z[ii, , drop = FALSE], data$Y[ii], means[[g]], u)
    })
    names(vars) <- c('C', 'E')
    sf <- if (source == 'logit') u$hct_unrestricted_logit(z[training, , drop = FALSE], R[training]) else
      u$hct_unrestricted_calibrate(z[training, , drop = FALSE], R[training])
    zv <- z[validation, , drop = FALSE]
    for (g in c('T', 'C', 'E')) predictions[validation, paste0('m', g)] <- drop(zv %*% means[[g]]$beta)
    # Match the full assembled mixture when stratified-fold rounding changes
    # the training trial/external ratio. This uses counts, never held-out Y.
    prior_factor <- (sum(R) / sum(1 - R)) / (sum(R[training]) / sum(1 - R[training]))
    linear <- drop(zv %*% sf$eta) + log(prior_factor)
    predictions[validation, 'p'] <- plogis(linear)
    predictions[validation, 'odds'] <- exp(linear)
    for (g in c('C', 'E')) predictions[validation, paste0('v', g)] <- exp(drop(zv %*% vars[[g]]$beta))
    for (g in c('C', 'E')) fit_diagnostics[[length(fit_diagnostics) + 1L]] <-
      data.frame(fold = fold, source = g, as.list(vars[[g]]$diagnostics), check.names = FALSE)
    warnings <- c(warnings, sf$warnings, unlist(lapply(vars, '[[', 'warnings')))
    models[[as.character(fold)]] <- list(means = means, variances = vars, source = sf, source_prior_factor = prior_factor)
  }
  if (any(!is.finite(predictions)) || any(predictions[, 'p'] <= 0 | predictions[, 'p'] >= 1) ||
    any(predictions[, c('odds', 'vC', 'vE')] <= 0))
    stop('Nonfinite or nonpositive nuisance predictions; no clipping.')
  list(
    predictions = as.data.frame(predictions), training_diagnostics = do.call(rbind, fit_diagnostics),
    models = models, warnings = unique(warnings)
  )
}

cgm_est_score <- function(data, predictions, e = .5) {
  R <- data$R
  T <- R * data$A
  C <- R * (1 - data$A)
  E <- 1 - R
  a <- predictions
  delta <- a$mC - a$mE
  ratio <- a$vC / a$vE
  lambda <- a$p * (1 - e) / (a$p * (1 - e) + (1 - a$p) * ratio)
  rT <- data$Y - a$mT
  rC <- data$Y - a$mC
  rE <- data$Y - a$mE
  raw_fixed <- R * (a$mT - lambda * a$mC - (1 - lambda) * a$mE) +
    T / e * rT - C * lambda / (1 - e) * rC - E * (1 - lambda) * a$odds * rE
  multiplier <- delta * lambda * (1 - lambda)
  correction_source <- -multiplier * (R - a$p) / (1 - a$p)
  correction_internal_variance <- multiplier * C / (1 - e) * (rC^2 / a$vC - 1)
  correction_external_variance <- -multiplier * E * a$odds * (rE^2 / a$vE - 1)
  raw <- raw_fixed + correction_source + correction_internal_variance + correction_external_variance
  estimate <- sum(raw) / sum(R)
  phi <- (raw - R * estimate) / mean(R)
  if (any(!is.finite(c(raw, phi, lambda))) || any(lambda <= 0 | lambda >= 1) || abs(mean(phi)) > 1e-10)
    stop('Full law-dependent EIF target equation failed.')
  list(
    estimate = estimate, phi = phi, lambda = lambda, ratio = ratio,
    corrections = c(
      fixed_score_estimate = sum(raw_fixed) / sum(R),
      source_correction = sum(correction_source) / sum(R),
      internal_variance_correction = sum(correction_internal_variance) / sum(R),
      external_variance_correction = sum(correction_external_variance) / sum(R)
    ),
    max_absolute_phi = max(abs(phi)), EIF_mean = mean(phi)
  )
}
