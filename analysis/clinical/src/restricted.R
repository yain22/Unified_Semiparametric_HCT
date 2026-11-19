# Exact function extracts; see source_manifest.json. Sourcing performs no I/O.

hct_eif_formula <- function(rhs, response) {
  if (!inherits(rhs, 'formula') || length(rhs) != 2L)
    stop('Nuisance terms must be one-sided formulas, e.g. ~ X1 + X2.')
  as.formula(call('~', as.name(response), rhs[[2L]]), env = environment(rhs))
}

# Fits m1 on trial-treated outcomes, the COMMON m0 on both control sources,
# and p on all source indicators. Source-specific outcome regressions are not
# substituted for the common-control model required by this r=1 EIF.

hct_eif_nuisance <- function(data, outcome_terms, source_terms = outcome_terms,
                             fold_id = NULL, clip = .01) {
  design <- hct_eif_validate(data)
  if (length(clip) != 1L || !is.finite(clip) || clip <= 0 || clip >= .5)
    stop('clip must be strictly between zero and one half.')
  yf <- hct_eif_formula(outcome_terms, 'Y')
  sf <- hct_eif_formula(source_terms, 'R')
  # Avoid accidentally regressing on the response, treatment, or source.
  if (any(c('Y', 'R', 'A') %in% c(all.vars(outcome_terms), all.vars(source_terms))))
    stop('Nuisance terms must contain only baseline covariates.')
  if (is.null(fold_id)) {
    fold <- rep(1L, nrow(data))
    crossfit <- FALSE
  } else {
    if (length(fold_id) != nrow(data) || anyNA(fold_id) ||
      !is.numeric(fold_id) || any(!is.finite(fold_id)) ||
      length(unique(fold_id)) < 2L) stop('Invalid fold_id.')
    fold <- match(fold_id, sort(unique(fold_id)))
    crossfit <- TRUE
    if (any(table(design$group, fold) == 0L))
      stop('Every fold must contain all three T/C/E groups.')
  }
  m1 <- m0 <- praw <- rep(NA_real_, nrow(data))
  details <- list()
  for (k in sort(unique(fold))) {
    valid <- which(fold == k)
    train <- if (crossfit) which(fold != k) else seq_len(nrow(data))
    Ttrain <- train[data$R[train] == 1 & data$A[train] == 1]
    Ctrain <- train[data$A[train] == 0]
    fit1 <- lm(yf, data = data[Ttrain, , drop = FALSE])
    fit0 <- lm(yf, data = data[Ctrain, , drop = FALSE])
    if (anyNA(coef(fit1)) || anyNA(coef(fit0)))
      stop('Rank-deficient outcome nuisance fit; simplify terms or increase sample size.')
    captured <- character()
    fitp <- withCallingHandlers(
      glm(sf,
        data = data[train, , drop = FALSE], family = binomial(),
        control = glm.control(maxit = 50L)
      ),
      warning = function(w) {
        captured <<- c(captured, conditionMessage(w))

        invokeRestart('muffleWarning')
      }
    )
    if (!isTRUE(fitp$converged) || anyNA(coef(fitp)))
      stop('Source-probability fit did not converge or was rank deficient.')
    m1[valid] <- predict(fit1, newdata = data[valid, , drop = FALSE])
    m0[valid] <- predict(fit0, newdata = data[valid, , drop = FALSE])
    praw[valid] <- predict(fitp, newdata = data[valid, , drop = FALSE], type = 'response')
    details[[length(details) + 1L]] <- data.frame(
      fold = k, n_train = length(train),
      n_valid = length(valid), n_train_T = length(Ttrain), n_train_controls = length(Ctrain),
      source_converged = fitp$converged, source_iterations = fitp$iter,
      source_warning = paste(unique(captured), collapse = ' | ')
    )
  }
  if (any(!is.finite(c(m1, m0, praw)))) stop('Nonfinite nuisance prediction.')
  p <- pmin(1 - clip, pmax(clip, praw))
  list(
    m1 = m1, m0 = m0, p = p, p_raw = praw, fold_id = if (crossfit) fold else NULL,
    diagnostics = list(
      cross_fitted = crossfit, folds = length(unique(fold)),
      clip = clip, clipped_count = sum(p != praw), clipped_fraction = mean(p != praw),
      raw_p_min = min(praw), raw_p_max = max(praw),
      max_control_weight = max(p / (1 - design$e * p)),
      fits = do.call(rbind, details)
    )
  )
}

# The empirical score uses the realized assembled mixture q and known e.
# Centering separately within the design's fixed sampling groups removes their
# random-count components. var(phi)/N is supplied only as an iid-mixture diagnostic.

hct_eif_from_predictions <- function(data, m1, m0, p, e = NULL,
                                     estimate = NULL, alpha = .05) {
  design <- hct_eif_validate(data, e)
  N <- design$N
  e <- design$e
  q <- design$q
  if (any(lengths(list(m1, m0, p)) != N) || any(!is.finite(c(m1, m0, p))) ||
    any(p <= 0 | p >= 1) || length(alpha) != 1L || !is.finite(alpha) ||
    alpha <= 0 || alpha >= 1) stop('Invalid nuisance predictions or alpha.')
  trial <- data$R == 1
  treated <- trial & data$A == 1
  control <- data$A == 0
  h <- p / (1 - e * p)
  residual_term <- treated / (q * e) * (data$Y - m1) - control * h / q * (data$Y - m0)
  plugin <- mean((m1 - m0)[trial])
  one_step <- plugin + mean(residual_term)
  if (is.null(estimate)) estimate <- one_step
  if (length(estimate) != 1L || !is.finite(estimate)) stop('Invalid estimate.')
  phi <- data$R / q * (m1 - m0 - estimate) + residual_term
  components <- do.call(rbind, lapply(levels(design$group), function(g) {
    values <- phi[design$group == g]
    data.frame(
      group = g, n = length(values), EIF_mean = mean(values),
      EIF_variance = var(values), variance_contribution = length(values) * var(values) / N^2
    )
  }))
  se <- sqrt(sum(components$variance_contribution))
  list(
    estimate = estimate, se = se, ci = estimate + qnorm(c(alpha / 2, 1 - alpha / 2)) * se,
    p_value = pnorm(estimate / se, lower.tail = FALSE), EIF = phi, EIF_mean = mean(phi),
    residual_score_mean = mean(residual_term), plugin = plugin, one_step = one_step,
    se_iid_mixture = sqrt(var(phi) / N), variance_components = components,
    predictions = data.frame(m1 = m1, m0 = m0, p = p, delta = m1 - m0),
    e = e, q = q, variance_regime = 'independent_T_C_E_samples_with_fixed_counts',
    validity = 'requires_identification_and_nuisance_conditions; working_model_under_drift'
  )
}
