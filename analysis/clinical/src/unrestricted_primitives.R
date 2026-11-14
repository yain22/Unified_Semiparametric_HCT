# Exact function extracts; see source_manifest.json. Sourcing performs no I/O.

hct_unrestricted_ols <- function(z, y, weights = rep(1, length(y))) {
  if (nrow(z) <= ncol(z)) stop('Outcome regression has no positive residual degrees of freedom.')
  fit <- lm.wfit(z, y, weights)
  if (fit$rank < ncol(z) || any(!is.finite(fit$coefficients)))
    stop('Rank-deficient outcome regression; no fallback.')
  Q <- qr.Q(fit$qr)[, seq_len(ncol(z)), drop = FALSE]
  h <- rowSums(Q^2)
  if (any(!is.finite(h)) || any(h >= 1 - 1e-10))
    stop('Outcome regression has unit leverage; HC2 unsupported.')
  list(beta = fit$coefficients, leverage = h, rank = fit$rank, df = nrow(z) - fit$rank)
}

hct_unrestricted_logit <- function(g, S, weights = rep(1, length(S))) {
  captured <- character()
  fit <- withCallingHandlers(glm.fit(g, S,
    weights = weights, family = binomial(),
    control = glm.control(epsilon = 1e-10, maxit = 100L)
  ), warning = function(w) {
    captured <<- c(captured, conditionMessage(w))
    invokeRestart('muffleWarning')
  })
  if (!isTRUE(fit$converged) || fit$rank < ncol(g) ||
    any(!is.finite(fit$coefficients)) || isTRUE(fit$boundary) ||
    any(grepl('numerically 0 or 1|did not converge', captured)))
    stop('Source logistic fit failed or reached numerical support boundary; no clipping.')
  list(eta = fit$coefficients, iterations = fit$iter, warnings = unique(captured))
}

# Damped Newton solves sum_i w_i g_i (E_i exp(g_i eta)-S_i)=0.
# Standardization is solely numerical; returned eta and equations use raw g.

hct_unrestricted_group_variance <- function(influence, group) {
  N <- length(influence)
  components <- do.call(rbind, lapply(levels(group), function(a) {
    v <- influence[group == a]
    data.frame(
      group = a, n = length(v), influence_mean = mean(v),
      influence_variance = var(v), variance_contribution = length(v) * var(v) / N^2
    )
  }))
  list(se = sqrt(sum(components$variance_contribution)), components = components)
}
