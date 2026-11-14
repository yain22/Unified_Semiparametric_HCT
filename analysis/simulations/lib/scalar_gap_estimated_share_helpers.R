# Positive control-variance fits from HC2-adjusted squared residuals.

est_share_import <- function(path, names) {
  env <- new.env(parent = .GlobalEnv)
  for (ex in parse(path)) if (is.call(ex) && identical(ex[[1]], as.name('<-')) &&
    is.symbol(ex[[2]]) && as.character(ex[[2]]) %in% names) eval(ex, envir = env)
  stopifnot(all(vapply(names, exists, logical(1), envir = env, inherits = FALSE)))
  env
}

est_share_load <- function(paper) {
  list(u = est_share_import(
    file.path(paper, 'lib/hct_unrestricted_comparators.R'),
    c('hct_unrestricted_ols', 'hct_unrestricted_logit', 'hct_unrestricted_group_variance')
  ))
}

est_share_variance <- function(g, residual, leverage) {
  response <- residual^2 / (1 - leverage)
  center <- colMeans(g)
  scale <- apply(g, 2, sd)
  center[1] <- 0
  scale[1] <- 1
  if (any(!is.finite(response)) || any(response < 0) || mean(response) <= 0 || any(scale <= 0)) stop('Invalid variance training data.')
  gs <- sweep(sweep(g, 2, center, '-'), 2, scale, '/')
  coefficients <- c(log(mean(response)), rep(0, ncol(g) - 1L))
  state <- function(b) {
    eta <- drop(gs %*% b)
    ratio <- response * exp(-eta)
    if (any(!is.finite(c(eta, ratio)))) return(NULL)
    list(objective = sum(eta + ratio), gradient = drop(crossprod(gs, 1 - ratio)), ratio = ratio)
  }
  for (iter in seq_len(100L)) {
    cur <- state(coefficients)
    if (is.null(cur)) stop('Nonfinite conditional variance objective.')
    if (max(abs(cur$gradient)) / nrow(g) < 1e-10) break
    step <- solve(crossprod(gs, gs * cur$ratio), cur$gradient)
    decrement <- sum(step * cur$gradient)
    rate <- 1
    repeat {
      nxt <- state(coefficients - rate * step)
      if (!is.null(nxt) && nxt$objective <= cur$objective - 1e-4 * rate * decrement + 1e-12 * max(1, abs(cur$objective))) break
      rate <- rate / 2
      if (rate < 2^-35) stop('Variance Newton line search failed; no fallback.')
    }
    coefficients <- coefficients - rate * step
  }
  cur <- state(coefficients)
  if (is.null(cur) || max(abs(cur$gradient)) / nrow(g) >= 1e-10) stop('Variance Newton fit did not converge.')
  beta <- coefficients / scale
  beta[1] <- beta[1] - sum(center[-1] * beta[-1])
  list(beta = beta, iterations = iter, equation_error = max(abs(cur$gradient)) / nrow(g))
}
