# HCT estimating equations, analytic Jacobian, and joint sandwich variance.
# All methods share the fitted source odds and source-specific variance models.
# R-matched uses their fitted variance ratio; R-unit fixes its working ratio at 1.
# HC2 training leverages are treated as design-fixed in the parameter Jacobian.
fg_load <- function(paper) {
  env <- new.env(parent = environment())
  sys.source(file.path(paper, 'lib/scalar_gap_estimated_share_helpers.R'), envir = env)
  list(u = env$est_share_load(paper)$u, variance = env$est_share_variance, paper = paper)
}
fg_prepare <- function(data, helpers, mean_basis = 'linear', variance_basis = 'constant') {
  stopifnot(mean_basis %in% c('linear', 'quadratic'), variance_basis %in% c('constant', 'linear', 'x1'))
  X <- as.matrix(data[paste0('X', 1:10)])
  z <- cbind('(Intercept)' = 1, X)
  if (mean_basis == 'quadratic') z <- cbind(z, X^2)
  g <- cbind('(Intercept)' = 1, X)
  v <- if (variance_basis == 'constant') matrix(1, nrow(data), 1, dimnames = list(NULL, '(Intercept)')) else if (variance_basis == 'x1') g[, 1:2, drop = FALSE] else g
  S <- data$S
  T <- S * data$A
  C <- S * (1 - data$A)
  E <- 1 - S
  P <- C + E
  group <- factor(ifelse(E == 1, 'E', ifelse(T == 1, 'T', 'C')), levels = c('T', 'C', 'E'))
  stopifnot(all(S %in% 0:1), all(data$A %in% 0:1), all(data$A[E == 1] == 0), all(table(group) > 1), all(is.finite(data$Y)))
  dims <- c(bT = ncol(z), bC = ncol(z), bE = ncol(z), eta = ncol(g), gC = ncol(v), gE = ncol(v), bP = ncol(z), theta = 1L)
  ends <- cumsum(dims)
  starts <- ends - dims + 1L
  ix <- lapply(seq_along(dims), function(j) seq.int(starts[j], ends[j]))
  names(ix) <- names(dims)
  parameters <- numeric(sum(dims))
  leverage <- numeric(nrow(data))
  fits <- list()
  for (n in c('T', 'C', 'E', 'P')) {
    indicator <- switch(n,
      T = T,
      C = C,
      E = E,
      P = P
    )
    ii <- which(indicator == 1)
    fit <- helpers$u$hct_unrestricted_ols(z[ii, , drop = FALSE], data$Y[ii])
    parameters[ix[[paste0('b', n)]]] <- fit$beta
    fits[[n]] <- fit
    if (n != 'P') leverage[ii] <- fit$leverage
  }
  source_fit <- helpers$u$hct_unrestricted_logit(g, S)
  parameters[ix$eta] <- source_fit$eta
  variance_fits <- list()
  for (n in c('C', 'E')) {
    ii <- which(group == n)
    residual <- data$Y[ii] - drop(z[ii, , drop = FALSE] %*% fits[[n]]$beta)
    fit <- helpers$variance(v[ii, , drop = FALSE], residual, leverage[ii])
    parameters[ix[[paste0('g', n)]]] <- fit$beta
    variance_fits[[n]] <- fit
  }
  list(
    data = data, z = z, g = g, v = v, ix = ix, parameters = parameters, leverage = leverage,
    S = S, T = T, C = C, E = E, P = P, group = group, e = sum(T) / sum(S), q = mean(S), N = nrow(data), nR = sum(S),
    mean_basis = mean_basis, variance_basis = variance_basis, source_warnings = source_fit$warnings,
    variance_fits = variance_fits
  )
}
fg_system <- function(d, parameters, method = 'U', jacobian = TRUE) {
  stopifnot(method %in% c('U', 'R-matched', 'R-unit'), length(parameters) == length(d$parameters))
  ix <- d$ix
  z <- d$z
  g <- d$g
  v <- d$v
  Y <- d$data$Y
  S <- d$S
  T <- d$T
  C <- d$C
  E <- d$E
  e <- d$e
  a <- 1 - e
  N <- d$N
  mT <- drop(z %*% parameters[ix$bT])
  mC <- drop(z %*% parameters[ix$bC])
  mE <- drop(z %*% parameters[ix$bE])
  mP <- drop(z %*% parameters[ix$bP])
  eta <- drop(g %*% parameters[ix$eta])
  w <- exp(eta)
  p <- plogis(eta)
  vC <- exp(drop(v %*% parameters[ix$gC]))
  vE <- exp(drop(v %*% parameters[ix$gE]))
  r <- if (method == 'R-unit') rep(1, N) else vC / vE
  L <- w * a / (w * a + r)
  H <- L * (1 - L)
  D <- mC - mE
  rT <- Y - mT
  rC <- Y - mC
  rE <- Y - mE
  rP <- Y - mP
  RC <- rC^2 / vC
  RE <- rE^2 / vE
  if (method == 'U') {
    B <- -S + C / a * (RC - 1) + E * w * (2 - RE)
    raw <- S * (mT - L * mC - (1 - L) * mE) + T / e * rT - C * L / a * rC - E * (1 - L) * w * rE + D * H * B
  } else {
    aC <- w / (w * a + r)
    bE <- w * r / (w * a + r)
    K <- C * aC + E * bE
    raw <- S * (mT - mP) + T / e * rT - K * rP
  }
  psi <- matrix(0, N, length(parameters))
  for (n in c('T', 'C', 'E', 'P')) {
    residual <- switch(n,
      T = rT,
      C = rC,
      E = rE,
      P = rP
    )
    psi[, ix[[paste0('b', n)]]] <- z * (d[[n]] * residual)
  }
  psi[, ix$eta] <- g * (S - p)
  psi[, ix$gC] <- v * (C * (RC / (1 - d$leverage) - 1))
  psi[, ix$gE] <- v * (E * (RE / (1 - d$leverage) - 1))
  psi[, ix$theta] <- raw - S * parameters[ix$theta]
  if (any(!is.finite(c(psi, L, vC, vE))) || any(L <= 0 | L >= 1)) stop('Invalid fitted score or nuisance prediction.')
  out <- list(
    equations = psi, mean = colMeans(psi), raw = raw, L = L, p = p, w = w, vC = vC, vE = vE,
    mT = mT, mC = mC, mE = mE, mP = mP
  )
  if (jacobian) {
    J <- matrix(0, length(parameters), length(parameters))
    for (n in c('T', 'C', 'E', 'P')) {
      ii <- ix[[paste0('b', n)]]
      J[ii, ii] <- -crossprod(z, z * d[[n]]) / N
    }
    J[ix$eta, ix$eta] <- -crossprod(g, g * (p * (1 - p))) / N
    J[ix$gC, ix$gC] <- -crossprod(v, v * (C * RC / (1 - d$leverage))) / N
    J[ix$gE, ix$gE] <- -crossprod(v, v * (E * RE / (1 - d$leverage))) / N
    J[ix$gC, ix$bC] <- -2 * crossprod(v, z * (C * rC / ((1 - d$leverage) * vC))) / N
    J[ix$gE, ix$bE] <- -2 * crossprod(v, z * (E * rE / ((1 - d$leverage) * vE))) / N
    derivative <- list(bT = S - T / e, bC = rep(0, N), bE = rep(0, N), bP = rep(0, N), eta = rep(0, N), gC = rep(0, N), gE = rep(0, N))
    if (method == 'U') {
      FL <- -S * D - C / a * rC + E * w * rE
      derivative$bC <- L * (-S + C / a) + H * B - 2 * D * H * C / a * rC / vC
      derivative$bE <- (1 - L) * (-S + E * w) - H * B + 2 * D * H * E * w * rE / vE
      derivative$eta <- H * FL - E * (1 - L) * w * rE + D * H * ((1 - 2 * L) * B + E * w * (2 - RE))
      derivative$gC <- -H * FL + D * H * (-(1 - 2 * L) * B - C / a * RC)
      derivative$gE <- H * FL + D * H * ((1 - 2 * L) * B + E * w * RE)
    } else {
      derivative$bP <- -S + K
      derivative$eta <- -K * (1 - L) * rP
      if (method == 'R-matched') {
        derivative$gC <- (C * aC * (1 - L) - E * bE * L) * rP
        derivative$gE <- -derivative$gC
      }
    }
    for (n in names(derivative)) {
      feature <- if (n %in% c('bT', 'bC', 'bE', 'bP')) z else if (n == 'eta') g else v
      J[ix$theta, ix[[n]]] <- colMeans(feature * derivative[[n]])
    }
    J[ix$theta, ix$theta] <- -d$q
    out$jacobian <- J
  }
  out
}
fg_state <- function(d, parameters, method, helpers, keep = FALSE) {
  theta <- d$ix$theta
  nu <- setdiff(seq_along(parameters), theta)
  parameters[theta] <- sum(fg_system(d, parameters, method, FALSE)$raw) / d$nR
  state <- fg_system(d, parameters, method)
  adjustment <- solve(t(state$jacobian[nu, nu, drop = FALSE]), state$jacobian[theta, nu])
  influence <- drop((state$equations[, theta] - state$equations[, nu, drop = FALSE] %*% adjustment) / d$q)
  working <- state$equations[, theta] / d$q
  se <- helpers$u$hct_unrestricted_group_variance(influence, d$group)$se
  working_se <- helpers$u$hct_unrestricted_group_variance(working, d$group)$se
  if (max(abs(state$mean)) > 1e-7 || any(!is.finite(c(se, working_se))) || min(se, working_se) <= 0) stop('Joint-M equation/variance check failed.')
  answer <- list(
    estimate = parameters[theta], se = se, working_se = working_se,
    influence = influence, working_influence = working, lambda = state$L, group = d$group, N = d$N, nR = d$nR,
    nuisance_equation_max = max(abs(state$mean[nu])), EIF_mean = mean(working), warnings = d$source_warnings,
    max_leverage = max(d$leverage), source_probability_range = range(state$p), variance_ratio_range = range(state$vC / state$vE),
    method = method, mean_basis = d$mean_basis, variance_basis = d$variance_basis
  )
  if (keep) {
    answer$prepared <- d
    answer$parameters <- parameters
    answer$state <- state
  }
  answer
}
fg_fit <- function(
  data, helpers, mean_basis = 'linear', variance_basis = 'constant', keep = FALSE,
  methods = c('U', 'R-matched', 'R-unit')
) {
  d <- fg_prepare(data, helpers, mean_basis, variance_basis)
  setNames(lapply(methods, function(method) fg_state(d, d$parameters, method, helpers, keep)), methods)
}
fg_path <- function(fit0, fit1) {
  stopifnot(identical(fit0$group, fit1$group), fit0$N == fit1$N)
  cv <- function(x, y) sum(vapply(levels(fit0$group), function(gr) {
    ii <- fit0$group == gr
    sum(ii) * cov(x[ii], y[ii]) / fit0$N^2
  }, numeric(1)))
  poly <- function(x, y) {
    s <- y - x
    c(v0 = cv(x, x), v1 = 2 * cv(x, s), v2 = cv(s, s))
  }
  list(
    intercept = fit0$estimate, slope = fit1$estimate - fit0$estimate,
    variance = poly(fit0$influence, fit1$influence), working_variance = poly(fit0$working_influence, fit1$working_influence),
    influence0 = fit0$influence, influence_slope = fit1$influence - fit0$influence,
    working0 = fit0$working_influence, working_slope = fit1$working_influence - fit0$working_influence
  )
}
