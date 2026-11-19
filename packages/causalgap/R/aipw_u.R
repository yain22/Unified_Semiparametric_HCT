aipw_u <- function(data, nuisance) {
  if (!is.data.frame(data) || !is.data.frame(nuisance) || nrow(data) == 0L ||
      nrow(data) != nrow(nuisance))
    stop("Supply two data frames with the same, positive number of rows.", call. = FALSE)
  obs <- c("R", "A", "Y")
  fits <- c("p", "e", "m_treated", "m_internal", "m_external", "v_internal", "v_external")
  if (!all(obs %in% names(data)) || !all(fits %in% names(nuisance)))
    stop("Missing required column in data or nuisance. See ?aipw_u.", call. = FALSE)
  cols <- c(as.list(data[obs]), as.list(nuisance[fits]))
  if (!all(vapply(cols, function(x) is.numeric(x) && is.null(dim(x)), logical(1))))
    stop("Each required column must be a numeric vector.", call. = FALSE)
  s <- data$R
  a <- data$A
  y <- data$Y
  if (!all(s %in% 0:1) || !all(a %in% 0:1) || any(a > s) || any(!is.finite(y)))
    stop("R and A must be binary, external records must have A=0, and Y must be finite.",
         call. = FALSE)
  if (sum(s) == 0)
    stop("At least one observation from the trial is required.", call. = FALSE)
  p <- nuisance$p
  e <- nuisance$e
  vi <- nuisance$v_internal
  ve <- nuisance$v_external
  shares <- share_probabilities(p, e, vi, ve)
  li <- shares$internal
  le <- shares$external
  if (any((s == 1 & p == 0) | (s == 0 & p == 1)))
    stop("An observation is inconsistent with its source probability.", call. = FALSE)

  both <- p > 0 & p < 1
  trial <- s == 1
  tt <- trial & a == 1
  ic <- trial & a == 0
  ec <- !trial & both
  mt <- nuisance$m_treated
  mi <- nuisance$m_internal
  me <- nuisance$m_external
  if (any(!is.finite(mt[trial])) || any(!is.finite(mi[trial | both])) ||
      any(!is.finite(me[both])))
    stop("A required conditional mean is missing or nonfinite.", call. = FALSE)

  q <- mean(s)
  terms <- matrix(0, nrow(data), 7L, dimnames = list(NULL, c(
    "standardization", "treated", "internal_control", "external_control",
    "source_probability", "internal_variance", "external_variance")))
  t_trial <- mt[trial] - mi[trial]
  terms[trial, 1L] <- t_trial
  shared_trial <- trial & both
  delta <- mi - me
  terms[shared_trial, 1L] <- terms[shared_trial, 1L] + le[shared_trial] * delta[shared_trial]
  terms[tt, 2L] <- (y[tt] - mt[tt]) / e[tt]
  r_ic <- y[ic] - mi[ic]
  terms[ic, 3L] <- -li[ic] * r_ic / (1 - e[ic])
  odds_ec <- p[ec] / (1 - p[ec])
  r_ec <- y[ec] - me[ec]
  terms[ec, 4L] <- -odds_ec * le[ec] * r_ec

  weight_derivative <- delta * li * le
  # Inside the trial (R-p)/(1-p) equals one, including close to p=1.
  terms[shared_trial, 5L] <- -weight_derivative[shared_trial]
  terms[ec, 5L] <- weight_derivative[ec] * odds_ec
  shared_ic <- ic & both
  terms[shared_ic, 6L] <- weight_derivative[shared_ic] / (1 - e[shared_ic]) *
    ((y[shared_ic] - mi[shared_ic])^2 / vi[shared_ic] - 1)
  terms[ec, 7L] <- -weight_derivative[ec] * odds_ec * (r_ec^2 / ve[ec] - 1)

  trial_score <- numeric(nrow(data))
  trial_score[trial] <- t_trial
  trial_score[tt] <- trial_score[tt] + terms[tt, 2L]
  trial_score[ic] <- trial_score[ic] - r_ic / (1 - e[ic])
  tce <- sum(terms) / sum(s)
  ate <- sum(trial_score) / sum(s)
  terms[, 1L] <- terms[, 1L] - s * tce
  terms <- terms / q
  phi <- data.frame(tce = rowSums(terms), trial = (trial_score - s * ate) / q)
  phi$gap <- phi$tce - phi$trial
  if (any(!is.finite(as.matrix(phi))))
    stop("Nonfinite score. Inspect the nuisance predictions and overlap.", call. = FALSE)
  list(estimate = c(tce = tce, trial = ate, gap = tce - ate),
       influence = phi, components = as.data.frame(terms), lambda = li)
}
