# qlogis(p) + log(1-e) + log(v_external/v_internal) is the log odds
# of using internal controls. Keep the two log variances separate.
share_probabilities <- function(p, e, v_internal, v_external) {
  args <- list(p, e, v_internal, v_external)
  if (!all(vapply(args, is.numeric, logical(1))) ||
      any(vapply(args, function(x) !is.null(dim(x)), logical(1)))) {
    stop("Inputs must be numeric vectors.", call. = FALSE)
  }
  sizes <- lengths(args)
  n <- max(sizes)
  if (n == 0L || !all(sizes %in% c(1L, n)))
    stop("Each input must have length one or the common vector length.", call. = FALSE)
  p <- rep_len(p, n)
  e <- rep_len(e, n)
  vi <- rep_len(v_internal, n)
  ve <- rep_len(v_external, n)
  if (any(!is.finite(p) | p < 0 | p > 1))
    stop("Source probabilities must be between zero and one.", call. = FALSE)
  trial <- p > 0
  if (any(!is.finite(e[trial]) | e[trial] <= 0 | e[trial] >= 1))
    stop("Treatment probability must be strictly between zero and one on trial support.",
         call. = FALSE)
  both <- p > 0 & p < 1
  if (any(!is.finite(vi[both]) | vi[both] <= 0 |
          !is.finite(ve[both]) | ve[both] <= 0))
    stop("Both conditional control variances must be finite and positive on shared support.",
         call. = FALSE)
  log_ratio <- stats::qlogis(p[both]) + log1p(-e[both]) + log(ve[both]) - log(vi[both])
  inside <- outside <- rep(NA_real_, n)
  inside[p == 1] <- 1
  outside[p == 1] <- 0
  inside[both] <- stats::plogis(log_ratio)
  outside[both] <- stats::plogis(-log_ratio)
  list(internal = inside, external = outside)
}

precision_share <- function(p, e, v_internal, v_external) {
  share_probabilities(p, e, v_internal, v_external)$internal
}
