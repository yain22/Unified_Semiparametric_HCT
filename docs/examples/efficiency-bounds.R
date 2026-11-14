# Reproduce the efficiency-bound explorer in base R.
# One covariate profile, iid mixture sampling and known trial allocation.
# Group means are 1 (treated), 0 (internal control) and -delta (external
# control). Outcomes are Gaussian at the illustrated law. The EIF allows
# the conditional outcome laws and the trial source probability to vary.
#
# Run: source("efficiency-bounds.R")
# Vary the inputs: efficiency_bounds(delta = 0.2, vE = 4, plan = "equal-total")
# nT, nI and nE are expected enrollment counts used to set mixture
# probabilities. B = E[phi^2] is per record. Divide by the design's own N
# for the study variance bound, then take its square root for the SE.

efficiency_bounds <- function(nT = 100, nI = 50, nE = 150,
                              vT = 1, vI = 1, vE = 1, delta = 0,
                              plan = c("retain", "equal-total")) {
  plan <- match.arg(plan)
  inputs <- list(nT, nI, nE, vT, vI, vE, delta)
  finite_scalar <- function(value) {
    is.numeric(value) && length(value) == 1L && is.finite(value)
  }
  stopifnot(all(vapply(inputs, finite_scalar, logical(1))),
            all(c(nT, nI, nE, vT, vI, vE) > 0))

  groups <- c("treated", "internal", "external")
  variances <- setNames(c(vT, vI, vE), groups)
  means <- setNames(c(1, 0, -delta), groups)
  hybrid_counts <- setNames(c(nT, nI, nE), groups)
  rct_counts <- hybrid_counts
  ect_counts <- hybrid_counts
  rct_counts["external"] <- 0
  ect_counts["internal"] <- 0
  if (plan == "equal-total") {
    rct_counts["internal"] <- nI + nE
    ect_counts["external"] <- nI + nE
  }
  designs <- list(unrestricted = hybrid_counts, transportable = hybrid_counts,
                  rct = rct_counts, ect = ect_counts)

  evaluate_design <- function(model) {
    counts <- designs[[model]]
    N <- sum(counts)
    positive <- counts > 0
    row <- data.frame(model = model, nT = counts[[1L]], nI = counts[[2L]],
                      nE = counts[[3L]], N = N, target = NA_real_,
                      B = NA_real_, variance = NA_real_, SE = NA_real_)
    # At different control means the same HCT law is outside this submodel.
    if (model == "transportable" && delta != 0) return(row)

    # Conditional on group, the EIF is quadratic in the standardized outcome.
    # Its square has degree at most four. These three normal quadrature nodes
    # and weights integrate all such polynomials exactly.
    obs <- expand.grid(group = groups[positive], node = 1:3,
                       stringsAsFactors = FALSE)
    z <- c(-sqrt(3), 0, sqrt(3))[obs$node]
    probability <- (counts[obs$group] / N) * c(1/6, 2/3, 1/6)[obs$node]
    residual <- sqrt(variances[obs$group]) * z
    treated <- obs$group == "treated"
    internal <- obs$group == "internal"
    external <- obs$group == "external"
    score <- numeric(nrow(obs))
    score[treated] <- residual[treated] / (counts[["treated"]] / N)
    if (model %in% c("rct", "ect")) {
      # With one covariate profile the trial-standardization term vanishes.
      control <- if (model == "rct") "internal" else "external"
      score[!treated] <- -residual[!treated] / (counts[[control]] / N)
      row$target <- means[["treated"]] - means[[control]]
    } else {
      q <- (nT + nI) / N
      e <- nT / (nT + nI)
      internal_precision <- nI / vI
      external_precision <- nE / vE
      lambda <- internal_precision / (internal_precision + external_precision)
      external_share <- external_precision / (internal_precision + external_precision)
      row$target <- means[["treated"]] -
        (lambda * means[["internal"]] + external_share * means[["external"]])
      score[internal] <- -lambda * residual[internal] / (q * (1 - e))
      score[external] <- -external_share * residual[external] / (1 - q)
      if (model == "unrestricted") {
        # Random source membership and the two estimated control variances
        # contribute even with a single covariate value when delta is nonzero.
        S <- as.integer(!external)
        correction <- -(S - q) / (1 - q)
        correction[internal] <- correction[internal] + (z[internal]^2 - 1) / (1 - e)
        correction[external] <- correction[external] - q * (z[external]^2 - 1) / (1 - q)
        score <- score + delta * lambda * external_share * correction / q
      }
    }
    stopifnot(abs(sum(probability) - 1) < 1e-12,
              abs(sum(probability * score)) < 1e-11 * max(1, max(abs(score))))
    row$B <- sum(probability * score^2)
    row$variance <- row$B / N
    row$SE <- sqrt(row$variance)
    row
  }
  out <- do.call(rbind, lapply(names(designs), evaluate_design))
  rownames(out) <- NULL
  out
}

# Retain each design's groups, then recruit 200 controls in every design.
retained <- efficiency_bounds()
equal_total <- efficiency_bounds(plan = "equal-total")
unequal_variances <- efficiency_bounds(vE = 4, plan = "equal-total")
discrepant <- efficiency_bounds(delta = 0.2)
matches <- function(actual, expected) {
  isTRUE(all.equal(actual, expected, tolerance = 1e-12))
}
stopifnot(
  matches(retained$N, c(300, 300, 150, 250)),
  matches(retained$target, rep(1, 4)),
  matches(retained$variance, c(0.015, 0.015, 0.03, 1/60)),
  matches(equal_total$N, rep(300, 4)),
  matches(equal_total$variance, rep(0.015, 4)),
  matches(unequal_variances$variance, c(0.01 + 1/87.5, 0.01 + 1/87.5, 0.015, 0.03)),
  matches(discrepant$target, c(1.15, NA, 1, 1.2)),
  matches(discrepant$B, c(4.528125, NA, 4.5, 25/6)),
  matches(discrepant$SE, sqrt(discrepant$variance))
)
examples <- list("Retained groups, equal means" = retained,
                 "Equal total enrollment, equal means" = equal_total,
                 "Equal total enrollment, external variance 4" = unequal_variances,
                 "Retained groups, discrepancy 0.2" = discrepant)
for (label in names(examples)) {
  cat("\n", label, "\n", sep = "")
  print(examples[[label]], row.names = FALSE, digits = 8)
}
