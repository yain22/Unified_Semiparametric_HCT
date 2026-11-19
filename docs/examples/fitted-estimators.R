# Six hypothetical observations from Figure S1. Run with base R.
# Set the external outcome shift to 1 to repeat the example with higher
# observed external-control outcomes. X is the first baseline covariate.
records <- data.frame(
  id = sprintf("P%02d", 1:6),
  R = c(1, 1, 1, 1, 0, 0),
  A = c(1, 1, 0, 0, 0, 0),
  X = c(0.20, 0.70, 0.30, 0.80, 0.25, 0.65),
  Y = c(8.4, 10.1, 7.2, 8.5, 6.9, 8.0)
)

worked_estimates <- function(external_shift = 0, data = records) {
  data$Y[data$R == 0] <- data$Y[data$R == 0] + external_shift
  with(data, {
    it <- R * A
    ic <- R * (1 - A)
    ec <- 1 - R
    n_trial <- sum(R)
    # Fixed illustrative predictions for the unrestricted and ECT scores.
    e <- 0.5
    p <- 0.55 + 0.15 * X
    m1 <- 8 + 2 * X
    m10 <- 7 + 1.5 * X
    m0 <- 6.5 + 1.5 * X
    v10 <- 1
    v0 <- 1.5
    rho <- v10 / v0
    d <- p * (1 - e) + (1 - p) * rho
    lambda <- p * (1 - e) / d
    odds <- p / (1 - p)
    delta <- m10 - m0
    correction <- delta * lambda * (1 - lambda)
    hct_u <- cbind(
      standardize = R * (m1 - lambda * m10 - (1 - lambda) * m0),
      treated = it / e * (Y - m1),
      internal = -ic * lambda / (1 - e) * (Y - m10),
      external = -ec * odds * (1 - lambda) * (Y - m0),
      source = -correction * (R - p) / (1 - p),
      internal_var = correction * ic / (1 - e) * ((Y - m10)^2 / v10 - 1),
      external_var = -correction * ec * odds * ((Y - m0)^2 / v0 - 1)
    )

    # Refit Y on X across both control sources for every displayed scenario.
    pooled <- lm(Y ~ X, data = data, subset = A == 0)
    mp <- predict(pooled, newdata = data)
    hct_r <- cbind(
      standardize = R * (m1 - mp),
      treated = it / e * (Y - m1),
      pooled = -p * (ic + ec * rho) / d * (Y - mp)
    )

    # The ECT observes the treated trial subset and the external controls.
    # Its source odds are e times the HCT source odds under this restriction.
    keep <- it == 1 | ec == 1
    ect <- cbind(treated = R * (Y - m0), external = -ec * e * odds * (Y - m0))[keep, ]

    # The A row of the least-squares operator gives each participant's weight.
    trial <- subset(data, R == 1)
    trial_contributions <- function(formula) {
      fit <- lm(formula, data = trial)
      D <- model.matrix(fit)
      weights <- solve(crossprod(D), t(D))["A", ]
      contribution <- weights * trial$Y
      stopifnot(isTRUE(all.equal(sum(contribution), unname(coef(fit)["A"]))))
      matrix(contribution, ncol = 1, dimnames = list(trial$id, "OLS weight x outcome"))
    }

    rownames(hct_u) <- rownames(hct_r) <- id
    rownames(ect) <- id[keep]
    components <- list(AIPW_U = hct_u, AIPW_R = hct_r, ECT = ect,
                       ANOVA = trial_contributions(Y ~ A), ANCOVA = trial_contributions(Y ~ A + X))
    denominator <- c(AIPW_U = n_trial, AIPW_R = n_trial, ECT = sum(it), ANOVA = 1, ANCOVA = 1)
    numerator <- vapply(components, sum, numeric(1))
    list(summary = data.frame(numerator, denominator, estimate = numerator / denominator),
         components = components)
  })
}

if (sys.nframe() == 0L) {
  for (shift in c(0, 1)) {
    cat("\nExternal outcome shift:", shift, "\n")
    fit <- worked_estimates(shift)
    print(round(fit$summary, 5))
    cat("AIPW-U contributions before dividing by the trial sample size:\n")
    print(round(fit$components$AIPW_U, 5))
  }
}
