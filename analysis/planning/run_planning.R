# Exact Gaussian planning calculations for the main and supplementary figures.
# Uses base R; no random numbers or fitted estimators.
args <- commandArgs(trailingOnly = TRUE)
if (any(args %in% c("--help", "-h"))) {
  cat(paste(
    "Usage: Rscript run_planning.R [--output-dir PATH] [--python-reference-dir PATH]",
    "",
    "Compute Table 1, planning curves, gap thresholds, and exact reference values.",
    "--output-dir            Output folder (default: results/planning-r).",
    "--python-reference-dir  Compare with a complete run_planning.py output folder.",
    "--help                  Show this message.",
    sep = "\n"
  ), "\n")
  quit(status = 0)
}
value <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args) || startsWith(args[i + 1], "--")) stop("Missing value for ", flag)
  args[i + 1]
}
unknown <- args[startsWith(args, "--") & !args %in% c("--output-dir", "--python-reference-dir")]
if (length(unknown)) stop("Unknown option: ", paste(unknown, collapse = ", "))
filearg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script <- normalizePath(sub("^--file=", "", filearg), mustWork = TRUE)
out <- value("--output-dir", file.path(dirname(script), "..", "..", "results", "planning-r"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out, mustWork = TRUE)
reference <- value("--python-reference-dir")
alpha <- .025
theta_star <- 2
r_squared <- .25
gamma <- .10
z <- qnorm(1 - alpha)
scenarios <- data.frame(
  scenario = c("I", "II", "III"), sigma = c(7.5, 10, 10),
  treated = c(120, 120, 240), rct_controls = c(120, 120, 240),
  hct_controls = c(60, 60, 120), external_controls = c(300, 300, 600)
)
prob <- function(gap, se) data.frame(
  coverage = pnorm(z - gap / se) - pnorm(-z - gap / se),
  null_rejection = pnorm(gap / se - z),
  efficacy_declaration = pnorm((theta_star + gap) / se - z)
)
for (directory in c("tables", "curves", "validation")) {
  dir.create(file.path(out, directory), recursive = TRUE, showWarnings = FALSE)
}
put <- function(x, name, directory = "tables") {
  write.csv(x, file.path(out, directory, name), row.names = FALSE, na = "")
}
panel_a <- panel_b <- thresholds <- common <- matched <- baseline <- list()
counter <- 0L
cindex <- 0L
for (i in seq_len(nrow(scenarios))) {
  s <- scenarios[i, ]
  v <- (1 - r_squared) * s$sigma^2
  nt <- s$treated
  nr <- s$rct_controls
  nc <- s$hct_controls
  ne <- s$external_controls
  share <- nc / (nc + ne)
  stopifnot(abs(share - 1 / 6) < 1e-14)
  variances <- c(
    "AIPW-U (HCT)" = v * (1 / nt + 1 / (nc + ne)),
    "AIPW (ECT)" = v * (1 / nt + 1 / ne), "ANOVA (RCT)" = s$sigma^2 * (1 / nt + 1 / nr),
    "ANCOVA (RCT)" = v * (1 / nt + 1 / nr)
  )
  stopifnot(abs(variances[1] - v * (1 / nt + share^2 / nc + (1 - share)^2 / ne)) < 1e-14)
  for (method in names(variances)) {
    counter <- counter + 1L
    variance <- unname(variances[method])
    se <- sqrt(variance)
    borrowing <- !grepl("RCT", method, fixed = TRUE)
    zero <- prob(0, se)
    stopifnot(abs(zero$coverage - .95) < 1e-14, abs(zero$null_rejection - alpha) < 1e-14)
    panel_a[[counter]] <- data.frame(
      scenario = s$scenario, sigma = s$sigma, method = method,
      n_treated = nt, n_internal = if (!borrowing) nr else if (method == "AIPW-U (HCT)") nc else 0,
      n_external = if (borrowing) ne else 0, standard_error = se, power = zero$efficacy_declaration,
      type_I = zero$null_rejection
    )
    endpoints <- prob(if (borrowing) c(-.5, 0, .5) else c(0, 0, 0), se)
    panel_b[[counter]] <- data.frame(
      scenario = s$scenario, method = method,
      coverage_abs_0_5 = endpoints$coverage[3],
      rejection_minus_0_5 = endpoints$null_rejection[1], efficacy_minus_0_5 = endpoints$efficacy_declaration[1],
      rejection_plus_0_5 = endpoints$null_rejection[3], efficacy_plus_0_5 = endpoints$efficacy_declaration[3]
    )
    tol <- if (borrowing) se * (z + qnorm(gamma)) else Inf
    if (borrowing) stopifnot(abs(prob(tol, se)$null_rejection - gamma) < 1e-14)
    external_share <- if (method == "AIPW-U (HCT)") 1 - share else 1
    thresholds[[counter]] <- data.frame(
      scenario = s$scenario, method = method,
      gamma = gamma, xi_tol = tol, within_plotted_range = tol <= .5,
      source_delta_tol = if (borrowing) tol / external_share else Inf
    )
    gap_axis <- seq(-.5, .5, length.out = 2001)
    actual_gap <- if (borrowing) gap_axis else rep(0, length(gap_axis))
    matched[[counter]] <- cbind(data.frame(
      scenario = s$scenario, method = method,
      plotted_gap = gap_axis, actual_gap = actual_gap, standard_error = se
    ), prob(actual_gap, se))
    effect <- seq(-2, 4, length.out = 1201)
    baseline[[counter]] <- data.frame(
      scenario = s$scenario, method = method,
      theta = effect, power = pnorm(effect / se - z)
    )
    if (method == "ANOVA (RCT)") next
    cindex <- cindex + 1L
    label <- switch(method,
      "AIPW-U (HCT)" = "HCT",
      "AIPW (ECT)" = "ECT",
      "ANCOVA (RCT)" = "RCT ANCOVA"
    )
    multiplier <- switch(label,
      HCT = 1 - share,
      ECT = 1,
      "RCT ANCOVA" = 0
    )
    delta <- seq(-.5, .5, length.out = 1001)
    gap <- multiplier * delta
    common[[cindex]] <- cbind(
      s[rep(1, length(delta)), ], data.frame(
        design = label,
        source_discrepancy = delta, sampling_variance = variance, standard_error = se,
        gap_multiplier = multiplier, causal_gap = gap
      ), prob(gap, se),
      data.frame(causal_mse = variance + gap^2, interval_width = 2 * z * se)
    )
  }
}
panel_a <- do.call(rbind, panel_a)
panel_b <- do.call(rbind, panel_b)
common <- do.call(rbind, common)
matched <- do.call(rbind, matched)
baseline <- do.call(rbind, baseline)
thresholds <- do.call(rbind, thresholds)
put(panel_a, "table1_panel_a.csv")
put(panel_b, "table1_panel_b.csv")
put(thresholds, "thresholds.csv")
put(common, "common_discrepancy.csv", "curves")
put(common[common$source_discrepancy %in% c(-.5, 0, .5), ], "common_discrepancy_endpoints.csv")
put(matched, "matched_gap.csv", "curves")
put(baseline, "baseline_power.csv", "curves")
crossings <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i) {
  s <- scenarios[i, ]
  v <- .75 * s$sigma^2
  rct_var <- v * (1 / s$treated + 1 / s$rct_controls)
  vars <- v * c(1 / s$treated + 1 / (s$hct_controls + s$external_controls), 1 / s$treated + 1 / s$external_controls)
  data.frame(
    scenario = s$scenario, design = c("HCT", "ECT"),
    absolute_delta_at_equal_mse_to_rct_ancova = sqrt(rct_var - vars) / c(5 / 6, 1)
  )
}))
put(crossings, "mse_crossings.csv")
# Main Figure 4 and Web Figures 8, 10, 11: unit-variance references.
exact_scalar <- function(label, nt, nc, ne, sigma = 1, grid = seq(-.3, .3, by = .025)) {
  lambda <- nc / (nc + ne)
  variance <- sigma^2 * (1 / nt + 1 / (nc + ne))
  se <- sqrt(variance)
  gap <- if (ne == 0) rep(0, length(grid)) else grid
  data.frame(
    method = label, xi_axis = grid, achieved_gap = gap,
    source_delta = if (ne == 0) rep(NA_real_, length(grid)) else gap / (1 - lambda),
    internal_share = lambda, NT = nt, NC = nc, N0 = ne, sigma = sigma, zero_gap_benchmark = ne == 0,
    variance = variance, coverage_theta = pnorm(z - gap / se) - pnorm(-z - gap / se),
    bias_theta = gap, mse_theta = variance + gap^2, width = 2 * z * se,
    null_rejection = pnorm(gap / se - z), power = pnorm((.4 + gap) / se - z)
  )
}
designs <- list(
  c("RCT: zero-gap reference", 100, 50, 0), c("ECT (100, 0)", 100, 0, 150),
  c("HCT (100, 50)", 100, 50, 150), c("HCT (200, 100)", 200, 100, 150), c("HCT (300, 150)", 300, 150, 150)
)
scalar_main <- do.call(rbind, lapply(designs, function(d)
  exact_scalar(d[1], as.numeric(d[2]), as.numeric(d[3]), as.numeric(d[4]), grid = seq(-.3, .3, length.out = 401))))
put(scalar_main, "scalar_design_sensitivity.csv", "curves")
signed <- scaling <- list()
for (ne in c(150, 300)) {
  signed[[as.character(ne)]] <- do.call(rbind, lapply(designs, function(d) {
    x <- exact_scalar(d[1], as.numeric(d[2]), as.numeric(d[3]), if (as.numeric(d[4]) == 0) 0 else ne)
    x$baseline_external_n <- ne
    x
  }))
  scaling[[as.character(ne)]] <- rbind(
    exact_scalar("Baseline", 100, 50, ne),
    exact_scalar("Higher variance", 100, 50, ne, 1.25),
    exact_scalar("All samples doubled", 200, 100, 2 * ne), exact_scalar("Trial doubled", 200, 100, ne)
  )
  scaling[[as.character(ne)]]$baseline_external_n <- ne
}
signed <- do.call(rbind, signed)
scaling <- do.call(rbind, scaling)
put(signed, "scalar_reference_signed.csv", "curves")
put(scaling, "scalar_reference_scaling.csv", "curves")
k <- c(1, 4, 16, 64)
se <- sqrt(.015 / k)
precision <- data.frame(k = k, se = se, half_width = z * se, xi_fixed = .15, xi_local = .15 / sqrt(k))
stopifnot(max(abs(precision$xi_local / precision$se - .15 / sqrt(.015))) < 1e-14)
put(precision, "precision_scaling.csv")
figure2 <- data.frame(
  delta = c(0, .03, 1.5), theta = 5, lambda = .4, TCE = 5 + .6 * c(0, .03, 1.5),
  gap = .6 * c(0, .03, 1.5), SE = .6 * sqrt(1 / 20 + .4^2 / 10 + .6^2 / 15)
)
figure2$gap_over_SE <- figure2$gap / figure2$SE
figure2$exact_coverage <- pnorm(z - figure2$gap_over_SE) - pnorm(-z - figure2$gap_over_SE)
put(figure2, "figure2_exact_values.csv")
comparisons <- list()
compare <- function(actual, path, keys, numeric_cols) {
  expected <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(nrow(actual) == nrow(expected))
  # Rows have the same documented scenario/design/grid order; check identifiers.
  stopifnot(all(vapply(keys, function(k) identical(as.character(actual[[k]]), as.character(expected[[k]])), logical(1))))
  errors <- vapply(numeric_cols, function(k) {
    stopifnot(identical(is.na(actual[[k]]), is.na(expected[[k]])))
    max(c(0, abs(actual[[k]] - expected[[k]])), na.rm = TRUE)
  }, numeric(1))
  stopifnot(all(is.finite(errors)), all(errors < 1e-12))
  relative_path <- paste(basename(dirname(path)), basename(path), sep = "/")
  data.frame(reference = relative_path, column = names(errors), max_abs_error = errors, passed = TRUE)
}
if (!is.null(reference)) {
  reference <- normalizePath(reference, mustWork = TRUE)
  curves <- file.path(reference, "curves")
  tables <- file.path(reference, "tables")
  comparisons[[1]] <- compare(
    matched, file.path(curves, "matched_gap.csv"), c("scenario", "method"),
    c("plotted_gap", "actual_gap", "standard_error", "coverage", "null_rejection", "efficacy_declaration")
  )
  comparisons[[2]] <- compare(baseline, file.path(curves, "baseline_power.csv"), c("scenario", "method"), c("theta", "power"))
  comparisons[[3]] <- compare(
    common, file.path(curves, "common_discrepancy.csv"), c("scenario", "design"),
    c(
      "source_discrepancy", "sampling_variance", "standard_error", "gap_multiplier", "causal_gap",
      "coverage", "null_rejection", "efficacy_declaration", "causal_mse", "interval_width"
    )
  )
  comparisons[[4]] <- compare(
    crossings, file.path(tables, "mse_crossings.csv"),
    c("scenario", "design"), "absolute_delta_at_equal_mse_to_rct_ancova"
  )
  scalar_columns <- c(
    "xi_axis", "achieved_gap", "source_delta", "internal_share", "NT", "NC", "N0",
    "sigma", "variance", "coverage_theta", "bias_theta", "mse_theta", "width", "null_rejection", "power"
  )
  comparisons[[5]] <- compare(scalar_main, file.path(curves, "scalar_design_sensitivity.csv"), "method", scalar_columns)
  comparisons[[6]] <- compare(signed, file.path(curves, "scalar_reference_signed.csv"), "method", scalar_columns)
  comparisons[[7]] <- compare(scaling, file.path(curves, "scalar_reference_scaling.csv"), "method", scalar_columns)
  comparisons[[8]] <- compare(precision, file.path(tables, "precision_scaling.csv"), "k", names(precision))
  comparisons[[9]] <- compare(figure2, file.path(tables, "targets_and_sampling.csv"), "delta", names(figure2))
  put(do.call(rbind, comparisons), "python_comparison.csv", "validation")
}
put(data.frame(check = c(
  "zero_gap_coverage", "zero_gap_rejection", "pooling_variance_identity",
  "threshold_equation", "python_comparison_requested"
), value = c(.95, .025, 1, 1, !is.null(reference))), "checks.csv", "validation")
writeLines(c(
  "Independent base-R known-input, fixed-weight Gaussian planning calculation.",
  "delta = internal-control mean minus external-control mean; causal gap = (1-lambda)*delta.",
  "Matched-gap axis: HCT delta=6*xi/5; ECT delta=xi; RCT actual gap zero.",
  "Common-discrepancy axis: HCT gap=5*delta/6; ECT gap=delta; RCT gap zero.",
  "RCT uses its own internal-control count, twice the HCT count; ECT has no internal controls.",
  "Reported CSV values are unrounded; format to three decimals only for publication.",
  "RCT positive-threshold tolerance is infinite because only gap zero is attainable.",
  "These calculations do not simulate fitted AIPW-U or its nuisance-estimation variance.",
  paste("R:", R.version.string)
), file.path(out, "validation", "provenance.txt"))
cat("Base-R planning calculations passed; output:", out, "\n")
