figure_font <- Sys.getenv("HCT_FIGURE_FONT", "sans")
# Local notation revision; numerical specifications and saved records unchanged.

# Analytical illustrations for the single-document HCT EIF manuscript.
# The surface figure adds seeded synthetic scatter, not Monte Carlo results.
# No manuscript fragments or LaTeX generation.
args0 <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args0), mustWork = TRUE))
cli <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% cli) {
  cat("Rscript SCRIPT.R [--output-dir PATH]\n")
  quit(status = 0)
}
if (length(cli) && (length(cli) != 2L || cli[[1]] != "--output-dir")) stop("Use --output-dir PATH")
output <- if (length(cli)) cli[[2]] else file.path(script_dir, "../../results/illustrations")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
fig_dir <- normalizePath(output, winslash = "/", mustWork = TRUE)
stopifnot(dir.exists(fig_dir), capabilities("cairo"))
dir.create(file.path(fig_dir, "notation_review", "unused"), recursive = TRUE, showWarnings = FALSE)

# Only affected included figures replace manuscript assets. Other deterministic
# examples are drawn to an audit folder, preserving untouched manuscript PDFs.
notation_targets <- c(
  "eif_information.pdf", "eif_source_discrepancy.pdf",
  "tce_weighting.pdf", "external-borrowing-anchor.pdf"
)
figure_output <- function(name) {
  if (name %in% notation_targets) file.path(fig_dir, name)
  else file.path(fig_dir, "notation_review", "unused", name)
}

r_values <- c(.25, 1, 4)
cols <- c("#1D4ED8", "#059669", "#B45309")
types <- c(3, 1, 2)
lambda <- function(p, e, r) p * (1 - e) / (p * (1 - e) + (1 - p) * r) # internal share
se_hct <- function(nE, r) sqrt(1 / 100 + 1 / (50 + nE * r))
# Check normalization and variance minimization on an interior grid.
g <- expand.grid(p = c(.1, .5, .9), e = c(.3, 2 / 3), r = r_values)
d <- with(g, p * (1 - e) + (1 - p) * r)
a <- g$p / d
b <- g$p * g$r / d
stopifnot(max(abs(g$p * (1 - g$e) * a + (1 - g$p) * b - g$p)) < 1e-12)
stopifnot(max(abs(g$p * (1 - g$e) * a^2 + (1 - g$p) / g$r * b^2 - g$p^2 / d)) < 1e-12)
stopifnot(abs(se_hct(150, 1) - sqrt(.015)) < 1e-12)
# Independent finite-mixture verification of the exact drift identity.
p0 <- .5
e0 <- 2 / 3
q0 <- .5
pb <- .6
eb <- .7
rb <- .8
mT <- 2
m0 <- 1
gT <- 2.3
g0 <- 1.2
db <- pb * (1 - eb) + (1 - pb) * rb
expectation <- (p0 * (gT - g0 - (mT - m0)) +
  p0 * e0 / eb * (mT - gT) -
  pb / db * (p0 * (1 - e0) + (1 - p0) * rb) * (m0 - g0)) / q0
drift <- (p0 * (eb - e0) / eb * (gT - mT) +
  (p0 * pb * (eb - e0) + (pb - p0) * rb) / db * (g0 - m0)) / q0
stopifnot(abs(expectation - drift) < 1e-12)
# Verify the closed-form targeting normal equations on fixed observations.
y <- c(1, 2, 4, 2, 3, 4, 1, 5)
it <- c(1, 1, 1, 0, 0, 0, 0, 0)
ic <- c(0, 0, 0, 1, 1, 0, 0, 0)
ie <- 1 - it - ic
pp <- seq(.25, .75, length.out = 8)
ee <- rep(2 / 3, 8)
rr <- seq(.5, 1.5, length.out = 8)
dd <- pp * (1 - ee) + (1 - pp) * rr
ht <- 1 / ee
h0 <- pp / dd
cc <- ic + ie * rr
initT <- rep(1.7, 8)
init0 <- rep(2.2, 8)
et <- sum(it * ht * (y - initT)) / sum(it * ht^2)
e0u <- sum(cc * h0 * (y - init0)) / sum(cc * h0^2)
stopifnot(abs(sum(it * ht * (y - initT - et * ht))) < 1e-12)
stopifnot(abs(sum(cc * h0 * (y - init0 - e0u * h0))) < 1e-12)
open_fig <- function(name, square = FALSE) {
  # The discrepancy panels have in-panel legends, so their axes can use the
  # bottom strip that the other figures reserve for a common legend.
  cairo_pdf(figure_output(name), width = 10.4, height = if (square) 5.53 else 4.9,
            family = figure_font, pointsize = 11)
  par(
    mfrow = c(1, 2), mar = c(4.3, 4.6, 2.3, .8),
    oma = c(if (square) .25 else 2.6, .1, .1, .1), pty = if (square) "s" else "m",
    las = 1, tcl = -.25, mgp = c(2.6, .7, 0), xaxs = "i", yaxs = "i"
  )
}
panel <- function(xlim, ylim, xlab, ylab, main, xaxt = "s") {
  plot(NA,
    xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, xaxt = xaxt
  )
  mtext(main, side = 3, line = .65, adj = 0, font = 2, cex = 1)
  grid(col = "#E5E7EB", lty = 1)
  box(col = "#555555")
}
open_fig("eif_information.pdf")
px <- seq(0, 1, length.out = 401)
panel(
  c(0, 1), c(0, 1), expression(p(x)), expression(lambda[r](x)),
  "(a) Internal control information"
)
for (j in seq_along(r_values)) lines(px, lambda(px, 2 / 3, r_values[j]), col = cols[j], lty = types[j], lwd = 2.2)
ne <- seq(0, 400, by = 2)
panel(
  c(0, 400), c(.095, .18), "Number of external controls", "Standard error",
  "(b) Precision at fixed trial size"
)
abline(h = sqrt(.03), lty = 4, col = "#444444", lwd = 1.6)
for (j in seq_along(r_values)) lines(ne, se_hct(ne, r_values[j]), col = cols[j], lty = types[j], lwd = 2.2)
par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
plot.new()
legend("bottom",
  legend = c("r = 0.25", "r = 1", "r = 4", "RCT Only"),
  col = c(cols, "#444444"), lty = c(types, 4), lwd = 2, bty = "n", horiz = TRUE, cex = .95, inset = .012
)
dev.off()
open_fig("eif_source_discrepancy.pdf", square = TRUE)
gap <- seq(-.5, .5, length.out = 401)
panel(
  c(-.5, .5), c(-.5, .5), expression(delta == m["1,0"] - m[0]), expression("Causal gap " * xi[rho]),
  "(a) Conditional source discrepancy"
)
abline(h = 0, v = 0, col = "#777777", lty = 3)
for (j in seq_along(r_values)) lines(gap, 150 * r_values[j] / (50 + 150 * r_values[j]) * gap,
  col = cols[j], lty = types[j], lwd = 2.2
)
legend("topleft",
  legend = expression(rho == 0.25, rho == 1, rho == 4),
  col = cols, lty = types, lwd = 2, bty = "n", cex = .85
)
panel(
  c(0, 400), c(0, .22), "Number of external controls", "Outcome units",
  "(b) Fixed discrepancy: \u03b4 = 0.2"
)
lines(ne, .2 * ne / (50 + ne), col = "#B45309", lwd = 2.2)
lines(ne, se_hct(ne, 1), col = "#1D4ED8", lwd = 2.2, lty = 2)
legend("topright",
  legend = c("Causal gap", "Standard error"),
  col = c("#B45309", "#1D4ED8"), lty = c(1, 2), lwd = 2, bty = "n", cex = .85
)
dev.off()
cat("Verified weight normalization, optimal variance, exact drift, and targeting equations.\n")
cat("Generated information and source-discrepancy figures; no simulation was run.\n")

# Normal-limit consequences of the standardized causal gap.
kappa <- seq(-3, 3, length.out = 601)
coverage <- function(k) pnorm(qnorm(.975) - k) - pnorm(-qnorm(.975) - k)
rejection <- function(k) 1 - pnorm(qnorm(.95) - k)
stopifnot(abs(coverage(0) - .95) < 1e-12, abs(rejection(0) - .05) < 1e-12)
stopifnot(max(abs(coverage(kappa) - coverage(-kappa))) < 1e-12)
qmix <- 9 / 14
px0 <- qmix * .5 / (qmix * .5 + (1 - qmix) * .1)
px1 <- qmix * .5 / (qmix * .5 + (1 - qmix) * .9)
stopifnot(max(abs(lambda(c(px0, px1), 2 / 3, 1) - c(.75, .25))) < 1e-12)
stopifnot(abs(mean(c(.25, .75) * c(.2, -.2)) + .05) < 1e-12)
stopifnot(rejection(1) > rejection(0), rejection(-1) < rejection(0))
open_fig("eif_gap_inference.pdf")
panel(
  c(-3, 3), c(0, 1), expression(kappa), expression("Limiting Cov. (" * theta * ")"),
  "(a) Nominal 95% causal interval"
)
abline(h = .95, lty = 3, col = "#555555", lwd = 1.5)
lines(kappa, coverage(kappa), col = "#1D4ED8", lwd = 2.2)
abline(v = 0, lty = 3, col = "#AAAAAA")
panel(
  c(-3, 3), c(0, 1), expression(kappa), "Limiting rejection probability",
  "(b) One-sided 5% superiority test"
)
abline(h = .05, lty = 3, col = "#555555", lwd = 1.5)
lines(kappa, rejection(kappa), col = "#B45309", lwd = 2.2)
abline(v = 0, lty = 3, col = "#AAAAAA")
dev.off()
cat("Verified heterogeneous-gap example and local coverage/rejection formulas.\n")

# TCE is the contrast built from a transported control surface; delta and xi
# follow Lee (2026): delta = m_{1,0} - m_0, xi = TCE - trial causal effect.
# The surface figure uses synthetic observations around stipulated mean curves.
# Remaining markers identify analytical quantities, not patient data.
col_t <- "#22B94C"
col_c <- "#2563B8"
col_e <- "#7427BA"
col_h <- "#D77A10"
# Pair the same covariates and residuals across the two transportability
# conditions. Only the external mean changes between them.
surface_points <- local({
  set.seed(20260926)
  lapply(c("T", "C", "E"), function(group)
    data.frame(x = sort(runif(18, .035, .965)), error = rnorm(18, sd = .6)))
})
draw_surface_points <- function(meanT, meanC, meanE, cex = .62) {
  for (k in 1:3) {
    z <- surface_points[[k]]
    mu <- list(meanT, meanC, meanE)[[k]]
    points(z$x, mu(z$x) + z$error,
      pch = c(17, 15, 16)[k],
      col = adjustcolor(c(col_t, col_c, col_e)[k], alpha.f = .50), cex = cex
    )
  }
}
open_tce <- function(name, height = 5.35, bottom = 3.7) {
  cairo_pdf(figure_output(name),
    width = 10.4, height = height,
    family = figure_font, pointsize = 11
  )
  par(
    mfrow = c(1, 2), mar = c(4.3, 4.6, 2.7, 1), oma = c(bottom, .1, .1, .1),
    las = 1, tcl = -.25, mgp = c(2.6, .7, 0), xaxs = "i", yaxs = "i"
  )
}
overlay_legend <- function(labels, colors, ltys = NULL, pchs = NULL, cex = .86, ncol = 2) {
  par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
  plot.new()
  legend("bottom",
    legend = labels, col = colors, lty = ltys, pch = pchs,
    lwd = 2, bty = "n", ncol = ncol, cex = cex, inset = .024, x.intersp = .75, y.intersp = 1.2
  )
}
theta_demo <- 5
delta_demo <- 1.5
share_demo <- .4 # internal-control share
stopifnot(
  abs(theta_demo + (1 - share_demo) * delta_demo - 5.9) < 1e-12,
  abs(theta_demo + delta_demo - 6.5) < 1e-12,
  abs((1 - share_demo) * delta_demo - .9) < 1e-12
)
open_tce("tce_surfaces.pdf", height = 5.6, bottom = 3.4)
xx <- seq(0, 1, length.out = 301)
meanC <- function(x) 18 - 2 * x
meanT <- function(x) meanC(x) + theta_demo
for (j in 1:2) {
  shift <- if (j == 1) 0 else delta_demo
  meanE <- function(x) meanC(x) - shift
  meanH <- function(x) meanC(x) - (1 - share_demo) * shift
  panel(
    c(0, 1), c(13, 24.5), "Baseline covariate x", "Outcome",
    if (j == 1) "(a) Mean transportability holds" else
      "(b) Mean transportability is violated"
  )
  draw_surface_points(meanT, meanC, meanE, cex = .72)
  lines(xx, meanT(xx), col = col_t, lwd = 2)
  lines(xx, meanE(xx), col = col_e, lwd = 2)
  lines(xx, meanC(xx), col = col_c, lwd = 2)
  lines(xx, meanH(xx), col = col_h, lwd = 2, lty = 2)
  # Constant treatment effect and constant source shift make these vertical
  # differences identical to their trial-population averages.
  arrows(.19, meanC(.19), .19, meanT(.19), code = 3, length = .07, lwd = 1.2)
  text(.22, meanC(.19) + 2.6, expression(theta == 5), adj = 0, cex = .95)
  if (j == 1) {
    text(.50, 19.2, expression(tilde(theta)[lambda] == tilde(theta)[lambda == 0] ~ ~"=" ~ ~theta), cex = .93)
    text(.52, 13.7, expression(delta(x) == 0 ~ ~"," ~ ~ xi[lambda] == 0), cex = .92)
  } else {
    arrows(.58, meanH(.58), .58, meanT(.58), code = 3, length = .07, lwd = 1.2)
    text(.61, meanH(.58) + 3.1, expression(tilde(theta)[lambda] == 5.9), adj = 0, cex = .92)
    arrows(.87, meanH(.87), .87, meanC(.87), code = 3, length = .06, lwd = 1.2)
    text(.83, meanC(.87) + 1.1, expression(xi[lambda] == 0.9), adj = 1, cex = .9)
    text(.055, 13.7, expression(delta(x) == 1.5 ~ ~"," ~ ~ tilde(theta)[lambda == 0] == 6.5),
      adj = 0, cex = .9
    )
  }
}
overlay_legend(
  c(
    "Trial treated mean", "Internal control mean (observed)",
    "External control mean (observed)", "Hybrid transported control"
  ),
  c(col_t, col_c, col_e, col_h),
  ltys = c(1, 1, 1, 2),
  pchs = c(17, 15, 16, NA), cex = .88
)
dev.off()

# A constant-share path and a heterogeneous example distinguish the two gaps.
stratum_delta <- c(.2, -.2)
stratum_share <- c(.75, .25) # internal-control shares
trial_mass <- c(.5, .5)
raw_contribution <- trial_mass * stratum_delta
weighted_contribution <- trial_mass * (1 - stratum_share) * stratum_delta
xiE <- sum(raw_contribution)
xiL <- sum(weighted_contribution)
stopifnot(
  max(abs(raw_contribution - c(.1, -.1))) < 1e-12,
  max(abs(weighted_contribution - c(.025, -.075))) < 1e-12,
  abs(xiE) < 1e-12, abs(xiL + .05) < 1e-12
)
# Population covariance uses trial masses, not the sample covariance divisor.
covariance <- sum(trial_mass * (stratum_share - sum(trial_mass * stratum_share)) *
  (stratum_delta - sum(trial_mass * stratum_delta)))
stopifnot(abs(xiL - ((1 - sum(trial_mass * stratum_share)) * xiE - covariance)) < 1e-12)
open_tce("tce_weighting.pdf", height = 5.8, bottom = 4.2)
lambda <- seq(0, 1, length.out = 301)
panel(
  c(0, 1), c(4.85, 6.7), expression("Constant internal share " * lambda),
  "Treatment-minus-control contrast", "(a) A family of transported targets"
)
lines(lambda, theta_demo + (1 - lambda) * delta_demo, col = col_h, lwd = 2.2)
abline(h = theta_demo, col = "#666666", lty = 3)
points(c(1, share_demo, 0), c(5, 5.9, 6.5),
  pch = c(15, 18, 16),
  col = c(col_c, col_h, col_e), cex = 1.1
)
text(.96, 5.13, expression(theta == 5), adj = 1, cex = .94)
text(.40, 6.09, expression(tilde(theta)[0.4] == 5.9), adj = .5, cex = .91)
text(.04, 6.60, expression(tilde(theta)[lambda == 0] == 6.5), adj = 0, cex = .91)
arrows(.4, 5, .4, 5.9, code = 3, length = .07, col = "#333333", lwd = 1.2)
text(.44, 5.45, expression(xi[0.4] == 0.9), adj = 0, cex = .89)

mat <- rbind(raw_contribution, weighted_contribution)
bp <- barplot(mat, beside = TRUE, plot = FALSE)
panel(c(0, max(bp) + .8), c(-.137, .137), "Trial covariate stratum", "Contribution to causal gap",
  "(b) Cancellation changes after weighting",
  xaxt = "n"
)
abline(h = 0, col = "#555555", lwd = 1)
for (j in 1:2) for (i in 1:2) {
  rect(bp[i, j] - .46, min(0, mat[i, j]), bp[i, j] + .46, max(0, mat[i, j]),
    col = c(col_e, col_h)[i], border = NA
  )
  text(bp[i, j], mat[i, j] + if (mat[i, j] > 0) .010 else -.010,
    sprintf("%+.3f", mat[i, j]),
    cex = .89
  )
}
axis(1, at = colMeans(bp), labels = c("X = 0", "X = 1"), tick = FALSE)
text(colMeans(bp), rep(.130, 2),
  labels = expression(
    delta(0) == +0.2 ~ ~";" ~ ~ lambda(0) == 0.75,
    delta(1) == -0.2 ~ ~";" ~ ~ lambda(1) == 0.25
  ), cex = .76
)
text(colMeans(bp), rep(-.129, 2),
  labels = c("Trial mass = 0.5", "Trial mass = 0.5"), cex = .78
)
mtext(expression(xi[lambda == 0] == 0 ~ ~" but " ~ ~ xi[lambda] == -0.05), side = 1, line = 4.05, cex = .88)
overlay_legend(c("Unweighted stratum contribution", "Weighted stratum contribution"),
  c(col_e, col_h),
  pchs = c(15, 15), cex = .88
)
dev.off()

# TCE-centered spans display population geometry, not Monte Carlo coverage.
# All three group sizes grow proportionally to k; internal share stays 0.25.
scales <- c(1, 4, 16, 64)
nT <- 100 * scales
nC <- 50 * scales
nE <- 150 * scales
share_count <- nC / (nC + nE)
se_count <- sqrt(1 / nT + 1 / (nC + nE))
delta_fixed <- rep(.2, length(scales))
delta_local <- .2 / sqrt(scales)
xi_fixed <- (1 - share_count) * delta_fixed
xi_local <- (1 - share_count) * delta_local
half_width <- qnorm(.975) * se_count
stopifnot(
  max(abs(share_count - .25)) < 1e-12,
  max(abs(se_count - sqrt(.015 / scales))) < 1e-12,
  max(abs(xi_fixed - .15)) < 1e-12,
  max(abs(xi_local - .15 / sqrt(scales))) < 1e-12,
  max(abs(xi_local / se_count - .15 / sqrt(.015))) < 1e-12
)
open_tce("tce_precision.pdf", height = 5.35, bottom = 3.2)
for (j in 1:2) {
  centers <- if (j == 1) xi_fixed else xi_local
  col0 <- if (j == 1) col_h else col_e
  plot(NA,
    xlim = c(-.15, .42), ylim = c(.5, 4.5), xlab = "Treatment-effect scale",
    ylab = "Sample-size multiplier k", yaxt = "n"
  )
  mtext(if (j == 1) "(a) Fixed source discrepancy" else
    "(b) Local source discrepancy", side = 3, line = 1.7, adj = 0, cex = 1, font = 2)
  axis(2, at = 4:1, labels = scales, las = 1)
  grid(col = "#E5E7EB", lty = 1)
  box(col = "#555555")
  abline(v = 0, col = col_c, lty = 3, lwd = 1.8)
  yy <- 4:1
  for (i in seq_along(scales)) {
    segments(centers[i] - half_width[i], yy[i], centers[i] + half_width[i], yy[i],
      col = col0, lwd = 2.1
    )
    segments(c(centers[i] - half_width[i], centers[i] + half_width[i]), yy[i] - .09,
      c(centers[i] - half_width[i], centers[i] + half_width[i]), yy[i] + .09,
      col = col0, lwd = 1.5
    )
  }
  lines(centers, yy, col = col0, lty = 2, lwd = 1.1)
  points(centers, yy, pch = 18, col = col0, cex = 1.05)
  text(centers + .025, yy + .24, labels = sprintf("TCE = %.5f", centers), cex = .77, adj = 0)
  mtext(
    if (j == 1) expression(delta[k] == 0.2 ~ ~"," ~ ~ xi[k] == 0.15) else
      expression(delta[k] == 0.2 / sqrt(k) ~ ~"," ~ ~ xi[k] == 0.15 / sqrt(k)),
    side = 3, line = .1, cex = .87
  )
}
overlay_legend(
  c(
    "Trial causal effect: theta = 0",
    "Fixed gap: TCE-centered span",
    "Local gap: TCE-centered span"
  ),
  c(col_c, col_h, col_e),
  ltys = c(3, 1, 1), pchs = c(NA, 18, 18), cex = .83
)
dev.off()
cat("Verified SIM gap signs, transported targets, heterogeneous contributions, and precision scaling.\n")
cat("Generated three deterministic TCE illustrations; markers are not patient observations.\n")
print(data.frame(
  k = scales, nT = nT, nC = nC, nE = nE, SE = se_count,
  fixed_TCE = xi_fixed, local_TCE = xi_local,
  half_width_95 = half_width
), row.names = FALSE)

# Direct ECT--HCT comparison: more internal randomized controls both attenuate
# external-source drift and supply observed outcomes for learning the gap.
# HCT designs use n_T=2*n_C and n_E=150; the ECT endpoint has n_T=100.
# The gap and discrepancy SE shown here do not depend on n_T.
# At n_C=0 the trial control mean, hence delta, is not identified without an
# additional transportability assumption. No finite discrepancy SE is plotted.
anchor_nC <- 0:200
anchor_nE <- 150
anchor_delta <- .2
anchor_gap <- function(nc, rho) anchor_delta * anchor_nE * rho / (nc + anchor_nE * rho)
anchor_se <- function(nc) {
  if (any(nc <= 0)) stop("Discrepancy standard error requires an internal control sample.")
  sqrt(1 / nc + 1 / anchor_nE)
}
anchor_positive <- 5:200
stopifnot(
  all(vapply(r_values, function(rho)
    all(diff(anchor_gap(anchor_nC, rho)) < 0), logical(1))),
  all(diff(anchor_se(anchor_positive)) < 0),
  max(abs(vapply(r_values, function(rho) anchor_gap(0, rho), numeric(1)) - .2)) < 1e-12,
  abs(anchor_gap(50, 1) - .15) < 1e-12,
  abs(anchor_gap(150, 1) - .10) < 1e-12,
  abs(anchor_se(50) - sqrt(.02 + .006666666666666667)) < 1e-12,
  abs(anchor_se(50) - .1632993161855452) < 1e-12
)
open_tce("tce_anchor_comparison.pdf", height = 5.55, bottom = 3.5)
panel(
  c(0, 200), c(0, .245), expression("Number of internal controls " * n["1,0"]),
  expression("Transported causal gap " * xi[rho]),
  "(a) Internal controls attenuate the gap"
)
abline(h = anchor_delta, lty = 4, lwd = 1.6, col = "#444444")
for (j in seq_along(r_values))
  lines(anchor_nC, anchor_gap(anchor_nC, r_values[j]),
    col = cols[j], lty = types[j], lwd = 2.2
  )
text(100, .232, expression(n["1,0"] == 0 * ": externally controlled single arm"), cex = .8)
points(c(50, 150), c(.15, .10), col = cols[2], pch = 16, cex = .9)
text(50, .15, expression(xi[rho == 1] == 0.15), pos = 3, offset = .55, cex = .85, col = cols[2])
text(150, .10, expression(xi[rho == 1] == 0.10), pos = 1, offset = .6, cex = .85, col = cols[2])

panel(
  c(0, 200), c(0, .55), expression("Number of internal controls " * n["1,0"]),
  expression("Standard error of " * hat(delta)),
  "(b) The gap becomes empirically assessable"
)
lines(anchor_positive, anchor_se(anchor_positive), col = col_c, lwd = 2.2)
text(112, .495, "No internal controls:\nsource discrepancy is unidentified", cex = .86)
points(50, anchor_se(50), col = col_c, pch = 16, cex = .9)
text(56, anchor_se(50) + .033, "nC = 50: SE = 0.1633", adj = 0, cex = .84, col = col_c)
text(120, .345, expression(hat(delta) == bar(Y)["1,0"] - bar(Y)[0]), cex = .98)
text(120, .285, expression(SE(hat(delta)) == sqrt(1 / n["1,0"] + 1 / 150)), cex = .98)
overlay_legend(
  expression(
    rho == 0.25, rho == 1, rho == 4,
    "External-only gap = 0.2"
  ),
  c(cols, "#444444"),
  ltys = c(types, 4), cex = .88, ncol = 4
)
dev.off()
cat("Verified internal-control gap attenuation, discrepancy SE, and ECT boundary.\n")
cat("Anchor comparison uses nT=2*nC for HCT and nE=150; total trial size increases with nC.\n")
print(data.frame(
  nT = c(100, 300), nC = c(50, 150), nE = anchor_nE, rho = 1,
  hybrid_gap = anchor_gap(c(50, 150), 1),
  discrepancy_SE = anchor_se(c(50, 150))
), row.names = FALSE)

# Conference-width variants use native 6.75-inch width so label sizes remain
# readable when included across two columns. The full-size variants are kept.
open_compact <- function(name, height = 3.7, bottom = 3.0) {
  cairo_pdf(figure_output(name),
    width = 6.75, height = height,
    family = figure_font, pointsize = 10
  )
  par(
    mfrow = c(1, 2), mar = c(3.55, 3.65, 2, .85), oma = c(bottom, 0, 0, 0),
    las = 1, tcl = -.22, mgp = c(2.25, .58, 0), xaxs = "i", yaxs = "i",
    cex.axis = .9, cex.lab = .94
  )
}
panel_compact <- function(xlim, ylim, xlab, ylab, main) {
  plot(NA,
    xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, bty = "n"
  )
  mtext(main, side = 3, line = .65, adj = 0, font = 2, cex = .88)
  grid(col = "#E5E7EB", lty = 1)
  box(bty = "o", col = "#555555", lwd = .9)
  stopifnot(abs(diff(par("pin"))) < 1e-10)
}
legend_compact <- function(labels, colors, ltys, pchs = NULL, ncol = 2) {
  par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), pty = "m", plt = c(0, 1, 0, 1), xpd = NA)
  plot.new()
  legend("bottom",
    legend = labels, col = colors, lty = ltys, pch = pchs,
    lwd = 1.7, bty = "n", ncol = ncol, text.width = if (ncol == 4) NA else NULL,
    cex = .9, inset = .012,
    x.intersp = .7, y.intersp = 1.1
  )
}
open_compact("tce_surfaces.pdf", height = 3.05, bottom = 1.8)
local({
  # One exact, covariate-adjusted estimator supplies the inferential scale
  # in every panel. For Z = Y + 2X, estimate mean(Z_T) - .4 mean(Z_C)
  # - .6 mean(Z_E), with independent N(0, .6^2) outcome residuals.
  # Fixed group sizes give internal share .4 (external .6) and SE .18. Known common
  # slope and residual variance make the displayed Gaussian CI coverage
  # exact for the trial causal effect, conditional on the covariates.
  n_surface <- c(T = 20L, C = 10L, E = 15L)
  sd_surface <- .6
  share_surface <- n_surface["C"] / sum(n_surface[c("C", "E")])
  se_surface <- sd_surface * sqrt(1 / n_surface["T"] +
    share_surface^2 / n_surface["C"] +
    (1 - share_surface)^2 / n_surface["E"])
  discrepancy <- c(0, .03, 1.5)
  causal_gap <- (1 - share_surface) * discrepancy
  standardized_gap <- causal_gap / se_surface
  ci_coverage <- pnorm(qnorm(.975) - standardized_gap) -
    pnorm(-qnorm(.975) - standardized_gap)
  stopifnot(
    abs(share_surface - .4) < 1e-12, abs(se_surface - .18) < 1e-12,
    max(abs(causal_gap - c(0, .018, .9))) < 1e-12,
    max(abs(standardized_gap - c(0, .1, 5))) < 1e-12,
    abs(ci_coverage[1] - .95) < 1e-12,
    ci_coverage[2] > .948, ci_coverage[2] < .95, ci_coverage[3] < .002
  )
  # Synthetic points share their X values and residuals across scenarios;
  # only the external mean changes. These are not Monte Carlo replicates.
  set.seed(20260926)
  observations <- lapply(n_surface, function(n)
    data.frame(x = sort(runif(n, .035, .965)), error = rnorm(n, sd = sd_surface)))
  meanC <- function(x) 18 - 2 * x
  meanT <- function(x) 23 - 2 * x
  xx <- seq(0, 1, length.out = 301)
  par(
    mfrow = c(1, 3), cex = 1, mai = c(.40, .55, .40, .12), oma = c(1.85, 0, 0, 0),
    pty = "s", mgp = c(1.8, .5, 0), cex.axis = .82, cex.lab = .84
  )
  for (j in seq_along(discrepancy)) {
    meanE <- function(x) meanC(x) - discrepancy[j]
    meanH <- function(x) meanC(x) - causal_gap[j]
    plot(NA,
      xlim = c(0, 1), ylim = c(11, 24.5), xlab = "Baseline covariate x",
      ylab = if (j == 1) "Outcome" else "", xaxt = "n", yaxt = "n"
    )
    abline(v = c(0, .5, 1), h = seq(14, 24, 2), col = "#E5E7EB", lty = 1)
    axis(1, at = c(0, .5, 1), labels = c("0", "0.5", "1"))
    axis(2, at = seq(14, 24, 2))
    box(bty = "o", col = "#555555", lwd = .9)
    mtext(c(
      "(a) Transportability holds", "(b) Transportability fails",
      "(c) Transportability fails"
    )[j], side = 3, line = 1.20, adj = 0, cex = .82, font = 2)
    mtext(c(
      "Zero causal gap", "Small gap relative to SE",
      "Large gap relative to SE"
    )[j], side = 3, line = .05, adj = 0, cex = .82)
    for (k in 1:3) {
      z <- observations[[k]]
      mu <- list(meanT, meanC, meanE)[[k]]
      points(z$x, mu(z$x) + z$error,
        pch = c(17, 15, 16)[k],
        col = adjustcolor(c(col_t, col_c, col_e)[k], alpha.f = .43), cex = .56
      )
    }
    lines(xx, meanT(xx), col = col_t, lwd = 1.65)
    lines(xx, meanE(xx), col = col_e, lwd = 1.65)
    lines(xx, meanC(xx), col = col_c, lwd = 1.65)
    lines(xx, meanH(xx), col = col_h, lwd = 1.65, lty = 2)
    effect_x <- if (j == 2) .90 else .17
    arrows(effect_x, meanC(effect_x), effect_x, meanT(effect_x),
      code = 3, length = .045, lwd = .9
    )
    if (j == 2) {
      text(.90, 22.45, expression(theta == 5), adj = 1, cex = .82)
      # An explicitly labeled, separate horizontal scale resolves the tiny
      # difference without changing the main outcome scale or population means.
      # At x=.5 the three source-relative mean offsets are -.030, -.018, 0.
      # The inset sits in empty outcome space, above every control observation.
      zoom_x <- .5
      zoom_offset <- c(meanE(zoom_x), meanH(zoom_x), meanC(zoom_x)) - meanC(zoom_x)
      stopifnot(max(abs(zoom_offset - c(-.03, -.018, 0))) < 1e-12)
      zoom_position <- .13 + (zoom_offset + .03) / .03 * (.69 - .13)
      rect(.025, 18.90, .765, 21.35, col = "white", border = "#777777", lwd = .7)
      text(.395, 20.82, "Zoom: offsets from internal control\n(x = 0.5)", cex = .58)
      segments(.13, 20.06, .69, 20.06, col = "#777777", lwd = .8)
      segments(zoom_position, 19.89, zoom_position, 20.23,
        col = c(col_e, col_h, col_c), lwd = 1.3
      )
      points(zoom_position, rep(20.06, 3),
        pch = c(16, 18, 15),
        col = c(col_e, col_h, col_c), cex = .65
      )
      text(zoom_position, rep(19.35, 3), labels = c("-0.03", "-0.018", "0"), cex = .80)
    } else {
      text(.23, meanC(.17) + 2.55, expression(theta == 5), adj = 0, cex = .82)
    }
    if (j == 3) {
      arrows(.88, meanH(.88), .88, meanC(.88), code = 3, length = .035, lwd = .9)
      text(.84, meanC(.88) + 2.2, expression(xi[lambda] == 0.9), adj = 1, cex = .82)
    }
    text(.5, 13.35, bquote(TCE == .(c("5", "5.018", "5.9")[j]) ~ ~";" ~ ~
      xi[lambda] == .(c("0", "0.018", "0.9")[j])), cex = .80)
    text(.5, 12.45, bquote(xi[lambda] / SE == .(c("0", "0.1", "5")[j])), cex = .80)
    text(.5, 11.55, sprintf("Causal coverage = %.2f%%", 100 * ci_coverage[j]), cex = .80)
    stopifnot(abs(diff(par("pin"))) < 1e-10)
  }
  par(
    fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0),
    pty = "m", plt = c(0, 1, 0, 1), xpd = NA
  )
  plot.new()
  legend("bottom",
    legend = c(
      "Trial treated", "Internal control", "External control",
      "Hybrid control"
    ),
    col = c(col_t, col_c, col_e, col_h), lty = c(1, 1, 1, 2), pch = c(17, 15, 16, NA),
    lwd = 1.65, bty = "n", horiz = TRUE, cex = .82, inset = .025,
    x.intersp = .65, seg.len = 1.8
  )
  print(data.frame(
    delta = discrepancy, TCE = 5 + causal_gap, gap = causal_gap,
    SE = as.numeric(se_surface), gap_over_SE = standardized_gap,
    exact_95_coverage = ci_coverage
  ), row.names = FALSE)
})
dev.off()

open_compact("external-borrowing-anchor.pdf", height = 3.05, bottom = 1.85)
# Reduce only in-panel annotations; preserve axis, title, and legend sizing.
anchor_annotation_scale <- .86
# Reset cex after mfrow: the automatic 3-panel reduction would make labels
# smaller than 8 pt at the native two-column width.
par(
  mfrow = c(1, 3), cex = 1, mai = c(.40, .55, .40, .12), oma = c(1.85, 0, 0, 0),
  pty = "s", mgp = c(1.8, .5, 0), cex.axis = .82, cex.lab = .84
)
panel_compact(
  c(0, 200), c(0, .245), expression("Internal controls " * n["1,0"]),
  expression("Causal gap " * xi[rho]), "(a) Gap attenuation"
)
abline(h = anchor_delta, lty = 4, lwd = 1.4, col = "#444444")
for (j in seq_along(r_values))
  lines(anchor_nC, anchor_gap(anchor_nC, r_values[j]),
    col = cols[j], lty = types[j], lwd = 1.8
  )
text(100, .228, "ECT: no internal controls", cex = .82 * anchor_annotation_scale)
points(c(50, 100), anchor_gap(c(50, 100), 1), col = cols[2], pch = 16, cex = .72)
text(50, .15, "0.150", pos = 3, offset = .45, cex = .85 * anchor_annotation_scale, col = cols[2])
text(100, .12, "0.120", pos = 1, offset = .55, cex = .85 * anchor_annotation_scale, col = cols[2])
panel_compact(
  c(0, 200), c(0, .55), expression("Internal controls " * n["1,0"]),
  expression("Discrepancy SE: " * SE(hat(delta))), "(b) Learning discrepancy"
)
lines(anchor_positive, anchor_se(anchor_positive), col = col_c, lwd = 1.8)
text(110, .493, "No internal controls:\ndiscrepancy unidentified", cex = .82 * anchor_annotation_scale)
points(50, anchor_se(50), col = col_c, pch = 16, cex = .72)
text(66, anchor_se(50) + .030, "0.1633", adj = 0, cex = .85 * anchor_annotation_scale, col = col_c)
text(115, .340, expression(hat(delta) == bar(Y)["1,0"] - bar(Y)[0]), cex = .86 * anchor_annotation_scale)
text(112, .271, expression(SE(hat(delta)) == sqrt(1 / n["1,0"] + 1 / 150)), cex = .84 * anchor_annotation_scale)

# Along the 2:1 enrollment path, uncertainty in the treatment-effect
# estimator differs from uncertainty in the discrepancy estimator above.
# All outcome variances equal one and rho=1. The curve uses n_C>0 only;
# an ECT with fixed n_T is a separate design, not this path's endpoint.
anchor_hybrid_se <- function(nc) {
  if (any(nc <= 0)) stop("Hybrid standard error requires positive n_C on the 2:1 path.")
  sqrt(1 / (2 * nc) + 1 / (nc + anchor_nE))
}
anchor_standardized <- function(nc) anchor_gap(nc, 1) / anchor_hybrid_se(nc)
anchor_standardized_nC <- seq(1, 200, length.out = 600)
anchor_peak <- anchor_nE / sqrt(3)
anchor_comparison_nC <- c(50, 100)
stopifnot(
  max(abs(anchor_gap(anchor_comparison_nC, 1) - c(.15, .12))) < 1e-12,
  max(abs(anchor_standardized(anchor_comparison_nC) -
    c(1.224744871391589, 1.264911064067352))) < 1e-12,
  abs(optimize(anchor_standardized, c(1, 200), maximum = TRUE)$maximum -
    anchor_peak) < 1e-4,
  all(diff(anchor_standardized(seq(1, anchor_peak, length.out = 100))) > 0),
  all(diff(anchor_standardized(seq(anchor_peak, 200, length.out = 100))) < 0)
)
panel_compact(
  c(0, 200), c(0, 1.55), expression("Internal controls " * n["1,0"]),
  expression("Gap / estimator SE: " * xi[rho == 1] / s[H]), "(c) Standardized gap"
)
lines(anchor_standardized_nC, anchor_standardized(anchor_standardized_nC),
  col = cols[2], lwd = 1.8
)
points(anchor_comparison_nC, anchor_standardized(anchor_comparison_nC),
  col = cols[2], pch = 16, cex = .72
)
text(108, 1.44, expression(rho == 1 ~ ~"," ~ ~ n[1] == 2 * n["1,0"]), cex = .85 * anchor_annotation_scale)
segments(c(50, 100), anchor_standardized(anchor_comparison_nC) - .05,
  c(45, 135), c(.86, 1.03),
  col = cols[2], lwd = .7
)
text(45, .77, "1.225", cex = .85 * anchor_annotation_scale, col = cols[2])
text(139, .95, "1.265", cex = .85 * anchor_annotation_scale, col = cols[2])
text(105, .28, expression(s[H]^2 == frac(1, 2 * n["1,0"]) + frac(1, n["1,0"] + 150)), cex = .84 * anchor_annotation_scale)
legend_compact(
  expression(
    rho == 0.25, rho == 1, rho == 4,
    "External-only gap = 0.2"
  ),
  c(cols, "#444444"),
  ltys = c(types, 4), ncol = 4
)
dev.off()
cat("Saved conceptual figures.\n")
