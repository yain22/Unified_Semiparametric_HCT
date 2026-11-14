# Oracle evaluation paired to a completed fitted baseline run.
# U keeps all seven unrestricted score terms. ECT keeps its existing plug-in
# EIF variance convention; fixed-group ECT variance is saved as a diagnostic.
script_path <- normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]), winslash = '/')
paper <- dirname(script_path)
source(file.path(paper, 'lib/cli_guard.R'))
args <- commandArgs(trailingOnly = TRUE)
opt <- function(n, d) {
  if (n == 'output' && any(grepl('^--output-dir=', args))) n <- 'output-dir'
  h <- grep(paste0('^--', n, '='), args, value = TRUE)
  if (length(h)) sub(paste0('^--', n, '='), '', h) else d
}
seed_base <- as.integer(opt('seed-base', '52461004'))
out <- normalizePath(opt('output', file.path(paper, '../../results/oracle')), winslash = '/', mustWork = FALSE)
source_dir <- normalizePath(opt('baseline', file.path(paper, '../../results/baseline')), winslash = '/', mustWork = TRUE)
out <- simulation_guard_output(out, paper, source_dir)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
library(data.table)
setDTthreads(1L)
B <- as.integer(opt('reps', '5000'))
stopifnot(B >= 2L, is.finite(seed_base), seed_base > 0, seed_base + B < .Machine$integer.max)
nT <- 167L
nC <- 83L
nE <- 500L
nR <- nT + nC
e <- nT / nR
aa <- 1 - e
w <- nR / nE
p <- nR / (nR + nE)
grid <- (-12:12) / 40
pairing <- fread(file.path(source_dir, 'sample_pairing.csv'))[N_trial == nR]
fitted <- fread(file.path(source_dir, 'summary.csv'))[N_trial == nR]
source_paths <- file.path(source_dir, c('summary.csv', 'sample_pairing.csv', 'calibration.csv', 'protocol.json'))
names(source_paths) <- basename(source_paths)
source_hashes <- vapply(source_paths, function(f) digest::digest(file = f, algo = 'sha256'), character(1))
gh <- function(n) {
  J <- matrix(0, n, n)
  J[cbind(1:(n - 1), 2:n)] <- sqrt(1:(n - 1))
  J <- J + t(J)
  z <- eigen(J, symmetric = TRUE)
  list(x = z$values, w = z$vectors[1, ]^2)
}
g <- gh(256L)
gcheck <- gh(128L)
calc <- function(z) {
  L <- plogis(log(nC / nE) - .8 * z$x)
  c(K0 = sum(z$w * (1 - L)), K1 = sum(z$w * (1 - L) * z$x))
}
K <- calc(g)
stopifnot(max(abs(K - calc(gcheck))) < 1e-12)
cal_saved <- fread(file.path(source_dir, 'calibration.csv'))[N_trial == nR]
stopifnot(abs(K['K0'] - cal_saved$K0) < 1e-12, abs(K['K1'] - cal_saved$K1) < 1e-12)

generate <- function(b) {
  seed <- seed_base + b
  set.seed(seed)
  xe <- matrix(rnorm(nE * 10L), ncol = 10L)
  epsE <- rnorm(nE)
  xt <- matrix(rnorm(nT * 10L), ncol = 10L)
  xc <- matrix(rnorm(nC * 10L), ncol = 10L)
  epsT <- rnorm(nT)
  epsC <- rnorm(nC)
  X <- rbind(xt, xc, xe)
  f <- .5 + (.6 / sqrt(5)) * rowSums(X[, 1:5]) + (.25 / sqrt(5)) * rowSums(X[, 1:5]^2)
  error <- c(epsT, exp(.2 * xc[, 1]) * epsC, exp(-.2 * xe[, 1]) * epsE)
  data <- setNames(as.data.frame(X), paste0('X', 1:10))
  data$S <- c(rep(1L, nR), rep(0L, nE))
  data$A <- c(rep(1L, nT), rep(0L, nC + nE))
  data$Y <- f + error - (1 - data$S) * X[, 1]
  T <- data$S * data$A
  C <- data$S * (1 - data$A)
  E <- 1 - data$S
  list(
    data = data, x = X[, 1], f = f, error = error, T = T, C = C, E = E, S = data$S,
    z = c(epsT, epsC, epsE), seed = seed, group = c(rep('T', nT), rep('C', nC), rep('E', nE))
  )
}
gcov <- function(a, b, d) sum(vapply(c('T', 'C', 'E'), function(k) {
  z <- d$group == k
  sum(z) * cov(a[z], b[z])
}, numeric(1))) / nR^2
direct_U <- function(d, delta) {
  L <- plogis(log(nC / nE) - .8 * d$x)
  H <- L * (1 - L)
  Y <- d$data$Y - d$E * (delta - d$x)
  mT <- mC <- d$f
  mE <- d$f - delta
  d$S * (mT - L * mC - (1 - L) * mE) + d$T / e * (Y - mT) - d$C * L / aa * (Y - mC) -
    d$E * (1 - L) * w * (Y - mE) - delta * H * (d$S - p) / (1 - p) +
    delta * H * d$C / aa * ((Y - mC)^2 / exp(.4 * d$x) - 1) -
    delta * H * d$E * w * ((Y - mE)^2 / exp(-.4 * d$x) - 1)
}
coeff <- checks <- list()
moments <- list()
start <- Sys.time()
for (b in seq_len(B)) {
  d <- generate(b)
  saved <- pairing[replicate == b]
  full_hash <- digest::digest(d$data, algo = 'sha256')
  ect_hash <- digest::digest(d$data[d$data$A == 1 | d$data$S == 0, ], algo = 'sha256')
  stopifnot(
    nrow(saved) == 1L, saved$seed == d$seed,
    saved$full_sample_sha256 == full_hash, saved$ECT_sample_sha256 == ect_hash
  )
  L <- plogis(log(nC / nE) - .8 * d$x)
  H <- L * (1 - L)
  # Affine complete U score in the discrepancy intercept.
  noise <- d$T / e * d$error - d$C * L / aa * d$error - d$E * (1 - L) * w * d$error
  effect <- d$S * (1 - L) - H * (d$S - p) / (1 - p) + H * d$C / aa * (d$z^2 - 1) - H * d$E * w * (d$z^2 - 1)
  raw0 <- noise + effect * (d$x - K['K1'] / K['K0'])
  raw1 <- effect / K['K0']
  coeff[[2 * b - 1L]] <- data.table(
    method = 'U', replicate = b, seed = d$seed,
    intercept = sum(raw0) / nR, slope = sum(raw1) / nR, v0 = gcov(raw0, raw0, d),
    v1 = 2 * gcov(raw0, raw1, d), v2 = gcov(raw1, raw1, d), fixed_group_v0 = gcov(raw0, raw0, d)
  )
  zT <- d$x[d$T == 1] + d$error[d$T == 1]
  rE <- d$error[d$E == 1]
  est0 <- mean(zT) - mean(rE)
  odds <- nT / nE
  plugin <- (sum((zT - est0)^2) + odds^2 * sum(rE^2)) / nT^2
  fixed <- var(zT) / nT + var(rE) / nE
  coeff[[2 * b]] <- data.table(
    method = 'AIPW (ECT)', replicate = b, seed = d$seed,
    intercept = est0, slope = 1, v0 = plugin, v1 = 0, v2 = 0, fixed_group_v0 = fixed
  )
  moments[[b]] <- data.table(
    replicate = b, source_mean = mean(d$S), lambda_trial_mean = mean(L[d$S == 1]),
    x1_trial_mean = mean(d$x[d$S == 1]), control_error2 = mean(d$error[d$C == 1]^2),
    external_error2 = mean(rE^2), paired_full_hash = full_hash, paired_ect_hash = ect_hash
  )
  if (b <= 10L) for (xi in c(-.3, 0, .15, .3)) {
    delta <- (xi - K['K1']) / K['K0'] + d$x
    direct <- direct_U(d, delta)
    predicted <- raw0 + xi * raw1
    # Independent direct ECT expression from observed outcomes and true mean.
    m0 <- d$f - d$x - xi
    Y <- d$data$Y - d$E * xi
    eraw <- d$T * (Y - m0) - d$E * odds * (Y - m0)
    eestimate <- sum(eraw) / nT
    phi_num <- eraw - d$T * eestimate
    checks[[length(checks) + 1L]] <- data.table(
      replicate = b, xi = xi,
      U_score_error = max(abs(direct - predicted)),
      U_variance_error = abs(gcov(direct, direct, d) - gcov(predicted, predicted, d)),
      ECT_estimate_error = abs(eestimate - (est0 + xi)),
      ECT_plugin_error = abs(sum(phi_num^2) / nT^2 - plugin)
    )
  }
}
co <- rbindlist(coeff)
ch <- rbindlist(checks)
mo <- rbindlist(moments)
stopifnot(nrow(co) == 2L * B, max(as.matrix(ch[, -c('replicate', 'xi')])) < 1e-10)
metric <- function(z, xi) {
  estimate <- z$intercept + xi * z$slope
  V <- z$v0 + xi * z$v1 + xi^2 * z$v2
  stopifnot(all(is.finite(V)), all(V > 0))
  half <- qnorm(.975) * sqrt(V)
  data.table(
    coverage_theta = as.integer(abs(estimate) <= half), bias_theta = estimate,
    mse_theta = estimate^2, width = 2 * half, coverage_TCE = as.integer(abs(estimate - xi) <= half),
    bias_TCE = estimate - xi, estimated_variance = V
  )
}
summaries <- list()
analytics <- list()
for (method_name in c('U', 'AIPW (ECT)')) for (xi in grid) {
  z <- co[method == method_name]
  val <- metric(z, xi)
  row <- list(
    mode = 'Oracle nuisance', N_trial = nR, nT = nT, nC = nC, nE = nE,
    method = method_name, xi = xi, xi_index = as.integer(round(xi * 40)), achieved_gap = xi, n_reps = B, n_failed = 0
  )
  for (nm in names(val)) {
    row[[nm]] <- mean(val[[nm]])
    row[[paste0(nm, '_mcse')]] <- sd(val[[nm]]) / sqrt(B)
  }
  row$empirical_variance <- var(val$bias_theta)
  summaries[[length(summaries) + 1L]] <- as.data.table(row)
  if (method_name == 'U') {
    x <- g$x
    L <- plogis(log(nC / nE) - .8 * x)
    H <- L * (1 - L)
    delta <- (xi - K['K1']) / K['K0'] + x
    cc <- cbind((1 - L)^2, (1 - L)^2, w * H) * delta
    dd <- cbind(0 * x, H / aa, -w * H) * delta
    bb <- cbind(rep(1 / e, length(x)), -L / aa * exp(.2 * x), -(1 - L) * w * exp(-.2 * x))
    av <- function(a) sum(g$w * a)
    analytic_mean <- sum(c(nT, nC, nE) * apply(cc, 2, av)) / nR
    analytic_var <- sum(c(nT, nC, nE) * vapply(1:3, function(j) av(cc[, j]^2) - av(cc[, j])^2 + av(bb[, j]^2) + 2 * av(dd[, j]^2), numeric(1))) / nR^2
  } else {
    analytic_mean <- xi
    analytic_var <- 2 / nT + exp(.08) / nE
  }
  analytics[[length(analytics) + 1L]] <- data.table(
    method = method_name, xi = xi,
    expected_bias = analytic_mean, expected_sampling_variance = analytic_var, expected_mse = analytic_mean^2 + analytic_var,
    bias_z = (mean(val$bias_theta) - analytic_mean) / (sd(val$bias_theta) / sqrt(B)),
    empirical_variance = var(val$bias_theta), mean_estimated_variance = mean(val$estimated_variance)
  )
}
oracle <- rbindlist(summaries)
analytical <- rbindlist(analytics)
stopifnot(max(abs(analytical$expected_bias - analytical$xi)) < 1e-12)
if (B >= 1000L) stopifnot(max(abs(analytical$bias_z)) < 5) # MC screening is not meaningful for tiny smoke runs.
methods <- c('U', 'AIPW (ECT)', 'ANOVA (RCT)', 'ANCOVA (RCT)')
estimated <- copy(fitted[method %in% methods])
estimated[, mode := 'Estimated nuisance']
estimated[, c('nT', 'nC', 'nE') := list(nT, nC, nE)]
trial <- copy(estimated[method %in% c('ANOVA (RCT)', 'ANCOVA (RCT)')])
trial[, mode := 'Oracle nuisance']
combined <- rbindlist(list(oracle, trial, estimated), fill = TRUE)
stopifnot(nrow(combined) == 200L, all(combined$n_reps == B))
fwrite(co, file.path(out, 'oracle_replicate_coefficients.csv.gz'), compress = 'gzip')
fwrite(oracle, file.path(out, 'oracle_summary.csv'))
fwrite(combined, file.path(out, 'summary.csv'))
fwrite(combined[abs(xi) %in% c(0, .15)], file.path(out, 'key_gaps.csv'))
fwrite(analytical, file.path(out, 'analytic_checks.csv'))
fwrite(ch, file.path(out, 'direct_score_checks.csv'))
fwrite(mo, file.path(out, 'generator_pairing_checks.csv'))
protocol <- list(
  counts = c(treated = nT, internal = nC, external = nE), repetitions = B, seed = paste0(seed_base, ' + replicate'),
  covariates = 'Independent N10(0,I10) in each group; exact saved fitted draw order and datasets verified by SHA256.',
  means = 'mT=mC=f=.5+.6/sqrt(5)*sum(X1:5)+.25/sqrt(5)*sum(X1:5^2); mE=f-a-X1.',
  variances = 'vT=1;vC=exp(.4X1);vE=exp(-.4X1)', K0 = unname(K['K0']), K1 = unname(K['K1']),
  calibration = 'HCT a=(xi-K1)/K0; ECT a=xi; trial-only has zero causal-gap.',
  oracle_U = 'Complete seven-term H_eta at true nuisances, with source and both variance derivative corrections. Fixed-count within-group score covariance.',
  oracle_ECT = 'True external outcome mean and source odds; same plug-in EIF variance convention as fitted ECT. Fixed-group variance also saved as diagnostic.',
  trial_reference = 'Saved fitted ANOVA/ANCOVA results reused identically in both rows; no external nuisance stage.',
  fitted_source = 'summary.csv',
  limitation = 'Oracle is true-input evaluation of the specified complete score, not an efficiency bound for a smaller known-nuisance model. Zero aggregate gap retains heterogeneous discrepancy.',
  validation = list(failures = 0, paired_datasets = B, max_direct_error = max(as.matrix(ch[, -c('replicate', 'xi')])), max_bias_z = max(abs(analytical$bias_z))),
  source_sha256 = as.list(source_hashes), elapsed_seconds = as.numeric(difftime(Sys.time(), start, units = 'secs'))
)
stopifnot(identical(source_hashes, vapply(source_paths, function(f) digest::digest(file = f, algo = 'sha256'), character(1))))
jsonlite::write_json(protocol, file.path(out, 'protocol.json'), auto_unbox = TRUE, pretty = TRUE, digits = 16)
print(combined[xi == 0, .(mode, method, coverage_theta, bias_theta, mse_theta, width, coverage_TCE)])
print(protocol$validation)
