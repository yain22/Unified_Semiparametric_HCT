# Fitted baseline simulations at the two trial enrollment levels.
script_path <- normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]), winslash = '/')
paper <- dirname(script_path)
source(file.path(paper, 'lib/cli_guard.R'))
root <- paper
args <- commandArgs(trailingOnly = TRUE)
opt <- function(n, d) {
  if (n == 'output' && any(grepl('^--output-dir=', args))) n <- 'output-dir'
  h <- grep(paste0('^--', n, '='), args, value = TRUE)
  if (length(h)) sub(paste0('^--', n, '='), '', h) else d
}
phase <- opt('phase', 'full')
reps <- as.integer(opt('reps', '5000'))
workers <- as.integer(opt('workers', '8'))
seed_base <- as.integer(opt('seed-base', '52461004'))
stopifnot(reps >= 2L, workers >= 1L, workers <= 8L, is.finite(seed_base), seed_base > 0, seed_base + 100000 + reps < .Machine$integer.max)
out <- normalizePath(opt('output', file.path(paper, '../../results/baseline')), winslash = '/', mustWork = FALSE)
out <- simulation_guard_output(out, paper)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
thread_vars <- c('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')
do.call(Sys.setenv, as.list(setNames(rep('1', length(thread_vars)), thread_vars)))
library(data.table)
setDTthreads(1L)
helper <- file.path(paper, 'lib/fitted_helpers.R')
ect_helper <- file.path(paper, 'lib/scalar_gap_ect_helpers.R')
trial_helper <- file.path(paper, 'lib/trial_baselines.R')
sim_source <- file.path(paper, 'lib/ect_sim_source.R')
source(helper)
helpers <- fg_load(paper)
source(ect_helper)
ect_helpers <- scalar_ect_load(paper, sim_source)
source(trial_helper)
adjustment <- as.formula(paste('~', paste(paste0('X', 1:10), collapse = '+')))
designs <- data.frame(N_trial = c(250L, 500L), NT = c(167L, 333L), NC = c(83L, 167L), N0 = c(500L, 1000L))
gh <- function(n) {
  J <- matrix(0, n, n)
  J[cbind(1:(n - 1), 2:n)] <- sqrt(1:(n - 1))
  J <- J + t(J)
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = e$vectors[1, ]^2)
}
gh128 <- gh(128)
gh256 <- gh(256)
calibration <- rbindlist(lapply(seq_len(nrow(designs)), function(i) {
  d <- designs[i, ]
  a <- d$NC / d$N_trial
  w <- d$N_trial / d$N0
  weight <- function(x) 1 - plogis(log(w * a) - .8 * x)
  calc <- function(g, m) sum(g$w * weight(g$x) * if (m == 0) 1 else g$x)
  K0 <- calc(gh256, 0)
  K1 <- calc(gh256, 1)
  int0 <- integrate(function(x) weight(x) * dnorm(x), -Inf, Inf, abs.tol = 1e-12, rel.tol = 1e-12)$value
  int1 <- integrate(function(x) weight(x) * x * dnorm(x), -Inf, Inf, abs.tol = 1e-12, rel.tol = 1e-12)$value
  error <- max(abs(c(K0 - calc(gh128, 0), K1 - calc(gh128, 1), K0 - int0, K1 - int1)))
  stopifnot(error < 1e-11)
  data.frame(d,
    K = K0, gap_offset = K1, K0 = K0, K1 = K1, K_128 = calc(gh128, 0), K_integrate = int0,
    quadrature_error = error, zero_gap_intercept = -K1 / K0,
    K_unit = 1 - d$NC / (d$NC + d$N0), K_ECT = 1
  )
}))
fwrite(calibration, file.path(out, 'calibration.csv'))
protocol <- list(
  stage = 'Fitted baseline simulation with additive source discrepancy',
  N_trial = c(250, 500), counts = as.data.frame(designs), repetitions_per_size = reps, seed_base = seed_base,
  seed = paste0(seed_base, '+100000*(N_trial==500)+replicate'), grid = (-12:12) / 40,
  X = 'All groups independent N10(0,I10); identical covariate laws',
  draw_order = 'external X matrix column-major; epsE; treated X matrix; control X matrix; epsT; epsC; assemble T,C,E',
  f = '.5+(.6/sqrt5)*sum(X1:5)+(.25/sqrt5)*sum(X1:5^2)',
  means = 'mT=mC=f; mE=f-(a+X1); theta0; fixed slope1 for every method and gap',
  variances = 'vT=1;vC=exp(.4X1);vE=exp(-.4X1)',
  true_share = 'lambda(X1)=w*(1-e)/(w*(1-e)+exp(.8X1)); w=nRCT/nE=.5',
  calibration = 'K0=E(1-lambda),K1=E[(1-lambda)X1]; a_U=a_Rmatched=(xi-K1)/K0; a_Runit=xi/(1-nC/(nC+nE)); a_ECT=xi. 256-node Gaussian quadrature checked at128 and integrate.',
  pairing = 'All procedures share base covariates and errors. U/R-matched share identical shifted outcomes. ECT/R-unit own amplitudes target the same xi but distinct external outcome laws. Trial-only target gap is zero.',
  methods = c('U', 'R-matched', 'R-unit', 'AIPW (ECT)', 'ANOVA (RCT)', 'ANCOVA (RCT)'),
  HCT_nuisances = 'Source-specific and pooled OLS quadratic21; source logistic11; separate logvariance intercept+X1 (2), own HC2 residuals; no clipping/fallback.',
  HCT_inference = 'Complete law-dependent U and restricted/common-control-score R; full analytic nuisance-stacked joint-M sandwich for all three, fixed T/C/E group-centered meat.',
  ECT = 'Original SIM DR quad_lm/glm and published plug-in EIF variance; only treated/external samples; fixed-group variance saved separately.',
  trial = 'Existing observed-outcome ANOVA/linear10X ANCOVA, HC2 robust standard errors',
  failures = 'Retain all failed cases and warnings; no incomplete summaries, no tail exclusions.',
  path = 'Base data include the fixed external -X1 discrepancy; fit intercept shifts a=0,1. Point/influence affine, variance quadratic. Direct refits check gap endpoints andzero on first2 reps/design.',
  decision = 'Repetitions and seed are explicitly configured above; published settings use 5000 and seed_base=52461004. Fixed gamma=.4 and additive slope1; no selection or tuning after outcomes. Zero aggregate gap does NOT impose mean transportability. Separate delta=0 reference is saved and not included in gap-grid curves.',
  created_before_fits = format(Sys.time(), tz = 'UTC', usetz = TRUE)
)
jsonlite::write_json(protocol, file.path(out, 'protocol.json'), auto_unbox = TRUE, pretty = TRUE, digits = 16)
generate <- function(n, b) {
  d <- designs[designs$N_trial == n, ]
  seed <- seed_base + 100000L * (n == 500L) + b
  set.seed(seed)
  xe <- matrix(rnorm(d$N0 * 10L), ncol = 10L)
  epsE <- rnorm(d$N0)
  xt <- matrix(rnorm(d$NT * 10L), ncol = 10L)
  xc <- matrix(rnorm(d$NC * 10L), ncol = 10L)
  epsT <- rnorm(d$NT)
  epsC <- rnorm(d$NC)
  X <- rbind(xt, xc, xe)
  f <- .5 + (.6 / sqrt(5)) * rowSums(X[, 1:5]) + (.25 / sqrt(5)) * rowSums(X[, 1:5]^2)
  error <- c(epsT, exp(.2 * xc[, 1]) * epsC, exp(-.2 * xe[, 1]) * epsE)
  data <- setNames(as.data.frame(X), paste0('X', 1:10))
  data$S <- c(rep(1L, d$NT + d$NC), rep(0L, d$N0))
  data$A <- c(rep(1L, d$NT), rep(0L, d$NC + d$N0))
  data$Y <- f + error - (1 - data$S) * X[, 1]
  list(data = data, h = rep(1, nrow(data)), seed = seed, design = d)
}
shift <- function(data, h, k) {
  data$Y <- data$Y - (1 - data$S) * k * h
  data
}
ect_path <- function(a, b) {
  z <- b$influence - a$influence
  N <- length(z)
  group <- factor(c(rep('T', a$n1), rep('E', a$n0)), levels = c('T', 'E'))
  pv <- function(x, y) sum(vapply(levels(group), function(g) sum(group == g) * cov(x[group == g], y[group == g]), numeric(1))) / N^2
  list(
    intercept = a$estimate, slope = b$estimate - a$estimate, influence0 = a$influence, influence_slope = z,
    variance = c(sum(a$influence^2), 2 * sum(a$influence * z), sum(z^2)) / N^2,
    working_variance = c(pv(a$influence, a$influence), 2 * pv(a$influence, z), pv(z, z))
  )
}
one_dataset <- function(job) {
  gen <- generate(job$N_trial, job$replicate)
  data <- gen$data
  h <- gen$h
  cal <- calibration[calibration$N_trial == job$N_trial, ]
  xi <- (-12:12) / 40
  f0 <- fg_fit(data, helpers, mean_basis = 'quadratic', variance_basis = 'x1')
  f1 <- fg_fit(shift(data, h, 1), helpers, mean_basis = 'quadratic', variance_basis = 'x1')
  ect0 <- scalar_ect_fit(data, ect_helpers$sim)
  ect1 <- scalar_ect_fit(shift(data, h, 1), ect_helpers$sim)
  trial <- hct_trial_wald_comparators(data, adjustment)
  paths <- lapply(names(f0), function(m) fg_path(f0[[m]], f1[[m]]))
  names(paths) <- names(f0)
  paths[['AIPW (ECT)']] <- ect_path(ect0, ect1)
  for (m in names(trial)) paths[[m]] <- list(
    intercept = trial[[m]]$estimate, slope = 0,
    variance = c(trial[[m]]$se^2, 0, 0), working_variance = c(trial[[m]]$se^2, 0, 0)
  )
  rows <- direct <- diagnostics <- list()
  for (m in names(paths)) {
    path <- paths[[m]]
    is_trial <- m %in% names(trial)
    K <- if (m == 'R-unit') cal$K_unit else if (m == 'AIPW (ECT)' || is_trial) 1 else cal$K
    offset <- if (m %in% c('U', 'R-matched')) cal$gap_offset else 0
    k <- if (is_trial) rep(0, 25) else (xi - offset) / K
    target <- if (is_trial) rep(0, 25) else xi
    achieved <- if (is_trial) rep(0, 25) else K * k + offset
    stopifnot(max(abs(achieved - target)) < 1e-12)
    estimate <- path$intercept + k * path$slope
    variance <- path$variance[1] + k * path$variance[2] + k^2 * path$variance[3]
    wvariance <- path$working_variance[1] + k * path$working_variance[2] + k^2 * path$working_variance[3]
    if (any(!is.finite(c(estimate, variance, wvariance))) || any(c(variance, wvariance) <= 0)) stop('Nonfinite or nonpositive result.')
    se <- sqrt(variance)
    ws <- sqrt(wvariance)
    lo <- estimate - qnorm(.975) * se
    hi <- estimate + qnorm(.975) * se
    wlo <- estimate - qnorm(.975) * ws
    whi <- estimate + qnorm(.975) * ws
    rows[[m]] <- data.frame(
      N_trial = job$N_trial, setting = 'additive_discrepancy', method = m, replicate = job$replicate, seed = gen$seed,
      xi_index = -12:12, xi = xi, amplitude = k, K = K, gap_offset = offset, theta = 0, target_TCE = target, achieved_gap = achieved,
      estimate = estimate, se = se, coverage_theta = as.integer(lo <= 0 & hi >= 0), bias_theta = estimate, mse_theta = estimate^2,
      width = hi - lo, coverage_TCE = as.integer(lo <= target & hi >= target), bias_TCE = estimate - target, mse_TCE = (estimate - target)^2,
      estimated_variance = variance, working_se = ws, working_width = whi - wlo,
      working_coverage_theta = as.integer(wlo <= 0 & whi >= 0), working_coverage_TCE = as.integer(wlo <= target & whi >= target)
    )
    if (m %in% names(f0)) {
      a <- f0[[m]]
      b <- f1[[m]]
      diagnostics[[m]] <- data.frame(
        N_trial = job$N_trial, method = m, replicate = job$replicate, seed = gen$seed,
        estimate0 = path$intercept, amplitude_slope = path$slope, xi_slope = path$slope / K,
        lambda_trial_mean = mean(a$lambda[data$S == 1]), lambda_min = min(a$lambda), lambda_max = max(a$lambda),
        share_shift_error = max(abs(a$lambda - b$lambda)), max_leverage = max(a$max_leverage, b$max_leverage),
        min_variance_ratio = min(a$variance_ratio_range), max_variance_ratio = max(a$variance_ratio_range),
        nuisance_equation_max = max(a$nuisance_equation_max, b$nuisance_equation_max)
      )
      if (max(abs(a$lambda - b$lambda)) > 1e-8) stop('Share changed under OR-basis external shift.')
    }
    if (job$replicate <= 2L) for (i in c(1L, 13L, 25L)) {
      dd <- shift(data, h, k[i])
      if (m %in% names(f0)) {
        actual <- fg_fit(dd, helpers, mean_basis = 'quadratic', variance_basis = 'x1')[[m]]
        IFerr <- max(abs(actual$influence - path$influence0 - k[i] * path$influence_slope))
      } else if (m == 'AIPW (ECT)') {
        actual <- scalar_ect_fit(dd, ect_helpers$sim)
        IFerr <- max(abs(actual$influence - path$influence0 - k[i] * path$influence_slope))
      } else {
        actual <- hct_trial_wald_comparators(dd, adjustment)[[m]]
        IFerr <- 0
      }
      direct[[length(direct) + 1L]] <- data.frame(
        N_trial = job$N_trial, method = m, replicate = job$replicate, xi = xi[i],
        estimate_error = abs(actual$estimate - estimate[i]), SE_error = abs(actual$se - se[i]), IF_error = IFerr
      )
    }
  }
  direct <- rbindlist(direct)
  if (nrow(direct) && max(direct$estimate_error, direct$SE_error, direct$IF_error) > 1e-7) stop('Direct-refit path check failed.')
  # A separate transportability reference: restore external means to f.
  reference_data <- data
  reference_data$Y <- data$Y + (1 - data$S) * data$X1
  reference_fits <- fg_fit(reference_data, helpers, mean_basis = 'quadratic', variance_basis = 'x1')
  reference_fits[['AIPW (ECT)']] <- scalar_ect_fit(reference_data, ect_helpers$sim)
  reference_fits <- c(reference_fits, trial)
  reference <- rbindlist(lapply(names(reference_fits), function(m) {
    a <- reference_fits[[m]]
    data.frame(
      N_trial = job$N_trial, method = m, replicate = job$replicate, seed = gen$seed,
      delta_reference = 'identically_zero', theta = 0, target_TCE = 0, estimate = a$estimate, se = a$se,
      coverage_theta = as.integer(abs(a$estimate) <= qnorm(.975) * a$se),
      bias_theta = a$estimate, mse_theta = a$estimate^2, width = 2 * qnorm(.975) * a$se,
      coverage_TCE = as.integer(abs(a$estimate) <= qnorm(.975) * a$se),
      bias_TCE = a$estimate, mse_TCE = a$estimate^2, estimated_variance = a$se^2
    )
  }))
  list(
    success = TRUE, rows = rbindlist(rows), reference = reference, diagnostics = rbindlist(diagnostics), direct = direct,
    pairing = data.frame(
      N_trial = job$N_trial, replicate = job$replicate, seed = gen$seed,
      full_sample_sha256 = digest::digest(data, algo = 'sha256'), ECT_sample_sha256 = digest::digest(data[data$A == 1 | data$S == 0, ], algo = 'sha256'),
      ect_clipped = ect0$external_clipped + ect0$treated_clipped
    ),
    warnings = unique(c(unlist(lapply(f0, function(x) x$warnings)), unlist(lapply(f1, function(x) x$warnings)), ect0$warnings, ect1$warnings))
  )
}
one_job <- function(job) {
  if (file.exists(job$file)) return(job$file)
  ans <- tryCatch(one_dataset(job), error = function(e) list(success = FALSE, error = conditionMessage(e)))
  ans$N_trial <- job$N_trial
  ans$replicate <- job$replicate
  ans$seed <- seed_base + 100000L * (job$N_trial == 500L) + job$replicate
  saveRDS(ans, paste0(job$file, '.partial'), compress = FALSE)
  if (!file.rename(paste0(job$file, '.partial'), job$file)) stop('Checkpoint rename failed.')
  job$file
}
script <- normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]), winslash = '/')
input_paths <- c(script, list.files(file.path(paper, 'lib'), full.names = TRUE, pattern = '[.]R$'))
names(input_paths) <- substring(input_paths, nchar(paper) + 2L)
hashes <- vapply(input_paths, function(p) digest::digest(file = p, algo = 'sha256'), character(1))
key <- digest::digest(list(hashes, designs, calibration, protocol$seed, protocol$grid, R.version.string), algo = 'sha256')
checkdir <- file.path(out, 'checkpoints', substr(key, 1, 20))
dir.create(checkdir, recursive = TRUE, showWarnings = FALSE)
jobs <- list()
for (n in c(250L, 500L)) for (b in seq_len(reps)) jobs[[length(jobs) + 1L]] <- list(
  N_trial = n, replicate = b, file = file.path(checkdir, sprintf('n%d_rep%04d.rds', n, b))
)
start <- Sys.time()
cl <- parallel::makePSOCKcluster(workers, outfile = file.path(out, 'workers.log'))
parallel::clusterExport(cl, c(
  'paper', 'root', 'helper', 'ect_helper', 'trial_helper', 'sim_source', 'thread_vars', 'designs', 'calibration',
  'adjustment', 'seed_base', 'generate', 'shift', 'ect_path', 'one_dataset', 'one_job'
), envir = environment())
invisible(parallel::clusterEvalQ(cl, {
  do.call(Sys.setenv, as.list(setNames(rep('1', length(thread_vars)), thread_vars)))
  library(data.table)
  setDTthreads(1L)
  source(helper)
  helpers <- fg_load(paper)
  source(ect_helper)
  ect_helpers <- scalar_ect_load(paper, sim_source)
  source(trial_helper)
  if (requireNamespace('RhpcBLASctl', quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
  }
  NULL
}))
cat('Additive structural validation:', length(jobs), 'base samples;6 methods;', workers, 'workers.\n')
flush.console()
files <- tryCatch(parallel::parLapplyLB(cl, jobs, one_job), finally = parallel::stopCluster(cl))
ans <- lapply(unlist(files), readRDS)
failures <- rbindlist(lapply(ans, function(x) if (!x$success) data.frame(N_trial = x$N_trial, replicate = x$replicate, seed = x$seed, error = x$error) else NULL))
if (!ncol(failures)) failures <- data.table(N_trial = integer(), replicate = integer(), seed = integer(), error = character())
fwrite(failures, file.path(out, 'failures.csv'))
good <- Filter(function(x) x$success, ans)
raw <- rbindlist(lapply(good, '[[', 'rows'))
fwrite(raw, file.path(out, 'replicates.csv.gz'), compress = 'gzip')
reference <- rbindlist(lapply(good, '[[', 'reference'))
fwrite(reference, file.path(out, 'transportability_reference.csv.gz'), compress = 'gzip')
reference_metrics <- c('coverage_theta', 'bias_theta', 'mse_theta', 'width', 'se', 'coverage_TCE', 'bias_TCE', 'mse_TCE', 'estimated_variance')
reference_summary <- reference[,
  {
    r <- list(n_reps = .N)
    for (m in reference_metrics) {
      r[[m]] <- mean(.SD[[m]])
      r[[paste0(m, '_mcse')]] <- sd(.SD[[m]]) / sqrt(.N)
    }
    r$empirical_variance <- var(.SD$bias_theta)
    r
  },
  by = .(N_trial, method),
  .SDcols = reference_metrics
]
fwrite(reference_summary, file.path(out, 'transportability_reference_summary.csv'))

diagnostics <- rbindlist(lapply(good, '[[', 'diagnostics'))
fwrite(diagnostics, file.path(out, 'fit_diagnostics.csv'))
direct <- rbindlist(lapply(good, '[[', 'direct'))
fwrite(direct, file.path(out, 'direct_validation.csv'))
pairing <- rbindlist(lapply(good, '[[', 'pairing'))
fwrite(pairing, file.path(out, 'sample_pairing.csv'))
warnings <- rbindlist(lapply(good, function(x) if (length(x$warnings)) data.frame(N_trial = x$N_trial, replicate = x$replicate, warning = x$warnings) else NULL))
if (!ncol(warnings)) warnings <- data.table(N_trial = integer(), replicate = integer(), warning = character())
fwrite(warnings, file.path(out, 'warnings.csv'))
if (nrow(failures)) stop('Failures preserved: no incomplete full summary.')
stopifnot(nrow(raw) == 2L * reps * 6L * 25L, !anyDuplicated(raw[, .(N_trial, method, replicate, xi_index)]))
metrics <- c(
  'coverage_theta', 'bias_theta', 'mse_theta', 'width', 'se', 'coverage_TCE', 'bias_TCE', 'mse_TCE', 'estimated_variance',
  'working_se', 'working_width', 'working_coverage_theta', 'working_coverage_TCE'
)
summary <- raw[,
  {
    r <- list(n_reps = .N, n_failed = 0L)
    for (m in metrics) {
      r[[m]] <- mean(.SD[[m]])
      r[[paste0(m, '_mcse')]] <- sd(.SD[[m]]) / sqrt(.N)
    }
    r$empirical_variance <- var(.SD$bias_theta)
    r
  },
  by = .(N_trial, setting, method, xi_index, xi, amplitude, K, gap_offset, target_TCE, achieved_gap),
  .SDcols = metrics
]
fwrite(summary, file.path(out, 'summary.csv'))
fwrite(summary[xi_index %in% c(-12L, -6L, 0L, 6L, 12L)], file.path(out, 'key_gaps.csv'))
pairs <- list()
for (comparator in c('R-matched', 'R-unit')) {
  cols <- c('N_trial', 'replicate', 'xi_index', 'xi', 'seed', 'amplitude', metrics)
  j <- merge(raw[method == 'U', ..cols], raw[method == comparator, ..cols], by = c('N_trial', 'replicate', 'xi_index'), suffixes = c('_U', '_R'))
  stopifnot(nrow(j) == 2L * reps * 25L, all(j$seed_U == j$seed_R))
  if (comparator == 'R-matched') stopifnot(max(abs(j$amplitude_U - j$amplitude_R)) == 0)
  for (m in metrics) j[, (paste0(m, '_difference')) := get(paste0(m, '_U')) - get(paste0(m, '_R'))]
  dc <- paste0(metrics, '_difference')
  q <- j[,
    {
      r <- list(n_pairs = .N)
      for (m in dc) {
        r[[m]] <- mean(.SD[[m]])
        r[[paste0(m, '_mcse')]] <- sd(.SD[[m]]) / sqrt(.N)
      }
      r
    },
    by = .(N_trial, xi_index, xi = xi_U),
    .SDcols = dc
  ]
  q[, comparator := comparator]
  q[, method := 'U']
  pairs[[comparator]] <- q
}
paired <- rbindlist(pairs)
fwrite(paired, file.path(out, 'paired_differences.csv'))
fwrite(paired[xi_index %in% c(-12L, -6L, 0L, 6L, 12L)], file.path(out, 'paired_key_gaps.csv'))
stopifnot(identical(hashes, vapply(input_paths, function(p) digest::digest(file = p, algo = 'sha256'), character(1))))
validation <- list(
  success = TRUE, phase = phase, repetitions_per_size = reps, rows = nrow(raw), failures = nrow(failures),
  warnings = nrow(warnings), checkpoint_key = key, source_hashes = as.list(hashes), source_files_unchanged = TRUE,
  calibration_max_error = max(calibration$quadrature_error), reference_rows = nrow(reference), zero_gap_is_not_transportability = TRUE, direct_estimate_max_error = max(direct$estimate_error),
  direct_SE_max_error = max(direct$SE_error), direct_IF_max_error = max(direct$IF_error),
  U_Rmatched_identical_amplitudes = TRUE, ECT_any_clipping = sum(pairing$ect_clipped),
  elapsed_seconds = as.numeric(difftime(Sys.time(), start, units = 'secs'))
)
jsonlite::write_json(validation, file.path(out, 'validation.json'), auto_unbox = TRUE, pretty = TRUE, digits = 16)
print(summary[xi_index %in% c(0L, 6L), .(N_trial, method, xi, coverage_theta, bias_theta, mse_theta, width, coverage_TCE, bias_TCE)])
print(validation[c('success', 'rows', 'failures', 'elapsed_seconds')])
