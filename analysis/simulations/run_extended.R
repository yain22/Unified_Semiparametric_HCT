# Four scenarios varying external sample size, covariates, variances, or discrepancy.
script_path <- normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]), winslash = '/')
paper <- dirname(script_path)
source(file.path(paper, 'lib/cli_guard.R'))
args <- commandArgs(trailingOnly = TRUE)
opt <- function(n, d) {
  if (n == 'output' && any(grepl('^--output-dir=', args))) n <- 'output-dir'
  h <- grep(paste0('^--', n, '='), args, value = TRUE)
  if (length(h)) sub(paste0('^--', n, '='), '', h) else d
}
phase <- opt('phase', 'full')
reps <- as.integer(opt('reps', '5000'))
workers <- as.integer(opt('workers', '8'))
enrollment <- as.integer(opt('enrollment', '250'))
stopifnot(phase %in% c('pilot', 'full'), reps >= 2L, workers >= 1L, workers <= 8L,
          enrollment %in% c(250L, 500L))
out <- normalizePath(opt('output', file.path(paper, '../../results/extended')), winslash = '/', mustWork = FALSE)
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
designs <- data.table(
  id = 1:4,
  scenario = c('fewer_external', 'source_shift', 'stronger_heteroskedasticity', 'quadratic_discrepancy'),
  nT = if (enrollment == 250L) 167L else 333L,
  nE = enrollment * c(1L, 2L, 2L, 2L),
  source_shift = c(0, .5, 0, 0), variance_slope = c(.4, .4, .7, .4), quadratic_discrepancy = c(0, 0, 0, .5)
)
designs[, `:=`(nC = enrollment - nT, nR = enrollment)]
grid <- (-12:12) / 40
seed_base <- as.integer(opt('seed-base', if (phase == 'full') '74461005' else '73461005'))
seed_offset <- if (enrollment == 500L) 500000L else 0L
stopifnot(is.finite(seed_base), seed_base > 0, as.double(seed_base) + seed_offset + 400000 + reps < .Machine$integer.max)
gh <- function(n) {
  J <- matrix(0, n, n)
  J[cbind(1:(n - 1), 2:n)] <- sqrt(1:(n - 1))
  J <- J + t(J)
  z <- eigen(J, symmetric = TRUE)
  list(x = z$values, w = z$vectors[1, ]^2)
}
g128 <- gh(128L)
g256 <- gh(256L)
discrepancy_shape <- function(x, d) x + d$quadratic_discrepancy * (x^2 - 1)
true_odds <- function(x, d, ect = FALSE) {
  numerator <- if (ect) d$nT else d$nR
  exp(log(numerator / d$nE) - d$source_shift * x + d$source_shift^2 / 2)
}
true_share <- function(x, d) plogis(log(true_odds(x, d)) + log(d$nC / d$nR) - 2 * d$variance_slope * x)
calibration <- rbindlist(lapply(seq_len(nrow(designs)), function(i) {
  d <- designs[i]
  f <- function(g) {
    ext <- 1 - true_share(g$x, d)
    c(K0 = sum(g$w * ext), K1 = sum(g$w * ext * discrepancy_shape(g$x, d)))
  }
  a <- f(g256)
  b <- f(g128)
  c0 <- integrate(function(x) (1 - true_share(x, d)) * dnorm(x), -12, 12, abs.tol = 1e-12)$value
  c1 <- integrate(function(x) (1 - true_share(x, d)) * discrepancy_shape(x, d) * dnorm(x), -12, 12, abs.tol = 1e-12)$value
  error <- max(abs(c(a - b, a - c(c0, c1))))
  stopifnot(error < 1e-10)
  data.table(d,
    K0 = unname(a[1]), K1 = unname(a[2]), quadrature_error = error,
    zero_gap_intercept = -a[2] / a[1], ECT_shape_mean = 0
  )
}))
fwrite(calibration, file.path(out, 'calibration.csv'))
protocol <- list(
  stage = phase, repetitions_per_scenario = reps, scenarios = as.data.frame(designs), grid = grid,
  seed_base = seed_base, enrollment_seed_offset = seed_offset,
  seed = 'seed_base + enrollment_seed_offset + 100000 * scenario_id + replicate',
  original_baseline = list(nT = designs$nT[1], nC = designs$nC[1], nE = 2L * enrollment,
                           source_shift = 0, variance_slope = .4, quadratic_discrepancy = 0),
  factors = 'Exactly one factor differs from the original baseline in each scenario; settings fixed before fits.',
  covariates = 'Trial X~N10(0,I); external X~N10(mu*e1,I). Independent fixed-count groups.',
  means = 'mT=mC=f=.5+.6/sqrt(5)*sum(X1:5)+.25/sqrt(5)*sum(X1:5^2); mE=f-a-h(X1); h=x+q*(x^2-1). theta=0.',
  errors = 'Conditionally independent Gaussian; vT=1,vC=exp(gamma*X1),vE=exp(-gamma*X1).',
  source = 'omega_H=(nR/nE)*exp(-mu*x1+mu^2/2); omega_E=(nT/nE)*exp(-mu*x1+mu^2/2); p=omega/(1+omega).',
  share = 'lambda=omega_H*(1-e)/(omega_H*(1-e)+exp(2*gamma*x1)); e=nT/nR known.',
  gap = 'K0=E_trial(1-lambda),K1=E_trial[(1-lambda)h]; HCT a=(xi-K1)/K0; ECT a=xi because E_trial h=0; RCT actual gap zero.',
  fitting = 'Baseline fitting: separate quadratic21 OLS means; logistic11 source; separate log-link quasi-Gamma (1,X1) control variances with HC2 residual responses; correctly specified in all scenarios.',
  fitted_U = 'Complete seven-term unrestricted score and full analytic joint nuisance-target sandwich; fixed-count T/C/E centered meat; no clipping or fallback.',
  oracle_U = 'Same complete score at true source-specific means, design-specific source odds, and variances, including source and BOTH variance derivative terms; fixed-group score variance.',
  ECT = 'Only treated and external samples; published DR quadratic outcome/logistic source fit and plug-in EIF variance; fitted [.01,.99] probability guard retained and activations reported; oracle true odds are unguarded.',
  trial = 'Same observed trial outcomes, ANOVA/linear10X ANCOVA, HC2 SEs, normal 95% intervals; exactly reused in both modes.',
  pairing = 'Oracle and estimated modes receive identical generated covariates/errors and same HCT/ECT discrepancy calibration; mode pairing retained for differences; one sample hash per replicate.',
  evaluation = 'Genuine Monte Carlo draws; exact affine intercept-shift estimates/influences and quadratic variance paths evaluated on 25-point grid; direct refits independently check representative gaps.',
  failures = 'All failures and warnings saved; no silent exclusion, replacement or incomplete full summary.',
  MCSE = 'Sample SD of replicate contributions divided by sqrt(B), including paired mode differences.',
  interpretation = 'Oracle complete-score evaluation is not an efficiency bound for a model with known nuisance functions. Zero aggregate gap retains nonconstant source discrepancy.',
  created_before_fits = format(Sys.time(), tz = 'UTC', usetz = TRUE)
)
jsonlite::write_json(protocol, file.path(out, 'protocol.json'), auto_unbox = TRUE, pretty = TRUE, digits = 16)

generate <- function(d, b) {
  seed <- seed_base + seed_offset + 100000L * d$id + b
  set.seed(seed)
  xe <- matrix(rnorm(d$nE * 10L), ncol = 10L)
  xe[, 1] <- xe[, 1] + d$source_shift
  epsE <- rnorm(d$nE)
  xt <- matrix(rnorm(d$nT * 10L), ncol = 10L)
  xc <- matrix(rnorm(d$nC * 10L), ncol = 10L)
  epsT <- rnorm(d$nT)
  epsC <- rnorm(d$nC)
  X <- rbind(xt, xc, xe)
  x <- X[, 1]
  f <- .5 + (.6 / sqrt(5)) * rowSums(X[, 1:5]) + (.25 / sqrt(5)) * rowSums(X[, 1:5]^2)
  vC <- exp(d$variance_slope * x)
  vE <- exp(-d$variance_slope * x)
  T <- c(rep(1, d$nT), rep(0, d$nC + d$nE))
  C <- c(rep(0, d$nT), rep(1, d$nC), rep(0, d$nE))
  E <- 1 - T - C
  R <- T + C
  residual <- c(epsT, exp(d$variance_slope * xc[, 1] / 2) * epsC, exp(-d$variance_slope * xe[, 1] / 2) * epsE)
  data <- setNames(as.data.frame(X), paste0('X', 1:10))
  data$R <- R
  data$A <- T
  h <- discrepancy_shape(x, d)
  data$Y <- f + residual - E * h
  list(
    data = data, x = x, f = f, h = h, residual = residual, vC = vC, vE = vE, T = T, C = C, E = E, R = R,
    group = factor(ifelse(E == 1, 'E', ifelse(T == 1, 'T', 'C')), levels = c('T', 'C', 'E')), seed = seed
  )
}
shift <- function(data, a) {
  data$Y <- data$Y - (1 - data$R) * a
  data
}
gcov <- function(x, y, group, denominator) sum(vapply(levels(group), function(k) {
  ii <- group == k
  sum(ii) * cov(x[ii], y[ii])
}, numeric(1))) / denominator^2
oracle_U_raw <- function(g, d, a) {
  omega <- true_odds(g$x, d)
  p <- omega / (1 + omega)
  L <- true_share(g$x, d)
  H <- L * (1 - L)
  delta <- a + g$h
  e <- d$nT / d$nR
  ac <- 1 - e
  g$R * (1 - L) * delta + g$T / e * g$residual - g$C * L / ac * g$residual - g$E * (1 - L) * omega * g$residual -
    delta * H * (g$R - p) / (1 - p) + delta * H * g$C / ac * (g$residual^2 / g$vC - 1) -
    delta * H * g$E * omega * (g$residual^2 / g$vE - 1)
}
oracle_paths <- function(g, d) {
  raw0 <- oracle_U_raw(g, d, 0)
  raw1 <- oracle_U_raw(g, d, 1)
  slope <- raw1 - raw0
  U <- list(
    intercept = sum(raw0) / d$nR, slope = sum(slope) / d$nR,
    variance = c(gcov(raw0, raw0, g$group, d$nR), 2 * gcov(raw0, slope, g$group, d$nR), gcov(slope, slope, g$group, d$nR))
  )
  odds <- true_odds(g$x, d, TRUE)
  rawE <- g$T * (g$h + g$residual) - g$E * odds * g$residual
  estE <- sum(rawE) / d$nT
  phiE <- rawE - g$T * estE
  ECT <- list(intercept = estE, slope = 1, variance = c(sum(phiE^2) / d$nT^2, 0, 0))
  list(U = U, 'AIPW (ECT)' = ECT)
}
coefficient_row <- function(path, mode, method, d, b, seed) data.table(
  scenario = d$scenario,
  scenario_id = d$id, mode = mode, method = method, replicate = b, seed = seed,
  intercept = unname(path$intercept), slope = unname(path$slope),
  v0 = unname(path$variance[1]), v1 = unname(path$variance[2]), v2 = unname(path$variance[3])
)
one_dataset <- function(job) {
  d <- designs[id == job$id]
  cal <- calibration[id == job$id]
  b <- job$replicate
  g <- generate(d, b)
  data <- g$data
  f0 <- fg_fit(data, helpers, 'quadratic', 'x1', keep = b <= 2, methods = 'U')$U
  f1 <- fg_fit(shift(data, 1), helpers, 'quadratic', 'x1', methods = 'U')$U
  upath <- fg_path(f0, f1)
  ect <- scalar_ect_fit(data, ect_helpers$sim)
  paths <- list(U = upath, 'AIPW (ECT)' = list(intercept = ect$estimate, slope = 1, variance = c(ect$se^2, 0, 0)))
  trial <- hct_trial_wald_comparators(data, adjustment)
  oracle <- oracle_paths(g, d)
  for (m in names(trial)) {
    tp <- list(intercept = trial[[m]]$estimate, slope = 0, variance = c(trial[[m]]$se^2, 0, 0))
    paths[[m]] <- tp
    oracle[[m]] <- tp
  }
  coefficients <- rbindlist(c(
    lapply(names(paths), function(m) coefficient_row(paths[[m]], 'Estimated nuisance', m, d, b, g$seed)),
    lapply(names(oracle), function(m) coefficient_row(oracle[[m]], 'Oracle nuisance', m, d, b, g$seed))
  ))
  direct <- list()
  if (b <= 2L) {
    # Verify the independently written oracle expression against the same
    # complete-score system at true coefficients, including shifted p(x).
    prep <- f0$prepared
    pa <- prep$parameters
    ix <- prep$ix
    base <- c(.5, rep(.6 / sqrt(5), 5), rep(0, 5), rep(.25 / sqrt(5), 5), rep(0, 5))
    pa[ix$bT] <- base
    pa[ix$bC] <- base
    pa[ix$bE] <- base
    pa[ix$bE[1]] <- pa[ix$bE[1]] + d$quadratic_discrepancy
    pa[ix$bE[2]] <- pa[ix$bE[2]] - 1
    pa[ix$bE[12]] <- pa[ix$bE[12]] - d$quadratic_discrepancy
    pa[ix$eta] <- c(log(d$nR / d$nE) + d$source_shift^2 / 2, -d$source_shift, rep(0, 9))
    pa[ix$gC] <- c(0, d$variance_slope)
    pa[ix$gE] <- c(0, -d$variance_slope)
    oracle_error <- max(abs(fg_system(prep, pa, 'U', FALSE)$raw - oracle_U_raw(g, d, 0)))
    for (xi in c(-.3, 0, .3)) for (m in c('U', 'AIPW (ECT)')) {
      a <- if (m == 'U') (xi - cal$K1) / cal$K0 else xi
      pth <- paths[[m]]
      dd <- shift(data, a)
      actual <- if (m == 'U') fg_fit(dd, helpers, 'quadratic', 'x1', methods = 'U')$U else scalar_ect_fit(dd, ect_helpers$sim)
      predicted <- pth$intercept + a * pth$slope
      se <- sqrt(pth$variance[1] + a * pth$variance[2] + a^2 * pth$variance[3])
      IFerror <- if (m == 'U') max(abs(actual$influence - pth$influence0 - a * pth$influence_slope)) else max(abs(actual$influence - ect$influence))
      oraw <- oracle_U_raw(g, d, a)
      op <- oracle$U
      direct[[length(direct) + 1L]] <- data.table(
        scenario = d$scenario, replicate = b, method = m, xi = xi,
        estimate_error = abs(actual$estimate - predicted), SE_error = abs(actual$se - se), IF_error = IFerror,
        oracle_system_error = oracle_error,
        oracle_path_estimate_error = abs(sum(oraw) / d$nR - (op$intercept + a * op$slope)),
        oracle_path_variance_error = abs(gcov(oraw, oraw, g$group, d$nR) - (op$variance[1] + a * op$variance[2] + a^2 * op$variance[3]))
      )
    }
    direct <- rbindlist(direct)
    stopifnot(max(as.matrix(direct[, -c('scenario', 'replicate', 'method', 'xi')])) < 1e-7)
  }
  stopifnot(max(abs(f0$lambda - f1$lambda)) < 1e-8)
  list(
    success = TRUE, coefficients = coefficients, direct = rbindlist(list(direct), fill = TRUE),
    pairing = data.table(
      scenario = d$scenario, replicate = b, seed = g$seed,
      full_sample_sha256 = digest::digest(data, algo = 'sha256'), ECT_sample_sha256 = digest::digest(data[data$A == 1 | data$R == 0, ], algo = 'sha256')
    ),
    diagnostics = data.table(
      scenario = d$scenario, replicate = b, share_shift_error = max(abs(f0$lambda - f1$lambda)),
      nuisance_equation_max = max(f0$nuisance_equation_max, f1$nuisance_equation_max),
      lambda_trial_mean = mean(f0$lambda[data$R == 1]), lambda_min = min(f0$lambda), lambda_max = max(f0$lambda),
      p_min = f0$source_probability_range[1], p_max = f0$source_probability_range[2],
      ECT_clipped = ect$external_clipped + ect$treated_clipped, ECT_ESS = ect$external_ess,
      external_x1_mean = mean(g$x[g$E == 1]), trial_x1_mean = mean(g$x[g$R == 1])
    ),
    warnings = unique(c(f0$warnings, f1$warnings, ect$warnings))
  )
}
one_job <- function(job) {
  if (file.exists(job$file)) return(job$file)
  ans <- tryCatch(one_dataset(job), error = function(e) list(success = FALSE, error = conditionMessage(e)))
  ans$id <- job$id
  ans$replicate <- job$replicate
  ans$seed <- seed_base + seed_offset + 100000L * job$id + job$replicate
  saveRDS(ans, paste0(job$file, '.partial'), compress = FALSE)
  stopifnot(file.rename(paste0(job$file, '.partial'), job$file))
  job$file
}
script <- normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]), winslash = '/')
input_paths <- c(script, list.files(file.path(paper, 'lib'), full.names = TRUE, pattern = '[.]R$'))
names(input_paths) <- substring(input_paths, nchar(paper) + 2L)
hashes <- vapply(input_paths, function(p) digest::digest(file = p, algo = 'sha256'), character(1))
key <- digest::digest(list(hashes, designs, calibration, seed_base, R.version.string), algo = 'sha256')
checkdir <- file.path(out, 'checkpoints', substr(key, 1, 20))
dir.create(checkdir, recursive = TRUE, showWarnings = FALSE)
jobs <- list()
for (id in designs$id) for (b in seq_len(reps)) jobs[[length(jobs) + 1L]] <- list(id = id, replicate = b, file = file.path(checkdir, sprintf('s%d_rep%05d.rds', id, b)))
start <- Sys.time()
cl <- parallel::makePSOCKcluster(workers, outfile = file.path(out, 'workers.log'))
parallel::clusterExport(cl, c(
  'paper', 'helper', 'ect_helper', 'trial_helper', 'sim_source', 'thread_vars', 'designs', 'calibration',
  'adjustment', 'seed_base', 'seed_offset', 'discrepancy_shape', 'true_odds', 'true_share', 'generate', 'shift', 'gcov', 'oracle_U_raw', 'oracle_paths',
  'coefficient_row', 'one_dataset', 'one_job'
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
cat('Four extensions:', length(jobs), 'new paired samples;', workers, 'workers.\n')
flush.console()
files <- tryCatch(parallel::parLapplyLB(cl, jobs, one_job), finally = parallel::stopCluster(cl))
ans <- lapply(unlist(files), readRDS)
failures <- rbindlist(lapply(ans, function(x) if (!x$success) data.table(scenario_id = x$id, replicate = x$replicate, seed = x$seed, error = x$error) else NULL))
if (!ncol(failures)) failures <- data.table(scenario_id = integer(), replicate = integer(), seed = integer(), error = character())
fwrite(failures, file.path(out, 'failures.csv'))
good <- Filter(function(x) x$success, ans)
co <- rbindlist(lapply(good, '[[', 'coefficients'))
fwrite(co, file.path(out, 'replicate_coefficients.csv.gz'), compress = 'gzip')
diagnostics <- rbindlist(lapply(good, '[[', 'diagnostics'))
fwrite(diagnostics, file.path(out, 'fit_diagnostics.csv'))
direct <- rbindlist(lapply(good, '[[', 'direct'))
fwrite(direct, file.path(out, 'direct_refit_checks.csv'))
pairing <- rbindlist(lapply(good, '[[', 'pairing'))
fwrite(pairing, file.path(out, 'sample_pairing.csv'))
warnings <- rbindlist(lapply(good, function(x) if (length(x$warnings)) data.table(scenario_id = x$id, replicate = x$replicate, warning = x$warnings) else NULL))
if (!ncol(warnings)) warnings <- data.table(scenario_id = integer(), replicate = integer(), warning = character())
fwrite(warnings, file.path(out, 'warnings.csv'))
if (nrow(failures)) stop('Failures saved; no incomplete summary generated.')
stopifnot(nrow(co) == 4L * reps * 8L, !anyDuplicated(co[, .(scenario, mode, method, replicate)]))
rm(ans, good)
gc()
summaries <- paired <- analytic <- list()
metrics <- c('coverage_theta', 'bias_theta', 'mse_theta', 'width', 'coverage_TCE', 'bias_TCE', 'estimated_variance')
for (scenario_name in designs$scenario) {
  d <- designs[scenario == scenario_name]
  cal <- calibration[scenario == scenario_name]
  rawrows <- list()
  for (mode_name in c('Oracle nuisance', 'Estimated nuisance')) for (method_name in c('U', 'AIPW (ECT)', 'ANOVA (RCT)', 'ANCOVA (RCT)')) for (xi in grid) {
    z <- co[scenario == scenario_name & mode == mode_name & method == method_name]
    trial <- grepl('RCT', method_name)
    a <- if (trial) 0 else if (method_name == 'U') (xi - cal$K1) / cal$K0 else xi
    target <- if (trial) 0 else xi
    estimate <- z$intercept + a * z$slope
    V <- z$v0 + a * z$v1 + a^2 * z$v2
    stopifnot(length(estimate) == reps, all(is.finite(c(estimate, V))), all(V > 0))
    half <- qnorm(.975) * sqrt(V)
    rr <- data.table(
      scenario = scenario_name, mode = mode_name, method = method_name, replicate = z$replicate, seed = z$seed,
      xi = xi, xi_index = as.integer(round(xi * 40)), amplitude = a, achieved_gap = target,
      coverage_theta = as.integer(abs(estimate) <= half), bias_theta = estimate, mse_theta = estimate^2, width = 2 * half,
      coverage_TCE = as.integer(abs(estimate - target) <= half), bias_TCE = estimate - target, estimated_variance = V
    )
    rawrows[[length(rawrows) + 1L]] <- rr
    r <- list(
      scenario = scenario_name, mode = mode_name, method = method_name, xi = xi, xi_index = as.integer(round(xi * 40)),
      amplitude = a, achieved_gap = target, nT = d$nT, nC = d$nC, nE = d$nE, N_trial = d$nR, n_reps = reps, n_failed = 0L
    )
    for (m in metrics) {
      r[[m]] <- mean(rr[[m]])
      r[[paste0(m, '_mcse')]] <- sd(rr[[m]]) / sqrt(reps)
    }
    r$empirical_variance <- var(estimate)
    summaries[[length(summaries) + 1L]] <- as.data.table(r)
    if (mode_name == 'Oracle nuisance' && !trial) {
      if (method_name == 'U') {
        group_means <- group_vars <- numeric(3)
        for (j in 1:3) {
          x <- g256$x + if (j == 3) d$source_shift else 0
          w <- true_odds(x, d)
          L <- true_share(x, d)
          H <- L * (1 - L)
          delta <- a + discrepancy_shape(x, d)
          cmean <- if (j == 3) delta * H * w else delta * (1 - L)^2
          lin <- switch(j,
            rep(d$nR / d$nT, length(x)),
            -L / (d$nC / d$nR) * exp(d$variance_slope * x / 2),
            -(1 - L) * w * exp(-d$variance_slope * x / 2)
          )
          quad <- switch(j,
            rep(0, length(x)),
            delta * H / (d$nC / d$nR),
            -delta * H * w
          )
          group_means[j] <- sum(g256$w * cmean)
          group_vars[j] <- sum(g256$w * (cmean^2 + lin^2 + 2 * quad^2)) - group_means[j]^2
        }
        expected <- sum(c(d$nT, d$nC, d$nE) * group_means) / d$nR
        avar <- sum(c(d$nT, d$nC, d$nE) * group_vars) / d$nR^2
      } else {
        x <- g256$x + d$source_shift
        w <- true_odds(x, d, TRUE)
        expected <- xi
        avar <- (2 + 2 * d$quadratic_discrepancy^2) / d$nT + d$nE * sum(g256$w * w^2 * exp(-d$variance_slope * x)) / d$nT^2
      }
      stopifnot(abs(expected - xi) < 1e-10)
      analytic[[length(analytic) + 1L]] <- data.table(
        scenario = scenario_name, method = method_name, xi = xi,
        expected_mean = expected, expected_variance = avar, empirical_variance = var(estimate),
        mean_estimated_variance = mean(V), bias_z = (mean(estimate) - expected) / (sd(estimate) / sqrt(reps))
      )
    }
  }
  raw <- rbindlist(rawrows)
  fwrite(raw, file.path(out, paste0(scenario_name, '_replicates.csv.gz')), compress = 'gzip')
  keys <- c('scenario', 'method', 'replicate', 'xi_index', 'xi', 'seed')
  pp <- merge(raw[mode == 'Estimated nuisance'], raw[mode == 'Oracle nuisance'], by = keys, suffixes = c('_estimated', '_oracle'))
  stopifnot(nrow(pp) == reps * 4L * 25L, all(pp$achieved_gap_estimated == pp$achieved_gap_oracle))
  for (m in metrics) pp[, (m) := get(paste0(m, '_estimated')) - get(paste0(m, '_oracle'))]
  paired[[scenario_name]] <- pp[,
    {
      r <- list(n_pairs = .N)
      for (m in metrics) {
        r[[paste0(m, '_difference')]] <- mean(.SD[[m]])
        r[[paste0(m, '_difference_mcse')]] <- sd(.SD[[m]]) / sqrt(.N)
      }
      r
    },
    by = .(scenario, method, xi_index, xi),
    .SDcols = metrics
  ]
  rm(raw, rawrows, pp)
  gc()
}
summary <- rbindlist(summaries)
analytical <- rbindlist(analytic)
pairs <- rbindlist(paired)
fwrite(summary, file.path(out, 'summary.csv'))
fwrite(summary[xi_index %in% c(-12L, -6L, 0L, 6L, 12L)], file.path(out, 'key_gaps.csv'))
fwrite(pairs, file.path(out, 'paired_mode_differences.csv'))
fwrite(analytical, file.path(out, 'analytic_checks.csv'))
stopifnot(
  nrow(summary) == 4L * 2L * 4L * 25L, all(summary$n_reps == reps),
  identical(hashes, vapply(input_paths, function(p) digest::digest(file = p, algo = 'sha256'), character(1)))
)
trial_differences <- pairs[grepl('RCT', method)]
stopifnot(max(abs(as.matrix(trial_differences[, grep('_difference', names(trial_differences)), with = FALSE]))) == 0)
validation <- list(
  success = TRUE, phase = phase, repetitions_per_scenario = reps, scenarios = 4, paired_datasets = nrow(pairing),
  summary_rows = nrow(summary), failures = nrow(failures), warnings = nrow(warnings),
  source_files_unchanged = TRUE, source_hashes = as.list(hashes), checkpoint_key = key,
  calibration_max_error = max(calibration$quadrature_error),
  direct_estimate_max_error = max(direct$estimate_error), direct_SE_max_error = max(direct$SE_error),
  direct_IF_max_error = max(direct$IF_error), oracle_system_max_error = max(direct$oracle_system_error),
  oracle_path_estimate_max_error = max(direct$oracle_path_estimate_error), oracle_path_variance_max_error = max(direct$oracle_path_variance_error),
  max_oracle_bias_z = max(abs(analytical$bias_z)), ECT_total_probability_guards = sum(diagnostics$ECT_clipped),
  trial_reference_identical_between_modes = TRUE, elapsed_seconds = as.numeric(difftime(Sys.time(), start, units = 'secs'))
)
jsonlite::write_json(validation, file.path(out, 'validation.json'), auto_unbox = TRUE, pretty = TRUE, digits = 16)
print(summary[xi_index == 0, .(scenario, mode, method, coverage_theta, bias_theta, mse_theta, width, coverage_TCE)])
print(validation)
