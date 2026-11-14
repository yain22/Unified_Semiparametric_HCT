# Independent deterministic checks. Runs from any directory, writes only --output.
script <- normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]), winslash = '/')
here <- dirname(script)
source(file.path(here, 'lib/cli_guard.R'))
args <- commandArgs(trailingOnly = TRUE)
opt <- function(n, d) {
  if (n == 'output' && any(grepl('^--output-dir=', args))) n <- 'output-dir'
  h <- grep(paste0('^--', n, '='), args, value = TRUE)
  if (length(h)) sub(paste0('^--', n, '='), '', h) else d
}
out <- simulation_guard_output(opt('output', file.path(here, '../../results/validation')), here)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
library(data.table)
setDTthreads(1L)
source(file.path(here, 'lib/fitted_helpers.R'))
helpers <- fg_load(here)
source(file.path(here, 'lib/scalar_gap_ect_helpers.R'))
eh <- scalar_ect_load(here, file.path(here, 'lib/ect_sim_source.R'))
source(file.path(here, 'lib/trial_baselines.R'))
functions_only <- function(path) {
  env <- new.env(parent = .GlobalEnv)
  for (ex in parse(path)) if (is.call(ex) && identical(ex[[1]], as.name('<-')) && is.call(ex[[3]]) && identical(ex[[3]][[1]], as.name('function'))) eval(ex, envir = env)
  env
}
bg <- functions_only(file.path(here, 'run_baseline.R'))
bg$seed_base <- 52461004L
bg$designs <- data.frame(N_trial = c(250L, 500L), NT = c(167L, 333L), NC = c(83L, 167L), N0 = c(500L, 1000L))
eg <- functions_only(file.path(here, 'run_extended.R'))
eg$seed_base <- 74461005L
eg$seed_offset <- 0L  # The four extended validation samples use trial enrollment 250.
ed <- data.table(id = 1:4, scenario = c('fewer_external', 'source_shift', 'stronger_heteroskedasticity', 'quadratic_discrepancy'), nT = 167L, nC = 83L, nE = c(250L, 500L, 500L, 500L), source_shift = c(0, .5, 0, 0), variance_slope = c(.4, .4, .7, .4), quadratic_discrepancy = c(0, 0, 0, .5))
ed[, nR := nT + nC]
samples <- list(baseline250 = bg$generate(250, 1)$data, baseline500 = bg$generate(500, 1)$data)
for (j in 1:4) samples[[ed$scenario[j]]] <- eg$generate(ed[j], 1)$data
checks <- list()
failures <- list()
add <- function(name, value, tol = 1e-7) {
  checks[[length(checks) + 1L]] <<- data.table(check = name, error = as.numeric(value), tolerance = tol, pass = is.finite(value) && value <= tol)
  stopifnot(is.finite(value), value <= tol)
}
for (nm in names(samples)) {
  dat <- samples[[nm]]
  fit <- fg_fit(dat, helpers, 'quadratic', 'x1', keep = TRUE)$U
  d <- fit$prepared
  par <- fit$parameters
  # Complete seven-term expression, with no helper algebra reused.
  s <- fit$state
  L <- s$L
  H <- L * (1 - L)
  D <- s$mC - s$mE
  Y <- dat$Y
  a <- 1 - d$e
  direct <- d$S * (s$mT - L * s$mC - (1 - L) * s$mE) + d$T / d$e * (Y - s$mT) - d$C * L / a * (Y - s$mC) - d$E * (1 - L) * s$w * (Y - s$mE) - D * H * (d$S - s$p) / (1 - s$p) + D * H * d$C / a * ((Y - s$mC)^2 / s$vC - 1) - D * H * d$E * s$w * ((Y - s$mE)^2 / s$vE - 1)
  add(paste(nm, 'seven_term_score'), max(abs(direct - s$raw)), 1e-11)
  # All Jacobian columns at a perturbed (nonroot) parameter value.
  pert <- par + .017 * sin(seq_along(par))
  analytical <- fg_system(d, pert, 'U')$jacobian
  numericJ <- sapply(seq_along(par), function(k) {
    h <- 1e-5 * (1 + abs(pert[k]))
    up <- dn <- pert
    up[k] <- up[k] + h
    dn[k] <- dn[k] - h
    (fg_system(d, up, 'U', FALSE)$mean - fg_system(d, dn, 'U', FALSE)$mean) / (2 * h)
  })
  add(paste(nm, 'full_jacobian'), max(abs(numericJ - analytical)))
  # Full inverse of joint equations; independent of the block-elimination helper.
  allIF <- -s$equations %*% t(solve(s$jacobian))
  u <- allIF[, d$ix$theta]
  add(paste(nm, 'joint_influence'), max(abs(u - fit$influence)), 1e-10)
  meat <- matrix(0, ncol(s$equations), ncol(s$equations))
  for (gr in levels(d$group)) {
    z <- s$equations[d$group == gr, , drop = FALSE]
    z <- sweep(z, 2, colMeans(z), '-')
    meat <- meat + crossprod(z) * nrow(z) / (nrow(z) - 1)
  }
  bread <- solve(s$jacobian)
  V <- bread %*% meat %*% t(bread) / d$N^2
  add(paste(nm, 'group_centered_sandwich'), abs(sqrt(V[d$ix$theta, d$ix$theta]) - fit$se), 1e-11)
  # Independent HC2 matrix sandwich in the trial-only design.
  td <- dat[dat$S == 1, ]
  X <- model.matrix(~ A + ., td[, c('A', paste0('X', 1:10))])
  b <- solve(crossprod(X))
  r <- td$Y - X %*% (b %*% crossprod(X, td$Y))
  h <- rowSums((X %*% b) * X)
  Vt <- b %*% crossprod(X, X * as.numeric(r^2 / (1 - h))) %*% b
  tf <- hct_trial_wald_fit(dat, as.formula(paste('~', paste(paste0('X', 1:10), collapse = '+'))))
  add(paste(nm, 'trial_HC2'), abs(sqrt(Vt['A', 'A']) - tf$se), 1e-11)
  ect <- scalar_ect_fit(dat, eh$sim)
  tamper <- dat
  ii <- dat$S == 1 & dat$A == 0
  tamper$Y[ii] <- tamper$Y[ii] + 1000
  tamper$X1[ii] <- tamper$X1[ii] - 100
  unaffected <- scalar_ect_fit(tamper, eh$sim)
  add(paste(nm, 'ECT_excludes_internal'), max(abs(c(ect$estimate - unaffected$estimate, ect$se - unaffected$se))), 0)
  add(paste(nm, 'ECT_plugin'), abs(ect$se^2 - sum(ect$influence^2) / length(ect$influence)^2), 1e-12)
  # Failures must surface, never fall back to a different estimator.
  bad <- dat
  bad$X2 <- bad$X1
  failed <- inherits(try(fg_fit(bad, helpers, 'quadratic', 'x1'), silent = TRUE), 'try-error')
  failures[[length(failures) + 1L]] <- data.table(sample = nm, test = 'rank_deficient_mean', correctly_failed = failed)
  stopifnot(failed)
}
fwrite(rbindlist(checks), file.path(out, 'independent_checks.csv'))
fwrite(rbindlist(failures), file.path(out, 'failure_checks.csv'))
writeLines(capture.output(sessionInfo()), file.path(out, 'sessionInfo.txt'))
jsonlite::write_json(list(success = TRUE, samples = names(samples), checks = length(checks), maximum_jacobian_error = max(rbindlist(checks)[grepl('full_jacobian', check), error]), failure_tests = length(failures), scope = 'Independent analytic-score, numerical-Jacobian, full-inverse sandwich, HC2, retained ECT sample and rank-deficiency checks. Small deterministic checks; not a 5000-repetition rerun.'), file.path(out, 'validation.json'), auto_unbox = TRUE, pretty = TRUE, digits = 16)
print(rbindlist(checks))
