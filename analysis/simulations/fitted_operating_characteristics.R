# One-sided operating characteristics from the original fitted simulation draws.
# Run from any directory. No source results or manuscript files are overwritten.
# Optional: --baseline=DIR --quadratic=FILE --output=DIR --alternative=0.3

script <- normalizePath(sub('^--file=', '',
  grep('^--file=', commandArgs(FALSE), value=TRUE)[1]), winslash='/')
here <- dirname(script)
root <- normalizePath(file.path(here, '../..'), winslash='/')
args <- commandArgs(trailingOnly=TRUE)
option <- function(name, default) {
  hit <- grep(paste0('^--', name, '='), args, value=TRUE)
  if (length(hit) > 1L) stop('Repeated option: ', name)
  if (length(hit)) sub(paste0('^--', name, '='), '', hit) else default
}
if (any(!grepl('^--(baseline|quadratic|output|alternative)=', args)))
  stop('Options are --baseline=DIR --quadratic=FILE --output=DIR --alternative=NUMBER')
baseline <- normalizePath(option('baseline', file.path(root, 'results/baseline')),
  winslash='/', mustWork=TRUE)
quadratic_file <- normalizePath(option('quadratic', file.path(root,
  'analysis/comparators/quadratic_rct_replicates.csv.gz')), winslash='/', mustWork=TRUE)
output <- option('output', file.path(root, 'results/one-sided'))
dir.create(output, recursive=TRUE, showWarnings=FALSE)
output <- normalizePath(output, winslash='/', mustWork=TRUE)
if (identical(output, baseline) || identical(output, dirname(quadratic_file)))
  stop('Choose an output folder separate from the archived simulation inputs.')
theta_star <- as.numeric(option('alternative', '0.3'))
stopifnot(length(theta_star)==1L, is.finite(theta_star), theta_star>0)

library(data.table)
setDTthreads(1L)
thread_vars <- c('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')
do.call(Sys.setenv, as.list(setNames(rep('1', length(thread_vars)), thread_vars)))
paper <- here
source(file.path(paper, 'lib/fitted_helpers.R'))
helpers <- fg_load(paper)
source(file.path(paper, 'lib/scalar_gap_ect_helpers.R'))
ect_helpers <- scalar_ect_load(paper, file.path(paper, 'lib/ect_sim_source.R'))
source(file.path(paper, 'lib/trial_baselines.R'))

input_files <- c(file.path(baseline, c('replicates.csv.gz',
  'transportability_reference.csv.gz', 'protocol.json', 'validation.json',
  'calibration.csv', 'sample_pairing.csv')), quadratic_file,
  file.path(paper, 'run_baseline.R'),
  list.files(file.path(paper, 'lib'), pattern='[.]R$', full.names=TRUE))
hashes <- vapply(input_files, function(f) digest::digest(file=f, algo='sha256'), character(1))
protocol <- jsonlite::read_json(file.path(baseline, 'protocol.json'), simplifyVector=TRUE)
validation <- jsonlite::read_json(file.path(baseline, 'validation.json'), simplifyVector=TRUE)
stopifnot(isTRUE(validation$success), validation$failures==0,
  protocol$repetitions_per_size==5000L)
designs <- as.data.frame(protocol$counts)
seed_base <- as.integer(protocol$seed_base)
calibration <- fread(file.path(baseline, 'calibration.csv'))
pairing <- fread(file.path(baseline, 'sample_pairing.csv'))
methods <- c('U', 'AIPW (ECT)', 'ANOVA (RCT)', 'ANCOVA (RCT)', 'Quadratic ANCOVA (RCT)')
method_labels <- c('AIPW-U (HCT)', 'AIPW (ECT)', 'ANOVA (RCT)',
  'Linear ANCOVA (RCT)', 'Quadratic ANCOVA (RCT)')
names(method_labels) <- methods

columns <- c('N_trial', 'method', 'replicate', 'seed', 'xi_index', 'xi',
  'estimate', 'se', 'target_TCE')
raw <- fread(file.path(baseline, 'replicates.csv.gz'), select=columns)
raw <- raw[method %in% methods]
reference <- fread(file.path(baseline, 'transportability_reference.csv.gz'),
  select=setdiff(columns, c('xi_index', 'xi')))
reference <- reference[method %in% methods]
quadratic <- fread(quadratic_file,
  select=c('N_trial', 'method', 'replicate', 'seed', 'estimate', 'se'))
stopifnot(nrow(raw)==1000000L, nrow(reference)==40000L, nrow(quadratic)==10000L,
  !anyDuplicated(raw[, .(N_trial, method, replicate, xi_index)]),
  !anyDuplicated(reference[, .(N_trial, method, replicate)]),
  !anyDuplicated(quadratic[, .(N_trial, method, replicate)]))

# Import only the generator definition; sourcing the driver would rerun it.
for (expression in parse(file.path(paper, 'run_baseline.R'))) {
  if (is.call(expression) && identical(expression[[1]], as.name('<-')) &&
      identical(expression[[2]], as.name('generate')))
    eval(expression)
}
stopifnot(exists('generate', inherits=FALSE))
linear <- as.formula(paste('~', paste(paste0('X', 1:10), collapse='+')))
quadratic_basis <- as.formula(paste('~',
  paste(c(paste0('X', 1:10), paste0('I(X', 1:10, '^2)')), collapse='+')))
fit_methods <- function(data) {
  list('U'=fg_fit(data, helpers, mean_basis='quadratic',
    variance_basis='x1', methods='U')$U,
    'AIPW (ECT)'=scalar_ect_fit(data, ect_helpers$sim),
    'ANOVA (RCT)'=hct_trial_wald_fit(data, ~1),
    'ANCOVA (RCT)'=hct_trial_wald_fit(data, linear),
    'Quadratic ANCOVA (RCT)'=hct_trial_wald_fit(data, quadratic_basis))
}

# A constant treatment effect shifts treated outcomes only. OLS contains the
# treatment indicator/intercept; control fits, residuals, and source fits stay
# fixed. Check this directly, including the complete HCT sandwich, before
# using a translated estimate to evaluate the alternative.
checks <- list()
index <- 0L
for (n in c(250L, 500L)) for (b in c(1L, 2500L, 5000L)) {
  gen <- generate(n, b)
  saved_pair <- pairing[N_trial==n & replicate==b]
  stopifnot(nrow(saved_pair)==1L,
    digest::digest(gen$data, algo='sha256')==saved_pair$full_sample_sha256)
  cal <- calibration[N_trial==n]
  for (case in c('matched_gap', 'transportability')) {
    gaps <- if (case=='matched_gap') c(-.3, -.1, 0, .1, .3) else 0
    for (gap in gaps) {
      hct_data <- gen$data
      ect_data <- gen$data
      if (case=='matched_gap') {
        hct_data$Y <- hct_data$Y-(1-hct_data$S)*(gap-cal$K1)/cal$K0
        ect_data$Y <- ect_data$Y-(1-ect_data$S)*gap
      } else {
        hct_data$Y <- hct_data$Y+(1-hct_data$S)*hct_data$X1
        ect_data <- hct_data
      }
      null <- fit_methods(hct_data)
      null[['AIPW (ECT)']] <- scalar_ect_fit(ect_data, ect_helpers$sim)
      hct_data$Y <- hct_data$Y+theta_star*hct_data$S*hct_data$A
      ect_data$Y <- ect_data$Y+theta_star*ect_data$S*ect_data$A
      alternative <- fit_methods(hct_data)
      alternative[['AIPW (ECT)']] <- scalar_ect_fit(ect_data, ect_helpers$sim)
      for (method_name in methods) {
        if (method_name=='Quadratic ANCOVA (RCT)') {
          old <- quadratic[N_trial==n & replicate==b]
        } else if (case=='transportability') {
          old <- reference[N_trial==n & replicate==b & method==method_name]
        } else {
          old <- raw[N_trial==n & replicate==b & method==method_name & abs(xi-gap)<1e-10]
        }
        stopifnot(nrow(old)==1L)
        a <- null[[method_name]]
        z <- alternative[[method_name]]
        if_error <- if (!is.null(a$influence)) max(abs(z$influence-a$influence)) else NA_real_
        index <- index+1L
        checks[[index]] <- data.table(N_trial=n, replicate=b, seed=gen$seed,
          case=case, gap=gap, method=method_name, alternative=theta_star,
          archived_estimate_error=abs(a$estimate-old$estimate),
          archived_se_error=abs(a$se-old$se),
          translation_estimate_error=abs(z$estimate-a$estimate-theta_star),
          translation_se_error=abs(z$se-a$se),
          translation_influence_error=if_error)
      }
    }
  }
  cat('Verified original draws and treatment shift: nRCT=', n, ', replicate=', b, '\n', sep='')
  flush.console()
}
checks <- rbindlist(checks)
errors <- unlist(checks[, .(archived_estimate_error, archived_se_error,
  translation_estimate_error, translation_se_error, translation_influence_error)])
stopifnot(max(errors, na.rm=TRUE)<1e-8)
fwrite(checks, file.path(output, 'treatment-shift-validation.csv'))

# Use saved estimates and their own saved SEs, never two-sided coverage.
qgrid <- rbindlist(lapply(-12:12, function(k) {
  q <- copy(quadratic)
  q[, `:=`(xi_index=k, xi=k/40, target_TCE=0)]
  q
}))
qreference <- copy(quadratic)
qreference[, target_TCE:=0]
raw <- rbindlist(list(raw, qgrid), use.names=TRUE)
reference <- rbindlist(list(reference, qreference), use.names=TRUE)
raw[, case:='matched_gap']
reference[, `:=`(case='transportability', xi_index=0L, xi=0)]
draws <- rbindlist(list(raw, reference), use.names=TRUE)
stopifnot(all(is.finite(draws$estimate)), all(is.finite(draws$se)), all(draws$se>0),
  all(draws$seed==seed_base+100000L*(draws$N_trial==500L)+draws$replicate))
critical <- qnorm(.975)
draws[, `:=`(reject_null=as.integer(estimate-critical*se>=0),
  declare_alternative=as.integer(estimate+theta_star-critical*se>=0))]
summary <- draws[, {
  p0 <- mean(reject_null)
  p1 <- mean(declare_alternative)
  list(repetitions=.N, null_count=sum(reject_null), null_rejection=p0,
    null_mcse=sqrt(p0*(1-p0)/.N), alternative=theta_star,
    alternative_count=sum(declare_alternative), efficacy_declaration=p1,
    efficacy_mcse=sqrt(p1*(1-p1)/.N), actual_gap=unique(target_TCE),
    mean_se=mean(se))
}, by=.(case, N_trial, method, xi_index, xi)]
stopifnot(all(summary$repetitions==5000L))
summary[, method_order:=match(method, methods)]
setorder(summary, case, N_trial, method_order, xi_index)
summary[, method_order:=NULL]
fwrite(summary, file.path(output, 'fitted-operating-characteristics.csv'))
selected <- summary[case=='matched_gap' & xi_index %in% c(-4L, 0L, 4L)]
fwrite(selected, file.path(output, 'fitted-operating-characteristics-selected.csv'))
fwrite(summary[case=='transportability'],
  file.path(output, 'fitted-operating-characteristics-transportability.csv'))
fwrite(draws[(case=='matched_gap' & xi_index %in% c(-4L,0L,4L)) | case=='transportability',
  .(case, N_trial, method, replicate, seed, xi, actual_gap=target_TCE, estimate, se,
    reject_null, declare_alternative)],
  file.path(output, 'fitted-operating-characteristics-selected-draws.csv.gz'), compress='gzip')

cell <- function(p, mcse) sprintf('$%.4f\\,(%.4f)$', p, mcse)
table <- c('\\begin{table}[!htbp]', '\\centering',
  paste0('\\caption{Fitted-procedure one-sided operating characteristics from the original 5,000 paired repetitions at each trial size. Entries are probabilities with empirical binomial Monte Carlo standard errors in parentheses. The decision is $\\widehat\\theta-z_{.975}\\widehat s\\ge0$. HCT and ECT have the indicated matched causal-gaps; their external mean intercepts differ. The RCT gap remains zero, so its entries repeat across columns. Zero matched gap need not imply conditional mean transportability. All original nuisance fits and SE procedures are retained. The entry $1.0000$ denotes 5,000 declarations in 5,000 repetitions.}'),
  '\\label{tab:fitted-one-sided}',
  '\\renewcommand{\\baselinestretch}{1}\\fontsize{9.5}{11}\\selectfont',
  '\\setlength{\\tabcolsep}{3pt}', '\\renewcommand{\\arraystretch}{1.12}',
  '\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lrrr@{}}',
  '\\toprule',
  'Analysis & $\\xi=-.1$ & $\\xi=0$ & $\\xi=.1$\\\\',
  '\\midrule')
for (metric in c('null', 'alternative')) {
  heading <- if (metric=='null') 'A. Causal-null rejection ($\\theta=0$)' else
    paste0('B. Efficacy declaration ($\\theta_\\star=', theta_star, '$)')
  table <- c(table, paste0('\\multicolumn{4}{l}{\\textbf{', heading, '}}\\\\'))
  for (n in c(250L,500L)) {
    table <- c(table, paste0('\\multicolumn{4}{l}{Trial enrollment $n_{\\mathrm{RCT}}=', n, '$}\\\\'))
    for (m in methods) {
      z <- selected[N_trial==n & method==m][order(xi)]
      entries <- if (metric=='null') cell(z$null_rejection,z$null_mcse) else
        cell(z$efficacy_declaration,z$efficacy_mcse)
      table <- c(table, paste0(paste(c(method_labels[[m]], entries), collapse=' & '), '\\\\'))
    }
    table <- c(table, '\\addlinespace')
  }
  if (metric=='null') table <- c(table, '\\midrule')
}
table <- c(table, '\\bottomrule', '\\end{tabular*}', '\\end{table}')
writeLines(table, file.path(output, 'fitted-one-sided-table.tex'))

stopifnot(identical(hashes,
  vapply(input_files, function(f) digest::digest(file=f, algo='sha256'), character(1))))
report <- list(success=TRUE, source='original archived fitted-procedure draws',
  repetitions_per_size=5000,
  sample_sizes=designs, alternative=theta_star, critical_value=critical,
  decision='estimate - qnorm(.975)*se >= 0',
  alternative_evaluation='Add a constant theta_star to treated outcomes; exact location equivariance, verified by direct refits before translating archived estimates.',
  refit_seed_replicates=c(1,2500,5000), refit_gaps=c(-.3,-.1,0,.1,.3),
  refit_cases=c('matched_gap','transportability'), refit_rows=nrow(checks),
  maximum_archived_estimate_error=max(checks$archived_estimate_error),
  maximum_archived_se_error=max(checks$archived_se_error),
  maximum_translation_estimate_error=max(checks$translation_estimate_error),
  maximum_translation_se_error=max(checks$translation_se_error),
  maximum_translation_influence_error=max(checks$translation_influence_error,na.rm=TRUE),
  mcse='sqrt(p_hat*(1-p_hat)/5000); values are marginal Monte Carlo standard errors, not independent between methods or gaps.',
  conditional_transportability='Separate exact delta(x)=0 reference, not the heterogeneous-discrepancy grid at xi=0.',
  no_full_simulation_rerun=TRUE, unchanged_source_inputs=TRUE,
  source_sha256=as.list(setNames(unname(hashes), basename(names(hashes)))), R_version=R.version.string)
jsonlite::write_json(report,file.path(output,'fitted-operating-characteristics-validation.json'),
  pretty=TRUE,auto_unbox=TRUE,digits=16)
writeLines(capture.output(sessionInfo()), file.path(output, 'session-info.txt'))
print(selected[, .(N_trial,method,xi,null_rejection,null_mcse,
  efficacy_declaration,efficacy_mcse)])
print(summary[case=='transportability', .(N_trial,method,null_rejection,
  null_mcse,efficacy_declaration,efficacy_mcse)])
