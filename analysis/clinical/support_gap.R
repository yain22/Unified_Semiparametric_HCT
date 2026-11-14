# Six-month CGM estimates and matching-score gap with the registry eligibility rule.
# Rscript support_gap.R --clinical=<clinical code directory> --cohort=<authorized cache>
#   --registry-zip=<registry archive> --analysis=primary|subsequent|near-index
#   --bootstrap-reps=2000 --out=<aggregate output directory>
# In place of the cache, supply --trial-zip and --registry-zip to rebuild in memory.
# With --analysis=primary or near-index the registry archive is optional with a cache.
opts <- commandArgs(trailingOnly = TRUE)
opt <- function(key, default = NULL) {
  prefix <- paste0("--", key, "=")
  value <- opts[startsWith(opts, prefix)]
  if (length(value) > 1L) stop("Duplicate option: ", key)
  if (length(value)) substring(value, nchar(prefix) + 1L) else default
}
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
clinical <- opt("clinical", dirname(normalizePath(script)))
keys <- sub("=.*$", "", sub("^--", "", opts))
if (any(!grepl("^--[^=]+=", opts)) ||
    any(!keys %in% c("clinical", "out", "analysis", "bootstrap-reps", "cohort", "registry-zip", "trial-zip"))) {
  stop("Use --name=value options for the clinical entry point; see README.md.")
}
source(file.path(clinical, "run_clinical.R"))
h <- clinical_load(clinical)
analysis <- match.arg(opt("analysis", "primary"), c("primary", "subsequent", "near-index"))
out <- clinical_guard_output(opt("out", file.path(clinical, "../../results/cgm", analysis)), clinical)
b_text <- opt("bootstrap-reps", "0")
if (!grepl("^[0-9]+$", b_text)) stop("--bootstrap-reps must be 0 or an integer of at least 2.")
B <- suppressWarnings(as.integer(b_text))
stopifnot(is.finite(B), B == 0L || B >= 2L)
cohort <- opt("cohort")
registry_zip <- opt("registry-zip")
trial_zip <- opt("trial-zip")
if (!is.null(cohort) && !is.null(trial_zip)) stop("Choose a cohort cache or both source archives.")
if (analysis == "subsequent" && is.null(registry_zip)) stop("Supply --registry-zip for utilization status.")
if (any(!file.exists(c(cohort, registry_zip, trial_zip)))) stop("A supplied data file is missing.")
if (is.null(cohort)) {
  stopifnot(!is.null(trial_zip), !is.null(registry_zip))
  source(file.path(clinical, "cohort/cohort_6m.R"))
  z <- load_cgm_hct_6m_cohorts(file.path(clinical, "cohort"),
    spec = cgm_6m_spec(trial_zip, registry_zip),
    builder_path = file.path(clinical, "cohort/cohort_6m.R"))
} else {
  z <- readRDS(cohort)
}
trial <- z$trial
external <- z$external_primary
stopifnot(nrow(trial) == 442L, nrow(external) == 273L, sum(trial$A) == 229L,
          !anyDuplicated(trial$id), !anyDuplicated(external$id))
shared <- function(d) d$age >= 8 & d$duration >= 1 & d$hba1c0 < 10
stopifnot(all(shared(external)), !any(trial$age < 8 | trial$duration < 1))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_result <- function(x, name) {
  write.csv(x, file.path(out, paste0(name, ".csv")), row.names = FALSE)
}
if (B == 0L && any(file.exists(file.path(out, c("full_results.csv", "bootstrap_draws.csv"))))) {
  stop("Existing bootstrap results in this output directory; choose a new directory for a point-only run.")
}
if (analysis == "primary") {
  if (!is.null(z$trial_flow)) write_result(z$trial_flow, "trial_flow")
  if (!is.null(z$baseline_summary)) write_result(z$baseline_summary, "baseline_summary")
}
fields <- c("S", "A", "Y", cgm6m_baseline)
all_data <- rbind(trial[fields], external[fields])
fold <- h$cgm_folds(all_data, seed = 9302026L)
groups <- list(trial_CGM = trial[trial$A == 1, ], internal = trial[trial$A == 0, ],
               external = external)
source_check <- function(g, label) {
  # Trial source records measure days relative to randomization. Registry dates
  # are months relative to a common origin within participant.
  elapsed <- if (all(g$S == 1L)) {
    (g$followup6_days - g$baseline_sample_days) / (365.25/12)
  } else {
    g$followup_month - g$baseline_month
  }
  stopifnot(max(abs(elapsed - (g$followup_months - g$baseline_months))) < 1e-10)
  qu <- quantile(elapsed, c(0,.25,.5,.75,1), names = FALSE)
  data.frame(group = label, n = nrow(g), age_below_8 = sum(g$age < 8),
             duration_below_1 = sum(g$duration < 1), baseline_ge_10 = sum(g$hba1c0 >= 10),
             eligible = sum(shared(g)), lab_months_min = qu[1], q25 = qu[2],
             median = qu[3], q75 = qu[4], lab_months_max = qu[5],
             above_7 = sum(elapsed > 7), baseline_at_index = sum(g$baseline_months == 0))
}
write_result(do.call(rbind, Map(source_check,groups,names(groups))), "eligibility_and_timing")
ranges <- function(g,label) do.call(rbind,lapply(cgm6m_baseline,function(v) {
  b <- range(external[[v]])
  data.frame(group=label,variable=v,min=min(g[[v]]),max=max(g[[v]]),
             below_external_min=sum(g[[v]]<b[1]),above_external_max=sum(g[[v]]>b[2]))
}))
write_result(do.call(rbind, Map(ranges,groups,names(groups))), "observed_covariate_ranges")
keep <- rep(TRUE,nrow(external))
if (analysis == "subsequent") {
  stopifnot(!is.null(registry_zip))
  if (!is.null(z$registry_indexed)) {
    eligible <- z$registry_indexed
    include <- eligible$id %in% external$id
    stopifnot(nrow(eligible) == 18579L, sum(include) == 273L, sum(eligible$observed_hba1c_followup) == 274L)
    groups <- list(Included = eligible[include, , drop = FALSE], Excluded = eligible[!include, , drop = FALSE])
    write.csv(clinical_baseline(groups), file.path(out, 'registry_selection_baseline.csv'), row.names = FALSE)
    differences <- do.call(rbind, lapply(cgm6m_baseline, function(v) {
      a <- groups$Included[[v]]
      b <- groups$Excluded[[v]]
      data.frame(
        variable = v, difference = mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE),
        standardized_difference = (mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE)) / sqrt((var(a, na.rm = TRUE) + var(b, na.rm = TRUE)) / 2)
      )
    }))
    write.csv(differences, file.path(out, 'registry_selection_differences.csv'), row.names = FALSE)
  }
  status <- clinical_status(external, registry_zip)
  write_result(as.data.frame(table(status)), "subsequent_cgm_status")
  keep <- status != "Documented subsequent CGM use"
  stopifnot(sum(!keep) == 17L)
} else if (analysis == "near-index") {
  keep <- external$baseline_month == external$index_month
}
keep <- c(rep(TRUE,nrow(trial)), keep)
d <- all_data[keep, ]
fold <- fold[keep]
group <- h$cgm_group(d)
# Keep the original five bootstrap strata within each observed T/C/E group;
# the historical postbaseline sensitivity ordered its strata alphabetically.
strata <- split(seq_len(nrow(d)), interaction(
  if (analysis=="subsequent") as.character(group) else group, fold, drop=TRUE))
write_result(as.data.frame(table(group,fold)), "partition")
normal_interval <- function(estimate, se) estimate + c(-1,1)*qnorm(.975)*se
working_se <- function(psi, g) h$hct_unrestricted_group_variance(psi, g)$se
fit <- function(d, full = FALSE) {
  u <- shared(d)
  S <- d$S; T <- S*d$A; C <- S*(1-d$A); E <- 1-S
  R <- as.logical(S); take <- function(x) as.logical(x)
  N <- nrow(d); nR <- sum(R); group <- h$cgm_group(d)
  stopifnot(all(u[!R]), all(d$A[!R] == 0))
  design <- model.matrix(cgm6m_basis,d)
  means <- lapply(list(T, C, E), function(g) {
    rows <- take(g); h$hct_unrestricted_ols(design[rows,,drop=FALSE],d$Y[rows])
  })
  pred <- function(beta) drop(design %*% beta)
  mT <- pred(means[[1]]$beta); mC <- pred(means[[2]]$beta); mE <- numeric(N)
  mE[u] <- drop(design[u,,drop=FALSE]%*%means[[3]]$beta)
  variances <- lapply(2:3, function(j) {
    rows <- take(list(T,C,E)[[j]])
    h$cgm_est_variance(design[rows,,drop=FALSE],d$Y[rows],means[[j]],h)$beta
  })
  vC <- vE <- numeric(N)
  vC[u] <- exp(drop(design[u,,drop=FALSE]%*%variances[[1]]))
  vE[u] <- exp(drop(design[u,,drop=FALSE]%*%variances[[2]]))
  source_fit <- h$hct_unrestricted_logit(design[u,,drop=FALSE],S[u])
  p <- lambda <- rep(1,N);odds <- numeric(N)
  eta <- drop(design[u,,drop=FALSE]%*%source_fit$eta)
  p[u] <- plogis(eta);odds[u] <- exp(eta)
  lambda[u] <- (.5*odds[u])/(.5*odds[u]+vC[u]/vE[u])
  stopifnot(all(is.finite(odds[u])),all(p[u] > 0 & p[u] < 1),
            all(lambda[u] > 0 & lambda[u] < 1))
  trial_score <- S*(mT-mC)+2*T*(d$Y-mT)-2*C*(d$Y-mC)
  # The uncentered HCT and internal-only numerators coincide outside shared
  # eligibility. Form the gap score only where an external control can occur.
  a <- which(u)
  diff_mean <- mC[a]-mE[a]
  l <- lambda[a]; pp <- p[a]; od <- odds[a]
  rc <- d$Y[a]-mC[a]; re <- d$Y[a]-mE[a]
  gap <- numeric(N)
  gap[a] <- S[a]*(1-l)*diff_mean + 2*C[a]*(1-l)*rc - E[a]*(1-l)*od*re +
    diff_mean*l*(1-l)*(-(S[a]-pp)/(1-pp) + 2*C[a]*(rc^2/vC[a]-1) -
                       E[a]*od*(re^2/vE[a]-1))
  # Independent verification against the original unrestricted score on shared
  # profiles, reweighted by the overall trial denominator.
  score_shared <- h$cgm_est_score(d[a, ],
    data.frame(mT=mT[a],mC=mC[a],mE=mE[a],p=pp,odds=od,vC=vC[a],vE=vE[a]))
  hct <- (sum(S[a])*score_shared$estimate + sum(trial_score[!u]))/nR
  xi <- sum(gap)/nR; augmented <- sum(trial_score)/nR
  stopifnot(abs(hct-xi-augmented)<1e-10)
  # Common-control comparator on shared support, internal-control regression
  # outside it. There are no external observations beyond the eligibility cap.
  pooled <- !take(d$A) & u
  mP <- mC
  mP[u] <- drop(design[u,,drop=FALSE] %*%
    h$hct_unrestricted_ols(design[pooled,,drop=FALSE],d$Y[pooled])$beta)
  pr <- p
  pr[u] <- pmin(1-1e-6,pmax(1e-6,p[u]))
  rs <- S*(mT-mP)+2*T*(d$Y-mT)-(1-d$A)*pr/(1-.5*pr)*(d$Y-mP)
  restricted <- sum(rs)/nR
  # Under the ECT observation law restrict the treated standardization
  # population before fitting the source model and evaluating the contrast.
  eobs <- u & take(T+E); nde <- sum(eobs); nte <- sum(T[eobs])
  se <- S[eobs]; ze <- design[eobs,,drop=FALSE]
  logit_ect <- h$hct_unrestricted_logit(ze,se)
  ep <- plogis(drop(ze%*%logit_ect$eta)); eg <- pmin(.99,pmax(.01,ep))
  eo <- eg/(1-eg); observed_resid <- d$Y[eobs]-mE[eobs]
  eraw <- (se-(1-se)*eo)*observed_resid
  ect <- sum(eraw)/nte
  it <- h$hct_trial_wald_comparators(d,cgm6m_basis)
  est <- c("ANOVA (RCT)"=it[["ANOVA (RCT)"]]$estimate,"ANCOVA (RCT)"=it[["ANCOVA (RCT)"]]$estimate,
           "AIPW-R"=restricted,"AIPW-U (est)"=hct,"AIPW (ECT)"=ect,
           matching_gap=xi,trial_augmented=augmented,
           control_difference=mean(diff_mean[S[a]==1]),plugin_gap=sum(S[a]*(1-l)*diff_mean)/nR)
  diagnostics <- c(hct_ps_min=min(p[u]),hct_ps_max=max(p[u]),ect_ps_min=min(ep),
    ect_ps_max=max(ep),ect_lower_hits=sum(ep<.01),ect_upper_hits=sum(ep>.99),
    ect_external_hits=sum(!se & ep!=eg),lambda_min=min(lambda[R]),lambda_max=max(lambda[R]),
    lambda_mean=mean(lambda[R]))
  if (!full) return(c(est, diagnostics))
  se_method <- c(it[["ANOVA (RCT)"]]$se,it[["ANCOVA (RCT)"]]$se,
    working_se((rs-S*restricted)/mean(S),group),
    working_se((trial_score+gap-S*hct)/mean(S),group),
    sqrt(mean(((eraw-se*ect)/mean(se))^2)/nde))
  probabilities <- list(trial_treated=p[take(T)&u],internal=p[take(C)&u],external=p[take(E)],
                        ect_treated=ep[take(se)],ect_external=ep[!take(se)])
  weights <- list(hct_odds=odds[take(E)],hct_external=((1-lambda)*odds)[take(E)],ect_odds=eo[!take(se)])
  list(est=est, diagnostics=diagnostics, se=se_method, probs=probabilities, weights=weights)
}
point <- fit(d,full=TRUE)
summary_point <- data.frame(method=names(point$est),estimate=unname(point$est))
summary_point$analytic_se <- c(point$se,rep(NA_real_,nrow(summary_point)-length(point$se)))
summary_point$analytic_lower <- summary_point$estimate-qnorm(.975)*summary_point$analytic_se
summary_point$analytic_upper <- summary_point$estimate+qnorm(.975)*summary_point$analytic_se
write_result(summary_point,"point_estimates")
write_result(data.frame(quantity=names(point$diagnostics),value=unname(point$diagnostics)), "score_diagnostics")
summaries <- function(lst, concentration=FALSE) do.call(rbind,lapply(names(lst), function(label) {
  x <- lst[[label]];q <- unname(quantile(x,c(0,.05,.5,.95,1)))
  s <- data.frame(quantity=label,n=length(x),min=q[1],p05=q[2],median=q[3],p95=q[4],max=q[5])
  if(concentration) {
    w <- x/sum(x);s$effective <- 1/sum(w*w);s$max_share <- max(w)
    s$top_tenth_share <- sum(head(sort(w,decreasing=TRUE),ceiling(length(w)/10)))
  }
  s
}))
write_result(summaries(point$probs), "source_probabilities")
write_result(summaries(point$weights, TRUE), "external_weights")
if (B) {
  # Rows remain aligned with the primary archive. No identifiers, profiles,
  # fitted participant predictions, or record-level timing leave memory.
  draws <- lapply(seq_len(B),function(b) {
    set.seed(19302026L + 37L*b)
    ids <- unlist(lapply(strata,function(ii)sample(ii,length(ii),replace=TRUE)),use.names=FALSE)
    ids <- c(ids[d$S[ids]==1L],ids[d$S[ids]==0L])
    value <- fit(d[ids,,drop=FALSE])
    if(b %% 250L==0L)message("Completed ",b," of ",B," paired refits.")
    value
  })
  draw <- as.data.frame(do.call(rbind,draws),check.names=FALSE)
  stopifnot(nrow(draw)==B,all(is.finite(as.matrix(draw))),
    max(abs(draw[["AIPW-U (est)"]]-draw$matching_gap-draw$trial_augmented))<1e-10)
  write_result(data.frame(replicate=seq_len(B), draw,check.names=FALSE), "bootstrap_draws")
  summary_point$bootstrap_se <- vapply(draw[names(point$est)],sd,numeric(1))
  summary_point$bootstrap_lower <- summary_point$estimate-qnorm(.975)*summary_point$bootstrap_se
  summary_point$bootstrap_upper <- summary_point$estimate+qnorm(.975)*summary_point$bootstrap_se
  write_result(summary_point, "full_results")
}

