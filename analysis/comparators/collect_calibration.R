# Summarize each borrowing procedure against its own statistical target.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
here <- dirname(normalizePath(script))
args <- commandArgs(TRUE)
usage <- "Rscript collect_calibration.R [--output-dir DIR] [--simulation-dir DIR]"
if (identical(args, "--help")) { cat(usage, "\n"); quit(status = 0L) }
if (length(args) %% 2L) stop(usage, call. = FALSE)
options <- if (length(args)) args[seq(1L, length(args), 2L)] else character()
if (anyDuplicated(options) || any(!options %in% c("--output-dir", "--simulation-dir")))
  stop(usage, call. = FALSE)
value <- function(key, default = NULL) {
  index <- match(key, options)
  if (is.na(index)) default else args[2L * index]
}
out <- value("--output-dir", file.path(dirname(dirname(here)), "results", "calibration"))
input_root <- value("--simulation-dir")
if (is.null(input_root)) {
  files <- c(baseline = "inputs/baseline-summary.csv", transportability = "inputs/transportability-summary.csv",
             extended = "inputs/extended-summary.csv", oracle = "inputs/oracle-summary.csv")
  input_root <- here
} else {
  files <- c(baseline = "baseline/summary.csv", transportability = "baseline/transportability_reference_summary.csv",
             extended = "extended/summary.csv", oracle = "oracle/oracle_summary.csv")
}
if (dir.exists(out) && identical(normalizePath(out), here))
  stop("Choose an output directory separate from the reference tables.", call. = FALSE)

summarize_model <- function(setting) {
  x <- read.csv(file.path(input_root, files[[setting]]), fileEncoding = "UTF-8-BOM",
                stringsAsFactors = FALSE)
  x <- x[x$method %in% c("U", "AIPW (ECT)"), ]
  if (!nrow(x)) stop("No borrowing estimates in the ", setting, " summary.", call. = FALSE)
  mean_se <- x$width / (2 * qnorm(.975))
  if ("se" %in% names(x)) stopifnot(all(abs(mean_se - x$se) < 1e-12))
  mode <- if (setting == "oracle") "Oracle nuisance" else "Estimated nuisance"
  data.frame(setting = if ("scenario" %in% names(x)) x$scenario else setting,
             mode = if ("mode" %in% names(x)) x$mode else mode, N_trial = x$N_trial, method = x$method,
             delta_identically_zero = as.integer(setting == "transportability"),
             gap = if ("xi" %in% names(x)) x$xi else 0, repetitions = x$n_reps,
             bias_TCE = x$bias_TCE, bias_mcse = x$bias_TCE_mcse,
             empirical_sd = sqrt(x$empirical_variance), mean_se = mean_se,
             coverage_TCE = x$coverage_TCE, coverage_mcse = x$coverage_TCE_mcse)
}
calibration <- do.call(rbind, lapply(names(files), summarize_model))
reps <- unique(calibration$repetitions)
stopifnot(length(reps) == 1L, all(is.finite(reps)), reps > 0)

# The paper displays pointwise agreement separately from a zero average gap.
baseline <- expand.grid(gap = c(-.3, 0, .3), n = c(250, 500))
extra <- c("fewer_external", "source_shift", "stronger_heteroskedasticity", "quadratic_discrepancy")
selected <- rbind(
  data.frame(setting = "transportability", method = rep(c("U", "AIPW (ECT)"), each = 2),
             n = rep(c(250, 500), 2), gap = 0),
  data.frame(setting = "baseline", method = "U", n = baseline$n, gap = baseline$gap),
  data.frame(setting = "baseline", method = "AIPW (ECT)", n = c(250, 500), gap = 0),
  data.frame(setting = extra, method = "U", n = 250, gap = 0))
setting_labels <- c(transportability = "Pointwise $\\delta\\equiv0$", baseline = "Heterogeneous baseline",
                    fewer_external = "Fewer external controls", source_shift = "Shifted source covariates",
                    stronger_heteroskedasticity = "Stronger heteroskedasticity", quadratic_discrepancy = "Quadratic discrepancy")
caption <- paste0("\\caption{Calibration for the statistical target, based on ",
                  format(reps, big.mark = ",", scientific = FALSE, trim = TRUE),
                  " repetitions per row. Bias is measured against the statistical target of each procedure, ",
                  "SD is the empirical standard deviation, and SE is the mean reported standard error. ",
                  "The pointwise transportability check ($\\delta\\equiv0$) is separate from heterogeneous ",
                  "discrepancies with zero aggregate causal-gap. Additional scenarios are shown at zero ",
                  "aggregate gap.}\\label{tab:own-tce-calibration}")
lines <- c("\\begin{table}[!htbp]", "\\centering", caption, "\\small",
           "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}>{\\raggedright\\arraybackslash}p{3.5cm}lrrrrr@{}}",
           "\\toprule", "Setting & Procedure & $n_{\\mathrm{RCT}}$ & $\\xi_{\\lambda,\\delta}$ & Bias & SD / SE & Coverage\\\\",
           "\\midrule")
for (i in seq_len(nrow(selected))) {
  s <- selected[i, ]
  row <- calibration[calibration$setting == s$setting & calibration$method == s$method &
                       calibration$N_trial == s$n & abs(calibration$gap - s$gap) < 1e-12 &
                       calibration$mode == "Estimated nuisance", ]
  if (nrow(row) != 1L) stop("Cannot uniquely match table row ", i, ".", call. = FALSE)
  first <- i == 1L || s$setting != selected$setting[i - 1L]
  if (i > 1L && first) lines <- c(lines, "\\addlinespace[3pt]")
  label <- if (first) setting_labels[[s$setting]] else ""
  method <- if (s$method == "U") "AIPW-U" else "ECT AIPW"
  lines <- c(lines, sprintf("%s & %s & %d & $%.2f$ & $%.4f$ & $%.4f$ / $%.4f$ & $%.4f$ \\\\",
                            label, method, s$n, s$gap, row$bias_TCE, row$empirical_sd,
                            row$mean_se, row$coverage_TCE))
}
sources <- lapply(files, function(p) list(path = p, sha256 = digest::digest(file = file.path(input_root, p),
                                                                           algo = "sha256", serialize = FALSE)))
manifest <- list(sources = sources, rows = nrow(calibration),
                 mean_se = "Mean width/(2 z_.975), checked against saved mean SE where present",
                 distinction = paste("delta-identically-zero uses dedicated transportability_reference_summary;",
                                     "it is never inferred from zero aggregate gap."))
dir.create(out, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(out)) stop("Could not create output directory.", call. = FALSE)
write.csv(calibration, file.path(out, "own_tce_calibration.csv"), row.names = FALSE)
writeLines(c(lines, "\\bottomrule", "\\end{tabular*}", "\\end{table}"),
           file.path(out, "own-tce-calibration-table.tex"))
jsonlite::write_json(manifest, file.path(out, "calibration-provenance.json"), pretty = TRUE, auto_unbox = TRUE)
cat("Wrote", nrow(calibration), "calibration rows.\n")
