# Open HCT.Rproj, source this file, then run reproduce("tables") or reproduce("figures").
reproduce <- function(what = c("tables", "figures", "validate", "smoke", "refit"), workers = 1L) {
  what <- match.arg(what)
  if (!file.exists("HCT.Rproj")) stop("Open HCT.Rproj or set the working directory to its folder.")
  r <- file.path(R.home("bin"), "Rscript")
  run <- function(file, ...) {
    args <- c("--vanilla", file.path("analysis", file), ...)
    status <- system2(r, vapply(args, shQuote, character(1), USE.NAMES=FALSE))
    if (status != 0L) stop("Stopped at ", file, " (exit ", status, ").")
  }
  task <- function(...) run("paper_analysis.R", ...)
  if (what == "tables") return(invisible(run("paper_tables.R")))
  if (what == "validate") return(invisible(task("validate")))
  if (what == "figures") {
    out <- file.path("results", "figures")
    for (script in c("illustrate_hct_standardization_pictorial.R", "illustrate_hct_eif.R"))
      run(file.path("illustrations", script), "--output-dir", file.path(out, "illustrations"))
    run("planning/plot_planning.R", "--output-dir", file.path(out,"planning"))
    run("simulation_figures/plot_simulations.R", "--output-dir", file.path(out,"simulations"))
    run("clinical/render_clinical.R", "--out", file.path(out,"clinical"))
    run("clinical/render_descriptive.R", "--out", file.path(out, "clinical"))
    map <- read.csv(file.path("analysis", "paper-map.csv"), stringsAsFactors = FALSE)
    index <- map[startsWith(map$exhibit, "Figure "), c("exhibit", "figure_output", "script_or_source")]
    input <- file.path(out, index$figure_output)
    artwork <- endsWith(index$script_or_source, ".pdf")
    if (!all(file.copy(index$script_or_source[artwork], input[artwork], overwrite = TRUE)))
      stop("Could not copy the manuscript artwork listed in analysis/paper-map.csv.")
    if (anyDuplicated(index$exhibit) || anyDuplicated(index$figure_output) || any(!file.exists(input)))
      stop("Figure map has duplicate or missing outputs. Check analysis/paper-map.csv.")
    number <- sub("^Figure ", "", index$exhibit)
    prefix <- ifelse(startsWith(number, "S"), "S", "")
    index$pdf <- sprintf("Figure-%s%02d.pdf", prefix, as.integer(sub("^S", "", number)))
    if (!all(file.copy(input, file.path(out, index$pdf), overwrite = TRUE)))
      stop("Could not write numbered figure PDFs.")
    write.csv(index[c("exhibit", "pdf", "figure_output", "script_or_source")],
              file.path(out, "figure-index.csv"), row.names = FALSE)
    # Discard the two superseded generated PDFs when refreshing an earlier run.
    retired <- file.path(out, c("Figure-S12.pdf",
                                "simulations/figures/oracle_estimated_nuisance_comparison.pdf"))
    if (!all(file.remove(retired[file.exists(retired)])))
      stop("Could not remove the superseded figure outputs.")
    message("Wrote ", nrow(index), " numbered PDFs and figure-index.csv in ", out, ".")
    return(invisible(out))
  }
  if (length(workers)!=1L || !is.numeric(workers) || !is.finite(workers) ||
      workers != floor(workers) || workers < 1 || workers > 8) stop("workers must be 1 through 8.")
  if (what == "smoke") {
    reps <- "4"; out <- file.path("results","smoke")
  } else {
    reps <- "5000"; out <- "results"
  }
  # The oracle comparison must use the same samples as the fitted baseline.
  dir <- function(x) file.path(out,x)
  task("baseline", "--reps", reps, "--workers", as.character(workers), "--output-dir", dir("baseline"))
  task("oracle", "--reps", reps, "--baseline-dir", dir("baseline"), "--output-dir", dir("oracle"))
  task("extended", "--reps", reps, "--workers", as.character(workers), "--output-dir", dir("extended"))
  task("quadratic", "--reps", reps, "--baseline-dir", dir("baseline"), "--output-dir", dir("quadratic"))
  if (what == "refit") {
    run("comparators/collect_calibration.R", "--simulation-dir", out)
    run("simulations/fitted_operating_characteristics.R",
        paste0("--quadratic=", dir("quadratic/quadratic_rct_replicates.csv.gz")))
    run("paper_tables.R", "--simulation-dir", out, "--output-dir", dir("refitted-tables"))
    run("simulations/plot_publication.R", "--baseline", dir("baseline"), "--oracle", dir("oracle"),
        "--extended", dir("extended"), "--output-dir", dir("refitted-figures"))
  }
  invisible(out)
}
