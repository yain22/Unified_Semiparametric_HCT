# Draw five publication figures from complete retained or regenerated simulations.
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))))
source_dir <- file.path(here, "..", "simulation_figures")
source(file.path(source_dir, "plot_simulations.R"))
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat("Rscript plot_publication.R --output-dir DIR [--baseline DIR] [--oracle DIR] [--extended DIR]\n",
      "Refitted --oracle and --extended directories include their enrollment-500 subdirectories.\n")
} else {
  archive <- archived_inputs(source_dir)
  options <- plot_options(args, c(archive[c("baseline", "oracle", "extended")], list("output-dir" = "")))
  out <- options[["output-dir"]]
  if (!nzchar(out)) stop("Specify --output-dir for the five publication PDFs.")
  options$oracle_500 <- if (any(grepl("^--oracle(=|$)", args)))
    file.path(options$oracle, "enrollment-500") else archive$oracle_500
  options$extended_500 <- if (any(grepl("^--extended(=|$)", args)))
    file.path(options$extended, "enrollment-500") else archive$extended_500
  render_simulations(options[names(archive)], out, file.path(out, "validation"))
}
