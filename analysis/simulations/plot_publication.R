# Draw six publication figures from complete retained or regenerated simulations.
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))))
source_dir <- file.path(here, "..", "simulation_figures")
source(file.path(source_dir, "plot_simulations.R"))
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat("Rscript plot_publication.R --output-dir DIR [--baseline DIR] [--oracle DIR] [--extended DIR]\n")
} else {
  options <- plot_options(args, c(archived_inputs(source_dir), list("output-dir" = "")))
  out <- options[["output-dir"]]
  if (!nzchar(out)) stop("Specify --output-dir for the six publication PDFs.")
  render_simulations(options[c("baseline", "oracle", "extended")], out, file.path(out, "validation"))
}
