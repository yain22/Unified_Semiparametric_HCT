# Analyses for prospective causal-gap sensitivity across trial designs.
# Run from any directory: Rscript /path/to/analysis/paper_analysis.R --help

paper_help <- function() {
  cat(paste(c(
    "Causal-gap sensitivity analysis",
    "Usage: Rscript analysis/paper_analysis.R TASK [--option VALUE]",
    "",
    "Tasks:",
    "  validate      Short independent estimator and synthetic-data checks",
    "  planning      Exact Gaussian planning calculations (base R)",
    "  baseline      Fitted baseline simulation at both sample sizes",
    "  oracle        Paired oracle comparison at both baseline enrollment levels",
    "  extended      Four additional simulation scenarios",
    "  quadratic     Paired quadratic RCT ANCOVA comparator",
    "  cgm           CGM analysis with explicitly supplied local data",
    "  illustrations Conceptual manuscript figures (base R/Cairo)",
    "",
    "Options:",
    "  --output-dir PATH       Separate output directory; default results/TASK",
    "  --reps INTEGER          Simulation repetitions (full paper: 5000)",
    "  --workers INTEGER       Baseline/extended workers, 1-8 (default 1)",
    "  --seed-base INTEGER     Optional simulation seed override; default published seed",
    "  --baseline-dir PATH     Completed baseline for oracle or exact quadratic pairing",
    "  --cohort PATH           Existing authorized CGM cohort_data.rds",
    "  --trial-zip PATH        Original trial archive, alternative to --cohort",
    "  --registry-zip PATH     Original registry archive; needed to rebuild or select",
    "  --bootstrap-reps N      CGM full refits (default 0; paper uses 2000)",
    "  --analysis NAME         CGM primary, near-index or subsequent (default primary)",
    "  --selection true        Run primary and subsequent analyses together",
    "",
    "Examples:",
    "  Rscript analysis/paper_analysis.R validate",
    "  Rscript analysis/paper_analysis.R baseline --reps 4 --workers 1 --output-dir results/smoke/baseline",
    "  Rscript analysis/paper_analysis.R oracle --reps 4 --baseline-dir results/smoke/baseline --output-dir results/smoke/oracle",
    "  Rscript analysis/paper_analysis.R extended --reps 5000 --workers 4",
    "  Rscript analysis/paper_analysis.R cgm --cohort data/cohort_data.rds --bootstrap-reps 2000",
    "",
    "See README.md for table and figure commands.",
    "Short runs verify code operation; they do not reproduce published Monte Carlo precision.",
    "No data are downloaded. Patient-level files are never part of this code package."
  ), collapse = "\n"), "\n")
}

paper_main <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) != 1L) stop("Run this entry point with Rscript.")
  here <- dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE))
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) || args[[1]] %in% c("--help", "-h", "help")) {
    paper_help()
    return(invisible(NULL))
  }
  task <- args[[1]]
  args <- args[-1]
  opts <- list()
  while (length(args)) {
    if (!startsWith(args[[1]], "--")) stop("Expected an option; see --help.")
    if (grepl("=", args[[1]], fixed = TRUE)) {
      key <- sub("=.*$", "", substring(args[[1]], 3))
      val <- sub("^[^=]*=", "", args[[1]])
      args <- args[-1]
    } else {
      if (length(args) < 2L || startsWith(args[[2]], "--")) stop("Missing value for ", args[[1]])
      key <- substring(args[[1]], 3)
      val <- args[[2]]
      args <- args[-c(1, 2)]
    }
    if (key %in% names(opts)) stop("Duplicate option: ", key)
    opts[[key]] <- val
  }
  allowed <- switch(task,
    validate = character(),
    planning = character(),
    illustrations = character(),
    baseline = c("reps", "workers", "seed-base"),
    extended = c("reps", "workers", "seed-base"),
    oracle = c("reps", "seed-base", "baseline-dir"),
    quadratic = c("reps", "baseline-dir"),
    cgm = c("cohort", "trial-zip", "registry-zip", "bootstrap-reps", "selection", "analysis"),
    stop("Unknown task: ", task, "; see --help.")
  )
  invalid <- setdiff(names(opts), c("output-dir", allowed))
  if (length(invalid)) stop("Invalid option for ", task, ": ", paste(invalid, collapse = ", "))
  get <- function(key, default = NULL) if (is.null(opts[[key]])) default else opts[[key]]
  integer_option <- function(key, default, minimum = 1L, maximum = .Machine$integer.max) {
    value <- get(key, as.character(default))
    if (!grepl("^[0-9]+$", value)) stop(key, " must be an integer.")
    value <- suppressWarnings(as.integer(value))
    if (is.na(value) || value < minimum || value > maximum) stop(key, " is outside the allowed range.")
    value
  }
  output <- get("output-dir", file.path(here, "../results", task))
  resolve_path <- function(path) {
    path <- path.expand(path)
    if (file.exists(path) || dir.exists(path)) return(normalizePath(path, winslash = "/", mustWork = TRUE))
    parent <- dirname(path)
    leaf <- basename(path)
    if (identical(parent, path)) stop("Cannot resolve output path: ", path)
    resolved <- resolve_path(parent)
    if (leaf == "..") return(dirname(resolved))
    if (leaf == ".") return(resolved)
    file.path(resolved, leaf)
  }
  output <- resolve_path(output)
  # Accidental publication-directory writes are rejected by the entry point.
  protected <- c(
    file.path(here, "../assets"), file.path(here, "../sections"),
    file.path(here, "../supplement-fragments")
  )
  canonical <- function(x) {
    value <- resolve_path(x)
    if (.Platform$OS.type == "windows") tolower(value) else value
  }
  inside <- function(x, parent) canonical(x) == canonical(parent) || startsWith(canonical(x), paste0(canonical(parent), "/"))
  if (inside(output, here))
    stop("Choose an output directory outside analysis; the default is results/TASK.")
  if (canonical(output) == canonical(file.path(here, "..")))
    stop("Choose a separate results directory, not the manuscript root.")
  for (p in protected) if (canonical(output) == canonical(p) || startsWith(canonical(output), paste0(canonical(p), "/")))
    stop("Choose an output directory outside published assets/reference results.")
  if (task == "oracle" || (task == "quadratic" && !is.null(get("baseline-dir")))) {
    baseline <- get("baseline-dir")
    if (is.null(baseline) || !dir.exists(baseline)) stop("Supply --baseline-dir for a completed baseline run.")
    if (inside(output, baseline) || inside(baseline, output))
      stop("Output and baseline input directories must be separate and non-overlapping.")
  }
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  output <- normalizePath(output, winslash = "/", mustWork = TRUE)
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  if (!file.exists(rscript)) stop("Cannot locate Rscript in this R installation.")
  run <- function(relative, values = character()) {
    script <- file.path(here, relative)
    if (!file.exists(script)) stop("Missing analysis dependency: ", relative)
    cat("\nRunning ", relative, "\n", sep = "")
    flush.console()
    status <- system2(rscript, args = vapply(c(script, values), shQuote, character(1)))
    if (!identical(as.integer(status), 0L)) stop(relative, " failed (exit ", status, ").")
  }
  if (task == "validate") {
    run("simulations/validate.R", paste0("--output-dir=", file.path(output, "simulations")))
    run("clinical/validate_synthetic.R", paste0("--out=", file.path(output, "clinical")))
    run("planning/run_planning.R", c("--output-dir", file.path(output, "planning")))
    run("comparators/quadratic_rct_comparator.R", c("--repetitions", "4", "--output-dir", file.path(output, "quadratic")))
  } else if (task == "planning") {
    run("planning/run_planning.R", c("--output-dir", output))
  } else if (task %in% c("baseline", "extended")) {
    reps <- integer_option("reps", 5000L, 2L)
    workers <- integer_option("workers", 1L, 1L, 8L)
    seed <- integer_option("seed-base", if (task == "baseline") 52461004L else 74461005L)
    run(
      paste0("simulations/run_", task, ".R"),
      c(paste0("--reps=", reps), paste0("--workers=", workers), paste0("--seed-base=", seed), paste0("--output-dir=", output))
    )
  } else if (task == "oracle") {
    baseline <- get("baseline-dir")
    shared <- c(
      paste0("--baseline=", normalizePath(baseline, winslash = "/")),
      paste0("--reps=", integer_option("reps", 5000L, 2L)),
      paste0("--seed-base=", integer_option("seed-base", 52461004L))
    )
    run("simulations/run_oracle.R", c(shared, "--enrollment=250", paste0("--output-dir=", output)))
    run("simulations/run_oracle.R", c(shared, "--enrollment=500",
                                      paste0("--output-dir=", file.path(output, "enrollment-500"))))
  } else if (task == "quadratic") {
    baseline <- get("baseline-dir")
    run("comparators/quadratic_rct_comparator.R", c("--repetitions", as.character(integer_option("reps", 5000L, 1L, 5000L)),
        "--output-dir", output, if (!is.null(baseline)) c("--baseline-dir", baseline)))
  } else if (task == "cgm") {
    cohort <- get("cohort")
    trial <- get("trial-zip")
    registry <- get("registry-zip")
    if (is.null(cohort) && (is.null(trial) || is.null(registry))) stop("Supply --cohort or both --trial-zip and --registry-zip.")
    if (!is.null(cohort) && !is.null(trial)) stop("Use --cohort or source archives, not both.")
    selection <- get("selection", "false")
    if (!selection %in% c("true", "false")) stop("--selection must be true or false.")
    if (selection == "true" && is.null(registry)) stop("Selection diagnostics require --registry-zip.")
    scenario <- get("analysis", "primary")
    if (!scenario %in% c("primary", "subsequent", "near-index")) stop("Unknown CGM analysis.")
    if (selection == "true" && !is.null(get("analysis"))) stop("Use --analysis or --selection.")
    if (scenario == "subsequent" && is.null(registry)) stop("--registry-zip is required for subsequent use.")
    extra <- paste0("--bootstrap-reps=", integer_option("bootstrap-reps", 0L, 0L))
    for (key in c("cohort", "trial-zip", "registry-zip")) if (!is.null(get(key))) {
      if (!file.exists(get(key))) stop("Local input not found: ", key)
      extra <- c(extra, paste0("--", key, "=", normalizePath(get(key), winslash = "/")))
    }
    scenarios <- if (selection == "true") c("primary", "subsequent") else scenario
    for (name in scenarios) {
      dest <- if (length(scenarios) == 1L && !is.null(get("output-dir"))) output else file.path(output, name)
      run("clinical/support_gap.R", c(extra, paste0("--analysis=", name), paste0("--out=", dest)))
    }
  } else if (task == "illustrations") {
    for (name in c("eif", "target_pictorial", "anchor_pictorial", "standardization_pictorial"))
      run(paste0("illustrations/illustrate_hct_", name, ".R"), c("--output-dir", output))
  }
  writeLines(
    c(
      paste("Task:", task), paste("R:", R.version.string), 
      "Status: completed", "Published manuscript assets and reference results were not replaced."
    ),
    file.path(output, "run-completed.txt")
  )
  cat("\nCompleted ", task, ". Results: ", output, "\n", sep = "")
  invisible(output)
}

if (sys.nframe() == 0L) paper_main()
