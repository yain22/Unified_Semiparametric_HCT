# Six-month CGM analysis. Explicit local inputs; aggregate outputs only.
# Source this file for cohort and score helpers. The analysis entry point is support_gap.R.
cgm6m_baseline <- c('hba1c0', 'age', 'duration', 'male', 'pump')
cgm6m_basis <- ~ hba1c0 + age + duration + male + pump

clinical_guard_output <- function(out, package_dir) {
  canonical <- function(path) {
    # Resolve from the deepest existing ancestor so nonexistent output paths
    # and '..' components are checked before creating any directories.
    path <- path.expand(path)
    remaining <- character()
    while (!dir.exists(path)) {
      if (file.exists(path)) stop('Output path contains an existing file.')
      parent <- dirname(path)
      if (identical(parent, path)) stop('Cannot resolve output path.')
      remaining <- c(basename(path), remaining)
      path <- parent
    }
    path <- normalizePath(path, winslash = '/', mustWork = TRUE)
    for (part in remaining) {
      if (part == '.') next
      path <- if (part == '..') dirname(path) else file.path(path, part)
    }
    path <- normalizePath(path, winslash = '/', mustWork = FALSE)
    sub('/+$', '', path)
  }
  resolved <- canonical(out)
  compare <- function(x) if (.Platform$OS.type == 'windows') tolower(x) else x
  target <- compare(resolved)
  root <- compare(canonical(package_dir))
  protected <- vapply(
    file.path(package_dir, c('inputs', 'src', 'cohort')),
    function(x) compare(canonical(x)), character(1)
  )
  if (target == root || any(vapply(protected, function(p)
    target == p || startsWith(target, paste0(p, '/')), logical(1))))
    stop('Choose an aggregate output directory outside the clinical package root and its inputs/src/cohort directories.')
  resolved
}

clinical_load <- function(package_dir) {
  h <- new.env(parent = .GlobalEnv)
  for (f in c('unrestricted_primitives.R', 'restricted.R', 'trial_wald.R', 'folds.R', 'estimated_share.R'))
    sys.source(file.path(package_dir, 'src', f), envir = h)
  # Exact validation contract used by cgm6m_load; nominal e=.5 is a working
  # choice after selection, not a claim that retained counts are balanced.
  h$hct_eif_validate <- function(data, e = NULL) {
    if (!is.data.frame(data) || !all(c('Y', 'S', 'A') %in% names(data)) ||
      !is.numeric(data$Y) || any(!is.finite(data$Y)) || anyNA(data) ||
      !all(data$S %in% 0:1) || !all(data$A %in% 0:1) || any(data$A[data$S == 0] != 0))
      stop('Invalid complete six-month analysis data.')
    group <- factor(ifelse(data$S == 0, 'E', ifelse(data$A == 1, 'T', 'C')), levels = c('T', 'C', 'E'))
    counts <- table(group)
    if (any(counts < 2L)) stop('All groups need at least two records.')
    if (is.null(e)) e <- .5
    if (length(e) != 1L || !is.finite(e) || abs(e - .5) > 1e-12) stop('Original nominal assignment probability must be .5.')
    list(group = group, counts = counts, e = e, q = mean(data$S), N = nrow(data))
  }
  h
}

clinical_status <- function(external, registry_zip) {
  entry <- utils::unzip(registry_zip, list = TRUE)$Name
  entry <- entry[grepl('Registry - Longitudinal/Data Tables/VisitsLong\\.txt$', entry)]
  stopifnot(length(entry) == 1L)
  con <- unz(registry_zip, entry)
  visits <- read.delim(con,
    sep = '|', quote = '', comment.char = '', check.names = FALSE,
    colClasses = 'character', na.strings = c('', 'NA')
  )
  try(close(con), silent = TRUE)
  ix <- match(as.character(visits$PtID), external$source_id)
  delta <- suppressWarnings(as.numeric(visits$ExamDtMnC)) - external$index_month[ix]
  state <- tolower(trimws(visits$CGMUsed))
  valid <- !is.na(ix) & !is.na(delta) & delta > 0 & delta <= external$followup_months[ix] & state %in% c('yes', 'no')
  assessed <- unique(as.character(visits$PtID[which(valid)]))
  initiators <- unique(as.character(visits$PtID[which(valid & state == 'yes')]))
  ifelse(external$source_id %in% initiators, 'Documented subsequent CGM use',
    ifelse(external$source_id %in% assessed, 'Only nonuse recorded', 'Post-index status unassessed')
  )
}

clinical_baseline <- function(groups) {
  do.call(rbind, lapply(names(groups), function(g) do.call(rbind, lapply(cgm6m_baseline, function(v) {
    x <- groups[[g]][[v]]
    data.frame(
      group = g, variable = v, n_total = length(x), n_observed = sum(is.finite(x)), missing = sum(!is.finite(x)),
      mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE),
      count = if (v %in% c('male', 'pump')) sum(x == 1, na.rm = TRUE) else NA,
      percentage = if (v %in% c('male', 'pump')) 100 * mean(x, na.rm = TRUE) else NA
    )
  }))))
}


if (sys.nframe() == 0L) {
  script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  source(file.path(dirname(normalizePath(script)), "support_gap.R"))
}
