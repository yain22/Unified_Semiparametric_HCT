# Path checks run before any directory creation or output write.
simulation_resolve_path <- function(path) {
  stopifnot(length(path) == 1L, !is.na(path), nzchar(path))
  p <- path.expand(path)
  tail <- character()
  while (!file.exists(p) && !dir.exists(p)) {
    parent <- dirname(p)
    if (identical(parent, p)) stop('Cannot resolve output path: ', path)
    tail <- c(basename(p), tail)
    p <- parent
  }
  p <- normalizePath(p, winslash = '/', mustWork = TRUE)
  for (part in tail) {
    if (part == '..') p <- dirname(p) else if (part != '.') p <- file.path(p, part)
  }
  p <- gsub('\\\\', '/', p)
  sub('/+$', '', p)
}
simulation_path_within <- function(path, root) {
  if (.Platform$OS.type == 'windows') {
    path <- tolower(path)
    root <- tolower(root)
  }
  identical(path, root) || startsWith(path, paste0(root, '/'))
}
simulation_guard_output <- function(output, package, inputs = character()) {
  output <- simulation_resolve_path(output)
  package <- simulation_resolve_path(package)
  allowed <- simulation_resolve_path(file.path(package, 'output'))
  if ((simulation_path_within(output, package) || simulation_path_within(package, output)) &&
    !simulation_path_within(output, allowed))
    stop('Output must not overlap the simulation source package; its output/ subtree is allowed.')
  if (grepl('(^|/)(lib|inputs|reference|references)(/|$)', tolower(output)))
    stop('Output must not be inside a lib, inputs, reference, or references directory.')
  if (file.exists(output) && !dir.exists(output)) stop('Output path is an existing file.')
  for (input in inputs) {
    input <- simulation_resolve_path(input)
    if (simulation_path_within(output, input) || simulation_path_within(input, output))
      stop('Output must be disjoint from the baseline input directory (including ancestors and descendants).')
  }
  invisible(output)
}
