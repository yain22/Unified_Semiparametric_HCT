# Exact function extracts; see source_manifest.json. Sourcing performs no I/O.

cgm_group <- function(data) factor(ifelse(data$R == 0, 'E', ifelse(data$A == 1, 'T', 'C')), levels = c('T', 'C', 'E'))

cgm_folds <- function(data, seed = 9302026L, V = 5L) {
  set.seed(seed)
  group <- cgm_group(data)
  out <- integer(nrow(data))
  for (g in levels(group)) {
    ids <- which(group == g)
    out[ids] <- sample(rep(seq_len(V), length.out = length(ids)))
  }
  out
}
