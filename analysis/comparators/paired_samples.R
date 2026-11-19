# Paired baseline samples for the quadratic RCT comparison and gap analysis.
# A record's position and the order of its columns enter both fingerprints.
paired_sample <- function(n, b) {
  nT <- if (n == 250L) 167L else 333L
  nC <- n - nT
  nE <- 2L * n
  seed <- 52461004L + 100000L * (n == 500L) + b
  set.seed(seed)
  xe <- matrix(rnorm(nE * 10L), ncol = 10L)
  epsE <- rnorm(nE)
  xt <- matrix(rnorm(nT * 10L), ncol = 10L)
  xc <- matrix(rnorm(nC * 10L), ncol = 10L)
  epsT <- rnorm(nT)
  epsC <- rnorm(nC)
  X <- rbind(xt, xc, xe)
  f <- .5 + (.6 / sqrt(5)) * rowSums(X[, 1:5]) + (.25 / sqrt(5)) * rowSums(X[, 1:5]^2)
  error <- c(epsT, exp(.2 * xc[, 1]) * epsC, exp(-.2 * xe[, 1]) * epsE)
  data <- setNames(as.data.frame(X), paste0('X', 1:10))
  data$R <- c(rep(1L, n), rep(0L, nE))
  data$A <- c(rep(1L, nT), rep(0L, nC + nE))
  data$Y <- f + error - (1 - data$R) * X[, 1]
  list(data = data, seed = seed)
}

sample_values <- function(data) {
  stopifnot(identical(names(data), c(paste0("X", 1:10), "R", "A", "Y")),
            all(data$R %in% 0:1), all(data$A %in% 0:1))
  values <- as.double(as.matrix(data))
  stopifnot(all(is.finite(values)))
  values
}

# Pack all 13 columns in order, including all external outcomes. Binning on a
# binary grid avoids differences in the last few floating-point bits across
# math libraries. Integer-valued doubles represent the bins exactly, and the
# explicitly little-endian bytes do not depend on the R serialization version.
portable_sample_hash <- function(values, phase) {
  stopifnot(length(values) %% 13L == 0L, length(phase) == 1L, phase %in% 0:255)
  bins <- floor(values * 2^24 + phase / 256)
  stopifnot(all(abs(bins) < 2^53))
  payload <- c(length(values) / 13L, 13, phase, bins)
  digest::digest(writeBin(payload, raw(), size = 8L, endian = "little"),
                 serialize = FALSE, algo = "sha256")
}

# Historical raw hashes include the original column name. Only the temporary
# serialization copy uses that archived schema; analysis data retain R.
archived_sample_hash <- function(data) {
  names(data)[names(data) == "R"] <- "S"
  digest::digest(data, algo = "sha256")
}

check_paired_sample <- function(data, seed, reference, fingerprint, baseline = NULL) {
  stopifnot(nrow(reference) == 1L, nrow(fingerprint) == 1L, seed == reference$seed,
            reference$N_trial == fingerprint$N_trial, reference$replicate == fingerprint$replicate)
  raw_hash <- digest::digest(data, algo = "sha256")
  if (!is.null(baseline)) {
    stopifnot(nrow(baseline) == 1L, seed == baseline$seed,
              reference$N_trial == baseline$N_trial, reference$replicate == baseline$replicate)
    if (!baseline$full_sample_sha256 %in% c(raw_hash, archived_sample_hash(data)))
      stop("Exact pairing with the supplied baseline failed for trial size ", reference$N_trial,
           ", replicate ", reference$replicate, ".", call. = FALSE)
  }
  hash <- portable_sample_hash(sample_values(data), fingerprint$phase)
  if (hash != fingerprint$portable_sha256)
    stop("Full-sample fingerprint mismatch for trial size ", reference$N_trial,
         ", replicate ", reference$replicate, ".", call. = FALSE)
  list(raw = raw_hash, exact = archived_sample_hash(data) == reference$full_sample_sha256, portable = hash)
}

read_trial_fits <- function(here) {
  read.csv(file.path(here, "quadratic_rct_replicates.csv.gz"), check.names = FALSE)
}

# Compare both treatment coefficients and HC2 standard errors at the same seed.
check_paired_fits <- function(data, reference, quadratic_reference) {
  stopifnot(nrow(quadratic_reference) == 1L, nrow(reference) == 1L,
            quadratic_reference$N_trial == reference$N_trial,
            quadratic_reference$replicate == reference$replicate,
            quadratic_reference$seed == reference$seed)
  terms <- paste0("X", 1:10)
  linear <- hct_trial_wald_fit(data, reformulate(terms))
  quadratic <- hct_trial_wald_fit(data, reformulate(c(terms, paste0("I(", terms, "^2)"))))
  errors <- c(linear_estimate = abs(linear$estimate - reference$linear_estimate),
              linear_se = abs(linear$se - reference$linear_se),
              quadratic_estimate = abs(quadratic$estimate - quadratic_reference$estimate),
              quadratic_se = abs(quadratic$se - quadratic_reference$se))
  if (any(!is.finite(errors)) || any(errors >= 1e-10))
    stop("Archived trial estimates or HC2 SEs disagree for trial size ", reference$N_trial,
         ", replicate ", reference$replicate, ".", call. = FALSE)
  list(linear = linear, quadratic = quadratic, errors = errors)
}

# Construct reference hashes only when all underlying samples exactly match the
# original byte-level archive. A phase per sample keeps every original value
# at least 1e-12 away from a bin boundary; the grid width is 2^-24.
if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  if (length(args) != 1L || !grepl("[.]csv[.]gz$", args[1]) || file.exists(args[1]))
    stop("Usage: Rscript paired_samples.R NEW_FILE.csv.gz", call. = FALSE)
  script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  here <- dirname(normalizePath(script))
  reference <- read.csv(file.path(here, "inputs", "paired_reference.csv"))
  phase <- integer(nrow(reference))
  hash <- character(nrow(reference))
  for (j in seq_len(nrow(reference))) {
    row <- reference[j, ]
    d <- paired_sample(row$N_trial, row$replicate)
    stopifnot(row$seed == d$seed, identical(archived_sample_hash(d$data), row$full_sample_sha256))
    values <- sample_values(d$data)
    scaled <- values * 2^24
    # Widely separated neighboring candidates reduce the number of retries.
    for (offset in (0:255 * 73L + 128L) %% 256L) {
      bins <- scaled + offset / 256
      if (min(abs(bins - round(bins))) > 2^24 * 1e-12) break
    }
    stopifnot(min(abs(bins - round(bins))) > 2^24 * 1e-12)
    phase[j] <- offset
    hash[j] <- portable_sample_hash(values, offset)
  }
  index <- data.frame(N_trial = reference$N_trial, replicate = reference$replicate,
                      phase = phase, portable_sha256 = hash)
  con <- gzfile(args[1], "wt")
  tryCatch(write.csv(index, con, row.names = FALSE, quote = FALSE), finally = close(con))
  cat("Certified", nrow(index), "portable fingerprints against the archived raw hashes.\n")
}
