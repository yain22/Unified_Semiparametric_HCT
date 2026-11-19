# Figure S1: six invented observations illustrating the trial designs.
# Save this file beside observed-records.csv, then run Rscript observed-records.R.
# Or use source("docs/examples/observed-records.R", chdir = TRUE) in the project.
# When sourced with chdir = TRUE, R temporarily uses the script directory.
script_dir <- "."
if (sys.nframe() == 0L) {
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(arg)) script_dir <- dirname(sub("^--file=", "", arg[1L]))
}
records <- read.csv(file.path(script_dir, "observed-records.csv"))

# R = 1 denotes trial membership; A = 1 denotes treatment.
records$I1 <- with(records, R * A)
records$I10 <- with(records, R * (1L - A))
records$I0 <- 1L - records$R
stopifnot(all(rowSums(records[c("I1", "I10", "I0")]) == 1L))

# A common-mean restriction changes the hybrid model, not the rows observed.
# Removing the external or internal controls changes which groups supply data.
designs <- list(
  "Unrestricted HCT" = records,
  "Common-mean HCT" = records,
  "RCT" = subset(records, R == 1L),
  "ECT" = subset(records, I1 == 1L | I0 == 1L)
)
for (design in names(designs)) {
  cat("\n", design, "\n", sep = "")
  print(designs[[design]], row.names = FALSE)
}
