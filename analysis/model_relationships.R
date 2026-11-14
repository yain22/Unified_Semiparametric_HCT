# Observed-data models for hybrid, randomized, and externally controlled trials.
args <- commandArgs(TRUE)
if (any(args %in% c("--help", "-h"))) {
  cat("Usage: Rscript analysis/model_relationships.R [--output-dir PATH]\n")
  quit()
}
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg))))
out <- file.path(root, "results", "model-relationships")
if (length(args)) {
  if (length(args) != 2L || args[1L] != "--output-dir" || !nzchar(args[2L]))
    stop("Expected --output-dir PATH.")
  out <- args[2L]
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)

draw_models <- function() {
  library(grid)
  grid.newpage()
  ink <- "#172331"
  blue <- "#2667A8"
  green <- "#167565"
  orange <- "#AE5C17"
  gray <- "#566373"
  label <- function(x, y, text, size, colour = ink, bold = FALSE) {
    grid.text(text, x, y, gp = gpar(fontfamily = "sans", fontsize = size,
              col = colour, fontface = if (bold) "bold" else "plain", lineheight = .96))
  }
  box <- function(x, y, width, height, fill, border) {
    grid.roundrect(x, y, width, height, r = unit(.035, "inches"),
                   gp = gpar(fill = fill, col = border, lwd = 1))
  }
  ellipse <- function(x, y, rx, ry, fill, border) {
    angle <- seq(0, 2 * pi, length.out = 201L)
    grid.polygon(x + rx * cos(angle), y + ry * sin(angle),
                 gp = gpar(fill = fill, col = border, lwd = 1.35))
  }
  branch <- function(x, bend, colour) {
    grid.lines(c(x, bend, bend), c(.42, .42, .305),
               arrow = arrow(type = "closed", length = unit(.045, "inches")),
               gp = gpar(col = colour, fill = colour))
  }
  box(.5, .535, .975, .91, "#FCFCFD", gray)
  label(.5, .946, "Unified causal-gap sensitivity framework", 7.7, bold = TRUE)
  ellipse(.5, .646, .463, .251, "#EDF4FC", blue)
  label(.5, .806, "Unrestricted HCT\nobserved-data model", 7.1, blue, TRUE)
  label(.5, .719, expression(italic(M)["HCT,"~U]), 8.8, blue)
  ellipse(.5, .552, .365, .123, "#DFF1EB", green)
  label(.5, .619, "HCT under transportability", 6.8, green, TRUE)
  label(.5, .554, expression(italic(M)["HCT,trans"] %subset% italic(M)["HCT,"~U]), 7.9, green)
  label(.5, .491, expression(paste(m["1,0"](x) == m[0](x), "  ", symbol("\336"), "  ",
                                   xi[lambda*","*delta] == 0)), 7.2, green)
  branch(.312, .253, gray)
  branch(.683, .751, orange)
  label(.118, .364, "Remove external\ncontrol records", 5.4, gray)
  label(.884, .364, "Remove internal\ncontrol records", 5.4, orange)
  label(.349, .353, expression(lambda %==% 1), 7.8, gray)
  label(.653, .353, expression(lambda %==% 0), 7.8, orange)
  box(.253, .204, .432, .181, "#F1F3F6", gray)
  box(.751, .204, .432, .181, "#FFF3E8", orange)
  label(.253, .237, "RCT", 7.9, gray, TRUE)
  label(.751, .237, "ECT", 7.9, orange, TRUE)
  label(.253, .155, expression(italic(M)["RCT"]), 8.8, gray)
  label(.751, .155, expression(italic(M)["ECT"]), 8.8, orange)
  label(.5, .034, "Design reductions redefine the observation law", 6.1)
}
for (extension in c("pdf", "png")) {
  path <- file.path(out, paste0("observed_model_relationships.", extension))
  if (extension == "pdf") {
    cairo_pdf(path, width = 3.25, height = 2.12)
  } else {
    png(path, width = 3.25, height = 2.12, units = "in", res = 360, type = "cairo")
  }
  tryCatch(draw_models(), finally = dev.off())
}
cat("Wrote observed_model_relationships.pdf and .png\n")
