# No patient data required: full-score algebra and nontrivial variance checks.
args <- commandArgs(trailingOnly = TRUE)
script <- normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]), winslash = '/')
source(file.path(dirname(script), 'run_clinical.R'))
h <- clinical_load(dirname(script))
set.seed(7102026L)
n <- 300L
d <- data.frame(
  R = c(rep(1, 200), rep(0, 100)), A = c(rep(1, 100), rep(0, 200)),
  hba1c0 = rnorm(n), age = rnorm(n), duration = rnorm(n), male = rbinom(n, 1, .45), pump = rbinom(n, 1, .6)
)
d$Y <- .2 * d$age - .3 * d$A + .25 * (1 - d$R) + rnorm(n) * exp(.15 * d$hba1c0)
fit <- h$cgm_est_nuisance(d, h)
a <- fit$predictions
s <- h$cgm_est_score(d, a)
e <- .5
R <- d$R
T <- R * d$A
C <- R * (1 - d$A)
E <- 1 - R
lambda <- a$p * (1 - e) / (a$p * (1 - e) + (1 - a$p) * a$vC / a$vE)
delta <- a$mC - a$mE
# Independent assembly of each unrestricted-score term, with raw residuals.
fixed <- R * (a$mT - lambda * a$mC - (1 - lambda) * a$mE) + T / e * (d$Y - a$mT) -
  C * lambda / (1 - e) * (d$Y - a$mC) - E * (1 - lambda) * a$odds * (d$Y - a$mE)
cs <- -delta * lambda * (1 - lambda) * (R - a$p) / (1 - a$p)
cc <- delta * lambda * (1 - lambda) * C / (1 - e) * ((d$Y - a$mC)^2 / a$vC - 1)
ce <- -delta * lambda * (1 - lambda) * E * a$odds * ((d$Y - a$mE)^2 / a$vE - 1)
term_error <- max(abs(unname(s$corrections) - c(sum(fixed), sum(cs), sum(cc), sum(ce)) / sum(R)))
stopifnot(
  term_error < 1e-12, abs(s$estimate - sum(fixed + cs + cc + ce) / sum(R)) < 1e-12,
  abs(sum(s$corrections) - s$estimate) < 1e-12, abs(mean(s$phi)) < 1e-12
)
# At pointwise equal fitted control means all three derivative terms vanish.
same <- a
same$mE <- same$mC
s0 <- h$cgm_est_score(d, same)
stopifnot(max(abs(s0$corrections[-1])) == 0)
# Outcome-unit translation leaves the contrast and complete influence unchanged.
shift <- d
shift$Y <- d$Y + 7
ashift <- a
ashift[c('mT', 'mC', 'mE')] <- ashift[c('mT', 'mC', 'mE')] + 7
ss <- h$cgm_est_score(shift, ashift)
translation_error <- max(abs(c(ss$estimate - s$estimate, ss$phi - s$phi)))
stopifnot(translation_error < 1e-10)
# Group-centered variance must use n_g sample-var(phi_g)/N^2, not iid variance.
g <- h$cgm_group(d)
v <- h$hct_unrestricted_group_variance(s$phi, g)$se^2
manual_v <- sum(vapply(split(s$phi, g), function(x) sum((x - mean(x))^2) * length(x) / (length(x) - 1), numeric(1))) / n^2
stopifnot(abs(v - manual_v) < 1e-14)
# Numerical derivatives of the precision target verify all correction signs.
target <- function(p, vc, ve) a$mT - (p * (1 - e) / (p * (1 - e) + (1 - p) * vc / ve)) * a$mC -
  (1 - p * (1 - e) / (p * (1 - e) + (1 - p) * vc / ve)) * a$mE
eps <- 1e-6
fdp <- (target(a$p + eps, a$vC, a$vE) - target(a$p - eps, a$vC, a$vE)) / (2 * eps)
fdc <- (target(a$p, a$vC + eps, a$vE) - target(a$p, a$vC - eps, a$vE)) / (2 * eps)
fde <- (target(a$p, a$vC, a$vE + eps) - target(a$p, a$vC, a$vE - eps)) / (2 * eps)
derivative_error <- max(abs(c(
  fdp + delta * lambda * (1 - lambda) / (a$p * (1 - a$p)),
  fdc - delta * lambda * (1 - lambda) / a$vC, fde + delta * lambda * (1 - lambda) / a$vE
)))
stopifnot(derivative_error < 1e-7)
stopifnot(max(fit$training_diagnostics$variance_equation_error) < 1e-10)
out_arg <- grep('^--out=', args, value = TRUE)
result <- list(
  passed = TRUE, synthetic_n = n, term_error = term_error, translation_error = translation_error,
  derivative_error = derivative_error, group_variance_error = abs(v - manual_v),
  variance_equation_max_error = max(fit$training_diagnostics$variance_equation_error),
  tests = c('four-term score', 'zero discrepancy', 'outcome translation', 'fixed-group variance', 'share derivatives', 'variance solver')
)
out <- if (length(out_arg)) sub('^--out=', '', out_arg[1]) else file.path(dirname(script), '../../results/validation-clinical')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(result, file.path(out, 'synthetic_validation.json'), pretty = TRUE, auto_unbox = TRUE, digits = 16)
print(result)
