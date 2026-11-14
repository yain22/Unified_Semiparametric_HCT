# Exact finite version of the example in the EIF lab.
# Each row is equally likely; the 72 rows enumerate the law, not a random sample.
library(causalgap)

hybrid_example <- function(delta = 0.2, external_variance = 1) {
  z <- c(-sqrt(3), 0, 0, 0, 0, sqrt(3))
  x <- rep(0:1, each = 36)
  s <- rep(rep(c(1, 1, 1, 0, 0, 0), each = 6), 2)
  a <- rep(rep(c(1, 1, 0, 0, 0, 0), each = 6), 2)
  m <- 0.4 * x
  fits <- data.frame(p = rep(0.5, length(x)), e = 2/3,
                     m_treated = m + 1, m_internal = m, m_external = m - delta,
                     v_internal = 1, v_external = external_variance)
  y <- ifelse(a == 1, m + 1, ifelse(s == 1, m, m - delta)) +
    ifelse(s == 1, 1, sqrt(external_variance)) * rep(z, 12)
  aipw_u(data.frame(S = s, A = a, Y = y), fits)
}

fit <- hybrid_example()
print(fit$estimate)                          # tce = 1.15, trial = 1, gap = 0.15
print(unique(fit$lambda))                    # internal-control share = 0.25
print(mean(fit$influence$tce^2))             # variance bound per record = 4.528125
print(sqrt(mean(fit$influence$tce^2) / 300)) # asymptotic SE at N = 300

# At common control means, the unrestricted gap can still vary along paths
# that move the sources apart: its EIF need not vanish at a zero gap.
common <- hybrid_example(delta = 0)
print(common$estimate["gap"])                # zero
print(mean(common$influence$gap^2))          # 4.5

stopifnot(isTRUE(all.equal(unname(fit$estimate), c(1.15, 1, 0.15))),
          isTRUE(all.equal(mean(fit$influence$tce^2), 4.528125)),
          isTRUE(all.equal(mean(common$influence$gap^2), 4.5)))
