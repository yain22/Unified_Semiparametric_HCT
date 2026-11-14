library(causalgap)

stopifnot(all.equal(precision_share(.6, .5, 1, 1), 3/7),
          all.equal(precision_share(c(1,0), c(.5,NA_real_), NA_real_, NA_real_),
                    c(1, NA_real_)))

# Integrate the example in Section 4 over a finite residual distribution
# with second and fourth moments 1 and 3. There are two treated, one trial
# control, and three external blocks, yielding p=q=1/2 and e=2/3.
u <- c(-sqrt(3), 0, 0, 0, 0, sqrt(3))
s <- rep(c(1, 1, 1, 0, 0, 0), each = 6)
a <- rep(c(1, 1, 0, 0, 0, 0), each = 6)
n <- length(s)
group_mean <- ifelse(a == 1, 1, ifelse(s == 1, .2, 0))
d <- data.frame(S = s, A = a, Y = group_mean + u)
f <- data.frame(p = rep(.5,n), e = 2/3, m_treated = 1,
                m_internal = .2, m_external = 0, v_internal = 1, v_external = 1)
ans <- aipw_u(d, f)
stopifnot(all.equal(ans$estimate, c(tce=.95, trial=.8, gap=.15)),
          all.equal(rowSums(ans$components), ans$influence$tce),
          max(abs(colMeans(ans$influence))) < 1e-12,
          all.equal(mean(ans$influence$tce^2), 4.528125),
          all.equal(mean(rowSums(ans$components[, 1:4])^2), 4.5),
          all.equal(ans$lambda, rep(.25, n)))

# A regular tilt of the external outcome distribution changes its mean
# and variance. The full EIF reproduces its finite-difference derivative.
tilt <- function(z) z + .4 * (z^2 - 1)
path_target <- function(eps) {
  w <- exp(eps * tilt(u))
  w <- w / sum(w)
  m <- sum(w*u)
  v <- sum(w*(u-m)^2)
  lam <- precision_share(.5, 2/3, 1, v)
  1 - lam*.2 - (1-lam)*m
}
h <- 1e-5
deriv <- (path_target(h) - path_target(-h)) / (2*h)
score <- (1-s) * rep(tilt(u), 6)
stopifnot(abs(deriv - mean(ans$influence$tce * score)) < 1e-8,
          abs(deriv - mean(ans$influence$gap * score)) < 1e-8)

# Trial-only and external-only strata may omit undefined regressions.
f$p <- s
f$e[s == 0] <- NA_real_
f$m_internal[s == 0] <- f$m_treated[s == 0] <- NA_real_
f$m_external <- f$v_internal <- f$v_external <- NA_real_
edge <- aipw_u(d, f)
stopifnot(abs(edge$estimate["gap"]) < 1e-12,
          all(abs(edge$influence$gap) < 1e-12),
          all(edge$components[s == 0, ] == 0))

f$p[1] <- 0
stopifnot(inherits(try(aipw_u(d, f), silent = TRUE), "try-error"))
