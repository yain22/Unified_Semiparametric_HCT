\name{aipw_u}
\alias{aipw_u}
\title{Evaluate the unrestricted hybrid-trial estimator and its influence function}
\usage{aipw_u(data, nuisance)}
\arguments{
  \item{data}{Data frame containing the numeric columns \code{R} (1 for trial,
  0 for external), \code{A} (1 for treated, 0 for control), and \code{Y}
  (observed outcome). External records must have \code{A = 0}.}
  \item{nuisance}{Data frame of predictions at each participant's baseline
  covariates in the same row order as \code{data}. Supply \code{p} (trial
  membership probability), \code{e} (known treatment allocation probability
  within the trial), conditional outcome means \code{m_treated},
  \code{m_internal}, \code{m_external}, and conditional control variances
  \code{v_internal}, \code{v_external}.}
}
\description{
Solve the estimating equation for the precision-weighted borrowing target and
the matched trial-only contrast. Evaluate the empirical influence functions at
the resulting estimates. All required nuisance functions are supplied by the
caller, which permits predictions from sample splitting or a specified design law.
}
\details{
The internal and external control means are combined using
\code{\link{precision_share}}. The conditional treatment contrast is then
averaged over trial covariates. The augmented estimator includes outcome
residuals and three corrections for variation in source probability and
the two control variances. These corrections account for how the precision
share changes with the observed-data law.
The gap is the borrowing target minus the trial contrast. Under trial
consistency, randomization and positivity the latter identifies the
trial treatment effect.

On trial-only profiles (\code{p = 1}), only trial outcomes and predictions
for the randomized comparison are used. On external-only profiles
(\code{p = 0}), all scores are zero. On shared profiles the two control
variances must be finite and strictly positive. The treated mean is
required on trial records and the two control means on shared profiles.
The internal mean is also required on trial-only profiles. Predictions
outside these requirements may be missing.

For independent identically distributed mixture sampling and regular nuisance
estimation, the sample variance of an influence column divided by the sample
size estimates the corresponding estimator variance. For independent groups
with counts fixed by design, use within-group influence variances and the
sampling fractions appropriate to that design.
}
\value{
A list with \code{estimate} (named numeric vector containing \code{tce},
\code{trial} and \code{gap}), \code{influence} (one row per observation and
one column per target), \code{components} (seven signed, centered
contributions that sum across each row to the borrowing-target influence
value), and \code{lambda} (internal share at each profile).
}
\examples{
# Toy hybrid trial with a common covariate profile, two participants per group.
d <- data.frame(R = c(1, 1, 1, 1, 0, 0), A = c(1, 1, 0, 0, 0, 0),
                Y = c(3, 5, 2, 4, 1, 3))
f <- data.frame(p = rep(2/3, 6), e = .5,
                m_treated = 4, m_internal = 3, m_external = 2,
                v_internal = 1, v_external = 1)
ans <- aipw_u(d, f)
ans$estimate
ans$components
}
