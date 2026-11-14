\name{precision_share}
\alias{precision_share}
\title{Internal share of the precision-weighted control mean}
\usage{precision_share(p, e, v_internal, v_external)}
\arguments{
  \item{p}{Trial membership probability given baseline covariates.}
  \item{e}{Treatment allocation probability in the trial given baseline covariates.}
  \item{v_internal, v_external}{Conditional outcome variances in the two control groups.}
}
\description{Calculate the internal share separately at each covariate profile.}
\details{
For a shared profile the internal control precision is proportional to
\eqn{p(1-e)/v_i} and the external control precision to \eqn{(1-p)/v_e}.
The share is the first precision divided by their sum. Variances must be
strictly positive on shared support. Treatment probability must be strictly
between zero and one on trial support. Each argument may be a scalar or a
vector of the common length.
}
\value{A numeric vector of internal shares. Trial-only profiles have share
one. External-only profiles have no trial target and return \code{NA_real_}.
Inputs outside the required support may be missing.}
\examples{
precision_share(p = c(.6, 1, 0), e = c(.5, .5, NA),
                v_internal = c(1, NA, NA), v_external = c(1, NA, NA))
}
\seealso{\code{\link{aipw_u}}}
