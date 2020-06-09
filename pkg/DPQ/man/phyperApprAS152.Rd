\name{phyperApprAS152}
\alias{phyperApprAS152}
\title{Normal Approximation to cumulative Hyperbolic Distribution -- AS 152}
\description{
  Compute the normal approximation (via \code{\link{pnorm}(.)} from AS 152
  to the cumulative hyperbolic distribution function \code{\link{phyper}()}.
}
\usage{
phyperApprAS152(q, m, n, k)
}
\arguments{
  \item{q}{vector of quantiles representing the number of white balls
    drawn without replacement from an urn which contains both black and
    white balls.}
  \item{m}{the number of white balls in the urn.}
  \item{n}{the number of black balls in the urn.}
  \item{k}{the number of balls drawn from the urn, hence must be in
    \eqn{0,1,\dots, m+n}.}
}
%% \details{
%% }
\value{
  a \code{\link{numeric}} vector of the same length (etc) as \code{q}.
}
\references{
  Lund, Richard E. (1980)
  Algorithm AS 152: Cumulative Hypergeometric Probabilities.
  \emph{Journal of the Royal Statistical Society. Series C (Applied Statistics)}, \bold{29}(2), 221--223.
  \doi{10.2307/2986315}

  Shea, B. (1989)
  Remark AS R77: A Remark on Algorithm AS 152: Cumulative Hypergeometric Probabilities.
  \emph{JRSS C (Applied Statistics)}, \bold{38}(1), 199--204.
  \doi{10.2307/2347696}

  Berger, R. (1991)
  Algorithm AS R86: A Remark on Algorithm AS 152: Cumulative Hypergeometric Probabilities.
  \emph{JRSS C (Applied Statistics)}, \bold{40}(2), 374--375.
  \doi{10.2307/2347606}
}
\author{Martin Maechler, 19 Apr 1999}
\note{% ~/src/ApStat/152-chyper.f ~/src/ApStat/152-chyper-1.c
  I have Fortran (and C code translated from Fortran) which says
 \preformatted{
   ALGORITHM AS R77  APPL. STATIST. (1989), VOL.38, NO.1
   Replaces AS 59 and AS 152
   Incorporates AS R86 from vol.40(2)
 }
}

\seealso{
  \code{\link{phyper}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (q, m, n, k)
{
    kk <- n
    nn <- m
    mm <- m + n
    ll <- q
    mean <- kk * nn/mm
    sig <- sqrt(mean * (mm - nn)/mm * (mm - kk)/(mm - 1))
    pnorm(ll + 1/2, mean = mean, sd = sig)
  }
}
\keyword{distribution}
