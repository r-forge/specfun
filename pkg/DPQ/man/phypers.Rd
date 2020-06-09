\name{phypers}
\alias{phypers}
\title{The Four (4) Symmetric phyper() calls.}
\description{
  Compute the four (4) symmetric \code{\link{phyper}()} calls which
  mathematically would be identical but in practice typically slightl
  differ numerically.
}
\usage{
phypers(m, n, k, q = .suppHyper(m, n, k))
}
\arguments{
  \item{m}{the number of white balls in the urn.}
  \item{n}{the number of black balls in the urn.}
  \item{k}{the number of balls drawn from the urn, hence must be in
    \eqn{0,1,\dots, m+n}.}
  \item{q}{vector of quantiles representing the number of white balls
    drawn without replacement from an urn which contains both black and
    white balls.  The default ......... FIXME
  }
}
%% \details{
%% }
\value{
  a list with components
  \item{q}{Description of 'comp1'}
  \item{phyp}{a numeric \code{\link{matrix}} of 4 columns with the 4
    different calls to \code{\link{phyper}()} which are theoretically
    equivalent because of mathematical symmetry.}
}
\references{
  Johnson et al
}
\author{Martin Maechler}
\seealso{
 \code{\link{phyper}}.
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (m, n, k, q = .suppHyper(m, n, k))
{
    N <- m + n
    pm <- cbind(ph = phyper(q, m, n, k), p2 = phyper(q, k, N -
        k, m), Ip2 = phyper(m - 1 - q, N - k, k, m, lower.tail = FALSE),
        Ip1 = phyper(k - 1 - q, n, m, k, lower.tail = FALSE))
    stopifnot(all.equal(pm[, 1], pm[, 2]), all.equal(pm[, 2],
        pm[, 3]), all.equal(pm[, 3], pm[, 4]))
    list(q = q, phyp = pm)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
