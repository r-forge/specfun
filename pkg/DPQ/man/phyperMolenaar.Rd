\name{phyperMolenaar}
\title{Molenaar's Normal Approximations to the Hypergeometric Distribution}
\alias{phyper1molenaar}
\alias{phyper2molenaar}
\description{
  Compute Molenaar's two normal approximations to the (cumulative
  hypergeometric distribution \code{\link{phyper}()}.
}
\usage{
phyper1molenaar(q, m, n, k)
phyper2molenaar(q, m, n, k)
}
\arguments{
  \item{q}{(vector of) the number of white balls drawn without replacement
    from an urn which contains both black and white balls.}
  \item{m}{the number of white balls in the urn.}
  \item{n}{the number of black balls in the urn.}
  \item{k}{the number of balls drawn from the urn, hence in \eqn{0,1,\dots,m+n}.}
}
\details{
  Both approximations are from page 261 of Johnson, Kotz & Kemp (1992).
  \code{phyper1molenaar} is formula \eqn{(6.91)}, and
  \code{phyper2molenaar} is formula \eqn{(6.92)}.
}
\value{
 a \code{\link{numeric}} vector, with the length the maximum of the
 lengths of \code{q, m, n, k}.
}
\references{
  Johnson, Kotz & Kemp (1992): p.261
}
\author{Martin Maechler}
\seealso{
  \code{\link{phyper}},
  \code{\link{pnorm}}.
}
\examples{
 ## TODO -- maybe see  ../tests/hyper-dist-ex.R
}
\keyword{distribution}
