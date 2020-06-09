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
  \item{q}{ . }
  \item{m}{ . }
  \item{n}{ . }
  \item{k}{ . }
}
\details{
  Both approximations are from page 261 of J Johnson, Kotz & Kemp (1992).
  \code{phyper1molenaar} is formula \eqn{(6.91)}, and
  \code{phyper2molenaar} is formula \eqn{(6.92)}.
}
\value{
 . . .
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
 ## TODO
}
\keyword{distribution}
