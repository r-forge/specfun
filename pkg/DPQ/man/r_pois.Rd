\name{r_pois}
\alias{r_pois}
\alias{r_pois_expr}
\alias{plRpois}
\title{Compute the i-th term of the Poisson Distribution Series}% FIXME?
\description{
Purpose: r_v(i) :=  (v^i / i!) / e_{i-1}(v), where
     e_n(x) := 1 + x + x^2/2! + .... + x^n/n! (n-th partial of exp(x))

As function of i :
 o  Can this be put in a simple formula?
 o  When is it maximal?
 o  When does r_{lambda}(i) become smaller than (f+2i-x)/x = a + b*i ?

}
\usage{
r_pois(i, lambda)
r_pois_expr  # the R expression() for the asymptotic branch of r_pois()

plRpois(lambda, iset = 1:(2*lambda), do.main = TRUE,
        log = 'xy', cex = 0.4, col = c("red","blue"))
}
\arguments{
  \item{i}{integer ..}
  \item{lambda}{non-negative number ... }
  \item{iset}{ ..... }
  \item{do.main}{logical specifying if a main \code{\link{title}} should be
        drawn (via \code{\link{mtext}(..)}).}
  \item{log}{string specifying if (and where) logarithmic scales should be
        used, see \code{\link{plot.default}()}.}
  \item{cex}{character expansion factor.}
  \item{col}{color}
}
\details{

  ...........


  \code{plRpois()} simply produces a \dQuote{nice} plot of \code{r_pois(ii, *)}
  vs \code{ii}.
}
\value{
  .........
}
\references{
  ...
}
\author{Martin Maechler, 20 Jan 2004}
%% \note{

%% }

\seealso{
  \code{\link{dpois}()}.
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (i, lambda)
{
    if (lambda < -log(2) * .Machine$double.min.exp)
        dpois(i, lambda)/ppois(i - 1, lambda)
    else exp(i * log(lambda) - lgamma(i + 1) + -lambda - ppois(i -
        1, lambda, log.p = TRUE))
  }

}% ex.
\keyword{distribution}
