\name{qpoisR}
\alias{qpoisR}
\title{Pure R Implementation of R's qpois()  with Tuning Parameters}
\description{%-->>> ../R/qpoisR.R <<<<<<<<<
  A pure R implementation, including many tuning parameter arguments, of
  \R's own Mathlib C code algorithm, but with more flexibility.

  It is using \code{\link{Vectorize}(qpoisR1, *)} where the hidden
  \code{qpoisR1} works for numbers (aka \sQuote{scalar}, length one)
  arguments only, the same as the C code.
}
\usage{
qpoisR(p, lambda, lower.tail = TRUE, log.p = FALSE,
       yLarge = 4096, # was hard wired to 1e5
       incF = 1/64,   # was hard wired to .001
       iShrink = 8,   # was hard wired to 100
       relTol = 1e-15,# was hard wired to 1e-15
       pfEps.n = 8,   # was hard wired to 64: "fuzz to ensure left continuity"
       pfEps.L = 2,   # was hard wired to 64:   "   "   ..
       fpf = 4, # *MUST* be >= 1 (did not exist previously)
       trace = 0)
}
\arguments{
  \item{p, lambda, lower.tail, log.p}{\code{\link{qpois}()} standard
    argument, see its help page.}
  \item{yLarge}{
  }
  \item{incF}{
}
  \item{iShrink}{
}
  \item{relTol}{
}
  \item{pfEps.n}{
}
  \item{pfEps.L}{
}
  \item{fpf}{
}
  \item{trace}{logical (or integer) specifying if (and how much) output
    should be produced from the algorithm.}
}
%% \details{
%%
%% }
\value{
  a numeric vector like \code{p} recycled to the common lengths of \code{p}
  and \code{lambda}.
}
%% \references{

%% }
\author{Martin Maechler}
\seealso{
  \code{\link{qpois}}.
}
\examples{
x <- 10*(15:25)
Pp <- ppois(x, lambda = 100, lower.tail = FALSE)  # no cancellation
qPp <- qpois(Pp, lambda = 100, lower.tail=FALSE)
table(x == qPp) # all TRUE ?
if(getRversion() >= "4.2") stopifnot(x == qPp) # R-devel
qpRp <- qpoisR(Pp, lambda = 100, lower.tail=FALSE)
all.equal(x, qpRp, tol = 0)
stopifnot(all.equal(x, qpRp, tol = 1e-15))
}
\keyword{distribution}
