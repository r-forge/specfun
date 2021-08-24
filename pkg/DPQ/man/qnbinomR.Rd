\name{qnbinomR}
\alias{qnbinomR}
\title{Pure R Implementation of R's qnbinom()  with Tuning Parameters}
\description{%-->>> ../R/qnbinomR.R <<<<<<<<<
  A pure R implementation, including many tuning parameter arguments, of
  \R's own Rmathlib C code algorithm, but with more flexibility.

  It is using \code{\link{Vectorize}(qnbinomR1, *)} where the hidden
  \code{qnbinomR1} works for numbers (aka \sQuote{scalar}, length one)
  arguments only, the same as the C code.
}
\usage{
qnbinomR(p, size, prob, mu, lower.tail = TRUE, log.p = FALSE,
         yLarge = 4096, # was hard wired to 1e5
         incF = 1/64,   # was hard wired to .001
         iShrink = 8,   # was hard wired to 100
         relTol = 1e-15,# was hard wired to 1e-15
         pfEps.n = 8,   # was hard wired to 64: "fuzz to ensure left continuity"
         pfEps.L = 2,   # was hard wired to 64:   "   "   ..
         fpf = 4, # *MUST* be >= 1 (did not exist previously)
         trace = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p, size, prob, mu, lower.tail, log.p}{\code{\link{qnbinom}()} standard
    argument, see its help page.}%  ~/R/D/r-devel/R/src/library/stats/man/NegBinomial.Rd
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
  a numeric vector like \code{p} recycled to the common lengths of
  \code{p}, \code{size}, and either \code{prob} or \code{mu}.
}
%% \references{

%% }
\author{Martin Maechler}
\seealso{
  \code{\link{qnbinom}}, \code{\link{qpois}}.
}
\examples{
set.seed(12)
x10 <- rnbinom(500, mu = 4,       size = 10) ; p10 <- pnbinom(x10, mu=4,       size=10)
x1c <- rnbinom(500, prob = 31/32, size = 100); p1c <- pnbinom(x1c, prob=31/32, size=100)
stopifnot(exprs = {
    x10 == qnbinom (p10, mu=4, size=10)
    x10 == qnbinomR(p10, mu=4, size=10)
    x1c == qnbinom (p1c, prob=31/32, size=100)
    x1c == qnbinomR(p1c, prob=31/32, size=100)
})
}
\keyword{distribution}
