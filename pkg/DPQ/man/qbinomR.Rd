\name{qbinomR}
\alias{qbinomR}
\title{Pure R Implementation of R's qbinom()  with Tuning Parameters}
\description{%-->>> ../R/qbinomR.R <<<<<<<<<
  A pure R implementation, including many tuning parameter arguments, of
  \R's own Mathlib C code algorithm, but with more flexibility.

  It is using \code{\link{Vectorize}(qbinomR1, *)} where the hidden
  \code{qbinomR1} works for numbers (aka \sQuote{scalar}, length one)
  arguments only, the same as the C code.
}
\usage{
qbinomR(p, size, prob, lower.tail = TRUE, log.p = FALSE,
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
  \item{p, size, prob, lower.tail, log.p}{\code{\link{qbinom}()} standard
    argument, see its help page.}%  ~/R/D/r-devel/R/src/library/stats/man/Binomial.Rd
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
  \code{p}, \code{size}, and \code{prob}.
}
%% \references{

%% }
\author{Martin Maechler}
\seealso{
  \code{\link{qbinom}}, \code{\link{qpois}}.
}
\examples{
set.seed(12)
pr <- (0:16)/16 # supposedly recycled
x10 <- rbinom(500, prob=pr, size =  10); p10 <- pbinom(x10, prob=pr, size= 10)
x1c <- rbinom(500, prob=pr, size = 100); p1c <- pbinom(x1c, prob=pr, size=100)
## stopifnot(exprs = {
table( x10  == (qp10  <- qbinom (p10, prob=pr, size= 10) ))
table( qp10 == (qp10R <- qbinomR(p10, prob=pr, size= 10) )); summary(warnings()) # 30 x NaN
table( x1c  == (qp1c  <- qbinom (p1c, prob=pr, size=100) ))
table( qp1c == (qp1cR <- qbinomR(p1c, prob=pr, size=100) )); summary(warnings()) # 30 x NaN
## })
}
\keyword{distribution}
