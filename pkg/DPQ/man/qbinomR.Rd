\name{qbinomR}
\alias{qbinomR}
\title{Pure R Implementation of R's qbinom()  with Tuning Parameters}
\description{%-->>> ../R/qbinomR.R <<<<<<<<<
  A pure R implementation, including many tuning parameter arguments, of
  \R's own Rmathlib C code algorithm, but with more flexibility.

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
  %% In R, since {r80271 | maechler | 2021-05-08}, the joint C code is in <Rsrc/src/nmath/qDiscrete_search.h
  \item{yLarge}{when \eqn{y >= y_L, y_L =} \code{yLarge}, the binary root
    finding search is made  \dQuote{cleverer}, taking larger increments,
    determined by \code{incF} and \code{iShrink}:}
  \item{incF}{a positive \dQuote{increment factor} (originally hardwired to
    0.001), used only when \code{y >= yLarge}; defines the initial increment in
    the search algorithm as \code{incr <- floor(incF * y)}.}
  \item{iShrink}{a positive increment shrinking factor, used only when
    \code{y >= yLarge} to define the new increment from the old one as
    \code{incr <- max(1, floor(incr/iShrink))} where the LHS was hardired
    original to \code{(incr/100)}.}
  \item{relTol}{relative tolerance, \eqn{ > 0}; the search terminates when
    the (integer!) increment is less than \code{relTol * y} or the previous
    increment was not larger than 1.}
  \item{pfEps.n}{fuzz factor to ensure left continuity in the \bold{n}ormal
    case \code{log.p=FALSE}; used to be hardwired to 64 (in \R up to 2021-05-08).}
  \item{pfEps.L}{fuzz factor to ensure left continuity in case
    \code{log.p=TRUE}; used to be hardwired to 64 (in \R up to 2021-05-08).}
  \item{fpf}{factor \eqn{f \ge 1}{f >= 1} for the \bold{n}ormal upper tail case
    (\code{log.p=FALSE, lower.tail=FALSE}): \code{p} is only \dQuote{fuzz-corrected},
    i.e., multiplied by \eqn{1+e} when \code{1 - p > fpf*e} for
    \code{e <- pfEps.n * c_e} and \eqn{c_e = 2^{-52}}, the \code{.Machine$double_epsilon}.
  }
  \item{trace}{logical (or integer) specifying if (and how much) output
    should be produced from the algorithm.}
}
\details{as mentioned on \code{\link{qbinom}} help page, % ~/R/D/r-devel/R/src/library/stats/man/Binomial.Rd
  \code{qbinom} uses the Cornish--Fisher Expansion to include a skewness
  correction to a normal approximation, thus defining \code{y := Fn(p, size, prob, ..)}.

  The following (root finding) binary search is tweaked by the
  \code{yLarge, ..., fpf} arguments.
}
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
