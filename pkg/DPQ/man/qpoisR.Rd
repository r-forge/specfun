\name{qpoisR}
\alias{qpoisR}
\title{Pure R Implementation of R's qpois()  with Tuning Parameters}
\description{%-->>> ../R/qpoisR.R <<<<<<<<<
  A pure R implementation, including many tuning parameter arguments, of
  \R's own Rmathlib C code algorithm, but with more flexibility.

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
  \item{yLarge}{a positive number; in \R up to 2021, was internally
    hardwired to \code{yLarge = 1e5}:  Uses more careful search for
    \eqn{y \ge y_L}{y >= yLarge}, where \eqn{y} is the initial approximate result,
    derived from a Cornish-Fisher expansiion.
    %% we could give it:
    %% sigma <- sqrt(lambda)
    %% gamma <- 1/sigma
    %% z <- qnorm(p, lower.tail=lower.tail, log.p=log.p)
    %% y <- round(mu + sigma * (z + gamma * (z*z - 1) / 6))
  }
  \item{incF}{a positive \dQuote{increment factor} (originally hardwired to
    0.001), used only when \code{y >= yLarge}; defines the initial increment in
    the search algorithm as \code{incr <- floor(incF * y)}.}
  \item{iShrink}{a positive increment shrinking factor, used only when
    \code{y >= yLarge} to define the new increment from the old one as
    \code{incr <- max(1, floor(incr/iShrink))} where the LHS was hardired
    original to \code{(incr/100)}.}
  \item{relTol}{originally hard wired to 1e-15, defines the convergence
    tolerance for the search iterations when \code{y >= yLarge}; the
    iterations stop when (new) \code{incr <= y * relTol}.}
  \item{pfEps.n, pfEps.L}{positive factors defining \dQuote{fuzz to ensure
      left continuity}, both originally hardwired to \code{64}.
    originally, the fuzz adjustment was
    \preformatted{p <- p * (1 - 64 *.Machine$double.eps)}
    Now,  \code{pfEps.L} is used \code{if(log.p)} is true and \code{pfEps.n}
    is used otherwise ("n"ormal case), and the adjustments also depend on
    \code{lower.tail}, and also on \code{fpf} : }
  \item{fpf}{a number larger than \code{1}, together with \code{pfEps.n}
    determines the fuzz-adjustment to \code{p} in the case
    \code{(lower=tail=FALSE, log.p=FALSE)}:
    with \code{e <- pfEps.n * .Machine$double.eps}, the adjustment
    \code{p <- p * (1 + e)} is made \emph{iff} \code{1 - p > fpf*e}.}
  \item{trace}{logical (or integer) specifying if (and how much) output
    should be produced from the algorithm.}
}
\details{
  The defaults and exact meaning of the algorithmic tuning arguments from
  \code{yLarge} to \code{fpf} were experimentally determined are subject to change.
}
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
## future: if(getRversion() >= "4.2") stopifnot(x == qPp) # R-devel
qpRp <- qpoisR(Pp, lambda = 100, lower.tail=FALSE)
all.equal(x, qpRp, tol = 0)
stopifnot(all.equal(x, qpRp, tol = 1e-15))
}
\keyword{distribution}
