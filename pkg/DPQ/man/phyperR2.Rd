\name{phyperR2}
\title{Pure R version of R's C level phyper()}
\alias{phyperR2}
\alias{pdhyper}
\description{
  Use pure \R functions to compute (less efficiently and usually even less
  accurately) hypergeometric (point) probabilities with the same
  "Welinder"-algorithm as \R's C level code has been doing since 2004.
}
\usage{
phyperR2(q, m, n, k, lower.tail = TRUE, log.p = FALSE, \dots)
pdhyper (q, m, n, k, log.p = FALSE,
         epsC = .Machine$double.eps, verbose = getOption("verbose"))
}
\arguments{
  \item{q}{vector of quantiles representing the number of white balls
    drawn without replacement from an urn which contains both black and
    white balls.
  }
  \item{m}{the number of white balls in the urn.}
  \item{n}{the number of black balls in the urn.}
  \item{k}{the number of balls drawn from the urn, hence must be in \eqn{0,1,\dots, m+n}.}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
  \item{log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{\dots}{further arguments, passed to \code{pdhyper()}.}
  \item{epsC}{a non-negative number, the computer epsilon to be used;
    effectively a relative convergence tolerance for the \code{while()}
    loop in \code{pdhyper()}.}
  \item{verbose}{logical indicating if the \code{pdhyper()} calls,
    typically one per \code{phyperR2()} call, should show how many terms
    have been computed and summed up.}
}
\note{
  For now, all arguments of these functions must be of length \bold{one}.
}
%% \details{

%% }
\value{
  a number, typically very close to the corresponding
  \code{phyper(q, m,n,k, ..)}.
}
\references{
  E-mail / bug report from Morten Welinder to R-devel / R-bugs
 \cr Subject: [Rd] phyper accuracy and efficiency (PR#6772)
 \cr Date: Thu, 15 Apr 2004 18:06:37 (CEST)
 ......

 The current version has very serious cancellation issues.  For example,
 if you ask for a small right-tail you are likely to get total cancellation.
 For example,  \code{phyper(59, 150, 150, 60, FALSE, FALSE)} gives 6.372680161e-14.
 The right answer is \code{dhyper(0, 150, 150, 60, FALSE)} which is 5.111204798e-22.

 \code{phyper} is also really slow for large arguments.

 Therefore, I suggest using the code below. This is a sniplet from Gnumeric ...
 The code isn't perfect.  In fact, if  x*(NR+NB)  is close to	n*NR,
 then this code can take a while.  Not longer than the old code, though.

 Thanks to Ian Smith for ideas.
}
\author{Martin Maechler, based on \R's C code originally provided by Morton
  Welinder, from the Gnumeric project.}
\seealso{
  \code{\link{phyper}}
}
\examples{
## same example as phyper()
m <- 10; n <- 7; k <- 8
vapply(0:9, phyperR2, 0.1, m=m, n=n, k=k)  ==  phyper(0:9, m,n,k)
## see *all* TRUE  (64b FC30)

## to see number of steps:
vapply(0:9, phyperR2, 0.1, m=m, n=n, k=k,  verbose=TRUE)


## Larger arguments:
k <- 100 ; x <- .suppHyper(k,k,k)
ph  <- phyper(x, k,k,k)
ph2 <- vapply(x, phyperR2, 0.1, m=k, n=k, k=k)
cbind(x, ph, ph2, rE = 1-ph2/ph)
stopifnot(abs(1 -ph2/ph) < 8e-16) # 64bit FC30: see -2.22e-16 <= rE <= 3.33e-16
}

\keyword{ distribution }
