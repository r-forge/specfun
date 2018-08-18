\name{DPQmpfr-package}
\alias{DPQmpfr-package}
\alias{DPQmpfr}
\docType{package}
\title{
\packageTitle{DPQmpfr}
}
\description{
\packageDescription{DPQmpfr}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{DPQmpfr}
\packageIndices{DPQmpfr}
%% ~~ An overview of how to use the package, including the most important ~~
%% ~~ functions ~~
}
\author{
\packageAuthor{DPQmpfr}

Maintainer: \packageMaintainer{DPQmpfr}
}
%% \references{
%% ~~ Literature or other references for background information ~~
%% }
\seealso{
  Packages \code{\link[DPQ:DPQ-package]{DPQ}},
           \code{\link[DPQ:DPQ-package]{Rmpfr}} are both used by this package.
}
\examples{
## An example how mpfr-numbers  "just work" with reasonable R functions:
.srch <- search() ; doAtt <- is.na(match("Rmpfr:package", .srch))
if(doAtt) require(Rmpfr)
nu.s <- 2^seqMpfr(mpfr(-30, 64), mpfr(100, 64), by = 1/mpfr(4, 64))
b0 <- DPQ::b_chi(nu.s)
b1 <- DPQ::b_chi(nu.s, one.minus=TRUE)
stopifnot(inherits(b0,"mpfr"), inherits(b1, "mpfr"),
          b0+b1 == 1,  diff(log(b1)) < 0)
plot(nu.s,          log(b1),  type="l", log="x")
plot(nu.s[-1], diff(log(b1)), type="l", log="x")
if(doAtt) # detach the package(s) we've attached above
  for(pkg in setdiff(search(), .srch)) detach(pkg, character.only=TRUE)
}
\keyword{ distribution }
\keyword{ math }
\keyword{ package }
