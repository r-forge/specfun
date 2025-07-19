\name{ss}
\alias{ss}
\alias{pnchisq}
\title{Noncentral Chisquare Terms \code{"mpfr"}-ized}
\description{
  Compute the terms of the series of the noncentral chisquare distribution
  function, using \code{"mpfr"}-arithmetic unless all three \code{x, df, ncp} are
  numeric, in which case package \CRANpkg{DPQ}'s \code{\link[DPQ]{ss}()} is called.
}
\usage{
ss(x, df, ncp, i.max = 10000, minPrec = 128L,
   useLv = !(expMin < -lambda && 1/lambda < expMax))

pnchisq(x, df, ncp = 0, lower.tail = TRUE, log.p = FALSE, i.max = 10000,
        ssr = ss(x = x, df = df, ncp = ncp, i.max = i.max))
}
\arguments{
  \item{x}{numeric vector (of \sQuote{quantiles}, i.e., abscissa values).}
  \item{df}{degrees of freedom \eqn{> 0}, maybe non-integer.}
  \item{ncp}{non-centrality parameter \eqn{\delta}{delta}; ....}
  \item{lower.tail, log.p}{\code{\link{logical}}s, ..}
  \item{i.max}{maximal number of terms to compute and return.}
  \item{minPrec}{minimal \code{prec} to be used in coercions to
    \code{"\link[Rmpfr]{mpfr}"} of those arguments among \code{x, df, ncp}
    which inherit not from \code{"mpfr"} already.}
  \item{useLv}{\code{\link{logical}} indicating if logarithmic scale should
    be used for \eqn{\lambda} computations.}
  \item{ssr}{the result of calling \code{ss(..)}.}
}
\details{
  Using mpfr-arithmetic (provided by package \CRANpkg{Rmpfr}) where
  underflow to zero is very rare, typically all \code{i.max} terms are
  non-zero and hence returned.
}
\value{
  a \code{\link{list}} with components
  \item{s}{the (non-zero) terms of the series, a (numeric or \code{"mpfr"})
    vector of length \eqn{\le}{<=} \code{i.max}.}
  \item{i1}{location (in \code{s[]}) of the first change from 0 to
    positive, typically \code{i1 == 1}, currently}
  \item{max}{(first) location of the maximal value in the series (i.e.,
    \code{\link{which.max}(s)}).}
}
%% \references{

%% }
\author{Martin Maechler}
%% \note{
%% }

\seealso{
  \CRANpkg{DPQ}'s \code{\link[DPQ]{ss}()} % both in ../../DPQ/man/pnchisqAppr.Rd
  and \code{\link[DPQ]{pnchisq_ss}()}.
}
\examples{
x <- 2000 # try both
x <- 640
xm <- Rmpfr::mpfr(x, 128)
ssm <- ss(xm, df=512, ncp=128) # currently returns all i.max = 10'000 terms ..
str(ssm)
Rmpfr::roundMpfr(head(ssm $ s), 32) # 1.6038108905e-28 4.0563310465e-26 5.1109287414e-24 ...
Rmpfr::roundMpfr(tail(ssm $ s), 32) # ... 1.2146182008e-6176 1.1844156032e-6177
plot(Rmpfr::asNumeric(log10(ssm$s)), type = "l")
with(ssm, lines(max, log10(s[max]), type = "h", col=2)); abline(h=0, v=0, lty=3)
str(Dss <- DPQ::ss(x, df=512, ncp=128) )
(nD <- length(Dss$s))
stopifnot(Rmpfr::all.equal(ssm$s[1:nD], Dss$s, tolerance = 1e-14))
          Rmpfr::all.equal(ssm$s[1:nD], Dss$s, tol=0) # Mean rel. diff.: 1.0597e-15

x <- 2000 ; xm <- Rmpfr::mpfr(x, 128)
 pn  <- pchisq (x , df=512, ncp=128) # prints as 1
(pnM <- pnchisq(xm, df=512, ncp=128))
stopifnot(pn <= 1, pnM < 1, 1-pnM < 1e-35)
Rmpfr::asNumeric(1- pn/pnM) # 2.5336e-13 currently
}
\keyword{distribution}
\keyword{math}
