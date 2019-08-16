\name{qbetaAppr}
\title{Compute (Approximate) Quantiles of the Beta Distribution}
\alias{qbetaAppr.1}
\alias{qbetaAppr.3}
\alias{qbetaAppr.2}
\alias{qbetaAppr.4}
\alias{qbetaAppr}
\alias{qbeta.R}
\description{
  Compute quantiles (inverse distribution values) for the beta distribution,
  using diverse approximations.
}
\usage{% R code : ../R/beta-fns.R
qbetaAppr.1(a, p, q, y = qnormUappr(a))
qbetaAppr.2(a, p, q, lower.tail=TRUE, log.p=FALSE, logbeta = lbeta(p,q))
qbetaAppr.3(a, p, q, lower.tail=TRUE, log.p=FALSE, logbeta = lbeta(p,q))
qbetaAppr.4(a, p, q, y = qnormUappr(a),
            verbose = getOption("verbose"))

qbetaAppr  (a, p, q, y = qnormUappr(a), logbeta= lbeta(p,q),
            verbose = getOption("verbose") && length(a) == 1)

qbeta.R    (alpha, p, q,
            lower.tail = TRUE, log.p = FALSE,
	    logbeta = lbeta(p,q),
	    low.bnd = 3e-308, up.bnd = 1-2.22e-16,
            method = c("AS109", "Newton-log"),
            tol.outer = 1e-15,
	    f.acu = function(a,p,q) max(1e-300, 10^(-13- 2.5/pp^2 - .5/a^2)),
	    fpu = .Machine$ double.xmin,
	    qnormU.fun = function(u, lu) qnormUappr(p=u, lp=lu)
          , R.pre.2014 = FALSE
	  , verbose = getOption("verbose")
          , non.finite.report = verbose
           )
}
\arguments{%--> ~/R/D/r-devel/R/src/library/stats/man/Beta.Rd
  %% ===>  qbeta(p, shape1, shape2, ncp = 0, lower.tail = TRUE, log.p = FALSE)
  %% but we use the "old" naming scheme here :
  \item{a, alpha}{vector of probabilities (otherwise, e.g., in
    \code{\link{qbeta}()}, called \code{p}).}
  \item{p, q}{the two shape parameters of the beta distribution; otherwise, e.g., in
    \code{\link{qbeta}()}, called \code{shape1} and \code{shape2}.}
  \item{y}{an approximation to \eqn{\Phi^{-1}(1-\alpha)} (aka
    \eqn{z_{1-\alpha}}) where \eqn{\Phi(x)} is the standard normal
    cumulative probability function and \eqn{\Phi{-1}(x)} its inverse,
    i.e., \R's \code{\link{qnorm}(x)}.}
  \item{lower.tail, log.p}{logical, see, e.g., \code{\link{qchisq}()}; must
    have length 1.}
  \item{logbeta}{must be \code{\link{lbeta}(p,q)}; mainly an option to pass
    a value already computed.}
  \item{verbose}{logical or integer indicating if and how much
    \dQuote{monitoring} information should be produced by the algorithm.}
  \item{low.bnd, up.bnd}{lower and upper bounds for ...TODO...}
  \item{method}{a string specifying the approximation method to be used.}
  \item{tol.outer}{the \dQuote{outer loop} convergence tolerance; the
    default \code{1e-15} has been hardwired in \R's \code{\link{qbeta}()}.}
  \item{f.acu}{a \code{\link{function}} with arguments \code{(a,p,q)} ...TODO...}
  \item{fpu}{a very small positive number.}
  \item{qnormU.fun}{a \code{\link{function}} with arguments \code{(u,lu)}
    to compute \dQuote{the same} as \code{\link{qnormUappr}()}, the upper
      standard normal quantile.}
  \item{R.pre.2014}{a \code{\link{logical}} ... TODO ...}
  \item{non.finite.report}{\code{\link{logical}} indicating if during the
    \dQuote{outer loop} refining iterations, if \code{y} becomes non finite
    and the iterations have to stop, it should be reported (before the
    current best value is returned).}
}
%% \details{
%% }
\value{
  ...
}
%% \references{
%%
%% }
\author{The R Core Team for the C version in \R's sources;
  Martin Maechler for the \R port.}
\seealso{
  \code{\link{qbeta}}.
}
\examples{
 qbeta.R(0.6, 2, 3) # 0.4445
 qbeta.R(0.6, 2, 3) - qbeta(0.6, 2,3) # almost 0

 qbetaRV <- Vectorize(qbeta.R, "alpha") # now can use
 curve(qbetaRV(x, 1.5, 2.5))
 curve(qbeta  (x, 1.5, 2.5), add=TRUE, lwd = 3, col = adjustcolor("red", 1/2))

 ## an example of disagreement (and doubt, as borderline, close to underflow):
 qbeta.R(0.5078, .01, 5) # ->  2.77558e-15    # but
 qbeta  (0.5078, .01, 5) # -> 1.776357e-15  now gives 4.651188e-31 !!!
 qbeta  (0.5078, .01, 5, ncp=0) #      also     gives 4.651188e-31
}
\keyword{distribution}
