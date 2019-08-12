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
qbetaAppr.1(a, p, q, y = qnormUappr(a, lower.tail=FALSE))
qbetaAppr.2(a, p, q, lower.tail=TRUE, log.p=FALSE, logbeta = lbeta(p,q))
qbetaAppr.3(a, p, q, lower.tail=TRUE, log.p=FALSE, logbeta = lbeta(p,q))
qbetaAppr.4(a, p, q, y = qnormUappr(a, lower.tail=FALSE),
            verbose = getOption("verbose"))

qbetaAppr  (a, p, q, y = qnormUappr(a, lower.tail=FALSE), logbeta= lbeta(p,q),
            verbose = getOption("verbose") && length(a) == 1)

qbeta.R    (alpha, p, q,
            lower.tail = TRUE, log.p = FALSE,
	    logbeta = lbeta(p,q),
	    low.bnd = 3e-308, up.bnd = 1-2.22e-16,
            method = c("AS109", "Newton-log"),
	    f.acu = function(a,p,q) max(1e-300, 10^(-13- 2.5/pp^2 - .5/a^2)),
	    fpu = .Machine$ double.xmin,
	    qnormU.fun = function(u, lu) qnormUappr(p=u, lp=lu, lower.tail=FALSE),
            R.pre.2014 = FALSE,
	    verbose = getOption("verbose"))
}

%% Undocumented arguments in documentation object 'qbetaAppr'
%%   ‘y’ ‘logbeta’ ‘low.bnd’ ‘up.bnd’ ‘method’ ‘f.acu’ ‘fpu’ ‘qnormU.fun’
%%   ‘R.pre.2014’
%% Documented arguments not in \usage in documentation object 'qbetaAppr':
%%   ‘tol’ ‘maxit’

\arguments{%--> ~/R/D/r-devel/R/src/library/stats/man/Beta.Rd
  %% ===>  qbeta(p, shape1, shape2, ncp = 0, lower.tail = TRUE, log.p = FALSE)
  %% but we use the "old" naming scheme here :
  \item{a, alpha}{vector of probabilities (otherwise, e.g., in
    \code{\link{qbeta}()}, called \code{p}).}
  \item{p, q}{the two shape parameters of the beta distribution; otherwise, e.g., in
    \code{\link{qbeta}()}, called \code{shape1} and \code{shape2}).}
  \item{lower.tail, log.p}{logical, see, e.g., \code{\link{qchisq}()}; must
    have length 1.}
%%% TODO: \item{method}{a string specifying the approximation method to be used.}
%%% -----  INSTEAD of all these functions ??
  \item{tol}{non-negative number, the convergence tolerance}
  %% \item{g}{??? FIXME eliminate ?????? FIXME ?????}
  \item{maxit}{the maximal number of iterations}
  \item{verbose}{logical indicating if the algorithm should produce
    \dQuote{monitoring} information.}
}
%% \details{
%% }
\value{
  ...
}
%% \references{
%%
%% }
\author{Martin Maechler}
\seealso{
  \code{\link{qbeta}}.
}
\examples{
  ## TODO
}
\keyword{distribution}
