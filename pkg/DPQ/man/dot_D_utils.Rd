\name{dot-D-utils}
\title{Distribution Utilities "dpq"}
\alias{.D_0}
\alias{.D_1}

\alias{.D_Clog}
\alias{.D_Cval}
\alias{.D_exp}

\alias{.D_LExp}
\alias{.D_log}
\alias{.D_Lval}
\alias{.D_qIv}
\alias{.D_val}

\alias{.DT_0}
\alias{.DT_1}

\alias{.DT_Cexp}
\alias{.DT_CIv}
\alias{.DT_Clog}
\alias{.DT_Cval}
\alias{.DT_exp}
\alias{.DT_log}
\alias{.DT_Log}
\alias{.DT_qIv}
\alias{.DT_val}
%\alias{.D_fexp}
%\alias{.D_rtxp}
\description{
  Utility functions for "dpq"-computations, parelling those in R's own
  C source \file{<Rsource>/src/nmath/dpq.h},
  (\dQuote{dpq} := \bold{d}ensity--\bold{p}robability--\bold{q}uantile).
}
\usage{
.D_0(log.p) # prob/density == 0   (for log.p=FALSE)
.D_1(log.p) # prob         == 1     "       "

.DT_0(lower.tail, log.p) # == 0  when (lower.tail=TRUE, log.p=FALSE)
.DT_1(lower.tail, log.p) # == 1  when     "                "

.D_Lval(p, lower.tail) # p    {L}ower
.D_Cval(p, lower.tail) # 1-p  {C}omplementary

.D_val (x, log.p)  #     x  in pF(x,..)
.D_qIv (p, log.p)  #     p  in qF(p,..)
.D_exp (x, log.p)  # exp(x)        unless log.p where it's  x
.D_log (p, log.p)  #     p           "      "     "    "   log(p)
.D_Clog(p, log.p)  #   1-p           "      "     "    "   log(1-p) == log1p(-)

.D_LExp(x, log.p)  ## [log](1 - exp(-x))     == log1p(- .D_qIv(x))) even more stable

.DT_val (x, lower.tail, log.p) # := .D_val(.D_Lval(x, lower.tail), log.p) ==    x  in pF
.DT_Cval(x, lower.tail, log.p) # := .D_val(.D_Cval(x, lower.tail), log.p) ==  1-x  in pF

.DT_qIv (p, lower.tail, log.p) # := .D_Lval(.D_qIv(p))	==    p	 in qF
.DT_CIv (p, lower.tail, log.p) # := .D_Cval(.D_qIv(p))	==  1-p  in qF

.DT_exp (x, lower.tail, log.p) #  exp( x )
.DT_Cexp(x, lower.tail, log.p) #  exp(1-x)

.DT_log (p, lower.tail, log.p) # log ( p )  in qF
.DT_Clog(p, lower.tail, log.p) # log (1-p)  in qF
.DT_Log (p, lower.tail)        # log ( p )  in qF(p,..,log.p=TRUE)

% currently unused and not exported:
% .D_fexp(f, x, log.p)
% .D_rtxp(rf, x, log.p)%
}
\arguments{
  \item{x}{numeric vector.}
  \item{p}{(log) probability--like numeric vector.}
  \item{lower.tail}{logical; if true, probabilities are \eqn{P[X \le x]},
    otherwise upper tail probabilities, \eqn{P[X > x]}.}
  \item{log.p}{logical; if true, probabilities \eqn{p} are given as
    \eqn{\log(p)} in argument \code{p}.}
}
%% \details{

%% }
\value{
  Typically a numeric vector \dQuote{as} \code{x} or \code{p}, respectively.
}
%% \references{
%% }
\author{Martin Maechler}

\seealso{
  \code{\link{log1mexp}()} which is called from \code{.D_LExp()} and \code{.DT_Log()}.
}
\examples{
FT <- c(FALSE, TRUE)
stopifnot(exprs = {
    .D_0(log.p = FALSE) ==    (0)
    .D_0(log.p = TRUE ) == log(0)
    identical(c(1,0), vapply(FT, .D_1, double(1)))
})

## all such functions in package DPQ:
eDPQ <- as.environment("package:DPQ")
ls.str(envir=eDPQ, pattern = "^[.]D", all.names=TRUE)
(nD <- local({ n <- names(eDPQ); n[startsWith(n, ".D")] }))
trimW <- function(ch) sub(" +$","", sub("^ +","", ch))
writeLines(vapply(sort(nD), function(nm) {
    B <- deparse(eDPQ[[nm]])
    sprintf("\%31s := \%s", trimW(sub("function ", nm, B[[1]])),
            paste(trimW(B[-1]), collapse=" "))
                  }, ""))

do.lowlog <- function(Fn, ...) {
    stopifnot(is.function(Fn),
              all(c("lower.tail", "log.p") \%in\% names(formals(Fn)))) %
    FT <- c(FALSE, TRUE) ; cFT <- c("F", "T")
    L <- lapply(FT, function(lo) sapply(FT, function(lg) Fn(..., lower.tail=lo, log.p=lg)))
    r <- simplify2array(L)
    `dimnames<-`(r, c(rep(list(NULL), length(dim(r)) - 2L),
                      list(log.p = cFT, lower.tail = cFT)))
}
do.lowlog(.DT_0)
do.lowlog(.DT_1)
do.lowlog(.DT_exp, x = 1/4) ; do.lowlog(.DT_exp, x = 3/4)
do.lowlog(.DT_val, x = 1/4) ; do.lowlog(.DT_val, x = 3/4)
do.lowlog(.DT_Cexp, x = 1/4) ; do.lowlog(.DT_Cexp, x = 3/4)
do.lowlog(.DT_Cval, x = 1/4) ; do.lowlog(.DT_Cval, x = 3/4)
do.lowlog(.DT_Clog, p = (1:3)/4) # w/ warn
do.lowlog(.DT_log,  p = (1:3)/4) # w/ warn
do.lowlog(.DT_qIv,  p = (1:3)/4)

## unfinished: FIXME, the above is *not* really checking
stopifnot(exprs = {

})

}
\keyword{distribution}
\keyword{utilities}
