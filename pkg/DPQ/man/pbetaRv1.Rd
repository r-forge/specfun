\name{pbetaRv1}
\alias{pbetaRv1}
\title{Pure R Implementation of Old pbeta()}
\description{
  \code{pbetaRv1()} is an implementation of the original
  (\dQuote{version 1} \code{\link{pbeta}()} function in \R (versions <=
  2.2.x), before we started using TOMS 708 \command{bratio()} instead, see
  that help page also for references.

  \code{pbetaRv1()} is basically a manual translation from C to \R of the
  underlying \code{pbeta_raw()} C function, see in \R's source tree at
  \url{https://svn.r-project.org/R/branches/R-2-2-patches/src/nmath/pbeta.c}

  For consistency within \R, we are using \R's argument names
  \code{(q, shape1, shape2)} instead of C code's
  \code{(x, pin,    qin   )}.

  It is only for the \emph{central} beta distribution.
}
\usage{
pbetaRv1(q, shape1, shape2, lower.tail = TRUE,
         eps = 0.5 * .Machine$double.eps,
         sml = .Machine$double.xmin,
         verbose = 0)
}
\arguments{
  \item{q, shape1, shape2}{non-negative numbers, q in \eqn{[0,1]}, see
    \code{\link{pbeta}}.}
  \item{lower.tail}{indicating if \eqn{F(q; *)} should be returned or
    the upper tail probability \eqn{1 - F(q)}.}
  \item{eps}{the tolerance used to determine congerence.  \code{eps} has
    been hard coded in C code to \code{0.5 * .Machine$double.eps} which is
    equal to \eqn{2^{-53}} or \code{1.110223e-16}.}
  \item{sml}{the smallest positive number on the typical platform.  The
    default \code{.Machine$double.xmin} is hard coded in the C code (as
    \code{DBL_MIN}), and this is equal to \eqn{2^{-1022}} or
    \code{2.225074e-308} on all current platforms.}
  \item{verbose}{integer indicating the amount of verbosity of
    diagnostic output, \code{0} means no output, \code{1} more, etc.}
}
%% \details{
%% %%  ~~ If necessary, more details than the description above ~~
%% }
\value{
  a number.
}
\note{The C code contains\cr \emph{
 This routine is a translation into C of a Fortran subroutine
 by W. Fullerton of Los Alamos Scientific Laboratory.}
}
\references{
  (From the C code:)

  Nancy E. Bosten and E.L. Battiste (1974).
  Remark on Algorithm 179 (S14): Incomplete Beta Ratio.
  \emph{Communications of the ACM}, \bold{17}(3), 156--7.% ACM := Association for Computing Machinery
  %% MM : ~/save/papers/Numerics/
}
\author{Martin Maechler}
%% \note{
%% }
\seealso{
  \code{\link{pbeta}}.
}
\examples{
all.equal(pbetaRv1(1/4, 2, 3),
          pbeta   (1/4, 2, 3))
set.seed(101)

N <- 1000
x <- sample.int(7, N, replace=TRUE) / 8
a <-   rlnorm(N)
b <- 5*rlnorm(N)
pbt <- pbeta(x, a, b)
for(i in 1:N) {
   stopifnot(all.equal(pbetaRv1(x[i], a[i], b[i]), pbt[i]))
   cat(".", if(i \%\% 20 == 0) paste0(i, "\n"))
}

}
\keyword{math}
\keyword{distribution}
