\name{gam1M}
\title{Compute  1/Gamma(x+1) - 1  Accurately}
\alias{gam1M}
\description{

   FIXME: "R's own" double prec version is now in package DPQ:
   	       e.g. ~/R/Pkgs/DPQ/man/gam1.Rd

   FIXME2: R-only implementation is in
                   ~/R/Pkgs/DPQ/TODO_R_versions_gam1_etc.R

%%___________ FIXME:  no interval restriction -- this is for  gam1() in DPQ ___________
  Computes \eqn{1/\Gamma(a+1) - 1} accurately in \eqn{[-0.5, 1.5]} for numeric argument \code{a};
  For \code{"mpfr"} numbers, the precision is increased intermediately such
  that \eqn{a+1} should not lose precision.
}
\usage{
gam1M(a, usePr = NULL)
}
\arguments{
  \item{a}{a numeric or numeric-alike, typically inheriting from \code{class} \code{"mpfr"}.}
  \item{usePr}{the precision to use; the default, \code{NULL}, means to use
    a default which depends on \code{a}, specifically \code{\link[Rmpfr]{getPrec}(a)}.}
}
\details{ %% MM: This is all correct --- but 100% irrelevant !!!
  %%  It would be the approximation at a = -1 which is outside [-1/2, 3/2] !!!!
  \url{https://dlmf.nist.gov/} states the well-know Taylor series for
  \deqn{\frac{1}{\Gamma(z)} = \sum_{k=1}^\infty c_k z^k}{%
    1 / Gamma(z) = sum(k=1..Inf; c_k z^k)}
  with \eqn{c_1 = 1},  \eqn{c_2 = \gamma}, (Euler's gamma, \eqn{\gamma = 0.5772...}), with
  recursion \eqn{c_k = (\gamma c_{k-1} - \zeta(2) c_{k-2} ... +(-1)^k \zeta(k-1) c_1) /(k-1)}.

  Hence,
  \deqn{\frac{1}{\Gamma(z+1)}    = z+1 + \sum_{k=2}^\infty c_k (z+1)^k,}{%
      1/Gamma(z+1)     = z+1 +  sum (k= 2..Inf | c_k (z+1)^k),}
  \deqn{\frac{1}{\Gamma(z+1)} -1 = z   + \gamma*(z+1)^2 + \sum_{k=3}^\infty c_k (z+1)^k.}{%
      1/Gamma(z+1) -1  = z   + \gamma*(z+1)^2 +  sum(k = 3..Inf | c_k (z+1)^k).}
  Consequently, for \eqn{\zeta_k := \zeta(k)},
  \eqn{c_3 = (\gamma^2 - \zeta_2)/2},
  \eqn{c_4 = \gamma^3/6 - \gamma \zeta_2/2 + \zeta_3/3}.
\preformatted{
  require(Rmpfr) # Const(), mpfr(), zeta()
  gam <- Const("gamma", 128)
  z <- zeta(mpfr(1:7, 128))
  (c3 <- (gam^2 -z[2])/2)                       # -0.655878071520253881077019515145
  (c4 <- (gam*c3 - z[2]*c2 + z[3])/3)           # -0.04200263503409523552900393488
  (c4 <- gam*(gam^2/6 - z[2]/2) + z[3]/3)
  (c5 <- (gam*c4 - z[2]*c3 + z[3]*c2 - z[4])/4) # 0.1665386113822914895017007951
  (c5 <- (gam^4/6 - gam^2*z[2] + z[2]^2/2 + gam*z[3]*4/3 - z[4])/4)
}

}
\value{
  a numeric-alike vector like \code{a}.
}
\references{
  TOMS 708, see \code{\link{pbeta}}
}
\author{Martin Maechler building on C code of TOMS 708}
%% \note{

%% }
\seealso{
 \code{\link{gamma}}.
}
\examples{
##' naive direct formula:
g1 <- function(u) 1/gamma(u+1) - 1

%% TODO: move this to DPQ (only on help page?) -- once it has gam1() {on "identical" C code}

##' @title gam1() from TOMS 708 -- translated to R (*and* vectorized)
##' @author Martin Maechler
gam1R <- function(a, chk=TRUE) { ##  == 1/gamma(a+1) - 1  -- accurately  ONLY for  -0.5 <= a <= 1.5
    if(!length(a)) return(a)
    ## otherwise:
    if(chk) stopifnot(-0.5 <= a, a <= 1.5) # if not, the computation below is non-sense!
    d  <- a - 0.5
    ## t := if(a > 1/2)  a-1  else  a  ==> t in [-0.5, 0.5]  <==>  |t| <= 0.5
    R <- t <- a
    dP <- d > 0
    t[dP] <- d[dP] - 0.5
    if(any(N <- (t < 0.))) { ## L30: */
        r <- c(-.422784335098468, -.771330383816272,
               -.244757765222226, .118378989872749, 9.30357293360349e-4,
               -.0118290993445146, .00223047661158249, 2.66505979058923e-4,
               -1.32674909766242e-4)
        s1 <- .273076135303957
        s2 <- .0559398236957378
        t_ <- t[N]
        top  <- (((((((r[9] * t_ + r[8]) * t_ + r[7]) * t_ + r[6]) * t_ + r[5]) * t_ + r[4]
        ) * t_ + r[3]) * t_ + r[2]) * t_ + r[1]
        bot <- (s2 * t_ + s1) * t_ + 1.
        w <- top / bot
        ## if (d > 0.) :
        if(length(iP <- which(dP[N])))
            R[N &  dP] <- (t_ * w)[iP] / a[N & dP]
        ## else d <= 0 :
        if(length(iN <- which(!dP[N])))
            R[N & !dP] <- a[N & !dP] * (w[iN] + 0.5 + 0.5)
    }
    if(any(Z <- (t == 0))) ## L10: a in {0, 1}
        R[Z] <- 0.
    if(any(P <- t > 0)) { ## t > 0;  L20: */
        p <- c( .577215664901533, -.409078193005776,
               -.230975380857675, .0597275330452234, .0076696818164949,
               -.00514889771323592, 5.89597428611429e-4 )
        q <- c(1., .427569613095214, .158451672430138, .0261132021441447, .00423244297896961)
        t <- t[P]
        top <- (((((p[7] * t + p[6])*t + p[5])*t + p[4])*t + p[3])*t + p[2])*t + p[1]
        bot <- (((q[5] * t + q[4]) * t + q[3]) * t + q[2]) * t + 1.
        w <- top / bot
        ## if (d > 0.) ## L21: */
        if(length(iP <- which(dP[P])))
            R[P &  dP] <- t[iP] / a[P &  dP] * (w[iP] - 0.5 - 0.5)
        ## else d <= 0 :
        if(length(iN <- which(!dP[P])))
            R[P & !dP] <- a[P & !dP] * w[iN]
    }
    R
} ## gam1R()

u <- seq(-.5, 1.5, by=1/16); set.seed(1); u <- sample(u) # permuted (to check logic)
%% u <- seq(-.5, 1.5, by=1/4) ## for testing
g11   <- vapply(u, gam1R, 1) # [-.5, 1.5]  == the interval for which the above gam1() was made
gam1. <- gam1R(u)
cbind(u, gam1., D = sfsmisc::relErrV(gam1., g1(u)))[order(u),]
                               # looks "too good", as we are not close (but different) to {0, 1}
stopifnot( identical(g11, gam1.) )
           all.equal(g1(u), gam1., tolerance = 0) # 6.7e-16  ("too good", see above)
stopifnot( all.equal(g1(u), gam1.) )

## Comparison using Rmpfr; slightly extending [-.5, 1.5] interval (and getting much closer to {0,1})
u <- seq(-0.525, 1.525, length.out = 2001)
uM <- Rmpfr::mpfr(u, 128)
gam1M. <- gam1M(uM)
relE <- Rmpfr::asNumeric(sfsmisc::relErrV(gam1M., gam1R(u, chk=FALSE)))
rbind(rErr = summary(relE),
     `|rE|` = summary(abs(relE)))
##            Min.    1st Qu.    Median      Mean   3rd Qu.     Max.
## rErr -3.280e-15 -3.466e-16 1.869e-17 1.526e-16 4.282e-16 1.96e-14
## |rE|  1.343e-19  2.363e-16 3.861e-16 6.014e-16 6.372e-16 1.96e-14
stopifnot(max(abs(relE)) < 1e-13)

relEtit <- expression("Relative Error of " ~~ gam1(u) \%~\%{} == frac(1, Gamma(u+1)) - 1) #%
plot(relE ~ u, type="l", ylim = c(-1,1) * 2.5e-15, main = relEtit)
grid(lty = 3); abline(v = c(-.5, 1.5), col = adjustcolor(4, 1/2), lty=2, lwd=2)

## what about the direct formula -- how bad is it really ?
relED <- Rmpfr::asNumeric(sfsmisc::relErrV(gam1M., g1(u)))

plot(relE ~ u, type="l", ylim = c(-1,1) * 1e-14, main = relEtit)
lines(relED ~ u, col = adjustcolor(2, 1/2), lwd = 2); abline(v = (-1:3)/2, lty=2, col="orange3")
mtext("comparing with direct formula   1/gamma(u+1) - 1", col=2, cex=3/4)
legend("top", c("gam1R(u)", "1/gamma(u+1) - 1"), col = 1:2, lwd=1:2, bty="n")
## direct is clearly *worse* , but not catastrophical
}
\keyword{math}
