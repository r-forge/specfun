\name{lsum}
\alias{lsum}
\title{Properly Compute the Logarithm of a Sum}
\usage{
lsum(lx, l.off = max(lx))
}
\arguments{
  \item{lx}{n-vector of values log(x_1),..,log(x_n).}
  \item{l.off}{the offset to substract and re-add; ideally in
    the order of the maximum of each column.}
}
\value{
  \deqn{
    log(x_1 + .. + x_n) = log(sum(x)) = log(sum(exp(log(x)))) =
    = log(exp(log(x_max))*sum(exp(log(x)-log(x_max)))) =
    = log(x_max) + log(sum(exp(log(x)-log(x_max))))) =
    = lx.max + log(sum(exp(lx-lx.max)))
  }
}
\description{
  Simple vector version of \code{copula:::lsum()} (CRAN package
  \CRANpkg{copula}). % ~/R/Pkgs/copula/R/special-func.R

  Properly compute \eqn{\log(x_1 + \ldots + x_n)}{log(x1 + .. + xn)}.
  for given \eqn{log(x_1),..,log(x_n)}.  Here, \eqn{x_i > 0}{xi > 0} for all \eqn{i}.
}
\seealso{
  \code{\link{lssum}()} which computes a sum in log scale
  with specified (typically alternating) signs.
}
\author{
  Originally, via paired programming: Marius Hofert and Martin Maechler.
}
\examples{
## The "naive" version :
lsum0 <- function(lx) log(sum(exp(lx)))

lx1 <- 10*(-80:70) # is easy
lx2 <- 600:750     # lsum0() not ok [could work with rescaling]
lx3 <- -(750:900)  # lsum0() = -Inf - not good enough
m3 <- cbind(lx1,lx2,lx3)
lx6 <- lx5 <- lx4 <- lx3
lx4[149:151] <- -Inf ## = log(0)
lx5[150] <- Inf
lx6[1] <- NA_real_
m6 <- cbind(m3,lx4,lx5,lx6)
stopifnot(exprs = {
  all.equal(lsum(lx1), lsum0(lx1))
  all.equal((ls1 <- lsum(lx1)),  700.000045400960403, tol=8e-16)
  all.equal((ls2 <- lsum(lx2)),  750.458675145387133, tol=8e-16)
  all.equal((ls3 <- lsum(lx3)), -749.541324854612867, tol=8e-16)
  ## identical: matrix-version <==> vector versions
  identical(lsum(lx4), ls3)
  identical(lsum(lx4), lsum(head(lx4, -3))) # the last three were -Inf
  identical(lsum(lx5), Inf)
  identical(lsum(lx6), lx6[1])
  identical((lm3 <- apply(m3, 2, lsum)), c(lx1=ls1, lx2=ls2, lx3=ls3))
  identical(apply(m6, 2, lsum), c(lm3, lx4=ls3, lx5=Inf, lx6=lx6[1]))
})
}
