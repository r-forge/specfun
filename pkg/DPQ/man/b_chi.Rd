\name{b_chi}
\title{Compute \eqn{E[\chi_\nu] / \sqrt{\nu}}{E[chi_nu]/sqrt(nu)}
  useful for t- and chi-Distributions}
\alias{b_chi}
\alias{b_chiAsymp}
\alias{lb_chi00}
\alias{lb_chi0}
\alias{lb_chiAsymp}
\alias{c_dt}
\alias{c_dtAsymp}
\alias{c_pt}
\description{
  \deqn{b_\chi(\nu) := E[\chi(\nu)] / \sqrt{\nu} = \frac{\sqrt{2/\nu}\Gamma((\nu+1)/2)}{\Gamma(\nu/2)},}{%
    b_chi(nu) := E[chi(nu)] / sqrt(nu) = sqrt(2/nu) Gamma((nu+1)/2) / Gamma(nu/2),}
  where \eqn{\chi(\nu)} denotes a chi-distributed random variable, i.e.,
  the square of a chi-squared variable, and \eqn{\Gamma(z)}{Gamma(z)} is
  the Gamma function, \code{\link{gamma}()} in \R.

  This is a relatively important auxiliary function when computing with
  non-central t distribution functions and approximations, specifically see
  Johnson et al.(1994), p.520, after (31.26a).

  Its logarithm,
  \deqn{lb_\chi(\nu) := log\bigl(\frac{\sqrt{2/\nu}\Gamma((\nu+1)/2)}{\Gamma(\nu/2)}\bigr),}{%
    lb_chi(nu) := log(sqrt(2/nu) Gamma((nu+1)/2) / Gamma(nu/2)),}
  is even easier to compute via \code{\link{lgamma}} and \code{\link{log}},
  and I have used Maple to derive an asymptotic expansion in
  \eqn{\frac{1}{\nu}} as well.

  Note that \eqn{lb_\chi(\nu)}{lb_chi(nu)} also appears in the formula
  for the t-density (\code{\link{dt}}) and distribution (tail) functions.
}
\usage{
b_chi     (nu, one.minus = FALSE, c1 = 341, c2 = 1000)
b_chiAsymp(nu, order = 2, one.minus = FALSE)
#lb_chi    (nu, ......) # not yet
lb_chiAsymp(nu, order)

c_dt(nu)       # warning("FIXME: current c_dt() is poor -- base it on lb_chi(nu) !")
c_dtAsymp(nu)  # deprecated in favour of lb_chi(nu)
c_pt(nu)       # warning("use better c_dt()") %---> FIXME deprecate even stronger ?
}

\arguments{
  \item{nu}{non-negative numeric vector of degrees of freedom.}
  \item{one.minus}{logical indicating if \eqn{1 - b()} should be
    returned instead of \eqn{b()}.}
  \item{c1, c2}{boundaries for different approximation intervals used:
    \cr for \code{0  < nu <= c1}, internal \code{b1()} is used,
    \cr for \code{c1 < nu <= c2}, internal \code{b2()} is used, and
    \cr for \code{c2 < nu}, the \code{b_chiAsymp()} function is used,
    (and you can use that explicitly, also for smaller \code{nu}).

    FIXME: \code{c1} and \code{c2} were defined when the only asymptotic
    expansion known to me was the \code{order = 2} one.
    A future version of \code{b_chi} will \emph{very likely} use
    \code{b_chiAsymp(*, order)} for higher orders, and the c1 and c2
    arguments will change, possibly be abolished.
  }
  \item{order}{the polynomial order in \eqn{\frac{1}{\nu}}{1/nu} of the
    asymptotic expansion of \eqn{b_\chi(\nu)} for \eqn{\nu\to\infty}{nu -> Inf}.

    The default, \code{order = 2} corresponds to the order you can get
    out of the Abramowitz and Stegun (6.1.47) formula.
    Higher order expansions were derived using \command{Maple} by Martin Maechler
    in 2002, see below, but implemented in \code{b_chiAsymp()} only in 2018.}
}

\value{
  a numeric vector of the same length as \code{nu}.
}
\details{
  One can see that \code{b_chi()} has the properties of a CDF of a
  continuous positive random variable:  It grows monotonely from
  \eqn{b_\chi(0) = 0} to (asymptotically) one.  Specifically,
  for large \code{nu}, \code{b_chi(nu) = b_chiAsymp(nu)} and
  \deqn{1 - b_\chi(\nu) \sim \frac{1}{4\nu}.}{1 - b_chi(nu) ~= 1/(4 nu).}

  More accurately, derived from Abramowitz and Stegun, 6.1.47 (p.257)
  for a= 1/2, b=0,
  \deqn{\Gamma(z + 1/2) / \Gamma(z) \sim \sqrt(z)*(1 - 1/(8z) + 1/(128 z^2) + O(1/z^3)),}{%
    \Gamma(z + 1/2) / Gamma(z) ~= sqrt(z)*(1 - 1/(8z) + 1/(128 z^2) + O(1/z^3)),}
  and applied for \eqn{b_\chi(\nu)} with \eqn{z = \nu/2}, we get
  \deqn{b_\chi(\nu) \sim 1 - (1/(4\nu) * (1 - 1/(8\nu)) + O(\nu^{-3})),}
  which has been implemented in \code{b_chiAsymp(*, order=2)} in 1999.

  Even more accurately, Martin Maechler, used Maple to derive an
  asymptotic expansion up to order 15, here reported up to order 5,
  namely with \eqn{r := \frac{1}{4\nu}}{r := 1/(4nu)},
  \deqn{b_\chi(\nu) = c_\chi(r) = 1 - r + \frac{1}{2}r^2 +
    \frac{5}{2}r^3 - \frac{21}{8}r^4 - \frac{399}{8}r^5 + O(r^6).}{%
    b_chi(nu) = c_chi(r) = 1 - r + 1/2 r^2 + 5/2 r^3 - 21/8 r^4 - 399/8 r^5 + O(r^6).}
}
\references{
  Johnson, Kotz,  Balakrishnan (1995)
  \emph{Continuous Univariate Distributions},
  Vol 2, 2nd Edition; Wiley.\cr
  Formula on page 520, after (31.26a)

  Abramowitz, M. and Stegun, I. A. (1972)
  \emph{Handbook of Mathematical Functions}. New York: Dover.
  \url{https://en.wikipedia.org/wiki/Abramowitz_and_Stegun} provides
  links to the full text which is in public domain.
}
\author{Martin Maechler}
\seealso{
  The t-distribution (base \R) page \code{\link{pt}}.
}
\examples{
curve(b_chi, 0, 20); abline(h=0:1, v=0, lty=3)
r <- curve(b_chi, 1e-10, 1e5, log="x")
with(r, lines(x, b_chi(x, one.minus=TRUE), col = 2))

## Zoom in to c1-region
rc1 <- curve(b_chi, 340.5, 341.5, n=1001)# nothing to see
e <- 1e-3; curve(b_chi, 341-e, 341+e, n=1001) # nothing
e <- 1e-5; curve(b_chi, 341-e, 341+e, n=1001) # see noise, but no jump
e <- 1e-7; curve(b_chi, 341-e, 341+e, n=1001) # see float "granularity"+"jump"

## Zoom in to c2-region
rc2 <- curve(b_chi, 999.5, 1001.5, n=1001) # nothing visible
e <- 1e-3; curve(b_chi, 1000-e, 1000+e, n=1001) # clear small jump
c2 <- 1500
e <- 1e-3; curve(b_chi(x,c2=c2), c2-e, c2+e, n=1001)# still
## - - - -
c2 <- 3000
e <- 1e-3; curve(b_chi(x,c2=c2), c2-e, c2+e, n=1001)# ok asymp clearly better!!
curve(b_chiAsymp, add=TRUE, col=adjustcolor("red", 1/3), lwd=3)
if(requireNamespace("Rmpfr")) {
 xm <- Rmpfr::seqMpfr(c2-e, c2+e, length.out=1000)

}
## - - - -
c2 <- 4000
e <- 1e-3; curve(b_chi(x,c2=c2), c2-e, c2+e, n=1001)# ok asymp clearly better!!
curve(b_chiAsymp, add=TRUE, col=adjustcolor("red", 1/3), lwd=3)

grCol <- adjustcolor("forest green", 1/2)
curve(b_chi,                    1/2, 1e11, log="x")
curve(b_chiAsymp, add = TRUE, col = grCol, lwd = 3)
## 1-b(nu) ~= 1/(4 nu) a power function <==> linear in log-log scale:
curve(b_chi(x, one.minus=TRUE), 1/2, 1e11, log="xy")
curve(b_chiAsymp(x, one.minus=TRUE), add = TRUE, col = grCol, lwd = 3)

%% FIXME: and see an example in package 'DPQmpfr'
}
\keyword{math}
\keyword{distribution}
