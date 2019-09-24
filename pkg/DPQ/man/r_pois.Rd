\name{r_pois}
\alias{r_pois}
\alias{r_pois_expr}
\alias{plRpois}% ==>> ../R/pnchisq.R <<==
\title{Compute Relative Size of i-th term of Poisson Distribution Series}
\description{
  Compute \deqn{r_\lambda(i) := (\lambda^i / i!) / e_{i-1}(\lambda),}{r_v(i) := (v^i / i!) / e_{i-1}(v),}
  where \eqn{\lambda =}{v =}\code{lambda},  and
  \deqn{e_n(x) := 1 + x + x^2/2! + .... + x^n/n! } is the \eqn{n}-th
  partial sum of \eqn{\exp(x) = e^x}.

  Questions: As function of \eqn{i} \itemize{
    \item Can this be put in a simple formula, or at least be well
    approximated for large \eqn{\lambda}{v} and/or large \eqn{i}?
    \item For which \eqn{i} (\eqn{ := i_m(\lambda)}{ := i_m(v)}) is it maximal?
    \item When does \eqn{r_{\lambda}(i)} become smaller than (f+2i-x)/x = a + b*i ?
    %% FIXME Q 3: what is 'f' /  'a', 'b' ?????
  }
  NB: This is relevant in computations for non-central chi-squared (and
  similar non-central distribution functions) defined as weighted sum with
 \dQuote{Poisson weights}.
}
\usage{
r_pois(i, lambda)
r_pois_expr  # the R expression() for the asymptotic branch of r_pois()

plRpois(lambda, iset = 1:(2*lambda), do.main = TRUE,
        log = 'xy', type = "o", cex = 0.4, col = c("red","blue"),
        do.eaxis = TRUE, sub10 = "10")
}
\arguments{
  \item{i}{integer ..}
  \item{lambda}{non-negative number ... }
  \item{iset}{ ..... }
  \item{do.main}{\code{\link{logical}} specifying if a main
    \code{\link{title}} should be drawn via (\code{main = r_pois_expr}).}
  \item{type}{type of (line) plot, see \code{\link{lines}}.}
  \item{log}{string specifying if (and where) logarithmic scales should be
    used, see \code{\link{plot.default}()}.}
  \item{cex}{character expansion factor.}
  \item{col}{colors for the two curves.}
  \item{do.eaxis}{\code{\link{logical}} specifying if
        \code{\link[sfsmisc]{eaxis}()} (package \CRANpkg{sfsmisc}) should
        be used.}
  \item{sub10}{argument for \code{\link[sfsmisc]{eaxis}()} (with a
        different default than the original).}
}
\details{
  \code{r_pois()} is related to our series expansions and approximations
  for the non-central chisquared;
  in particular
  ........... % TODO


  \code{plRpois()} simply produces a \dQuote{nice} plot of \code{r_pois(ii, *)}
  vs \code{ii}.
}
\value{
  \describe{
    \item{\code{r_pois()}}{returns a numeric vector
      \eqn{r_\lambda(i)}{r_v(i)} values.}
    \item{\code{r_pois_expr()}}{an \code{\link{expression}}.}
  }
}
%% \references{
%%   TODO
%% }
\author{Martin Maechler, 20 Jan 2004}
%% \note{

%% }
\seealso{
  \code{\link{dpois}()}.
}
\examples{%%>> more [incl other plots] in ==>> ../tests/chisq-nonc-ex.R <<==
plRpois(12)
plRpois(120)
}% ex.
\keyword{distribution}
