\name{qbBaha2017}
\alias{qbBaha2017}
\docType{data}
\title{Accurate qbeta() values from Baharev et al (2017)'s Program}
\description{
  Compuate "accurate" \code{\link{qbeta}()} values from
  Baharev et al (2017)'s Program.
}
\usage{
data("qbBaha2017")
}
\format{
  FIXME: Their published table only shows 6 digits,
         but running their (32-bit statically linked) Linux executable
         \file{mindiffver} (from their github repos, see "source")
	 with their own \file{input.txt} gives \emph{12}
         digits accuracy, which we should be able to increase even more,
	 %% a  b  x_0  lambda_0  alpha  beta  eps_x  eps_lambda
	 %%  ~/R/MM/NUMERICS/dpq-functions/beta-gamma-etc/mindiffver/README.md
	 %%  ~/R/MM/NUMERICS/dpq-functions/beta-gamma-etc/Baharev_mindiffver_output.txt
	 see \url{https://github.com/baharev/mindiffver/blob/master/README.md}

 A numeric matrix, \eqn{9 \times 22}{9 * 22} with guaranteed accuracy
 \code{\link{qbeta}(0.95, a,b)} values, for
 \eqn{a = 0.5, 1, 1.5, 2, 2.5, 3, 5, 10, 25} and
 \eqn{b = }
 with \code{\link{str}()} \preformatted{
   num [1:9, 1:22] 0.902 0.95 0.966 0.975 0.98 ...
   - attr(*, "dimnames")=List of 2
    ..$ a: chr [1:9] "0.5" "1" "1.5" "2" ...
    ..$ b: chr [1:22] "1" "2" "3" "4" ...}
}
\details{
  MM constructed this data as follows  (TODO: say more..):
  \preformatted{
        ff <- "~/R/MM/NUMERICS/dpq-functions/beta-gamma-etc/Baharev_et_al-2017_table3.txt"
        qbB2017 <- t( data.matrix(read.table(ff)) )
        dimnames(qbB2017) <- dimnames(qbet)
        saveRDS(qbB2017, "..../qbBaha2017.rds")
    }
}

\source{% <<< FIXME: can have much better accuracy numbers !!
  This matrix comprises all entries of Table 3, p. 776 of
  \cr % ~/R/MM/NUMERICS/dpq-functions/beta-gamma-etc/Baharev*
  Baharev, A., Schichl, H. and Rév, E. (2017)
  Computing the noncentral-F distribution and the power of the F-test with
  guaranteed accuracy;
  \emph{Comput. Stat.} \bold{32}(2), 763--779.
  \doi{10.1007/s00180-016-0701-3}

  The paper mentions the first author's \file{github} repos where source code and
  executables are available from:
  \url{https://github.com/baharev/mindiffver/}
}
\examples{
data(qbBaha2017)
str(qbBaha2017)
str(ab <- lapply(dimnames(qbBaha2017), as.numeric))
stopifnot(ab$a == c((1:6)/2, 5, 10, 25),
          ab$b == c(1:15, 10*c(2:5, 10, 25, 50)))
matplot(ab$b, t(qbBaha2017)[,9:1], type="l", log = "x", xlab = "b",
        ylab = "qbeta(.95, a,b)",
        main = "Guaranteed accuracy 95\% percentiles of Beta distribution")
legend("right", paste("a = ", format(ab$a)),
       lty=1:5, col=1:6, bty="n")

## Relative error of R's qbeta() -- given that the table only shows 6
## digits, there is *no* relevant error: R's qbeta() is accurate enough:
x.ab <- do.call(expand.grid, ab)
matplot(ab$b, 1 - t(qbeta(0.95, x.ab$a, x.ab$b) / qbBaha2017),
        main = "rel.error of R's qbeta() -- w/ 6 digits, it is negligible",
        ylab = "1 - qbeta()/'true'",
        type = "l", log="x", xlab="b")
abline(h=0, col=adjustcolor("gray", 1/2))
}
\keyword{datasets}
