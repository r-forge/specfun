\name{pt_Witkovsky_Tab1}
\alias{pt_Witkovsky_Tab1}
\docType{data}
\title{Viktor Witosky's Table_1  pt() Examples}
\description{
  A data frame with 17  \code{pt()} examples from Witosky (2013)'s \sQuote{Table 1}.
  We provide the results for the FOSS Softwares
}
\usage{
data(pt_Witkovsky_Tab1)
}
\format{
  A data frame with 17 observations on the following 7 \code{\link{numeric}} variables.
  \describe{
    \item{\code{x}}{the abscissa, called \code{q} in \code{\link{pt}()}.}
    \item{\code{nu}}{the \emph{positive} degrees of freedom, called \code{df} in \code{\link{pt}()}.}
    \item{\code{delta}}{the noncentrality parameter, called \code{ncp} in \code{\link{pt}()}.}

    \item{\code{true_pnt}}{\dQuote{true} values (computed via higher precision, see Witkovsky(2013)).}
    \item{\code{NCTCDFVW}}{the \code{pt()} values computed with Witkovsky's
      matlab implementation.  Confirmed by using \command{octave} (on Fedora 40 Linux).
      These correspond to our \R (package \CRANpkg{DPQ}) \code{\link{pntVW13}()} values.}
    \item{\code{Boost}}{computed via the Boost C++ library; reported by Witkovsky.}
    \item{\code{R_3.3.0}}{computed by R version 3.3.0; confirmed to be
      identical using R 4.4.1}
  }
}
%% \details{

%% }
\source{% ../data/pt_Witkovsky_Tab1.tab --> ../Misc/pnt-Gil_etal-2023/Witkovsky_2013/
  The table was extracted (by MM) from the result of
  \command{pdftotext --layout <*>.pdf} from the publication.

  Note that the \file{arXiv} pre-publication has very slightly differing
  numbers in its \verb{R} column, e.g., first entry ending in \code{00200}
  instead of  \code{00111}
}
\references{
  \enc{Witkovský}{Witkovsky}, Viktor (2013)
  A Note on Computing Extreme Tail Probabilities of the Noncentral T Distribution with Large
  Noncentrality Parameter,
  \emph{Acta Universitatis Palackianae Olomucensis, Facultas Rerum Naturalium, Mathematica}
  %% GST23 have  \emph{Acta Univ Palack Olomuc Fac Rerum Natur Math.}
  \bold{52}(2), 131--143. % on arXiv --> ../Misc/pnt-Gil_etal-2023/Witkovsky_2013/arXiv-1306.5294v2/
}
\examples{
data(pt_Witkovsky_Tab1)
stopifnot(is.data.frame(d.W <- pt_Witkovsky_Tab1), # shorter
          nrow(d.W) >= 17)
mW <- as.matrix(d.W); row.names(mW) <- NULL # more efficient
colnames(mW)[1:3] #  "x" "nu" "delta"
## use 'R pt() - compatible' names:
(n3 <- names(formals(pt)[1:3]))# "q" "df" "ncp"
colnames(mW)[1:3] <- n3
ptR <- apply(mW[, 1:3], 1, \(a3) unname(do.call(pt, as.list(a3))))
cNm <- paste0("R_", with(R.version, paste(major, minor, sep=".")))
mW <- cbind(mW, `colnames<-`(cbind(ptR), cNm),
            relErr = sfsmisc::relErrV(mW[,"true_pnt"], ptR))
mW
## is current R better than R 3.3.0?  -- or even "the same"?
all.equal(ptR, mW[,"R_3.3.0"])              # still true in R 4.4.1
all.equal(ptR, mW[,"R_3.3.0"], tol = 1e-14) # (ditto)
table(ptR == mW[,"R_3.3.0"]) # {see only 4 (out of 17) *exactly* equal ??}
% 2^53* (d_R330 <- sfsmisc::relErrV(ptR, mW[,"R_3.3.0"])) %
}
\keyword{datasets}
