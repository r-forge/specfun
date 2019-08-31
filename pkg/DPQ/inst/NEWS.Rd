% Check from R:
%  news(db = tools:::.build_news_db_from_package_NEWS_Rd("~/R/Pkgs/DPQ/inst/NEWS.Rd"))
\name{NEWS}
\title{DPQ News}
\encoding{UTF-8}

\section{Version 0.3-3 (2019-08-31, svn r93)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Renamed (and generalized / merged) many functions to have less "." in names.

      \item New \code{pnbetaAS310()} function interfacing to my
      \emph{corrected} C version of 'ASA 310' (2007).

      \item New \code{algdiv()} function interfacing to the 'TOMS 708' version of
      our \code{logQab()}.

      \item New \code{pl2curves()} which generalizes (somewhat) previous
      function \code{p.dnchiB()}.

      \item Made \code{newton()} more flexible with new \code{xMin} and
      \code{xMax} arguments which are notably useful for \code{q*()}
      (quantile function) computations.

      Correspondingly replaced previous \code{qchisq2()} and
      \code{qchisqA()} with new \code{qchisqN()} (\sQuote{N}ewton).

      \item new \code{pnchisqRC()}, a version of R's C level non-central
      chi-squared, with additional options.

      \item \dQuote{new} \code{logspace.add()} and \code{logspace.sub()}.
    }
  }
}

\section{Version 0.3-0 [2018-08-28]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Move many of my up to 15 years old DPQ computation utilities
      into a package, to become public, \dQuote{bloggable}, etc.
    }
  }
}
