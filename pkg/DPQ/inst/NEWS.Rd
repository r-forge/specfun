% Check from R:
%  news(db = tools:::.build_news_db_from_package_NEWS_Rd("~/R/Pkgs/DPQ/inst/NEWS.Rd"))
\name{NEWS}
\title{DPQ News}
\encoding{UTF-8}

\section{Version 0.3-3 [2019-08-16]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item Renamed many functions to have less "." in names.

      \item New \code{pnbetaAS310()} function interfacing to my
      \emph{corrected} C version of 'ASA 310' (2007).

      \item New \code{algdiv()} function interfacing to TOMS 708 version of
      our \code{logQab}.
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
