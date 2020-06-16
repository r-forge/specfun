% Check from R:
%  news(db = tools:::.build_news_db_from_package_NEWS_Rd("~/R/Pkgs/DPQmpfr/inst/NEWS.Rd"))
\name{NEWS}
\title{DPQmpfr News}
\encoding{UTF-8}

\section{Version 0.2-0 (2018-08-16, svn r40)}{
  \subsection{FEATURES}{
    \itemize{
      \item Split off from original \CRANpkg{DPQ} package.  This
      \pkg{DPQmpfr} package will explicitly depend on both CRAN R packages
      \CRANpkg{gmp} and \CRANpkg{Rmpfr}.

      \item Started with a MPFR-ized version of \CRANpkg{DPQ}'s
      \code{dnt()} function for the non-central t-distribution density.
    }
  }
}


\section{Version 0.3-0 (2020-06-12, svn r...)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{dhyperQ()}, \code{phyperQ()}, and \code{phyperQall()} for
      \emph{exact} hypergeometric probabilities (as \dQuote{big rational}s,
      i.e., \CRANpkg{gmp} package's \code{"bigq"}).
    }
  }
}

