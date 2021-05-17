% Check from R:
%  news(db = tools:::.build_news_db_from_package_NEWS_Rd("~/R/Pkgs/DPQmpfr/inst/NEWS.Rd"))
\name{NEWS}
\title{DPQmpfr News}
\encoding{UTF-8}

\section{Version 0.3-1 (2021-05-17, svn r174)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{stirlerrM(n)} and \code{stirlerrSer(n,k)} providing
      Stirling's formula approximation errors for bigrational ("bigz") or
      "mpfr" numbers.
    }
  }
}

\section{Version 0.3-0 (2020-11-23, svn r156)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{dhyperQ()}, \code{phyperQ()}, and \code{phyperQall()} for
      \emph{exact} hypergeometric probabilities (as \dQuote{big rational}s,
      i.e., \CRANpkg{gmp} package's \code{"bigq"}).

      \item Finally export (and add to svn!) \code{dbetaD94()},
      \code{pbetaD94()}, and \code{qbetaD94()}.

      Add \code{log_scale} option to both \code{pbeta*()} and \code{qbeta*()}.

      \item \dQuote{mpfr-ize} more, notably \code{pnorm*()} from
      \CRANpkg{DPQ} (0.4-2). % ../R/dpq-mpfrized.R
    }
  }
}

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


