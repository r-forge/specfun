% Check from R:
%  news(db = tools:::.build_news_db_from_package_NEWS_Rd("~/R/Pkgs/DPQmpfr/inst/NEWS.Rd"))
\name{NEWS}
\title{DPQmpfr News}
\encoding{UTF-8}

% NB: svn revision numbering within all my R-forge `specfun` pkgs --- ..../R-forge/specfun/
\section{Version 0.3-3 (2024-08-19, svn r327)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{stirlerrSer(n, k)} now allows \code{k <= 22}.

      \item add some \verb{\link[<pkg>]{..}} as now required.

      \item \code{dtWVm()} is an \dQuote{mpfr-ized} version of \CRANpkg{DPQ}'s
      \code{\link[DPQ]{dtWV}()}.

      \item \code{gam1M()}, \code{lgamma1pM()}, and \code{algdivM()}
      each generalize the corresponding function in package \CRANpkg{DPQ}
      for the use with \CRANpkg{Rmpfr}'s \code{"mpfr"} numbers, notably
      avoiding cancellation by adding additional precision bits.

      \item New \code{pbeta_bpser()} using the 3 auxiliary functions for an
      mpfr-ized version of \R's \code{pbeta()}.

      \item More non-central t probability computations .. unfinished currently %% FIXME FIXME
    }
  }
  \subsection{Misc}{
    \itemize{
      \item internal \code{dpoisEr()}, gets arg \code{log2.min}, back to previous (2020)
      value, such that "correct underflow" in R does not give rel.err of -1.
    }
  }
  \subsection{DEPRECATED AND DEFUNCT}{
    \itemize{
      \item \code{dntJKBm()} is deprecated, as \CRANpkg{DPQ}'s
      \code{\link[DPQ]{dntJKBf}()} works already perfectly \dQuote{mpfr-ized}.
    }
  }
}

\section{Version 0.3-2 (2023-12-05, svn r281)}{
  \subsection{NEW FEATURES}{
    \itemize{
      %% not yet exported nor documented  {maybe go to DPQ vignette tools ?}
      \item utility functions for testing \code{dpois_raw()} etc, % ../R/dpoisEr.R
      notably for package \CRANpkg{DPQ}.

      \item new \dQuote{mpfr-ized} \code{pnormAsymp()} and \code{qnormAsymp()}

      \item new \code{ldexp(f, E)} utility.
    }
  }
  \subsection{Misc}{
    \itemize{
      \item URLs in \file{DESCRIPTION}
    }
  }
}

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
