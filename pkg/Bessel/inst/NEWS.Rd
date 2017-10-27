% Check from R:
%  news(db = tools:::.build_news_db_from_package_NEWS_Rd("~/R/Pkgs/Bessel/inst/NEWS.Rd"))
\name{NEWS}
\title{News for \R Package \pkg{Bessel}}
\encoding{UTF-8}

\section{Changes in version 0.6-0 (2017-10-27, svn r33)}{
  \subsection{New Features}{
    \itemize{
      \item deprecate \code{bI()} to \code{besselSi()} -- enable
      \code{log = TRUE} option and test it.

      \item new \code{besselJs()}, i.e., Taylor series for \eqn{J_n(x)};
      using \code{lssum()} for alternating exp/log sum.

      \item use C routine registration, adding safety and a bit of speed.

      \item improved documentation
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item fix \code{BesselI(<cmplx>, nu < 0, expo..=TRUE)} bug thanks
      to Hiroyuki Kawakatsu.

      \item better warning message (r29).
    }
  }
}

\section{Changes in version 0.5-5 (2013-12-10, svn r17)}{
  \subsection{Bug Fixes}{
    \itemize{
      \item Adapt to new standard: move vignettes to own
      \file{vignettes/} directory.
    }
  }
}

\section{Changes in version 0.5-4 (2012-03-21, svn r16)}{
  \subsection{New Features}{
    \itemize{
      \item Improved vignette
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item \code{BesselY(0, nu)} should not give an error message
	from \code{zbesy()}, but either -Inf (real 0) or something like
	\code{"Inf + NaN i"}.
      \item Cosmetic, mainly for \code{R CMD check}
    }
  }
}

\section{Changes in version 0.5-3 (2009-12-22 -- 2010-11-29, svn r8)}{
  \subsection{New Features}{
    \itemize{
      \item \code{besselI.nuAsym()} and  \code{besselK.nuAsym()}
      asymptotic versions for large \code{nu} and/or \code{x}
    }
  }
}

\section{ALL PREVIOUS changes (version 0.0-0)}{ \itemize{ \item
    See the  \file{ChangeLog} file in the sources.
}}

