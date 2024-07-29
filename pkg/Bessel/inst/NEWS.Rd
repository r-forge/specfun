% Check from R:
%  news(db = tools:::.build_news_db_from_package_NEWS_Rd("~/R/Pkgs/Bessel/inst/NEWS.Rd"))
\name{NEWS}
\title{News for \R Package \pkg{Bessel}}
\encoding{UTF-8}

\section{Changes in version 0.6-1 (2024-07-29, svn r323)}{
  \subsection{New Features}{
    \itemize{
      \item Using \code{cospi(nu)} instead of \code{cos(pi*nu)}, for
      negative \code{nu}, both in \code{BesselJ()} and \code{BesselY()} may
      increase accuracy.
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item Rd-crossrefs: replace \code{\\linkS4class\{mpfr\}} by
      \code{\\link[Rmpfr:mpfr-class]\{mpfr\}}.
    }
  }
}

\section{Changes in version 0.6-0 (2019-04-24, svn r67)}{
  \subsection{SIGNIFICANT USER-VISIBLE CHANGES}{
    \itemize{
      \item The exponential rescaling for the asymptotic version of the
      \eqn{I()} function has been changed (for non-real
      complex arguments) to be compatible to the scaling for
      \code{\link{BesselI}()}.
    }
  }
  \subsection{New Features}{
    \itemize{
      \item deprecate \code{bI()} to use \code{besselIs()} instead --
      enable \code{log = TRUE} option and test it.

      \item new \code{besselJs()}, i.e., Taylor series for \eqn{J_n(x)};
      using \code{lssum()} for alternating exp/log sum, also works with
      \CRANpkg{Rmpfr}'s \code{mpfr()} high precision numbers.

      \item use C routine registration, adding safety and a bit of speed.

      \item improved documentation for Bessel functions Ai,Bi, H, I,J,K,Y.

      \item \code{besselI.nuAsym()} and  \code{besselK.nuAsym()} now allow maximal
      \code{k.max = 5} (instead of \code{4}).

      \item new \code{besselKasym()} (9.7.2 of A.&S. p.378).

      \item All the TOMS 644 routines get an optional argument
      \code{verbose} indicating the verbosity level of output notably from
      C code.  The output itself may change and is \emph{not} part of the API.
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item fix \code{BesselI(<cmplx>, nu < 0, expo..=TRUE)} bug thanks
      to Hiroyuki Kawakatsu.

      \item better warning message (r29).

      \item the TOMS 644 routines Bessel[HIYJK] and Airy[AB] now all
      return \code{NaN} with a warning instead of stopping with an error
      in case of low level error \code{ierr == 4}.

      \item \code{BesselK(0,*)} is zero now, and \code{BesselH(m,
      z=0,..)} returns Inf/NaN combinations.% possibly to be improved
      Both no longer end with an error (from \code{ierr == 1}).

      \item \code{BesselK(<negative>, 3)} and \code{BesselJ(*)} now also
      return the correct imaginary part instead of just a numeric value.

      \item \code{besselIs(z, nu)} now also works for \code{complex} \code{z}.
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

