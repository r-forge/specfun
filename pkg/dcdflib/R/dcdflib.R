#### Automagically produced by the Perl program ./C2R (version 0.8, RCS($Revision: 1.11 $))
#### from inputfile `dcdflib.c'.

#### Probably do NOT EDIT (by hand) since it will be overwritten...

bgrat <- function(a, b, x, y, w, eps, ierr)
{
 .C("bgrat",
            a = as.double(a),
            b = as.double(b),
            x = as.double(x),
            y = as.double(y),
            w = as.double(w),
          eps = as.double(eps),
         ierr = as.integer(ierr)
  , PACKAGE = "dcdflib")
}
bratio <- function(a, b, x, y, w, w1, ierr)
{
 .C("bratio",
            a = as.double(a),
            b = as.double(b),
            x = as.double(x),
            y = as.double(y),
            w = as.double(w),
           w1 = as.double(w1),
         ierr = as.integer(ierr)
  , PACKAGE = "dcdflib")
}
cdfbet <- function(which, p, q, x, y, a, b, status, bound)
{
 .C("cdfbet",
        which = as.integer(which),
            p = as.double(p),
            q = as.double(q),
            x = as.double(x),
            y = as.double(y),
            a = as.double(a),
            b = as.double(b),
       status = as.integer(status),
        bound = as.double(bound)
  , PACKAGE = "dcdflib")
}
cdfbin <- function(which, p, q, s, xn, pr, ompr, status, bound)
{
 .C("cdfbin",
        which = as.integer(which),
            p = as.double(p),
            q = as.double(q),
            s = as.double(s),
           xn = as.double(xn),
           pr = as.double(pr),
         ompr = as.double(ompr),
       status = as.integer(status),
        bound = as.double(bound)
  , PACKAGE = "dcdflib")
}
cdfchi <- function(which, p, q, x, df, status, bound)
{
 .C("cdfchi",
        which = as.integer(which),
            p = as.double(p),
            q = as.double(q),
            x = as.double(x),
           df = as.double(df),
       status = as.integer(status),
        bound = as.double(bound)
  , PACKAGE = "dcdflib")
}
cdfchn <- function(which, p, q, x, df, pnonc, status, bound)
{
 .C("cdfchn",
        which = as.integer(which),
            p = as.double(p),
            q = as.double(q),
            x = as.double(x),
           df = as.double(df),
        pnonc = as.double(pnonc),
       status = as.integer(status),
        bound = as.double(bound)
  , PACKAGE = "dcdflib")
}
cdff <- function(which, p, q, f, dfn, dfd, status, bound)
{
 .C("cdff",
        which = as.integer(which),
            p = as.double(p),
            q = as.double(q),
            f = as.double(f),
          dfn = as.double(dfn),
          dfd = as.double(dfd),
       status = as.integer(status),
        bound = as.double(bound)
  , PACKAGE = "dcdflib")
}
cdffnc <- function(which, p, q, f, dfn, dfd, pnonc, status, bound)
{
 .C("cdffnc",
        which = as.integer(which),
            p = as.double(p),
            q = as.double(q),
            f = as.double(f),
          dfn = as.double(dfn),
          dfd = as.double(dfd),
        pnonc = as.double(pnonc),
       status = as.integer(status),
        bound = as.double(bound)
  , PACKAGE = "dcdflib")
}
cdfgam <- function(which, p, q, x, shape, scale, status, bound)
{
 .C("cdfgam",
        which = as.integer(which),
            p = as.double(p),
            q = as.double(q),
            x = as.double(x),
        shape = as.double(shape),
        scale = as.double(scale),
       status = as.integer(status),
        bound = as.double(bound)
  , PACKAGE = "dcdflib")
}
cdfnbn <- function(which, p, q, s, xn, pr, ompr, status, bound)
{
 .C("cdfnbn",
        which = as.integer(which),
            p = as.double(p),
            q = as.double(q),
            s = as.double(s),
           xn = as.double(xn),
           pr = as.double(pr),
         ompr = as.double(ompr),
       status = as.integer(status),
        bound = as.double(bound)
  , PACKAGE = "dcdflib")
}
cdfnor <- function(which, p, q, x, mean, sd, status, bound)
{
 .C("cdfnor",
        which = as.integer(which),
            p = as.double(p),
            q = as.double(q),
            x = as.double(x),
         mean = as.double(mean),
           sd = as.double(sd),
       status = as.integer(status),
        bound = as.double(bound)
  , PACKAGE = "dcdflib")
}
cdfpoi <- function(which, p, q, s, xlam, status, bound)
{
 .C("cdfpoi",
        which = as.integer(which),
            p = as.double(p),
            q = as.double(q),
            s = as.double(s),
         xlam = as.double(xlam),
       status = as.integer(status),
        bound = as.double(bound)
  , PACKAGE = "dcdflib")
}
cdft <- function(which, p, q, t, df, status, bound)
{
 .C("cdft",
        which = as.integer(which),
            p = as.double(p),
            q = as.double(q),
            t = as.double(t),
           df = as.double(df),
       status = as.integer(status),
        bound = as.double(bound)
  , PACKAGE = "dcdflib")
}
cdftnc <- function(which, p, q, t, df, pnonc, status, bound)
{
 .C("cdftnc",
        which = as.integer(which),
            p = as.double(p),
            q = as.double(q),
            t = as.double(t),
           df = as.double(df),
        pnonc = as.double(pnonc),
       status = as.integer(status),
        bound = as.double(bound)
  , PACKAGE = "dcdflib")
}
cumbet <- function(x, y, a, b, cum, ccum)
{
 .C("cumbet",
            x = as.double(x),
            y = as.double(y),
            a = as.double(a),
            b = as.double(b),
          cum = as.double(cum),
         ccum = as.double(ccum)
  , PACKAGE = "dcdflib")
}
cumbin <- function(s, xn, pr, ompr, cum, ccum)
{
 .C("cumbin",
            s = as.double(s),
           xn = as.double(xn),
           pr = as.double(pr),
         ompr = as.double(ompr),
          cum = as.double(cum),
         ccum = as.double(ccum)
  , PACKAGE = "dcdflib")
}
cumchi <- function(x, df, cum, ccum)
{
 .C("cumchi",
            x = as.double(x),
           df = as.double(df),
          cum = as.double(cum),
         ccum = as.double(ccum)
  , PACKAGE = "dcdflib")
}
cumchn <- function(x, df, pnonc, cum, ccum)
{
 .C("cumchn",
            x = as.double(x),
           df = as.double(df),
        pnonc = as.double(pnonc),
          cum = as.double(cum),
         ccum = as.double(ccum)
  , PACKAGE = "dcdflib")
}
cumf <- function(f, dfn, dfd, cum, ccum)
{
 .C("cumf",
            f = as.double(f),
          dfn = as.double(dfn),
          dfd = as.double(dfd),
          cum = as.double(cum),
         ccum = as.double(ccum)
  , PACKAGE = "dcdflib")
}
cumfnc <- function(f, dfn, dfd, pnonc, cum, ccum)
{
 .C("cumfnc",
            f = as.double(f),
          dfn = as.double(dfn),
          dfd = as.double(dfd),
        pnonc = as.double(pnonc),
          cum = as.double(cum),
         ccum = as.double(ccum)
  , PACKAGE = "dcdflib")
}
cumgam <- function(x, a, cum, ccum)
{
 .C("cumgam",
            x = as.double(x),
            a = as.double(a),
          cum = as.double(cum),
         ccum = as.double(ccum)
  , PACKAGE = "dcdflib")
}
cumnbn <- function(s, xn, pr, ompr, cum, ccum)
{
 .C("cumnbn",
            s = as.double(s),
           xn = as.double(xn),
           pr = as.double(pr),
         ompr = as.double(ompr),
          cum = as.double(cum),
         ccum = as.double(ccum)
  , PACKAGE = "dcdflib")
}
cumnor <- function(arg, result, ccum)
{
 .C("cumnor",
          arg = as.double(arg),
       result = as.double(result),
         ccum = as.double(ccum)
  , PACKAGE = "dcdflib")
}
cumpoi <- function(s, xlam, cum, ccum)
{
 .C("cumpoi",
            s = as.double(s),
         xlam = as.double(xlam),
          cum = as.double(cum),
         ccum = as.double(ccum)
  , PACKAGE = "dcdflib")
}
cumt <- function(t, df, cum, ccum)
{
 .C("cumt",
            t = as.double(t),
           df = as.double(df),
          cum = as.double(cum),
         ccum = as.double(ccum)
  , PACKAGE = "dcdflib")
}
cumtnc <- function(t, df, pnonc, cum, ccum)
{
 .C("cumtnc",
            t = as.double(t),
           df = as.double(df),
        pnonc = as.double(pnonc),
          cum = as.double(cum),
         ccum = as.double(ccum)
  , PACKAGE = "dcdflib")
}
dinvr <- function(status, x, fx, qleft, qhi)
{
 .C("dinvr",
       status = as.integer(status),
            x = as.double(x),
           fx = as.double(fx),
        qleft = as.logical(qleft),
          qhi = as.logical(qhi)
  , PACKAGE = "dcdflib")
}
dstinv <- function(zsmall, zbig, zabsst, zrelst, zstpmu, zabsto, zrelto)
{
 .C("dstinv",
       zsmall = as.double(zsmall),
         zbig = as.double(zbig),
       zabsst = as.double(zabsst),
       zrelst = as.double(zrelst),
       zstpmu = as.double(zstpmu),
       zabsto = as.double(zabsto),
       zrelto = as.double(zrelto)
  , PACKAGE = "dcdflib")
}
dzror <- function(status, x, fx, xlo, xhi, qleft, qhi)
{
 .C("dzror",
       status = as.integer(status),
            x = as.double(x),
           fx = as.double(fx),
          xlo = as.double(xlo),
          xhi = as.double(xhi),
        qleft = as.logical(qleft),
          qhi = as.logical(qhi)
  , PACKAGE = "dcdflib")
}
dstzr <- function(zxlo, zxhi, zabstl, zreltl)
{
 .C("dstzr",
         zxlo = as.double(zxlo),
         zxhi = as.double(zxhi),
       zabstl = as.double(zabstl),
       zreltl = as.double(zreltl)
  , PACKAGE = "dcdflib")
}
grat1 <- function(a, x, r, p, q, eps)
{
 .C("grat1",
            a = as.double(a),
            x = as.double(x),
            r = as.double(r),
            p = as.double(p),
            q = as.double(q),
          eps = as.double(eps)
  , PACKAGE = "dcdflib")
}
gratio <- function(a, x, ans, qans, ind)
{
 .C("gratio",
            a = as.double(a),
            x = as.double(x),
          ans = as.double(ans),
         qans = as.double(qans),
          ind = as.integer(ind)
  , PACKAGE = "dcdflib")
}
ftnstop <- function(msg)
{
 .C("ftnstop",
          msg = as.character(msg)
  , PACKAGE = "dcdflib")
}
