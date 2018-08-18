#### Automagically produced by the Perl program ./C2R (version 0.8, RCS($Revision: 1.11 $))
#### from inputfile `dcdflib.c'.

#### Probably do NOT EDIT (by hand) since it will be overwritten...

bgrat <- function(a, b, x, y, w, eps, ierr)
{
 .C(C_bgrat,
            a = as.double(a),
            b = as.double(b),
            x = as.double(x),
            y = as.double(y),
            w = as.double(w),
          eps = as.double(eps),
         ierr = as.integer(ierr)
  )
}
bratio <- function(a, b, x, y, w, w1, ierr)
{
 .C(C_bratio,
            a = as.double(a),
            b = as.double(b),
            x = as.double(x),
            y = as.double(y),
            w = as.double(w),
           w1 = as.double(w1),
         ierr = as.integer(ierr)
  )
}
cdfbet <- function(which, p, q, x, y, a, b, status, bound)
{
 .C(C_cdfbet,
        which = as.integer(which),
            p = as.double(p),
            q = as.double(q),
            x = as.double(x),
            y = as.double(y),
            a = as.double(a),
            b = as.double(b),
       status = as.integer(status),
        bound = as.double(bound)
  )
}
cdfbin <- function(which, p, q, s, xn, pr, ompr, status, bound)
{
 .C(C_cdfbin,
        which = as.integer(which),
            p = as.double(p),
            q = as.double(q),
            s = as.double(s),
           xn = as.double(xn),
           pr = as.double(pr),
         ompr = as.double(ompr),
       status = as.integer(status),
        bound = as.double(bound)
  )
}
cdfchi <- function(which, p, q, x, df, status, bound)
{
 .C(C_cdfchi,
        which = as.integer(which),
            p = as.double(p),
            q = as.double(q),
            x = as.double(x),
           df = as.double(df),
       status = as.integer(status),
        bound = as.double(bound)
  )
}
cdfchn <- function(which, p, q, x, df, pnonc, status, bound)
{
 .C(C_cdfchn,
        which = as.integer(which),
            p = as.double(p),
            q = as.double(q),
            x = as.double(x),
           df = as.double(df),
        pnonc = as.double(pnonc),
       status = as.integer(status),
        bound = as.double(bound)
  )
}
cdff <- function(which, p, q, f, dfn, dfd, status, bound)
{
 .C(C_cdff,
        which = as.integer(which),
            p = as.double(p),
            q = as.double(q),
            f = as.double(f),
          dfn = as.double(dfn),
          dfd = as.double(dfd),
       status = as.integer(status),
        bound = as.double(bound)
  )
}
cdffnc <- function(which, p, q, f, dfn, dfd, pnonc, status, bound)
{
 .C(C_cdffnc,
        which = as.integer(which),
            p = as.double(p),
            q = as.double(q),
            f = as.double(f),
          dfn = as.double(dfn),
          dfd = as.double(dfd),
        pnonc = as.double(pnonc),
       status = as.integer(status),
        bound = as.double(bound)
  )
}
cdfgam <- function(which, p, q, x, shape, scale, status, bound)
{
 .C(C_cdfgam,
        which = as.integer(which),
            p = as.double(p),
            q = as.double(q),
            x = as.double(x),
        shape = as.double(shape),
        scale = as.double(scale),
       status = as.integer(status),
        bound = as.double(bound)
  )
}
cdfnbn <- function(which, p, q, s, xn, pr, ompr, status, bound)
{
 .C(C_cdfnbn,
        which = as.integer(which),
            p = as.double(p),
            q = as.double(q),
            s = as.double(s),
           xn = as.double(xn),
           pr = as.double(pr),
         ompr = as.double(ompr),
       status = as.integer(status),
        bound = as.double(bound)
  )
}
cdfnor <- function(which, p, q, x, mean, sd, status, bound)
{
 .C(C_cdfnor,
        which = as.integer(which),
            p = as.double(p),
            q = as.double(q),
            x = as.double(x),
         mean = as.double(mean),
           sd = as.double(sd),
       status = as.integer(status),
        bound = as.double(bound)
  )
}
cdfpoi <- function(which, p, q, s, xlam, status, bound)
{
 .C(C_cdfpoi,
        which = as.integer(which),
            p = as.double(p),
            q = as.double(q),
            s = as.double(s),
         xlam = as.double(xlam),
       status = as.integer(status),
        bound = as.double(bound)
  )
}
cdft <- function(which, p, q, t, df, status, bound)
{
 .C(C_cdft,
        which = as.integer(which),
            p = as.double(p),
            q = as.double(q),
            t = as.double(t),
           df = as.double(df),
       status = as.integer(status),
        bound = as.double(bound)
  )
}
cdftnc <- function(which, p, q, t, df, pnonc, status, bound)
{
 .C(C_cdftnc,
        which = as.integer(which),
            p = as.double(p),
            q = as.double(q),
            t = as.double(t),
           df = as.double(df),
        pnonc = as.double(pnonc),
       status = as.integer(status),
        bound = as.double(bound)
  )
}
cumbet <- function(x, y, a, b, cum, ccum)
{
 .C(C_cumbet,
            x = as.double(x),
            y = as.double(y),
            a = as.double(a),
            b = as.double(b),
          cum = as.double(cum),
         ccum = as.double(ccum)
  )
}
cumbin <- function(s, xn, pr, ompr, cum, ccum)
{
 .C(C_cumbin,
            s = as.double(s),
           xn = as.double(xn),
           pr = as.double(pr),
         ompr = as.double(ompr),
          cum = as.double(cum),
         ccum = as.double(ccum)
  )
}
cumchi <- function(x, df, cum, ccum)
{
 .C(C_cumchi,
            x = as.double(x),
           df = as.double(df),
          cum = as.double(cum),
         ccum = as.double(ccum)
  )
}
cumchn <- function(x, df, pnonc, cum, ccum)
{
 .C(C_cumchn,
            x = as.double(x),
           df = as.double(df),
        pnonc = as.double(pnonc),
          cum = as.double(cum),
         ccum = as.double(ccum)
  )
}
cumf <- function(f, dfn, dfd, cum, ccum)
{
 .C(C_cumf,
            f = as.double(f),
          dfn = as.double(dfn),
          dfd = as.double(dfd),
          cum = as.double(cum),
         ccum = as.double(ccum)
  )
}
cumfnc <- function(f, dfn, dfd, pnonc, cum, ccum)
{
 .C(C_cumfnc,
            f = as.double(f),
          dfn = as.double(dfn),
          dfd = as.double(dfd),
        pnonc = as.double(pnonc),
          cum = as.double(cum),
         ccum = as.double(ccum)
  )
}
cumgam <- function(x, a, cum, ccum)
{
 .C(C_cumgam,
            x = as.double(x),
            a = as.double(a),
          cum = as.double(cum),
         ccum = as.double(ccum)
  )
}
cumnbn <- function(s, xn, pr, ompr, cum, ccum)
{
 .C(C_cumnbn,
            s = as.double(s),
           xn = as.double(xn),
           pr = as.double(pr),
         ompr = as.double(ompr),
          cum = as.double(cum),
         ccum = as.double(ccum)
  )
}
cumnor <- function(arg, result, ccum)
{
 .C(C_cumnor,
          arg = as.double(arg),
       result = as.double(result),
         ccum = as.double(ccum)
  )
}
cumpoi <- function(s, xlam, cum, ccum)
{
 .C(C_cumpoi,
            s = as.double(s),
         xlam = as.double(xlam),
          cum = as.double(cum),
         ccum = as.double(ccum)
  )
}
cumt <- function(t, df, cum, ccum)
{
 .C(C_cumt,
            t = as.double(t),
           df = as.double(df),
          cum = as.double(cum),
         ccum = as.double(ccum)
  )
}
cumtnc <- function(t, df, pnonc, cum, ccum)
{
 .C(C_cumtnc,
            t = as.double(t),
           df = as.double(df),
        pnonc = as.double(pnonc),
          cum = as.double(cum),
         ccum = as.double(ccum)
  )
}

grat1 <- function(a, x, r, p, q, eps)
{
 .C(C_grat1,
            a = as.double(a),
            x = as.double(x),
            r = as.double(r),
            p = as.double(p),
            q = as.double(q),
          eps = as.double(eps)
  )
}
gratio <- function(a, x, ans, qans, ind)
{
 .C(C_gratio,
            a = as.double(a),
            x = as.double(x),
          ans = as.double(ans),
         qans = as.double(qans),
          ind = as.integer(ind)
  )
}
