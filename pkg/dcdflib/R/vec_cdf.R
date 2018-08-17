#### Automagically produced by
###  cd ../src
'
 ./C2R vec_cdf.c | ./R2VR > vec_cdf.R
'

cdfbetV <- function(which, p, q, x, y, a, b, status, bound)
{
  len <- max(length(p), length(q), length(x), length(y), length(a), length(b), length(status), length(bound))
 .C("V_cdfbet",
        which = as.integer(which),
            p = rep(as.double(p), length = len),
            q = rep(as.double(q), length = len),
            x = rep(as.double(x), length = len),
            y = rep(as.double(y), length = len),
            a = rep(as.double(a), length = len),
            b = rep(as.double(b), length = len),
       status = rep(as.integer(status), length = len),
        bound = rep(as.double(bound), length = len),
          len = as.integer(len)
  , PACKAGE = "dcdflib")
}

cdfbinV <- function(which, p, q, s, xn, pr, ompr, status, bound)
{
  len <- max(length(p), length(q), length(s), length(xn), length(pr), length(ompr), length(status), length(bound))
 .C("V_cdfbin",
        which = as.integer(which),
            p = rep(as.double(p), length = len),
            q = rep(as.double(q), length = len),
            s = rep(as.double(s), length = len),
           xn = rep(as.double(xn), length = len),
           pr = rep(as.double(pr), length = len),
         ompr = rep(as.double(ompr), length = len),
       status = rep(as.integer(status), length = len),
        bound = rep(as.double(bound), length = len),
          len = as.integer(len)
  , PACKAGE = "dcdflib")
}

cdfchiV <- function(which, p, q, x, df, status, bound)
{
  len <- max(length(p), length(q), length(x), length(df), length(status), length(bound))
 .C("V_cdfchi",
        which = as.integer(which),
            p = rep(as.double(p), length = len),
            q = rep(as.double(q), length = len),
            x = rep(as.double(x), length = len),
           df = rep(as.double(df), length = len),
       status = rep(as.integer(status), length = len),
        bound = rep(as.double(bound), length = len),
          len = as.integer(len)
  , PACKAGE = "dcdflib")
}

cdfchnV <- function(which, p, q, x, df, pnonc, status, bound)
{
  len <- max(length(p), length(q), length(x), length(df), length(pnonc), length(status), length(bound))
 .C("V_cdfchn",
        which = as.integer(which),
            p = rep(as.double(p), length = len),
            q = rep(as.double(q), length = len),
            x = rep(as.double(x), length = len),
           df = rep(as.double(df), length = len),
        pnonc = rep(as.double(pnonc), length = len),
       status = rep(as.integer(status), length = len),
        bound = rep(as.double(bound), length = len),
          len = as.integer(len)
  , PACKAGE = "dcdflib")
}

cdffV <- function(which, p, q, f, dfn, dfd, status, bound)
{
  len <- max(length(p), length(q), length(f), length(dfn), length(dfd), length(status), length(bound))
 .C("V_cdff",
        which = as.integer(which),
            p = rep(as.double(p), length = len),
            q = rep(as.double(q), length = len),
            f = rep(as.double(f), length = len),
          dfn = rep(as.double(dfn), length = len),
          dfd = rep(as.double(dfd), length = len),
       status = rep(as.integer(status), length = len),
        bound = rep(as.double(bound), length = len),
          len = as.integer(len)
  , PACKAGE = "dcdflib")
}

cdffncV <- function(which, p, q, f, dfn, dfd, phonc, status, bound)
{
  len <- max(length(p), length(q), length(f), length(dfn), length(dfd), length(phonc), length(status), length(bound))
 .C("V_cdffnc",
        which = as.integer(which),
            p = rep(as.double(p), length = len),
            q = rep(as.double(q), length = len),
            f = rep(as.double(f), length = len),
          dfn = rep(as.double(dfn), length = len),
          dfd = rep(as.double(dfd), length = len),
        phonc = rep(as.double(phonc), length = len),
       status = rep(as.integer(status), length = len),
        bound = rep(as.double(bound), length = len),
          len = as.integer(len)
  , PACKAGE = "dcdflib")
}

cdfgamV <- function(which, p, q, x, shape, scale, status, bound)
{
  len <- max(length(p), length(q), length(x), length(shape), length(scale), length(status), length(bound))
 .C("V_cdfgam",
        which = as.integer(which),
            p = rep(as.double(p), length = len),
            q = rep(as.double(q), length = len),
            x = rep(as.double(x), length = len),
        shape = rep(as.double(shape), length = len),
        scale = rep(as.double(scale), length = len),
       status = rep(as.integer(status), length = len),
        bound = rep(as.double(bound), length = len),
          len = as.integer(len)
  , PACKAGE = "dcdflib")
}

cdfnbnV <- function(which, p, q, s, xn, pr, ompr, status, bound)
{
  len <- max(length(p), length(q), length(s), length(xn), length(pr), length(ompr), length(status), length(bound))
 .C("V_cdfnbn",
        which = as.integer(which),
            p = rep(as.double(p), length = len),
            q = rep(as.double(q), length = len),
            s = rep(as.double(s), length = len),
           xn = rep(as.double(xn), length = len),
           pr = rep(as.double(pr), length = len),
         ompr = rep(as.double(ompr), length = len),
       status = rep(as.integer(status), length = len),
        bound = rep(as.double(bound), length = len),
          len = as.integer(len)
  , PACKAGE = "dcdflib")
}

cdfnorV <- function(which, p, q, x, mean, sd, status, bound)
{
  len <- max(length(p), length(q), length(x), length(mean), length(sd), length(status), length(bound))
 .C("V_cdfnor",
        which = as.integer(which),
            p = rep(as.double(p), length = len),
            q = rep(as.double(q), length = len),
            x = rep(as.double(x), length = len),
         mean = rep(as.double(mean), length = len),
           sd = rep(as.double(sd), length = len),
       status = rep(as.integer(status), length = len),
        bound = rep(as.double(bound), length = len),
          len = as.integer(len)
  , PACKAGE = "dcdflib")
}

cdfpoiV <- function(which, p, q, s, xlam, status, bound)
{
  len <- max(length(p), length(q), length(s), length(xlam), length(status), length(bound))
 .C("V_cdfpoi",
        which = as.integer(which),
            p = rep(as.double(p), length = len),
            q = rep(as.double(q), length = len),
            s = rep(as.double(s), length = len),
         xlam = rep(as.double(xlam), length = len),
       status = rep(as.integer(status), length = len),
        bound = rep(as.double(bound), length = len),
          len = as.integer(len)
  , PACKAGE = "dcdflib")
}

cdftV <- function(which, p, q, t, df, status, bound)
{
  len <- max(length(p), length(q), length(t), length(df), length(status), length(bound))
 .C("V_cdft",
        which = as.integer(which),
            p = rep(as.double(p), length = len),
            q = rep(as.double(q), length = len),
            t = rep(as.double(t), length = len),
           df = rep(as.double(df), length = len),
       status = rep(as.integer(status), length = len),
        bound = rep(as.double(bound), length = len),
          len = as.integer(len)
  , PACKAGE = "dcdflib")
}

cdftncV <- function(which, p, q, t, df, pnonc, status, bound)
{
  len <- max(length(p), length(q), length(t), length(df), length(pnonc), length(status), length(bound))
 .C("V_cdftnc",
        which = as.integer(which),
            p = rep(as.double(p), length = len),
            q = rep(as.double(q), length = len),
            t = rep(as.double(t), length = len),
           df = rep(as.double(df), length = len),
        pnonc = rep(as.double(pnonc), length = len),
       status = rep(as.integer(status), length = len),
        bound = rep(as.double(bound), length = len),
          len = as.integer(len)
  , PACKAGE = "dcdflib")
}

