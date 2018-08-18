#### Automagically produced by the Perl program ./C2R (version 0.8, RCS($Revision: 1.11 $))
#### from inputfile `gaminv.c'.

#### Probably do NOT EDIT (by hand) since it will be overwritten...

gaminv <- function(a, x, x0, p, q, ierr)
{
 .C(C_gaminv,
            a = as.double(a),
            x = as.double(x),
           x0 = as.double(x0),
            p = as.double(p),
            q = as.double(q),
         ierr = as.integer(ierr)
  )
}
pni <- function(p, q, d, w, ierr)
{
 .C(C_pni,
            p = as.double(p),
            q = as.double(q),
            d = as.double(d),
            w = as.double(w),
         ierr = as.integer(ierr)
  )
}
