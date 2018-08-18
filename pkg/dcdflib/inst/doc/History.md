# Historical notes about the 'dcdflib' package:

## MM had this old package in my numerical experiments
	`~/R/MM/NUMERICS/dpq-functions/dcdflib`

### I had started a version of that to the C version,
   at `~/R/MM/NUMERICS/dpq-functions/dcdflib-c`  from 1997-1999 (!)
   but not used much.

- The original Fortran and C versions of the underlying library are still
   (2018-08) available from Netlib:

   http://www.netlib.org/random/  the NETLIB random section, distributes a
   tar file of the source code for the original CDFLIB library in C and FORTRAN
   which on 2018-08-17 'file' still says

>  last modified: Thu Apr 28 04:00:07 1994, max compression, from Unix

   and I keep a copy of these in  ~/src/BBrown/netlib-src/

	   59554 Apr 28  1994 dcdflib.c.tar.gz
	   63130 Apr 28  1994 dcdflib.f.tar.gz

- *But* do note that my own versions of DCDFLIB in `~/src/BBrown/` are
  slightly newer (1997 instead of 1994),

  *and* --- as John Burkardt's page

     https://people.sc.fsu.edu/~jburkardt/f_src/cdflib/cdflib.html

   also mentions --- that somewhat newer (only slightly updated) versions
   of DCDFLIB can now be "ordered"  "by registration" from

	https://biostatistics.mdanderson.org/SoftwareDownload/SingleSoftware.aspx?Software_Id=21

   which is the Unversity of Texas Md Anderson (Cancer) Center's "software kiosk".

- I am pretty positive I had used  f2c  myself of the Fortran DCDFLIB code,
  and worked with that generated "ANSI C" code.

- In the fortran code (i.e., before translation), I had removed the final
  'f' at the end of numeric constants, so they became double precision rather
  than "float" (single precision).



## Also partly used an S interface (at Statlib ??) by Donald H. MacQueen
	(mailto:macqueen1@llnl.gov) e.g. mentioned in `../../man/ptnc.Rd`
   with original S wrapper to the Fortran version of DCDFLIB.

