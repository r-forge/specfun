#### Automagically produced by the Perl program ./C2R (version 0.8, RCS($Revision: 1.11 $))
#### from inputfile `gratio.c'.

#### Probably do NOT EDIT (by hand) since it will be overwritten...

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
