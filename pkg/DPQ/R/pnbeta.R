pnbetaAppr2v1 <- ## pnbeta.appr2 <-
    function(x, a, b, ncp = 0, lower.tail=TRUE, log.p=FALSE)
{
  ## Purpose: "Approximation 2" of
  ## Chattamvelli, R. \& Shanmugam, F.  (1997), `Algorithm {AS 310}: ...
  ## JRSS C (Applied Statistics)}  46(1), 146--156; notably p.155
  ## ----------------------------------------------------------------------
  ## Arguments: as for pbeta(), but named more traditionally
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 16 Oct 2007, 11:36

                                        # λ <- ncp -- not liked by CRAN checks
    a2l <- 2*a + ncp
    c. <- (a2l/(2*a)) ^ (1/3)
    d <- ((b * x)/(a*(1-x))) ^ (1/3)
    ##        --      --- <==> allow to input "1-x" for x ~= 1
    ## sigma_1^2  and sigma_2 ^2 --> mu_j = 1 - sigma_j^2  (j = 1,2):
    s1.2 <- 4*(a + ncp) / (9 * a2l^2)
    s2.2 <- 1/(9 * b)
    mu.L <- c. * (1 - s1.2) - d   * (1 - s2.2)
    si.L <- sqrt(c.^2 * s1.2     + d^2 * s2.2)
    pnorm( - mu.L / si.L, lower.tail=lower.tail, log.p = log.p)
}

pnbetaAppr2 <- function(x, a, b, ncp = 0, lower.tail=TRUE, log.p=FALSE)
{
  ## Purpose: "Approximation 2" of
  ## Chattamvelli, R. \& Shanmugam, F.  (1997), `Algorithm {AS 310}: ...
  ## JRSS C (Applied Statistics)}  46(1), 146--156; notably p.155
  ## ----------------------------------------------------------------------
  ## Arguments: as for pbeta(), but named more traditionally
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 16 Aug 2018;
  ##         *Math* simplified: c == 1 (+ cancellation of "/a")
                                        # λ <- ncp
    a2l <- 2*a + ncp
    d <- ((2*b * x)/(a2l*(1-x))) ^ (1/3)
    ##        --      --- <==> allow to input "1-x" for x ~= 1
    ## sigma_1^2  and sigma_2 ^2 --> mu_j = 1 - sigma_j^2  (j = 1,2):
    s1.2 <- 4*(a + ncp) / (9 * a2l^2)
    s2.2 <- 1/(9 * b)
    mu.L <- (1 - s1.2) - d * (1 - s2.2)
    si.L <- sqrt(s1.2  + d^2 * s2.2)
    pnorm( - mu.L / si.L, lower.tail=lower.tail, log.p = log.p)
}
