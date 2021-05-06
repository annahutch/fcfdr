#' Function to perform cFDR for binary auxiliary covariates
#'
#' @param p p values for principal trait (vector of length n)
#' @param q binary auxillary data values (vector of length n)
#' @param chr chromosome each SNP resides (vector of length n)
#'
#' @importFrom Hmisc approxExtrap
#' 
#' @return dataframe of p, q and v values
#' @export
#'
binary_cfdr <- function(p, q, chr){

  unique_chr <- unique(chr)

  # split p and q into chromosomes
  p_res <- split(p, f = chr)
  q_res <- split(q, f = chr)
  minp=min(p)
  maxp=max(p)

  # prepare container for v
  v_res <- vector(mode = "list", length = length(unique_chr))

  ## prepare x for approxfun
  logx=seq(log10(minp),log10(maxp),length.out=1000)
  x=c(exp(logx),1)

  for(j in 1:length(unique_chr)){

    chrom <- unique_chr[j]

    p_loo <- p[-which(chr == chrom)] # df leave-one-chromosome-out
    q_loo <- q[-which(chr == chrom)] # df leave-one-chromosome-out
    ps <- p_res[[j]]
    qs <- q_res[[j]]

    q0 <- sum(q_loo == 1 & p_loo > 0.5)/sum(p_loo > 0.5)
    mult <- (sum(q_loo == 0 & p_loo > 1/2)/sum(q_loo == 1 & p_loo > 1/2))
    q0_sol=sapply(ps, function(p) max(sum(p_loo <= p & q_loo==0),1))
    q1_sol=sapply(ps, function(p) max(sum(p_loo <= p & q_loo==1),1))
    sol <- ifelse(qs==0,
                  mult*ps / q0_sol,
                  (1/mult)*ps / q1_sol)

    ## approx g0
    y=x/sapply(x, function(p) max(sum(p_loo <= p & q_loo==0),1))
    extr=Hmisc::approxExtrap(y,x,xout=unique(sol))
    invg0=approxfun(x=extr$x,y=pmax(pmin(extr$y,1),0),rule=2)
    ## invg0_spline=splinefun(x=y,y=x,method="hyman")

    ## approx g1
    y1=x/sapply(x, function(p) max(sum(p_loo <= p & q_loo==1),1))
    extr1=Hmisc::approxExtrap(y1,x,xout=unique(sol))
    invg1=approxfun(x=extr1$x,y=pmax(pmin(extr1$y,1),0),rule=2)

    p1=ifelse(qs==0,invg1(sol),ps)
    p0=ifelse(qs==1,invg0(sol),ps)
    v_res[[j]] <- p0*(1-q0) + p1*q0
  }

  data.frame(p = unsplit(p_res, f = chr), q = unsplit(q_res, f = chr), v = unsplit(v_res, f = chr))

}
