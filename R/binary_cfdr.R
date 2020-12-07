#' Function to perform cFDR for binary auxiliary covariates
#'
#' @param p p values for principal trait (vector of length n)
#' @param q binary auxillary data values (vector of length n)
#' @param chr chromosome each SNP resides (vector of length n)
#' @param p.thr p-value threshold for interpolation based on digits_highp and digits_lowp
#' @param digits_highp digits for rounding for p >= p.thr
#' @param digits_lowp digits for rounding for p < p.thr
#'
#' @return dataframe of p, q and v values
#' @export
#'
binary_cfdr <- function(p, q, chr, p.thr = 0.01, digits_highp = 1, digits_lowp = 2){
  
  unique_chr <- unique(chr)
  
  # split p and q into chromosomes
  p_res <- split(p, f = chr)
  q_res <- split(q, f = chr)
  
  # prepare container for v
  v_res <- vector(mode = "list", length = length(unique_chr))
  
  for(j in 1:length(unique_chr)){
    chrom <- unique_chr[j]
    p_loo <- p[-which(chr == chrom)] # df leave-one-chromosome-out
    q_loo <- q[-which(chr == chrom)] # df leave-one-chromosome-out
    ps <- c(p_res[[j]], 1, 1) # add p = 1 to prevent interpolating for p>ps[use]
    qs <- c(q_res[[j]], 0, 1)
    
    ## sample values to fit the function
    
    # sample less in p > p.thr, more in p < p.thr
    lps = ifelse(ps<p.thr, round(log(ps), digits_lowp), round(log(ps), digits_highp))
    use = c(ifelse(qs[-c(length(ps)-1, length(ps))]==0, !duplicated(lps * (qs==0)), !duplicated(lps * (qs!=0))), TRUE, TRUE) # add p = 1 to prevent interpolating for p>ps[use]
    
    q0 <- length(which(q_loo == 1 & p_loo > 0.5))/length(which(p_loo > 0.5))
    mult <- (length(which(q_loo == 0 & p_loo > 1/2))/length(which(q_loo == 1 & p_loo > 1/2)))
    
    q0_sol=sapply(ps, function(p) max(sum(p_loo <= p & q_loo==0),1))
    q1_sol=sapply(ps, function(p) max(sum(p_loo <= p & q_loo==1),1))
    
    sol <- ifelse(qs==0,
                  mult*ps / q0_sol,
                  (1/mult)*ps / q1_sol)
    
    sol.use=sol[use]
    
    f <- function(p, q, c) (p/max(1,length(which(p_loo <=  p & q_loo == q)))) - c
    tols=ifelse(ps<1e-4,ps,.Machine$double.eps^0.25)
    p1 <- ps
    for(i in which(qs==0 & use)){
      p1[i] <- min(1,uniroot(f, c(0,1), q = 1, c = sol[i], extendInt = "yes",
                             tol = tols[i],
                             maxiter = 5000)$root)
    }
    p0 <- ps
    for(i in which(qs==1)){
      p0[i] <- min(1,uniroot(f, c(0,1), q = 0, c = sol[i], extendInt = "yes",
                             tol = tols[i],
                             maxiter = 5000)$root)
    }
    
    
    f = approxfun(ps[use & qs==0],p1[use & qs==0])
    p1[!use & qs==0] = f(ps[!use & qs==0])
    f = approxfun(ps[use & qs==1],p0[use & qs==1])
    p0[!use & qs==1] = f(ps[!use & qs==1])
    
    v_res[[j]] <- p0*(1-q0) + p1*q0
    # remove the v-values for the artificially added p=1 points
    v_res[[j]] <- v_res[[j]][-c(length(ps)-1, length(ps))]
    
    ## fix kinks - fails - need to depend on qs too
    # ord=order(v_res[[j]])
    # revord=numeric(length(ord))
    # revord[ord]=seq_along(ord)
    # v_res[[j]] = cummax(v_res[[j]][ord])[revord]
  }
  data.frame(p = unsplit(p_res, f = chr), q = unsplit(q_res, f = chr), v = unlist(v_res))
}

#' Function to perform cFDR for binary auxiliary covariates
#'
#' @param p p values for principal trait (vector of length n)
#' @param q binary auxillary data values (vector of length n)
#' @param chr chromosome each SNP resides (vector of length n)
#'
#' @return dataframe of p, q and v values
#' @export
#'
faster_binary_cfdr <- function(p, q, chr,digits=2){

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
    extr=approxExtrap(y,x,xout=unique(sol))
    invg0=approxfun(x=extr$x,y=pmax(pmin(extr$y,1),0),rule=2)
    ## invg0_spline=splinefun(x=y,y=x,method="hyman")

    ## approx g1
    y1=x/sapply(x, function(p) max(sum(p_loo <= p & q_loo==1),1))
    extr1=approxExtrap(y1,x,xout=unique(sol))
    invg1=approxfun(x=extr1$x,y=pmax(pmin(extr1$y,1),0),rule=2)

    p1=ifelse(qs==0,invg1(sol),ps)
    p0=ifelse(qs==1,invg0(sol),ps)
    v_res[[j]] <- p0*(1-q0) + p1*q0
  }

  data.frame(p = unsplit(p_res, f = chr), q = unsplit(q_res, f = chr), v = unlist(v_res))

}
