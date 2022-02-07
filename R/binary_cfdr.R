#' Perform cFDR leveraging binary auxiliary covariates
#'
#' @param p p-values for principal trait (vector of length n)
#' @param q binary auxiliary data values (vector of length n)
#' @param group group membership of each SNP for leave-one-out procedure (vector of length n) (e.g. chromosome number or LD block)
#'
#' @importFrom Hmisc approxExtrap
#' 
#' @return data.frame of p, q and v values
#' 
#' @examples
#' 
#' # In this example, we generate some p-values (representing GWAS p-values)
#' # and some arbitrary auxiliary data values (e.g. representing functional genomic data).
#' # We use the parameters_in_locfdr() function to extract the parameters estimated by
#' # the locfdr function.
#' 
#' # generate p
#' set.seed(2)
#' n <- 1000
#' n1p <- 50 
#' zp <- c(rnorm(n1p, sd=5), rnorm(n-n1p, sd=1))
#' p <- 2*pnorm(-abs(zp))
#'
#' # generate q
#' q <- rbinom(n, 1, 0.1)
#'
#' group <- c(rep("A", n/2), rep("B", n/2)) 
#' 
#' binary_cfdr(p, q, group)
#' 
#' @export
#'
binary_cfdr <- function(p, q, group){

  unique_group <- unique(group)

  # split p and q into groups
  p_res <- split(p, f = group)
  q_res <- split(q, f = group)
  minp=min(p)
  maxp=max(p)

  # prepare container for v
  v_res <- vector(mode = "list", length = length(unique_group))

  ## prepare x for approxfun
  logx=seq(log10(minp),log10(maxp),length.out=1000)
  x=c(exp(logx),1)

  for(j in 1:length(unique_group)){

    this_group <- unique_group[j]

    p_loo <- p[-which(group == this_group)] # df leave-one-out
    q_loo <- q[-which(group == this_group)] # df leave-one-out
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

  p = unsplit(p_res, f = group)
  q = unsplit(q_res, f = group)
  v = unsplit(v_res, f = group)
  
  # correct plateau for low p in applications where p,q are very weakly correlated
  
  if(abs(cor(p,q))<0.01){
    
    # identify problematic points
    ind <- which( p < 0.01 & (v < 0.8*p | v > 1.2*p ))
    data_bad <- data.frame(p = p[ind], v = v[ind], q = q[ind])
    data_good <- data.frame(p = p[-ind], v = v[-ind], q = q[-ind])
    
    # find straight lines from data 
    lmout_q_0 <- lm(v~p-1, data = data_good[which(data_good$q==0),])
    lmout_q_1 <- lm(v~p-1, data = data_good[which(data_good$q==1),])
    
    # replace problematic points
    v[ind] <- ifelse(data_bad$q==0, predict(lmout_q_0, data.frame(p = data_bad$p)), predict(lmout_q_1, data.frame(p = data_bad$p)))
    
    if(ind > length(p)*0.5) warning("p,q have low correlation and >50% of v-values may be problematic - check results")
    
    
  }
  
  data.frame(p, q, v)

}
