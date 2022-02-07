#' @title Perform Flexible cFDR
#'
#' @description Performs Flexible cFDR for continuous auxiliary covariates
#'
#' @details If \code{maf} is specified, then the independent SNPs will be down-sampled to match the minor allele frequency distribution. 
#'
#' @param p p-values for principal trait (vector of length n)
#' @param q continuous auxiliary data values (vector of length n)
#' @param indep_index indices of independent SNPs
#' @param res_p number of grid points in x-direction (\code{p}) for KDE estimation
#' @param res_q number of grid points in y-direction (\code{q}) for KDE estimation
#' @param nxbin number of bins in x-direction (\code{p}) for hex-binning
#' @param gridp number of data points required in a KDE grid point for left-censoring
#' @param splinecorr logical value for whether spline correction should be implemented
#' @param dist_thr distance threshold for spline correction
#' @param locfdr_df \code{df} parameter in locfdr function
#' @param plot logical value for whether to produce plots to assess KDE fit
#' @param maf minor allele frequencies for SNPs to which \code{p} and \code{q} relate (optional and used to perform MAF matching)
#' @param check_indep_cor check that sign of the correlation between \code{p} and \code{q} is the same in the independent subset as in the whole
#' @param enforce_p_q_cor if \code{p} and \code{q} are negatively correlated, flip the sign on \code{q} values
#'
#' @importFrom locfdr locfdr
#' @importFrom spatstat.geom owin
#' @importFrom fields interp.surface
#' @importFrom polyCub polyCub
#' @importFrom hexbin hexbin
#' @importFrom bigsplines bigspline
#' @importFrom grDevices contourLines
#' @import stats
#'
#' @return List of length two: (1) data.frame of p-values, q-values and v-values (2) data.frame of auxiliary data (q_low used for left censoring, how many data-points were left censored and/or spline corrected)
#' 
#' @examples
#' \donttest{
#' # this is a long running example
#'  
#' # In this example, we generate some p-values (representing GWAS p-values)
#' # and some arbitrary auxiliary data values (e.g. representing functional genomic data).
#' # We use the flexible_cfdr() function to generate v-values using default parameter values.
#' 
#' # generate p
#' set.seed(1)
#' n <- 1000
#' n1p <- 50 
#' zp <- c(rnorm(n1p, sd=5), rnorm(n-n1p, sd=1))
#' p <- 2*pnorm(-abs(zp))
#'
#' # generate q
#' mixture_comp1 <- function(x) rnorm(x, mean = -0.5, sd = 0.5)
#' mixture_comp2 <- function(x) rnorm(x, mean = 2, sd = 1)
#' q <- c(mixture_comp1(n1p), mixture_comp2(n-n1p))
#'
#' n_indep <- n
#'
#' flexible_cfdr(p, q, indep_index = 1:n_indep)
#' }
#' 
#' @export
flexible_cfdr <- function(p, q, indep_index, res_p = 300, res_q = 500, nxbin = 1000, gridp = 50, splinecorr = TRUE, dist_thr = 0.5, locfdr_df = 10, plot = TRUE, maf = NULL, check_indep_cor = TRUE, enforce_p_q_cor = TRUE){

  # match MAF distribution of independent SNPs to that of whole
  if(!is.null(maf)) {

    if(length(maf) != length(p)) {
      stop("Mismatch in lengths of p and maf vectors") 
    }

    indep_index <- match_ind_maf(maf, indep_index)
  }

  # Suitable for auxiliary covariates other than p-values
  if(check_indep_cor) {
    if(sign(cor(p[indep_index], q[indep_index], method="spearman"))!= sign(cor(p, q, method="spearman"))) {
      stop('Correlation between p and q in whole dataset has a different sign to that in independent subset of SNPs')
      }
  }

  # ensure low q enriched for low p
  if(enforce_p_q_cor) {
    if(cor(p[indep_index], q[indep_index], method="spearman") < 0) {
      q <- -q
      }
  }

  zp = -qnorm(p/2) # convert p-vals to z-scores

  # define support for KDE (range of data +/- 10%)
  q_10 <- (max(q) - min(q)) * 0.1
  zp_10 <- (max(zp) - min(zp)) * 0.1
  lims <- c(0, max(zp) + zp_10, min(q) - q_10, max(q) + q_10) # c(xl, xu, yl, yu)

  # folded normal KDE only computed for independent SNPs (so BW computation isnt biased)
  p_ind <- p[indep_index]
  zp_ind <- zp[indep_index]
  q_ind <- q[indep_index]

  # bivariate density of zp and q
  kpq <- MASS::kde2d(c(zp_ind, -zp_ind), c(q_ind, q_ind), n = c(res_p, res_q), lims = lims)

  ### estimate P(Q<=q|H0)

  # find optimal mlests parameter values
  N = length(c(zp_ind, -zp_ind)); b = 4.3 * exp(-0.26 * log(N, 10)); med = median(c(zp_ind, -zp_ind))
  sc = diff(quantile(c(zp_ind, -zp_ind))[c(2,4)])/(2*qnorm(.75))
  mlests = locmle_local(c(zp_ind, -zp_ind), xlim=c(med, b*sc)) #locfdr:::locmle(c(zp_ind, -zp_ind), xlim=c(med, b*sc))
  names(mlests) = NULL
  
  # local FDR = P(H0|ZP=zp)
  #lfdr <- tryCatch(
  #  {
  #    locfdr(c(zp_ind, -zp_ind), bre = c(kpq$x[-length(kpq$x)] + diff(kpq$x)/2, kpq$x[length(kpq$x)] + diff(kpq$x)[length(diff(kpq$x))]/2, -c(kpq$x[-length(kpq$x)] + diff(kpq$x)/2, kpq$x[length(kpq$x)] + diff(kpq$x)[length(diff(kpq$x))]/2)), mlests = c(mlests[1], b*mlests[2]), plot = 0, df = locfdr_df)
  #  },
  #  warning=function(cond) {
  #    message("Warning from locfdr:")
  #    message(cond)
#      message("...")
#     message("Examine the fit to the data")
 #     message("See locfdr documentation here: https://cran.r-project.org/web/packages/locfdr/vignettes/locfdr-example.pdf")
 #     message("...")
 #     message("Alternatively, the fcfdr::parameters_in_locfdr function can be used to output the parameters used as input in locfdr::locfdr")
 #     invokeRestart("muffleWarning")
 #   }
 # )
  
  lfdr <- withCallingHandlers(
    
    locfdr(c(zp_ind, -zp_ind), bre = c(kpq$x[-length(kpq$x)] + diff(kpq$x)/2, kpq$x[length(kpq$x)] + diff(kpq$x)[length(diff(kpq$x))]/2, -c(kpq$x[-length(kpq$x)] + diff(kpq$x)/2, kpq$x[length(kpq$x)] + diff(kpq$x)[length(diff(kpq$x))]/2)), mlests = c(mlests[1], b*mlests[2]), plot = 1, df = locfdr_df),
    
    warning = function(cnd) {
      message("Warning from locfdr\nExamine the fit to the data\nSee locfdr documentation here: https://cran.r-project.org/web/packages/locfdr/vignettes/locfdr-example.pdf\nAlternatively, the fcfdr::parameters_in_locfdr function can be used to output the parameters used as input in locfdr::locfdr\n...")
    },
    locfdr(c(zp_ind, -zp_ind), bre = c(kpq$x[-length(kpq$x)] + diff(kpq$x)/2, kpq$x[length(kpq$x)] + diff(kpq$x)[length(diff(kpq$x))]/2, -c(kpq$x[-length(kpq$x)] + diff(kpq$x)/2, kpq$x[length(kpq$x)] + diff(kpq$x)[length(diff(kpq$x))]/2)), mlests = c(mlests[1], b*mlests[2]), plot = 1, df = locfdr_df)
  )
  

  # extract lfdr values for kpq$x values
  lfdr_vals <- lfdr$mat[res_p:length(lfdr$mat[,1]),][,2]
  lfdrmat <- matrix(lfdr_vals,nrow(kpq$z),ncol(kpq$z))
  Pr.pqh0 <- lfdrmat * kpq$z # prob(p,q,H0)
  Pr.q.h0 <-  colSums(Pr.pqh0)/sum(Pr.pqh0) # prob(q|H0) = prob(q,H0)/prob(H0)
  Pr.cq.h0 <- cumsum(Pr.q.h0) # prob (Q<q | H0)
  int_kqh0 <- t(outer(Pr.cq.h0, rep(1, res_p)))
  int_kqh0 <- int_kqh0/max(int_kqh0)

  ### estimate P(P<=p, Q<=q)

  cell_size <- (diff(range(lims[1], lims[2])) / res_p) * (diff(range(lims[3], lims[4])) / res_q)
  integ <- sum(kpq$z) * cell_size
  kpq_norm <- kpq$z/integ
  int_kpq <- t(apply(apply(kpq_norm[res_p:1,], 2, cumsum), 1, cumsum))[res_p:1, ]
  int_kpq <- int_kpq/max(int_kpq)

  # kgrid estimates P(P<=p, Q<=q)/P(Q<=q|H0)
  kgrid <- kpq
  kgrid$z <- exp(log(int_kpq)-log(int_kqh0))

  # avoid 0-0 errors
  kgrid$z[which(kgrid$z==0)] <- min(kgrid$z[which(kgrid$z>0)])
  
  if(plot == TRUE){
    hist(q_ind, freq = FALSE, xlab = "q", main = "Histogram of q with\nestimated density in red")
    lines(kpq$y, kpq_norm[1,], col =  "red")
    
    image(kpq, xlab = "Principal trait Z-scores", ylab = "q", main = "Estimated density from 2D KDE")
  }

  # cgrid is grid of cFDR values
  # estimated by p/kgrid
  cgrid <- vector(mode = "list", length = 3)
  cgrid$x = c(seq(0, 0.05, length.out = 101)[1:100], seq(0.05, lims[2], length.out = res_p))
  cgrid$y = seq(lims[3], lims[4], length.out = res_q)
  ptest = 2*pnorm(-cgrid$x)

  cfdrs <- matrix(rep(0, length(cgrid$y)*length(cgrid$x)), ncol = length(cgrid$y))

  for (i in 1:length(cgrid$y)) {
    xdenom=interp.surface(kgrid,cbind(cgrid$x, rep(cgrid$y[i], length(cgrid$x))))
    cfdrs[,i]=cummin(ptest/xdenom)
  }

  cgrid$z <- cfdrs

  # left-censoring
  q_grid <- seq(lims[3], lims[4], length.out = res_q)
  groups <- cut(q_ind, breaks = q_grid, include.lowest = TRUE, right = TRUE)
  # use left point of region where the number of points used for estimation in that bin is >gridp
  q_low <- q_grid[which(table(groups)>gridp)[1]]

  q_cens <- q
  q_cens[which(q_cens < q_low)] <- q_low

  # bivariate binning
  binned_res <- hexbin(zp, q_cens, xbins = nxbin, IDs = TRUE)

  bins <- cbind(binned_res@xcm, binned_res@ycm)

  ccut <- interp.surface(cgrid, bins) # cFDR values for binned data

  # L-curves are contours of cFDR curves
  cl <- lapply(ccut, function(l) grDevices::contourLines(x=cgrid, levels=l))

  cl_basic <- cl

  # rearrage L-curves so they are defined anticlockwise
  lengths <- lapply(cl, length)
  for(i in which(lengths>1)){
    first = unlist(lapply(cl[[i]], function(x) x$y[1]))
    cl[[i]] = cl[[i]][order(first)]
  }

  # remove contour coords in problematic low q regions
  for(i in which(lengths>1)){
      cl[[i]][[1]] <- NULL
    }

  # define newx and newy as joined segments
  for(i in 1:length(cl)){
    cl[[i]]$newx <- c()
    for(j in 1:length(cl[[i]])){
      cl[[i]]$newx <- c(cl[[i]]$newx, cl[[i]][[j]]$x)
      cl[[i]]$newy <- c(cl[[i]]$newy, cl[[i]][[j]]$y)
    }
  }

  Lx <- lapply(cl, function(x) 2*pnorm(-abs(x$newx)))
  Ly <- lapply(cl, function(x) x$newy)

  ### estimate P(P=p,Q=q|H0)=P(P=q|H0)*P(Q=q|H0)

  fph0 = function(x) dunif(x) # P|H0
  fqh0 <- approxfun(kpq$y, Pr.q.h0) # Q|H0

  fpqh0 <- function(s){
    fph0(s[,1])*fqh0(s[,2])
  }

  # normalise so integral is 1 over full region
  fullint <- polyCub(owin(poly=list(x = c(1, 0, 0, 1), y = c(lims[4], lims[4], lims[3], lims[3])), check = FALSE, fix = FALSE), fpqh0, method = "midpoint")

  fpqh0 <- function(s){
    fph0(s[,1])*fqh0(s[,2])/fullint
  }

  # integrate bivariate null over L region to get v-vals
  v_tmp <- mapply(function(X, Y){
    polyCub(owin(poly=list(x = c(X[length(X)], 0,0, c(X[1],X)), y = c(lims[4], lims[4], lims[3], c(lims[3],Y))), check = FALSE, fix  = FALSE), fpqh0, method = "midpoint")
  }, X = Lx, Y = Ly)

  # map back binned data to original data points
  v <- v_tmp[match(binned_res@cID, binned_res@cell)]

  v_b4spline <- v

  v[which(v==0)] <- 1e-100 # replace any 0 v-vals

  if(splinecorr == TRUE){ # spline correction

    spline_fit <- bigspline(x = q, y = log10(v/p), nknots = 5, rparm = NA)
    #pred_out <- predict.bigspline(spline_fit, newdata = seq(min(q), max(q), 0.05))

    distances <- abs(log10(v/p)-spline_fit$fitted.values)

    corrected_ind <- which(abs(log10(v/p)-spline_fit$fitted.values) > dist_thr)

    v[corrected_ind] <- pmin(10^spline_fit$fitted[corrected_ind]*p[corrected_ind], 1)

  }

  v[which(v>1)] <- 1 # fix bug where some v vals = 1 + 2.220446e-16
  
  # print warning if v-values have changed too much as something has likely gone wrong
  if( median(v) < 0.8*median(p) | median(v) > 1.2*median(p) ){
    
    warning('v-values very different to input p-values - check results (if q has a long tail then try left censoring)')
    
    }

  df <- data.frame(p, q, v)

  if(splinecorr == TRUE){
    
    return(list(df, data.frame(q_low = q_low, left_cens = length(which(q < q_low)), splinecorr = length(corrected_ind))))
    
  } else
    
    return(list(df, data.frame(q_low = q_low, left_cens = length(which(q < q_low)), splinecorr = NA)))
  
}

#' @title Function to downsample independent SNPs to match MAF distribution of whole set.
#'
#' @description Matches MAF distribution of independent set of SNPs to MAF distribution of whole set of SNPs to avoid MAF-based confounding.
#'
#' @details Must supply maf values from the whole data set, not just the independent SNPs.
#'
#' @param maf minor allele frequencies of (all) SNPs  
#' @param indep_index indices of independent SNPs
#'
#' @return indices of independent SNP in chosen in sample
match_ind_maf <- function(maf, indep_index) {
  breaks <- seq(0, 0.5, length=51)
  
  daf <- data.frame(indep_index = indep_index, maf = maf[indep_index])
  
  maf_interval <- as.character(cut(maf, breaks = breaks, include.lowest = T))
  
  daf$maf_interval <- maf_interval[indep_index]
  
  maf_interval_freq.whole <- table(maf_interval)
  
  maf_interval_freq.ind <- table(factor(daf$maf_interval,levels=unique(maf_interval)))[names( maf_interval_freq.whole)]
  
  maf_interval_freq.whole.relative <- maf_interval_freq.whole/sum(maf_interval_freq.whole)
  
  maf_interval_freq.ind.relative <- maf_interval_freq.ind/sum(maf_interval_freq.ind)
  
  w=which(maf_interval_freq.ind>0)
  
  max_sample_size <- floor(min(maf_interval_freq.ind[w]/maf_interval_freq.whole.relative[w]))
  
  scaled_interval_sample_sizes <- floor(maf_interval_freq.whole.relative*max_sample_size)
  
  indep_sample_index <- unlist(lapply(names(scaled_interval_sample_sizes), function(x) sample(subset(daf, maf_interval==x)$indep_index, size=scaled_interval_sample_sizes[x])))
  
  indep_sample_index
}

#' parameters_in_locfdr
#'
#' @param p p values for principal trait (vector of length n)
#' @param q continuous auxiliary data values (vector of length n)
#' @param indep_index indices of independent SNPs
#' @param res_p resolution for p
#' @param res_q resolution for q
#' @param maf minor allele frequencies for SNPs to which \code{p} and \code{q} relate (optional and used to perform MAF matching)
#' @param check_indep_cor check that sign of the correlation between \code{p} and \code{q} is the same in the independent subset as in the whole
#' @param enforce_p_q_cor if \code{p} and \code{q} are negatively correlated, flip the sign on \code{q} values
#'
#' @return list of values used as input into \code{locfdr::locfdr} function intrinsically in \code{flexible_cfdr}
#' 
#' @examples
#' 
#' # In this example, we generate some p-values (representing GWAS p-values)
#' # and some arbitrary auxiliary data values (e.g. representing functional genomic data).
#' # We use the parameters_in_locfdr() function to extract the parameters estimated by
#' # the locfdr function.
#' 
#' # generate p
#' set.seed(1)
#' n <- 1000
#' n1p <- 50 
#' zp <- c(rnorm(n1p, sd=5), rnorm(n-n1p, sd=1))
#' p <- 2*pnorm(-abs(zp))
#'
#' # generate q
#' mixture_comp1 <- function(x) rnorm(x, mean = -0.5, sd = 0.5)
#' mixture_comp2 <- function(x) rnorm(x, mean = 2, sd = 1)
#' q <- c(mixture_comp1(n1p), mixture_comp2(n-n1p))
#'
#' n_indep <- n
#'
#' parameters_in_locfdr(p, q, indep_index = 1:n_indep)
#' 
#' 
#' @export
#'
parameters_in_locfdr <- function(p, q, indep_index, res_p = 300, res_q = 500, maf = NULL, check_indep_cor = TRUE, enforce_p_q_cor = TRUE){
  
  # match MAF distribution of independent SNPs to that of whole
  if(!is.null(maf)) {
    
    if(length(maf) != length(p)) {
      stop("Mismatch in lengths of p and maf vectors") 
    }
    
    indep_index <- match_ind_maf(maf, indep_index)
  }
  
  # Suitable for auxiliary covariates other than p-values
  if(check_indep_cor) {
    if(sign(cor(p[indep_index], q[indep_index], method="spearman"))!= sign(cor(p, q, method="spearman"))) {
      stop('Correlation between p and q in whole dataset has a different sign to that in independent subset of SNPs')
    }
  }
  
  # ensure low q enriched for low p
  if(enforce_p_q_cor) {
    if(cor(p[indep_index], q[indep_index], method="spearman") < 0) {
      q <- -q
    }
  }
  
  zp = -qnorm(p/2) # convert p-vals to z-scores
  
  # define support for KDE (range of data +/- 10%)
  q_10 <- (max(q) - min(q)) * 0.1
  zp_10 <- (max(zp) - min(zp)) * 0.1
  lims <- c(0, max(zp) + zp_10, min(q) - q_10, max(q) + q_10) # c(xl, xu, yl, yu)
  
  # folded normal KDE only computed for independent SNPs (so BW computation isnt biased)
  p_ind <- p[indep_index]
  zp_ind <- zp[indep_index]
  q_ind <- q[indep_index]
  
  # bivariate density of zp and q
  kpq <- MASS::kde2d(c(zp_ind, -zp_ind), c(q_ind, q_ind), n = c(res_p, res_q), lims = lims)
  
  ### estimate P(Q<=q|H0)
  
  # find optimal mlests parameter values
  N = length(c(zp_ind, -zp_ind)); b = 4.3 * exp(-0.26 * log(N, 10)); med = median(c(zp_ind, -zp_ind))
  sc = diff(quantile(c(zp_ind, -zp_ind))[c(2,4)])/(2*qnorm(.75))
  mlests = locmle_local(c(zp_ind, -zp_ind), xlim=c(med, b*sc)) #locfdr:::locmle(c(zp_ind, -zp_ind), xlim=c(med, b*sc))
  names(mlests) = NULL
  
  return(list("zz" = c(zp_ind, -zp_ind), "bre" = c(kpq$x[-length(kpq$x)] + diff(kpq$x)/2, kpq$x[length(kpq$x)] + diff(kpq$x)[length(diff(kpq$x))]/2, -c(kpq$x[-length(kpq$x)] + diff(kpq$x)/2, kpq$x[length(kpq$x)] + diff(kpq$x)[length(diff(kpq$x))]/2)), "mlests" = c(mlests[1], b*mlests[2])))
}

