#' Function to perform cFDR for continuous auxiliary covariates
#'
#' @param p p values for principal trait (vector of length n)
#' @param q continuous auxillary data values (vector of length n)
#' @param indep_index indices of independent SNPs
#' @param nxbin number of bins in zp direction for hex-binning
#' @param res_p number of test points in x-direction (p)
#' @param res_q number of test points in y-direction (q)
#' @param gridp number of data points required in a KDE grid point for left-censoring
#' @param splinecorr logical value for whether spline correction should be implemented
#' @param dist_thr distance threshold for spline correction
#'
#' @rawNamespace import(dplyr, except = c(filter, lag))
#' @rawNamespace import(data.table, except = c(last, first, between, shift))
#' @rawNamespace import(MASS, except = c(select, area))
#' @import locfdr
#' @import spatstat
#' @import cfdr
#' @import fields
#' @import stats
#' @import polyCub
#' @import hexbin
#' @import bigsplines
<<<<<<< HEAD
#' @import grDevices
=======
#' @imports data.table
#' @imports MASS
#' @imports grDevices
#' 
#' @examples see vignette
>>>>>>> 9f4d11885d8c88bf477fe8dd9ee1195da98da514
#'
#' @return list of length two: (1) dataframe of p-values, q-values and v-values (2) dataframe of auxiliary data (q_low used for left censoring, how many data-points were left censored and/or spline corrected)
#' @export
flexible_cfdr <- function(p, q, indep_index, nxbin = 1000, res_p = 300, res_q = 500, gridp = 50, splinecorr = TRUE, dist_thr = 0.5){

  if( sign(cor(p[indep_index], q[indep_index], method="spearman"))!= sign(cor(p, q, method="spearman")) ) stop('Correlation between p and q in whole dataset has a different sign to that in independent subset of SNPs')

  # ensure low q enriched for low p
  if(cor(p[indep_index], q[indep_index], method="spearman") < 0) q <- -q

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
  mlests = locfdr:::locmle(c(zp_ind, -zp_ind), xlim=c(med, b*sc))
  names(mlests) = NULL

  # local FDR = P(H0|ZP=zp)
  lfdr <- locfdr(c(zp_ind, -zp_ind), bre = c(kpq$x[-length(kpq$x)] + diff(kpq$x)/2, kpq$x[length(kpq$x)] + diff(kpq$x)[length(diff(kpq$x))]/2, -c(kpq$x[-length(kpq$x)] + diff(kpq$x)/2, kpq$x[length(kpq$x)] + diff(kpq$x)[length(diff(kpq$x))]/2)), mlests = c(mlests[1], b*mlests[2]), plot = 0, df = 10)

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
    first = lapply(cl[[i]], function(x) x$y[1]) %>% unlist()
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

  df <- data.frame(p, q, v)

  return(list(df, data.frame(q_low = q_low, left_cens = length(which(q < q_low)), splinecorr = length(corrected_ind))))
}
