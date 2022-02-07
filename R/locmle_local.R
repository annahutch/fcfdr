#' @title locmle_local
#'
#' @description Verbatim copy of locmle from locfdr package
#'
#' @details Imported verbatim from locfdr package; reproduced here to avoid ::: use. For any documentation, see package locfdr on CRAN. Not exported
#'
#' @param z null
#' @param xlim null
#' @param Jmle null
#' @param d null
#' @param s null
#' @param ep null
#' @param sw null
#' @param Cov.in null
#' @return null
#' @noRd
locmle_local=function (z, xlim, Jmle = 35, d = 0, s = 1, ep = 1/1e+05, sw = 0, 
          Cov.in) 
{
  N = length(z)
  if (missing(xlim)) {
    if (N > 5e+05) 
      b = 1
    else b = 4.3 * exp(-0.26 * log(N, 10))
    xlim = c(median(z), b * diff(quantile(z)[c(2, 4)])/(2 * 
                                                          qnorm(0.75)))
  }
  aorig = xlim[1] - xlim[2]
  borig = xlim[1] + xlim[2]
  z0 = z[which(z >= aorig & z <= borig)]
  N0 = length(z0)
  Y = c(mean(z0), mean(z0^2))
  that = N0/N
  for (j in 1:Jmle) {
    bet = c(d/s^2, -1/(2 * s^2))
    aa = (aorig - d)/s
    bb = (borig - d)/s
    H0 = pnorm(bb) - pnorm(aa)
    fa = dnorm(aa)
    fb = dnorm(bb)
    H1 = fa - fb
    H2 = H0 + aa * fa - bb * fb
    H3 = (2 + aa^2) * fa - (2 + bb^2) * fb
    H4 = 3 * H0 + (3 * aa + aa^3) * fa - (3 * bb + bb^3) * 
      fb
    H = c(H0, H1, H2, H3, H4)
    r = d/s
    I = matrix(rep(0, 25), 5)
    for (i in 0:4) I[i + 1, 0:(i + 1)] = choose(i, 0:i)
    u1 = s^(0:4)
    II = pmax(row(I) - col(I), 0)
    II = r^II
    I = u1 * (I * II)
    E = as.vector(I %*% H)/H0
    E1 = E[2]
    E2 = E[3]
    E3 = E[4]
    E4 = E[5]
    mu = c(E1, E2)
    V = matrix(c(E2 - E1^2, E3 - E1 * E2, E3 - E1 * E2, E4 - 
                   E2^2), 2)
    bett = bet + solve(V, Y - mu)/(1 + 1/j^2)
    if (bett[2] > 0) 
      bett = bet + 0.1 * solve(V, Y - mu)/(1 + 1/j^2)
    if (is.na(bett[2])) 
      break
    else if (bett[2] >= 0) 
      break
    d = -bett[1]/(2 * bett[2])
    s = 1/sqrt(-2 * bett[2])
    if (sum((bett - bet)^2)^0.5 < ep) 
      break
  }
  if (is.na(bett[2])) {
    mle = rep(NA, 6)
    Cov.lfdr = NA
    Cor = matrix(NA, 3, 3)
  }
  else if (bett[2] >= 0) {
    mle = rep(NA, 6)
    Cov.lfdr = Cov.out = NA
    Cor = matrix(NA, 3, 3)
  }
  else {
    aa = (aorig - d)/s
    bb = (borig - d)/s
    H0 = pnorm(bb) - pnorm(aa)
    p0 = that/H0
    J = s^2 * matrix(c(1, 0, 2 * d, s), 2)
    JV = J %*% solve(V)
    JVJ = JV %*% t(J)
    mat2 = cbind(0, JVJ/N0)
    mat1 = c((p0 * H0 * (1 - p0 * H0))/N, 0, 0)
    mat = rbind(mat1, mat2)
    h = c(H1/H0, (H2 - H0)/H0)
    matt = c(1/H0, -(p0/s) * t(h))
    matt = rbind(matt, cbind(0, diag(2)))
    C = matt %*% (mat %*% t(matt))
    mle = c(p0, d, s, diag(C)^0.5)
    if (sw == 1) {
      sd = mle[4:6]
      Co = C/outer(sd, sd)
      dimnames(Co) = list(c("p0", "d", "s"), c("p0", "d", 
                                               "s"))
      Cor = Co[c(2, 3, 1), c(2, 3, 1)]
    }
    if (!missing(Cov.in)) {
      i0 = which(Cov.in$x > aa & Cov.in$x < bb)
      Cov.out = loccov_local(N, N0, p0, d, s, Cov.in$x, Cov.in$X, 
                       Cov.in$f, JV, Y, i0, H, h, Cov.in$sw)
    }
  }
  names(mle) = c("p0", "del0", "sig0", "sd.p0", "sd.del0", 
                 "sd.sig0")
  mle = mle[c(2, 3, 1, 5, 6, 4)]
  out = list(mle = mle)
  if (sw == 1) {
    Cor = list(Cor = Cor)
    out = c(out, Cor)
  }
  if (!missing(Cov.in)) {
    if (Cov.in$sw == 2) {
      pds. = list(pds. = Cov.out)
      out = c(out, pds.)
    }
    else if (Cov.in$sw == 3) {
      Ilfdr = list(Ilfdr = Cov.out)
      out = c(out, Ilfdr)
    }
    else {
      Cov.lfdr = list(Cov.lfdr = Cov.out)
      out = c(out, Cov.lfdr)
    }
  }
  if ((sw == 1) | !missing(Cov.in)) 
    return(out)
  else return(mle)
}


#' @title loccov_local
#'
#' @description Verbatim copy of loccov from locfdr package
#'
#' @details Imported verbatim from locfdr package; reproduced here to avoid ::: use. For any documentation, see package locfdr on CRAN. Not exported
#'
#' @param N NULL
#' @param N0 NULL
#' @param p0 NULL
#' @param d NULL
#' @param s NULL
#' @param x NULL
#' @param X NULL
#' @param f NULL
#' @param JV NULL
#' @param Y NULL
#' @param i0 NULL
#' @param H NULL
#' @param h NULL
#' @param sw NULL
#' @return NULL
#' @noRd
loccov_local=function (N, N0, p0, d, s, x, X, f, JV, Y, i0, H, h, sw) 
{
  M = rbind(1, x - Y[1], x^2 - Y[2])
  if (sw == 2) {
    K = length(x)
    K0 = length(i0)
    toprow = c(1 - N0/N, -t(h) %*% JV/s)
    botrow = cbind(0, JV/p0)
    mat = rbind(toprow, botrow)
    M0 = M[, i0]
    dpds.dy0 = mat %*% M0/N/H[1]
    dy0.dy = matrix(0, K0, K)
    dy0.dy[, i0] = diag(1, K0)
    dpds.dy = dpds.dy0 %*% dy0.dy
    rownames(dpds.dy) = c("p", "d", "s")
    return(dpds.dy)
  }
  else {
    xstd = (x - d)/s
    U = cbind(xstd - H[2]/H[1], xstd^2 - H[3]/H[1])
    M[, -i0] = 0
    dl0plus.dy = cbind(1 - N0/N, U %*% JV/s) %*% M/N/H[1]/p0
    G <- t(X) %*% (f * X)
    dl.dy = X %*% solve(G) %*% t(X)
    dlfdr.dy = dl0plus.dy - dl.dy
    if (sw == 3) 
      return(dlfdr.dy)
    else {
      Cov.lfdr = dlfdr.dy %*% (f * t(dlfdr.dy))
      return(Cov.lfdr)
    }
  }
}