#' Function to perform cFDR for binary auxiliary covariates
#'
#' @title binary_cfdr
#' @param p p values for principal trait (vector of length n)
#' @param q binary auxiliary data values (vector of length n)
#' @param chr SNP chromosomes (vector of length n)
#'
#' @return dataframe of p, q and v values
#' @export
#'
binary_cfdr <- function(p, q, chr){

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
    ps <- p_res[[j]]
    qs <- q_res[[j]]

    q0 <- length(which(q_loo == 1 & p_loo > 0.5))/length(which(p_loo > 0.5))

    mult <- (length(which(q_loo == 0 & p_loo > 1/2))/length(which(q_loo == 1 & p_loo > 1/2)))

    sol <- ifelse(qs==0,
                  mult*ps / (max(1, length(which(p_loo <=  ps & q_loo == 0)))),
                  (1/mult)*ps / (max(1, length(which(p_loo <=  ps & q_loo == 1)))))

    f <- function(p, q, c) (p/max(1,length(which(p_loo <=  ps & q_loo == q)))) - c

    p1 <- ps
    for(i in which(qs==0)){
      p1[i] <- ifelse(ps[i] < 1e-04,
                      min(1,uniroot(f, c(0, 1), q = 1, c = sol[i], extendInt = "yes", tol = ps[i], maxiter = 5000)$root),
                      min(1,uniroot(f, c(0, 1), q = 1, c = sol[i], extendInt = "yes", maxiter = 5000)$root))
    }

    p0 <- ps
    for(i in which(qs==1)){
      p0[i] <- ifelse(ps[i] < 1e-04,
                      min(1,uniroot(f, c(0, 1), q = 0, c = sol[i], extendInt = "yes", tol = ps[i], maxiter = 5000)$root),
                      min(1,uniroot(f, c(0, 1), q = 0, c = sol[i], extendInt = "yes", maxiter = 5000)$root))
    }

    v_res[[j]] <- p0*(1-q0) + p1*q0

  }

  data.frame(p = unsplit(p_res, f = chr), q = unsplit(q_res, f = chr), v = unlist(v_res))

}
