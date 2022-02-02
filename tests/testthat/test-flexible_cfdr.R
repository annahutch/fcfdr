test_that("flexible_cfdr throws an error with maf/p length mismatch", {
  set.seed(1)
  n <- 5000
  n1p <- 50 
  zp <- c(rnorm(n1p, sd=5), rnorm(n-n1p, sd=1))
  p <- 2*pnorm(-abs(zp))

  mixture_comp1 <- function(x) rnorm(x, mean = -0.5, sd = 0.5)
  mixture_comp2 <- function(x) rnorm(x, mean = 2, sd = 1)
  n <- length(p)
  z <- runif(n)

  q <- c(mixture_comp1(n1p), mixture_comp2(n-n1p))

  n_indep <- 2000

  maf <- runif(n-1000, min=0, max=0.5)
  
  expect_error(flexible_cfdr(p, q, indep_index = 1:n_indep, maf = maf), "Mismatch in lengths of p and maf vectors")
})

test_that("flexible_cfdr throws an error with flipped sign of principal-auxiliary covariate correlation between whole and independent (sub)sets", {
  set.seed(1)

  p <- numeric(5000)
  q <- numeric(5000)
  p[1:4000] <- runif(4000, min = 0, max=0.5)
  q[1:4000] <- 2*p[1:4000]

  p[4001:5000] <- runif(1000)
  # If X is uniform, then 1-X is negatively correlated with X
  q[4001:5000] <- 1-p[4001:5000]

  q <- -qnorm(q/2)
  
  expect_error(flexible_cfdr(p, q, indep_index = 4001:5000), "Correlation between p and q in whole dataset has a different sign to that in independent subset of SNPs")
})


test_that("match_ind_maf runs on package sample data", {
  load(file=system.file('testdata', 'T1D_df.RData', package='fcfdr', mustWork = T))

  set.seed(1)
  
  indep_index_sample <- match_ind_maf(maf = T1D_df$MAF, indep_index = which(T1D_df$LDAK_weight != 0))
  
  expect_equal(digest::digest(indep_index_sample), "0de885c64aff8809e0fc6e606a15d3be")
})


