test_that("flexible_cfdr runs on a smaller version of the vignette data", {
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

  res <- flexible_cfdr(p, q, indep_index = seq(1, n, 1))

  expect_equal(digest::digest(res),"b6f8cdd80fb1bdbc4e315a827698d3b9")
})

test_that("flexible_cfdr runs on a smaller version of the vignette data with MAF matching with no difference in MAF distribution", {
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

  maf <- runif(n, min=0, max=0.5)

  res <- flexible_cfdr(p, q, indep_index = 1:n_indep, maf = maf)
  
  expect_equal(digest::digest(res), "fe4a8d075a91a45ac7871ef67d833714")
})

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

test_that("flexible_cfdr runs on subset of package sample data when matching MAF and using RA p-value as the auxiliary covariate", {
  load(file=system.file('data', 'df.rda', package='fcfdr', mustWork = T))

  set.seed(42)

  res <- flexible_cfdr(df$p, q = -qnorm(df$RA_p/2), indep_index = which(df$ldak_weight != 0), maf = df$maf, check_indep_cor = FALSE, enforce_p_q_cor = FALSE)
  
  expect_equal(digest::digest(res), "f07705b80f5fc87fbed00872158046ad")
})

test_that("match_ind_maf runs on package sample data", {
  load(file=system.file('data', 'df.rda', package='fcfdr', mustWork = T))

  set.seed(1)

  indep_index_sample <- match_ind_maf(maf = df$maf, indep_index = which(df$ldak_weight != 0))

  expect_equal(digest::digest(indep_index_sample), "b5108cc5933a7c35f30aa718d30d93d2")
})

