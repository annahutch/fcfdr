test_that("flexible_cfdr runs to completion on a smaller version of the vignette data", {
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

  res <- flexible_cfdr(p, q, indep_index = seq(1, n, 1), match_maf=F)

  expect_equal(digest::digest(res),"b6f8cdd80fb1bdbc4e315a827698d3b9")
})

test_that("match_ind_maf runs on package sample data", {
  # Replace missing
  load(file=system.file('data', 'df.rda', package='fcfdr', mustWork = T))

  set.seed(1)

  indep_index_sample <- match_ind_maf(maf = df$maf, indep_index = which(df$ldak_weight != 0))

  expect_equal(digest::digest(indep_index_sample), "2f92b79b76d7584998fedfccf332e0a0")
})
