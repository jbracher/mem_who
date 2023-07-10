# Compute approximation of expected thresholds based on empirical mean vector and covariance matrix.

# set working directory
current_path = rstudioapi::getActiveDocumentContext()$path # get path of this file
setwd(dirname(current_path))

# analytical approximation of expected thresholds base on mean vector and covariance matrix.
approximate_expectations <- function(mu, Sigma, m, n){
  mean_ana <- mean(mu)
  var_ana <- m/(m*n - 1)*sum(mu^2 + diag(Sigma)) - (1/(n*(m*n - 1)) * sum(Sigma) + m/(n*(m*n - 1))*sum(mu)^2)
  sd_ana <- sqrt(var_ana)
  threshold_0.4 <- mean_ana - 0.25*sd_ana
  threshold_0.9 <- mean_ana + 1.28*sd_ana
  threshold_0.975 <- mean_ana + 1.96*sd_ana
  
  return(c(mean = mean_ana,
           var = var_ana,
           sd = sd_ana,
           medium_0.4 = threshold_0.4,
           high_0.9 = threshold_0.9,
           very_high_0.975 = threshold_0.975))
}

# analytical approximation of expected thresholds based on data in mem format from which mean
# and covariance matrix are estimated.
approximate_expectations_from_data0 <- function(dat, log = FALSE, m, n){
  
  # sort:
  for(i in 1:ncol(dat)) dat[, i] <- sort(dat[, i], decreasing = TRUE)
  
  sub <- dat[1:n, ]
  sub_t <- t(as.matrix(sub))
  if(log) sub_t <- log(sub_t)
  
  mu <- colMeans(sub_t)
  Sigma <- cov(sub_t)
  
  res <-approximate_expectations(mu = mu, Sigma = Sigma, m = m, n = n)
  
  if(log) res <- exp(res)
  
  return(res)
}

# analytical approximation of expected thresholds based on data in mem format from which mean
# and covariance matrix are estimated. Wrapper around approximate_expectations_from_data0
# to apply to several values of n and K.
approximate_expectations_from_data <- function(dat, log = FALSE, m, n){
  ret <- NULL
  for(i in seq_along(n)){
    res_temp <- approximate_expectations_from_data0(dat = dat, log = log, m = m[i], n = n[i])
    res_temp <- c(m = m[i], n = n[i], res_temp)
    if(is.null(ret)) ret <- res_temp else ret <- rbind(ret, res_temp)
  }
  rownames(ret) <- NULL
  
  return(ret)
}


# compute analytical approximations of expected thresholds for both data sets:

# load data:
dat_for_mem_fr <- read.csv("../Data/for_mem/ili_mem_fr.csv")
dat_for_mem_us <- read.csv("../Data/for_mem/ili_mem_us.csv")

approx_expectations_fr <- approximate_expectations_from_data(dat_for_mem_fr, log = TRUE, m = 5:15, n = round(30/(5:15)))
approx_expectations_us <- approximate_expectations_from_data(dat_for_mem_us, log = TRUE, m = 5:15, n = round(30/(5:15)))

approx_expectations1_fr <- approximate_expectations_from_data(dat_for_mem_fr, log = TRUE, m = 5:15, n = rep(1, 11))
approx_expectations1_us <- approximate_expectations_from_data(dat_for_mem_us, log = TRUE, m = 5:15, n = rep(1, 11))

approx_expectations_nolog_fr <- approximate_expectations_from_data(dat_for_mem_fr, log = FALSE, m = 5:15, n = round(30/(5:15)))
approx_expectations_nolog_us <- approximate_expectations_from_data(dat_for_mem_us, log = FALSE, m = 5:15, n = round(30/(5:15)))

approx_expectations1_nolog_fr <- approximate_expectations_from_data(dat_for_mem_fr, log = FALSE, m = 5:15, n = rep(1, 11))
approx_expectations1_nolog_us <- approximate_expectations_from_data(dat_for_mem_us, log = FALSE, m = 5:15, n = rep(1, 11))


# write out results:
save(approx_expectations_fr, approx_expectations_us,
     approx_expectations1_fr, approx_expectations1_us,
     approx_expectations_nolog_fr, approx_expectations_nolog_us,
     approx_expectations1_nolog_fr, approx_expectations1_nolog_us,
     file = "Results/approx_expectations.rda")
