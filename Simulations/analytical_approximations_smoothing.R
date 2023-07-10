# Compute approximation of expected thresholds for smoothing methods 
# based on empirical mean vector and covariance matrix.

# set working directory
current_path = rstudioapi::getActiveDocumentContext()$path # get path of this file
setwd(dirname(current_path))

# analytical approximation of expected thresholds base on mean vector and covariance matrix.
approximate_expectations_smoothing <- function(mu, Sigma, l, m){
  mean_ana <- mean(mu)
  var_ana <- mean(Sigma)
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
# and covariance matrix are estimated. Wraps around approximate_expectations_smoothing
approximate_expectations_from_data_smoothing0 <- function(dat, log = FALSE, l, m){
  dat_smooth <- data.frame(apply(dat, MARGIN = 2, FUN = rollmean, k = l))

  dat_peaks <- matrix(NA, ncol = ncol(dat), nrow = l)
  # sort:
  for(i in 1:ncol(dat)){
    ind_peak <- which.max(dat_smooth[, i])
    dat_peaks[, i] <- sort(dat[ind_peak + l - 1 - seq(0, by = 1, length.out = l), i], decreasing = TRUE)
  }
  
  dat_peaks_t <- t(as.matrix(dat_peaks))
  if(log) dat_peaks_t <- log(dat_peaks_t)
  
  mu <- colMeans(dat_peaks_t)
  Sigma <- cov(dat_peaks_t)
  
  res <- approximate_expectations_smoothing(mu = mu, Sigma = Sigma, l = l, m = m)
  
  if(log) res <- exp(res)
  
  return(res)
}

# analytical approximation of expected thresholds based on data in mem format from which mean
# and covariance matrix are estimated. Wrapper around approximate_expectations_from_data0
# to apply to several values of n and K.
approximate_expectations_from_data_smoothing <- function(dat, log = FALSE, l, m){
  ret <- NULL
  for(i in seq_along(l)){
    res_temp <- approximate_expectations_from_data_smoothing0(dat = dat, log = log, l = l[i], m = m[i])
    res_temp <- c(l = l[i], m = m[i], res_temp)
    if(is.null(ret)) ret <- res_temp else ret <- rbind(ret, res_temp)
  }
  rownames(ret) <- NULL
  
  return(ret)
}

# load data:
dat_for_mem_fr <- read.csv("../Data/for_mem/ili_mem_fr.csv")
dat_for_mem_us <- read.csv("../Data/for_mem/ili_mem_us.csv")


approx_expectations1_nolog_smoothed3_fr <- approximate_expectations_from_data_smoothing(dat_for_mem_fr, log = FALSE, l = rep(3, 11), m = 5:15)
approx_expectations1_nolog_smoothed3_us <- approximate_expectations_from_data_smoothing(dat_for_mem_us, log = FALSE, l = rep(3, 11), m = 5:15)

approx_expectations1_nolog_smoothed7_fr <- approximate_expectations_from_data_smoothing(dat_for_mem_fr, log = FALSE, l = rep(7, 11), m = 5:15)
approx_expectations1_nolog_smoothed7_us <- approximate_expectations_from_data_smoothing(dat_for_mem_us, log = FALSE, l = rep(7, 11), m = 5:15)


save(approx_expectations1_nolog_smoothed3_fr, # approx_expectations_nolog_es, approx_expectations_nolog_ch, approx_expectations_nolog_us,
     approx_expectations1_nolog_smoothed3_us,# approx_expectations1_nolog_es, approx_expectations1_nolog_ch, approx_expectations1_nolog_us,
     approx_expectations1_nolog_smoothed7_fr,
     approx_expectations1_nolog_smoothed7_fr,
     file = "Results/approx_expectations_smoothed.rda")
