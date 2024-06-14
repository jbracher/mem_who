############################################################
# Running permutation/simulation study for US

# extension to secular trends

# set working directory
current_path = rstudioapi::getActiveDocumentContext()$path # get path of this file
setwd(dirname(current_path))

library(mem)
library(ciTools)
# get functions
source("functions_simulations_mem.R")
source("mem_lm.R")

dat_for_mem <- read.csv("../Data/for_mem/ili_mem_us.csv")

# compute peaks of seasons:
peaks_test <- apply(dat_for_mem, 2, max)


# simulation settings:
range_i.seasons <- 5:15 # range to be explored for i.seasons
n_sim <- 500




##############################################
### Using one observation per season, yearly growth rate 0.03, trend accounted for
##############################################

# initialize matrices to store results (thresholds and exceedance proportions)
thresholds1_trend3in_us_0.4 <- matrix(NA, nrow = n_sim, ncol = length(range_i.seasons))
colnames(thresholds1_trend3in_us_0.4) <- paste0("i.seasons_", range_i.seasons)
thresholds1_trend3in_us_0.9 <- thresholds1_trend3in_us_0.975 <- thresholds1_trend3in_us_0.4

thresholds1_trend3in_us_0.4 <- matrix(NA, nrow = n_sim, ncol = length(range_i.seasons))
colnames(thresholds1_trend3in_us_0.4) <- paste0("i.seasons_", range_i.seasons)
thresholds1_trend3in_us_0.9 <- thresholds1_trend3in_us_0.975 <-
  exceedance1_trend3in_us_0.4 <- exceedance1_trend3in_us_0.9 <- exceedance1_trend3in_us_0.975 <-
  thresholds1_trend3in_us_0.4

for(k in 1:n_sim){
  
  set.seed(k)
  
  # generate bootstrapped time series:
  inds_train <- sample(1:ncol(dat_for_mem), size = max(range_i.seasons), replace = TRUE)
  dat_train <- dat_for_mem[, inds_train]
  
  # run through different numbers of past seasons used:
  for(a in seq_along(range_i.seasons)){
    
    # restrict to first a seasons (could also restrict to last n seasons, does not make a difference):
    dat_train_temp <- dat_train[, 1:range_i.seasons[a]]
    # scaling works such that (m+1)-th would remain unscaled
    scaling <- matrix(rep(1.03^((-ncol(dat_train_temp)):-1), each = nrow(dat_train_temp)),
                      ncol = ncol(dat_train_temp), nrow = nrow(dat_train_temp))
    dat_train_temp <- scaling*dat_train_temp
    
    # run memmodel
    memmodel_temp <- memmodel_lm(dat_train_temp, i.seasons = range_i.seasons[a],
                                 i.level.intensity = c(0.4, 0.9, 0.975),
                                 i.trafo = "log", 
                                 i.n.max = 1,
                                 include_trend = TRUE)
    
    # store thresholds:
    thresholds1_trend3in_us_0.4[k, a] <- memmodel_temp$intensity.thresholds[1]
    thresholds1_trend3in_us_0.9[k, a] <- memmodel_temp$intensity.thresholds[2]
    thresholds1_trend3in_us_0.975[k, a] <- memmodel_temp$intensity.thresholds[3]
    
    # compute and store exceedance proportions:
    exceedance1_trend3in_us_0.4[k, a] <- mean(peaks_test > thresholds1_trend3in_us_0.4[k, a])
    exceedance1_trend3in_us_0.9[k, a] <- mean(peaks_test > thresholds1_trend3in_us_0.9[k, a])
    exceedance1_trend3in_us_0.975[k, a] <- mean(peaks_test > thresholds1_trend3in_us_0.975[k, a])
  }
  print(k)
}

# save results:
save(thresholds1_trend3in_us_0.4, thresholds1_trend3in_us_0.9, thresholds1_trend3in_us_0.975,
     exceedance1_trend3in_us_0.4, exceedance1_trend3in_us_0.9, exceedance1_trend3in_us_0.975,
     file = "Results/results1_trend3in_us.rda")



##############################################
### Using one observation per season, yearly growth rate 0.03, trend not accounted for
##############################################

# initialize matrices to store results (thresholds and exceedance proportions)
thresholds1_trend3out_us_0.4 <- matrix(NA, nrow = n_sim, ncol = length(range_i.seasons))
colnames(thresholds1_trend3out_us_0.4) <- paste0("i.seasons_", range_i.seasons)
thresholds1_trend3out_us_0.9 <- thresholds1_trend3out_us_0.975 <- thresholds1_trend3out_us_0.4

thresholds1_trend3out_us_0.4 <- matrix(NA, nrow = n_sim, ncol = length(range_i.seasons))
colnames(thresholds1_trend3out_us_0.4) <- paste0("i.seasons_", range_i.seasons)
thresholds1_trend3out_us_0.9 <- thresholds1_trend3out_us_0.975 <-
  exceedance1_trend3out_us_0.4 <- exceedance1_trend3out_us_0.9 <- exceedance1_trend3out_us_0.975 <-
  thresholds1_trend3out_us_0.4

for(k in 1:n_sim){
  
  set.seed(k)
  
  # generate bootstrapped time series:
  inds_train <- sample(1:ncol(dat_for_mem), size = max(range_i.seasons), replace = TRUE)
  dat_train <- dat_for_mem[, inds_train]
  
  # run through different numbers of past seasons used:
  for(a in seq_along(range_i.seasons)){
    
    # restrict to first a seasons (could also restrict to last n seasons, does not make a difference):
    dat_train_temp <- dat_train[, 1:range_i.seasons[a]]
    # scaling works such that (m+1)-th would remain unscaled
    scaling <- matrix(rep(1.03^((-ncol(dat_train_temp)):-1), each = nrow(dat_train_temp)),
                      ncol = ncol(dat_train_temp), nrow = nrow(dat_train_temp))
    dat_train_temp <- scaling*dat_train_temp
    
    # run memmodel
    memmodel_temp <- memmodel_lm(dat_train_temp, i.seasons = range_i.seasons[a],
                                 i.level.intensity = c(0.4, 0.9, 0.975),
                                 i.trafo = "log", 
                                 i.n.max = 1,
                                 include_trend = FALSE)
    
    # store thresholds:
    thresholds1_trend3out_us_0.4[k, a] <- memmodel_temp$intensity.thresholds[1]
    thresholds1_trend3out_us_0.9[k, a] <- memmodel_temp$intensity.thresholds[2]
    thresholds1_trend3out_us_0.975[k, a] <- memmodel_temp$intensity.thresholds[3]
    
    # compute and store exceedance proportions:
    exceedance1_trend3out_us_0.4[k, a] <- mean(peaks_test > thresholds1_trend3out_us_0.4[k, a])
    exceedance1_trend3out_us_0.9[k, a] <- mean(peaks_test > thresholds1_trend3out_us_0.9[k, a])
    exceedance1_trend3out_us_0.975[k, a] <- mean(peaks_test > thresholds1_trend3out_us_0.975[k, a])
  }
  print(k)
}

# save results:
save(thresholds1_trend3out_us_0.4, thresholds1_trend3out_us_0.9, thresholds1_trend3out_us_0.975,
     exceedance1_trend3out_us_0.4, exceedance1_trend3out_us_0.9, exceedance1_trend3out_us_0.975,
     file = "Results/results1_trend3out_us.rda")



##############################################
### Using one observation per season, yearly growth rate -0.03, trend accounted for
##############################################

# initialize matrices to store results (thresholds and exceedance proportions)
thresholds1_trendm3in_us_0.4 <- matrix(NA, nrow = n_sim, ncol = length(range_i.seasons))
colnames(thresholds1_trendm3in_us_0.4) <- paste0("i.seasons_", range_i.seasons)
thresholds1_trendm3in_us_0.9 <- thresholds1_trendm3in_us_0.975 <- thresholds1_trendm3in_us_0.4

thresholds1_trendm3in_us_0.4 <- matrix(NA, nrow = n_sim, ncol = length(range_i.seasons))
colnames(thresholds1_trendm3in_us_0.4) <- paste0("i.seasons_", range_i.seasons)
thresholds1_trendm3in_us_0.9 <- thresholds1_trendm3in_us_0.975 <-
  exceedance1_trendm3in_us_0.4 <- exceedance1_trendm3in_us_0.9 <- exceedance1_trendm3in_us_0.975 <-
  thresholds1_trendm3in_us_0.4

for(k in 1:n_sim){
  
  set.seed(k)
  
  # generate bootstrapped time series:
  inds_train <- sample(1:ncol(dat_for_mem), size = max(range_i.seasons), replace = TRUE)
  dat_train <- dat_for_mem[, inds_train]
  
  # run through different numbers of past seasons used:
  for(a in seq_along(range_i.seasons)){
    
    # restrict to first a seasons (could also restrict to last n seasons, does not make a difference):
    dat_train_temp <- dat_train[, 1:range_i.seasons[a]]
    # scaling works such that (m+1)-th would remain unscaled
    scaling <- matrix(rep(0.97^((-ncol(dat_train_temp)):-1), each = nrow(dat_train_temp)),
                      ncol = ncol(dat_train_temp), nrow = nrow(dat_train_temp))
    dat_train_temp <- scaling*dat_train_temp
    
    # run memmodel
    memmodel_temp <- memmodel_lm(dat_train_temp, i.seasons = range_i.seasons[a],
                                 i.level.intensity = c(0.4, 0.9, 0.975),
                                 i.trafo = "log", 
                                 i.n.max = 1,
                                 include_trend = TRUE)
    
    # store thresholds:
    thresholds1_trendm3in_us_0.4[k, a] <- memmodel_temp$intensity.thresholds[1]
    thresholds1_trendm3in_us_0.9[k, a] <- memmodel_temp$intensity.thresholds[2]
    thresholds1_trendm3in_us_0.975[k, a] <- memmodel_temp$intensity.thresholds[3]
    
    # compute and store exceedance proportions:
    exceedance1_trendm3in_us_0.4[k, a] <- mean(peaks_test > thresholds1_trendm3in_us_0.4[k, a])
    exceedance1_trendm3in_us_0.9[k, a] <- mean(peaks_test > thresholds1_trendm3in_us_0.9[k, a])
    exceedance1_trendm3in_us_0.975[k, a] <- mean(peaks_test > thresholds1_trendm3in_us_0.975[k, a])
  }
  print(k)
}

# save results:
save(thresholds1_trendm3in_us_0.4, thresholds1_trendm3in_us_0.9, thresholds1_trendm3in_us_0.975,
     exceedance1_trendm3in_us_0.4, exceedance1_trendm3in_us_0.9, exceedance1_trendm3in_us_0.975,
     file = "Results/results1_trendm3in_us.rda")

##############################################
### Using one observation per season, yearly growth rate -0.03, trend not accounted for
##############################################

# initialize matrices to store results (thresholds and exceedance proportions)
thresholds1_trendm3out_us_0.4 <- matrix(NA, nrow = n_sim, ncol = length(range_i.seasons))
colnames(thresholds1_trendm3out_us_0.4) <- paste0("i.seasons_", range_i.seasons)
thresholds1_trendm3out_us_0.9 <- thresholds1_trendm3out_us_0.975 <- thresholds1_trendm3out_us_0.4

thresholds1_trendm3out_us_0.4 <- matrix(NA, nrow = n_sim, ncol = length(range_i.seasons))
colnames(thresholds1_trendm3out_us_0.4) <- paste0("i.seasons_", range_i.seasons)
thresholds1_trendm3out_us_0.9 <- thresholds1_trendm3out_us_0.975 <-
  exceedance1_trendm3out_us_0.4 <- exceedance1_trendm3out_us_0.9 <- exceedance1_trendm3out_us_0.975 <-
  thresholds1_trendm3out_us_0.4

for(k in 1:n_sim){
  
  set.seed(k)
  
  # generate bootstrapped time series:
  inds_train <- sample(1:ncol(dat_for_mem), size = max(range_i.seasons), replace = TRUE)
  dat_train <- dat_for_mem[, inds_train]
  
  # run through different numbers of past seasons used:
  for(a in seq_along(range_i.seasons)){
    
    # restrict to first a seasons (could also restrict to last n seasons, does not make a difference):
    dat_train_temp <- dat_train[, 1:range_i.seasons[a]]
    # scaling works such that (m+1)-th would remain unscaled
    scaling <- matrix(rep(0.97^((-ncol(dat_train_temp)):-1), each = nrow(dat_train_temp)),
                      ncol = ncol(dat_train_temp), nrow = nrow(dat_train_temp))
    dat_train_temp <- scaling*dat_train_temp
    
    # run memmodel
    memmodel_temp <- memmodel_lm(dat_train_temp, i.seasons = range_i.seasons[a],
                                 i.level.intensity = c(0.4, 0.9, 0.975),
                                 i.trafo = "log", 
                                 i.n.max = 1,
                                 include_trend = FALSE)
    
    # store thresholds:
    thresholds1_trendm3out_us_0.4[k, a] <- memmodel_temp$intensity.thresholds[1]
    thresholds1_trendm3out_us_0.9[k, a] <- memmodel_temp$intensity.thresholds[2]
    thresholds1_trendm3out_us_0.975[k, a] <- memmodel_temp$intensity.thresholds[3]
    
    # compute and store exceedance proportions:
    exceedance1_trendm3out_us_0.4[k, a] <- mean(peaks_test > thresholds1_trendm3out_us_0.4[k, a])
    exceedance1_trendm3out_us_0.9[k, a] <- mean(peaks_test > thresholds1_trendm3out_us_0.9[k, a])
    exceedance1_trendm3out_us_0.975[k, a] <- mean(peaks_test > thresholds1_trendm3out_us_0.975[k, a])
  }
  print(k)
}

# save results:
save(thresholds1_trendm3out_us_0.4, thresholds1_trendm3out_us_0.9, thresholds1_trendm3out_us_0.975,
     exceedance1_trendm3out_us_0.4, exceedance1_trendm3out_us_0.9, exceedance1_trendm3out_us_0.975,
     file = "Results/results1_trendm3out_us.rda")



##############################################
### Using one observation per season, yearly growth rate 0.00 (i.e., no trend), but trend accounted for
##############################################

# initialize matrices to store results (thresholds and exceedance proportions)
thresholds1_trend0in_us_0.4 <- matrix(NA, nrow = n_sim, ncol = length(range_i.seasons))
colnames(thresholds1_trend0in_us_0.4) <- paste0("i.seasons_", range_i.seasons)
thresholds1_trend0in_us_0.9 <- thresholds1_trend0in_us_0.975 <- thresholds1_trend0in_us_0.4

thresholds1_trend0in_us_0.4 <- matrix(NA, nrow = n_sim, ncol = length(range_i.seasons))
colnames(thresholds1_trend0in_us_0.4) <- paste0("i.seasons_", range_i.seasons)
thresholds1_trend0in_us_0.9 <- thresholds1_trend0in_us_0.975 <-
  exceedance1_trend0in_us_0.4 <- exceedance1_trend0in_us_0.9 <- exceedance1_trend0in_us_0.975 <-
  thresholds1_trend0in_us_0.4

for(k in 1:n_sim){
  
  set.seed(k)
  
  # generate bootstrapped time series:
  inds_train <- sample(1:ncol(dat_for_mem), size = max(range_i.seasons), replace = TRUE)
  dat_train <- dat_for_mem[, inds_train]
  
  # run through different numbers of past seasons used:
  for(a in seq_along(range_i.seasons)){
    
    # restrict to first a seasons (could also restrict to last n seasons, does not make a difference):
    dat_train_temp <- dat_train[, 1:range_i.seasons[a]]
    # scaling works such that (m+1)-th would remain unscaled
    scaling <- matrix(rep(1.00^((-ncol(dat_train_temp)):-1), each = nrow(dat_train_temp)),
                      ncol = ncol(dat_train_temp), nrow = nrow(dat_train_temp))
    dat_train_temp <- scaling*dat_train_temp
    
    # run memmodel
    memmodel_temp <- memmodel_lm(dat_train_temp, i.seasons = range_i.seasons[a],
                                 i.level.intensity = c(0.4, 0.9, 0.975),
                                 i.trafo = "log", 
                                 i.n.max = 1,
                                 include_trend = TRUE)
    
    # store thresholds:
    thresholds1_trend0in_us_0.4[k, a] <- memmodel_temp$intensity.thresholds[1]
    thresholds1_trend0in_us_0.9[k, a] <- memmodel_temp$intensity.thresholds[2]
    thresholds1_trend0in_us_0.975[k, a] <- memmodel_temp$intensity.thresholds[3]
    
    # compute and store exceedance proportions:
    exceedance1_trend0in_us_0.4[k, a] <- mean(peaks_test > thresholds1_trend0in_us_0.4[k, a])
    exceedance1_trend0in_us_0.9[k, a] <- mean(peaks_test > thresholds1_trend0in_us_0.9[k, a])
    exceedance1_trend0in_us_0.975[k, a] <- mean(peaks_test > thresholds1_trend0in_us_0.975[k, a])
  }
  print(k)
}

# save results:
save(thresholds1_trend0in_us_0.4, thresholds1_trend0in_us_0.9, thresholds1_trend0in_us_0.975,
     exceedance1_trend0in_us_0.4, exceedance1_trend0in_us_0.9, exceedance1_trend0in_us_0.975,
     file = "Results/results1_trend0in_us.rda")


##############################################
### Using one observation per season, yearly growth rate 0.00 (i.e., no trend), and trend not accounted for
##############################################

# initialize matrices to store results (thresholds and exceedance proportions)
thresholds1_trend0out_us_0.4 <- matrix(NA, nrow = n_sim, ncol = length(range_i.seasons))
colnames(thresholds1_trend0out_us_0.4) <- paste0("i.seasons_", range_i.seasons)
thresholds1_trend0out_us_0.9 <- thresholds1_trend0out_us_0.975 <- thresholds1_trend0out_us_0.4

thresholds1_trend0out_us_0.4 <- matrix(NA, nrow = n_sim, ncol = length(range_i.seasons))
colnames(thresholds1_trend0out_us_0.4) <- paste0("i.seasons_", range_i.seasons)
thresholds1_trend0out_us_0.9 <- thresholds1_trend0out_us_0.975 <-
  exceedance1_trend0out_us_0.4 <- exceedance1_trend0out_us_0.9 <- exceedance1_trend0out_us_0.975 <-
  thresholds1_trend0out_us_0.4

for(k in 1:n_sim){
  
  set.seed(k)
  
  # generate bootstrapped time series:
  inds_train <- sample(1:ncol(dat_for_mem), size = max(range_i.seasons), replace = TRUE)
  dat_train <- dat_for_mem[, inds_train]
  
  # run through different numbers of past seasons used:
  for(a in seq_along(range_i.seasons)){
    
    # restrict to first a seasons (could also restrict to last n seasons, does not make a difference):
    dat_train_temp <- dat_train[, 1:range_i.seasons[a]]
    # scaling works such that (m+1)-th would remain unscaled
    scaling <- matrix(rep(1.00^((-ncol(dat_train_temp)):-1), each = nrow(dat_train_temp)),
                      ncol = ncol(dat_train_temp), nrow = nrow(dat_train_temp))
    dat_train_temp <- scaling*dat_train_temp
    
    # run memmodel
    memmodel_temp <- memmodel_lm(dat_train_temp, i.seasons = range_i.seasons[a],
                                 i.level.intensity = c(0.4, 0.9, 0.975),
                                 i.trafo = "log", 
                                 i.n.max = 1,
                                 include_trend = FALSE)
    
    # store thresholds:
    thresholds1_trend0out_us_0.4[k, a] <- memmodel_temp$intensity.thresholds[1]
    thresholds1_trend0out_us_0.9[k, a] <- memmodel_temp$intensity.thresholds[2]
    thresholds1_trend0out_us_0.975[k, a] <- memmodel_temp$intensity.thresholds[3]
    
    # compute and store exceedance proportions:
    exceedance1_trend0out_us_0.4[k, a] <- mean(peaks_test > thresholds1_trend0out_us_0.4[k, a])
    exceedance1_trend0out_us_0.9[k, a] <- mean(peaks_test > thresholds1_trend0out_us_0.9[k, a])
    exceedance1_trend0out_us_0.975[k, a] <- mean(peaks_test > thresholds1_trend0out_us_0.975[k, a])
  }
  print(k)
}

# save results:
save(thresholds1_trend0out_us_0.4, thresholds1_trend0out_us_0.9, thresholds1_trend0out_us_0.975,
     exceedance1_trend0out_us_0.4, exceedance1_trend0out_us_0.9, exceedance1_trend0out_us_0.975,
     file = "Results/results1_trend0out_us.rda")


##############################################
### Using one observation per season, yearly growth rate 0.07, trend accounted for
##############################################

# initialize matrices to store results (thresholds and exceedance proportions)
thresholds1_trend7in_us_0.4 <- matrix(NA, nrow = n_sim, ncol = length(range_i.seasons))
colnames(thresholds1_trend7in_us_0.4) <- paste0("i.seasons_", range_i.seasons)
thresholds1_trend7in_us_0.9 <- thresholds1_trend7in_us_0.975 <- thresholds1_trend7in_us_0.4

thresholds1_trend7in_us_0.4 <- matrix(NA, nrow = n_sim, ncol = length(range_i.seasons))
colnames(thresholds1_trend7in_us_0.4) <- paste0("i.seasons_", range_i.seasons)
thresholds1_trend7in_us_0.9 <- thresholds1_trend7in_us_0.975 <-
  exceedance1_trend7in_us_0.4 <- exceedance1_trend7in_us_0.9 <- exceedance1_trend7in_us_0.975 <-
  thresholds1_trend7in_us_0.4

for(k in 1:n_sim){
  
  set.seed(k)
  
  # generate bootstrapped time series:
  inds_train <- sample(1:ncol(dat_for_mem), size = max(range_i.seasons), replace = TRUE)
  dat_train <- dat_for_mem[, inds_train]
  
  # run through different numbers of past seasons used:
  for(a in seq_along(range_i.seasons)){
    
    # restrict to first a seasons (could also restrict to last n seasons, does not make a difference):
    dat_train_temp <- dat_train[, 1:range_i.seasons[a]]
    # scaling works such that (m+1)-th would remain unscaled
    scaling <- matrix(rep(1.07^((-ncol(dat_train_temp)):-1), each = nrow(dat_train_temp)),
                      ncol = ncol(dat_train_temp), nrow = nrow(dat_train_temp))
    dat_train_temp <- scaling*dat_train_temp
    
    # run memmodel
    memmodel_temp <- memmodel_lm(dat_train_temp, i.seasons = range_i.seasons[a],
                                 i.level.intensity = c(0.4, 0.9, 0.975),
                                 i.trafo = "log", 
                                 i.n.max = 1,
                                 include_trend = TRUE)
    
    # store thresholds:
    thresholds1_trend7in_us_0.4[k, a] <- memmodel_temp$intensity.thresholds[1]
    thresholds1_trend7in_us_0.9[k, a] <- memmodel_temp$intensity.thresholds[2]
    thresholds1_trend7in_us_0.975[k, a] <- memmodel_temp$intensity.thresholds[3]
    
    # compute and store exceedance proportions:
    exceedance1_trend7in_us_0.4[k, a] <- mean(peaks_test > thresholds1_trend7in_us_0.4[k, a])
    exceedance1_trend7in_us_0.9[k, a] <- mean(peaks_test > thresholds1_trend7in_us_0.9[k, a])
    exceedance1_trend7in_us_0.975[k, a] <- mean(peaks_test > thresholds1_trend7in_us_0.975[k, a])
  }
  print(k)
}

# save results:
save(thresholds1_trend7in_us_0.4, thresholds1_trend7in_us_0.9, thresholds1_trend7in_us_0.975,
     exceedance1_trend7in_us_0.4, exceedance1_trend7in_us_0.9, exceedance1_trend7in_us_0.975,
     file = "Results/results1_trend7in_us.rda")

##############################################
### Using one observation per season, yearly growth rate 0.07, trend not accounted for
##############################################

# initialize matrices to store results (thresholds and exceedance proportions)
thresholds1_trend7out_us_0.4 <- matrix(NA, nrow = n_sim, ncol = length(range_i.seasons))
colnames(thresholds1_trend7out_us_0.4) <- paste0("i.seasons_", range_i.seasons)
thresholds1_trend7out_us_0.9 <- thresholds1_trend7out_us_0.975 <- thresholds1_trend7out_us_0.4

thresholds1_trend7out_us_0.4 <- matrix(NA, nrow = n_sim, ncol = length(range_i.seasons))
colnames(thresholds1_trend7out_us_0.4) <- paste0("i.seasons_", range_i.seasons)
thresholds1_trend7out_us_0.9 <- thresholds1_trend7out_us_0.975 <-
  exceedance1_trend7out_us_0.4 <- exceedance1_trend7out_us_0.9 <- exceedance1_trend7out_us_0.975 <-
  thresholds1_trend7out_us_0.4

for(k in 1:n_sim){
  
  set.seed(k)
  
  # generate bootstrapped time series:
  inds_train <- sample(1:ncol(dat_for_mem), size = max(range_i.seasons), replace = TRUE)
  dat_train <- dat_for_mem[, inds_train]
  
  # run through different numbers of past seasons used:
  for(a in seq_along(range_i.seasons)){
    
    # restrict to first a seasons (could also restrict to last n seasons, does not make a difference):
    dat_train_temp <- dat_train[, 1:range_i.seasons[a]]
    # scaling works such that (m+1)-th would remain unscaled
    scaling <- matrix(rep(1.07^((-ncol(dat_train_temp)):-1), each = nrow(dat_train_temp)),
                      ncol = ncol(dat_train_temp), nrow = nrow(dat_train_temp))
    dat_train_temp <- scaling*dat_train_temp
    
    # run memmodel
    memmodel_temp <- memmodel_lm(dat_train_temp, i.seasons = range_i.seasons[a],
                                 i.level.intensity = c(0.4, 0.9, 0.975),
                                 i.trafo = "log", 
                                 i.n.max = 1,
                                 include_trend = FALSE)
    
    # store thresholds:
    thresholds1_trend7out_us_0.4[k, a] <- memmodel_temp$intensity.thresholds[1]
    thresholds1_trend7out_us_0.9[k, a] <- memmodel_temp$intensity.thresholds[2]
    thresholds1_trend7out_us_0.975[k, a] <- memmodel_temp$intensity.thresholds[3]
    
    # compute and store exceedance proportions:
    exceedance1_trend7out_us_0.4[k, a] <- mean(peaks_test > thresholds1_trend7out_us_0.4[k, a])
    exceedance1_trend7out_us_0.9[k, a] <- mean(peaks_test > thresholds1_trend7out_us_0.9[k, a])
    exceedance1_trend7out_us_0.975[k, a] <- mean(peaks_test > thresholds1_trend7out_us_0.975[k, a])
  }
  print(k)
}

# save results:
save(thresholds1_trend7out_us_0.4, thresholds1_trend7out_us_0.9, thresholds1_trend7out_us_0.975,
     exceedance1_trend7out_us_0.4, exceedance1_trend7out_us_0.9, exceedance1_trend7out_us_0.975,
     file = "Results/results1_trend7out_us.rda")


##############################################
### Using one observation per season, yearly growth rate -0.07, trend accounted for
##############################################

# initialize matrices to store results (thresholds and exceedance proportions)
thresholds1_trendm7in_us_0.4 <- matrix(NA, nrow = n_sim, ncol = length(range_i.seasons))
colnames(thresholds1_trendm7in_us_0.4) <- paste0("i.seasons_", range_i.seasons)
thresholds1_trendm7in_us_0.9 <- thresholds1_trendm7in_us_0.975 <- thresholds1_trendm7in_us_0.4

thresholds1_trendm7in_us_0.4 <- matrix(NA, nrow = n_sim, ncol = length(range_i.seasons))
colnames(thresholds1_trendm7in_us_0.4) <- paste0("i.seasons_", range_i.seasons)
thresholds1_trendm7in_us_0.9 <- thresholds1_trendm7in_us_0.975 <-
  exceedance1_trendm7in_us_0.4 <- exceedance1_trendm7in_us_0.9 <- exceedance1_trendm7in_us_0.975 <-
  thresholds1_trendm7in_us_0.4

for(k in 1:n_sim){
  
  set.seed(k)
  
  # generate bootstrapped time series:
  inds_train <- sample(1:ncol(dat_for_mem), size = max(range_i.seasons), replace = TRUE)
  dat_train <- dat_for_mem[, inds_train]
  
  # run through different numbers of past seasons used:
  for(a in seq_along(range_i.seasons)){
    
    # restrict to first a seasons (could also restrict to last n seasons, does not make a difference):
    dat_train_temp <- dat_train[, 1:range_i.seasons[a]]
    # scaling works such that (m+1)-th would remain unscaled
    scaling <- matrix(rep(0.93^((-ncol(dat_train_temp)):-1), each = nrow(dat_train_temp)),
                      ncol = ncol(dat_train_temp), nrow = nrow(dat_train_temp))
    dat_train_temp <- scaling*dat_train_temp
    
    # run memmodel
    memmodel_temp <- memmodel_lm(dat_train_temp, i.seasons = range_i.seasons[a],
                                 i.level.intensity = c(0.4, 0.9, 0.975),
                                 i.trafo = "log", 
                                 i.n.max = 1,
                                 include_trend = TRUE)
    
    # store thresholds:
    thresholds1_trendm7in_us_0.4[k, a] <- memmodel_temp$intensity.thresholds[1]
    thresholds1_trendm7in_us_0.9[k, a] <- memmodel_temp$intensity.thresholds[2]
    thresholds1_trendm7in_us_0.975[k, a] <- memmodel_temp$intensity.thresholds[3]
    
    # compute and store exceedance proportions:
    exceedance1_trendm7in_us_0.4[k, a] <- mean(peaks_test > thresholds1_trendm7in_us_0.4[k, a])
    exceedance1_trendm7in_us_0.9[k, a] <- mean(peaks_test > thresholds1_trendm7in_us_0.9[k, a])
    exceedance1_trendm7in_us_0.975[k, a] <- mean(peaks_test > thresholds1_trendm7in_us_0.975[k, a])
  }
  print(k)
}

# save results:
save(thresholds1_trendm7in_us_0.4, thresholds1_trendm7in_us_0.9, thresholds1_trendm7in_us_0.975,
     exceedance1_trendm7in_us_0.4, exceedance1_trendm7in_us_0.9, exceedance1_trendm7in_us_0.975,
     file = "Results/results1_trendm7in_us.rda")

##############################################
### Using one observation per season, yearly growth rate -0.07, trend not accounted for
##############################################

# initialize matrices to store results (thresholds and exceedance proportions)
thresholds1_trendm7out_us_0.4 <- matrix(NA, nrow = n_sim, ncol = length(range_i.seasons))
colnames(thresholds1_trendm7out_us_0.4) <- paste0("i.seasons_", range_i.seasons)
thresholds1_trendm7out_us_0.9 <- thresholds1_trendm7out_us_0.975 <- thresholds1_trendm7out_us_0.4

thresholds1_trendm7out_us_0.4 <- matrix(NA, nrow = n_sim, ncol = length(range_i.seasons))
colnames(thresholds1_trendm7out_us_0.4) <- paste0("i.seasons_", range_i.seasons)
thresholds1_trendm7out_us_0.9 <- thresholds1_trendm7out_us_0.975 <-
  exceedance1_trendm7out_us_0.4 <- exceedance1_trendm7out_us_0.9 <- exceedance1_trendm7out_us_0.975 <-
  thresholds1_trendm7out_us_0.4

for(k in 1:n_sim){
  
  set.seed(k)
  
  # generate bootstrapped time series:
  inds_train <- sample(1:ncol(dat_for_mem), size = max(range_i.seasons), replace = TRUE)
  dat_train <- dat_for_mem[, inds_train]
  
  # run through different numbers of past seasons used:
  for(a in seq_along(range_i.seasons)){
    
    # restrict to first a seasons (could also restrict to last n seasons, does not make a difference):
    dat_train_temp <- dat_train[, 1:range_i.seasons[a]]
    # scaling works such that (m+1)-th would remain unscaled
    scaling <- matrix(rep(0.93^((-ncol(dat_train_temp)):-1), each = nrow(dat_train_temp)),
                      ncol = ncol(dat_train_temp), nrow = nrow(dat_train_temp))
    dat_train_temp <- scaling*dat_train_temp
    
    # run memmodel
    memmodel_temp <- memmodel_lm(dat_train_temp, i.seasons = range_i.seasons[a],
                                 i.level.intensity = c(0.4, 0.9, 0.975),
                                 i.trafo = "log", 
                                 i.n.max = 1,
                                 include_trend = FALSE)
    
    # store thresholds:
    thresholds1_trendm7out_us_0.4[k, a] <- memmodel_temp$intensity.thresholds[1]
    thresholds1_trendm7out_us_0.9[k, a] <- memmodel_temp$intensity.thresholds[2]
    thresholds1_trendm7out_us_0.975[k, a] <- memmodel_temp$intensity.thresholds[3]
    
    # compute and store exceedance proportions:
    exceedance1_trendm7out_us_0.4[k, a] <- mean(peaks_test > thresholds1_trendm7out_us_0.4[k, a])
    exceedance1_trendm7out_us_0.9[k, a] <- mean(peaks_test > thresholds1_trendm7out_us_0.9[k, a])
    exceedance1_trendm7out_us_0.975[k, a] <- mean(peaks_test > thresholds1_trendm7out_us_0.975[k, a])
  }
  print(k)
}

# save results:
save(thresholds1_trendm7out_us_0.4, thresholds1_trendm7out_us_0.9, thresholds1_trendm7out_us_0.975,
     exceedance1_trendm7out_us_0.4, exceedance1_trendm7out_us_0.9, exceedance1_trendm7out_us_0.975,
     file = "Results/results1_trendm7out_us.rda")
