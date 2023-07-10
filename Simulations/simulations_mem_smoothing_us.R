############################################################
# Running permutation/simulation study for France

# set working directory
current_path = rstudioapi::getActiveDocumentContext()$path # get path of this file
setwd(dirname(current_path))

library(mem)
library(zoo)
source("functions_simulations_mem.R")

# width of moving average
width_ma <- 7

dat_for_mem <- read.csv("../Data/for_mem/ili_mem_us.csv")
dat_for_mem_smoothed <- data.frame(apply(dat_for_mem, MARGIN = 2, FUN = rollmean, k = width_ma))

# compute peaks of seasons:
peaks_test <- apply(dat_for_mem, 2, max)
peaks_test_smoothed <- apply(dat_for_mem_smoothed, 2, max)


# simulation settings:
range_i.seasons <- 5:15 # range to be explored for i.seasons
n_sim <- 500


##############################################
### Using one observation per season
##############################################

# initialize matrices to store results (thresholds and exceedance proportions)
thresholds1_smoothed_us_0.4 <- matrix(NA, nrow = n_sim, ncol = length(range_i.seasons))
colnames(thresholds1_smoothed_us_0.4) <- paste0("i.seasons_", range_i.seasons)
thresholds1_smoothed_us_0.9 <- thresholds1_smoothed_us_0.975 <-
  exceedance1_raw_us_0.4 <- exceedance1_raw_us_0.9 <- exceedance1_raw_us_0.975 <-
  exceedance1_smoothed_us_0.4 <- exceedance1_smoothed_us_0.9 <- exceedance1_smoothed_us_0.975 <-
  thresholds1_smoothed_us_0.4

# iterate simulation runs:
for(k in 1:n_sim){
  
  set.seed(k)
  
  # generate bootstrapped time series:
  inds_train <- sample(1:ncol(dat_for_mem_smoothed), size = max(range_i.seasons), replace = TRUE)
  dat_train <- dat_for_mem_smoothed[, inds_train]
  
  # run through different numbers of past seasons used:
  for(a in seq_along(range_i.seasons)){
    
    # restrict to first a seasons (could also restrict to last n seasons, does not make a difference):
    dat_train_temp <- dat_train[, 1:range_i.seasons[a]]
    
    # run memmodel
    memmodel_temp <- memmodel(dat_train_temp, i.seasons = range_i.seasons[a],
                              i.level.intensity = c(0.4, 0.9, 0.975),
                              i.n.max = 1)
    
    # store thresholds:
    thresholds1_smoothed_us_0.4[k, a] <- memmodel_temp$intensity.thresholds[1]
    thresholds1_smoothed_us_0.9[k, a] <- memmodel_temp$intensity.thresholds[2]
    thresholds1_smoothed_us_0.975[k, a] <- memmodel_temp$intensity.thresholds[3]
    
    # compute and store exceedance proportions for raw peaks:
    exceedance1_raw_us_0.4[k, a] <- mean(peaks_test > thresholds1_smoothed_us_0.4[k, a])
    exceedance1_raw_us_0.9[k, a] <- mean(peaks_test > thresholds1_smoothed_us_0.9[k, a])
    exceedance1_raw_us_0.975[k, a] <- mean(peaks_test > thresholds1_smoothed_us_0.975[k, a])
    
    # compute and store exceedance proportions for smoothed peaks:
    exceedance1_smoothed_us_0.4[k, a] <- mean(peaks_test_smoothed > thresholds1_smoothed_us_0.4[k, a])
    exceedance1_smoothed_us_0.9[k, a] <- mean(peaks_test_smoothed > thresholds1_smoothed_us_0.9[k, a])
    exceedance1_smoothed_us_0.975[k, a] <- mean(peaks_test_smoothed > thresholds1_smoothed_us_0.975[k, a])
  }
  print(k)
}

# re-name objects before saving:
generic_names <- c("thresholds1_smoothed_us_0.4", "thresholds1_smoothed_us_0.9", "thresholds1_smoothed_us_0.975",
                   "exceedance1_raw_us_0.4", "exceedance1_raw_us_0.9", "exceedance1_raw_us_0.975",
                   "exceedance1_smoothed_us_0.4", "exceedance1_smoothed_us_0.9", "exceedance1_smoothed_us_0.975")
new_names <- gsub("_smoothed_", paste0("_smoothed", width_ma, "_"), generic_names)
new_names <- gsub("_raw_", paste0("_raw", width_ma, "_"), new_names)

for(i in seq_along(generic_names)){
  assign(new_names[i], get(generic_names[i]))
}

# save results:
save(list = new_names,
     file = paste0("Results/results1_smoothed", width_ma, "_us.rda"))
# load(paste0("Results/results1_smoothed", width_ma, "_us.rda"))




##############################################
### Using one observation per season and without log trafo
##############################################

# initialize matrices to store results (thresholds and exceedance proportions)
thresholds1_nolog_smoothed_us_0.4 <- matrix(NA, nrow = n_sim, ncol = length(range_i.seasons))
colnames(thresholds1_nolog_smoothed_us_0.4) <- paste0("i.seasons_", range_i.seasons)
thresholds1_nolog_smoothed_us_0.9 <- thresholds1_nolog_smoothed_us_0.975 <-
  exceedance1_nolog_smoothed_us_0.4 <- exceedance1_nolog_smoothed_us_0.9 <- exceedance1_nolog_smoothed_us_0.975 <-
  exceedance1_nolog_raw_us_0.4 <- exceedance1_nolog_raw_us_0.9 <- exceedance1_nolog_raw_us_0.975 <-
  thresholds1_nolog_smoothed_us_0.4

for(k in 1:n_sim){
  
  set.seed(k)
  
  # generate bootstrapped time series:
  inds_train <- sample(1:ncol(dat_for_mem_smoothed), size = max(range_i.seasons), replace = TRUE)
  dat_train <- dat_for_mem_smoothed[, inds_train]
  
  # run through different numbers of past seasons used:
  for(a in seq_along(range_i.seasons)){
    
    # restrict to first a seasons:
    dat_train_temp <- dat_train[, 1:range_i.seasons[a]]
    
    # run memmodel
    memmodel_temp <- memmodel(dat_train_temp, i.seasons = range_i.seasons[a],
                              i.level.intensity = c(0.4, 0.9, 0.975),
                              i.n.max = 1, i.type.intensity = 5)
    
    # store thresholds:
    thresholds1_nolog_smoothed_us_0.4[k, a] <- memmodel_temp$intensity.thresholds[1]
    thresholds1_nolog_smoothed_us_0.9[k, a] <- memmodel_temp$intensity.thresholds[2]
    thresholds1_nolog_smoothed_us_0.975[k, a] <- memmodel_temp$intensity.thresholds[3]
    
    # compute and store exceedance proportions for raw peaks:
    exceedance1_nolog_raw_us_0.4[k, a] <- mean(peaks_test > thresholds1_nolog_smoothed_us_0.4[k, a])
    exceedance1_nolog_raw_us_0.9[k, a] <- mean(peaks_test > thresholds1_nolog_smoothed_us_0.9[k, a])
    exceedance1_nolog_raw_us_0.975[k, a] <- mean(peaks_test > thresholds1_nolog_smoothed_us_0.975[k, a])
    
    # compute and store exceedance proportions for smoothed peaks:
    exceedance1_nolog_smoothed_us_0.4[k, a] <- mean(peaks_test_smoothed > thresholds1_nolog_smoothed_us_0.4[k, a])
    exceedance1_nolog_smoothed_us_0.9[k, a] <- mean(peaks_test_smoothed > thresholds1_nolog_smoothed_us_0.9[k, a])
    exceedance1_nolog_smoothed_us_0.975[k, a] <- mean(peaks_test_smoothed > thresholds1_nolog_smoothed_us_0.975[k, a])
  }
  print(k)
}

# re-name objects before saving:
generic_names <- c("thresholds1_nolog_smoothed_us_0.4", "thresholds1_nolog_smoothed_us_0.9", "thresholds1_nolog_smoothed_us_0.975",
                   "exceedance1_nolog_raw_us_0.4", "exceedance1_nolog_raw_us_0.9", "exceedance1_nolog_raw_us_0.975",
                   "exceedance1_nolog_smoothed_us_0.4", "exceedance1_nolog_smoothed_us_0.9", "exceedance1_nolog_smoothed_us_0.975")
new_names <- gsub("_smoothed_", paste0("_smoothed", width_ma, "_"), generic_names)
new_names <- gsub("_raw_", paste0("_raw", width_ma, "_"), new_names)

for(i in seq_along(generic_names)){
  assign(new_names[i], get(generic_names[i]))
}

# save results:
save(list = new_names,
     file = paste0("Results/results1_nolog_smoothed", width_ma, "_us.rda"))
# load(paste0("Results/results1_nolog_smoothed", width_ma, "_us.rda"))
