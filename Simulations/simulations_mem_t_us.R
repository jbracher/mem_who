############################################################
# Running permutation/simulation study for France

# Extension to t-distribution

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



#####################################################
### Run simulation using default settings of mem, but with t-distribution
#####################################################

# initialize matrices to store results (thresholds and exceedance proportions)
thresholds_t_us_0.4 <- matrix(NA, nrow = n_sim, ncol = length(range_i.seasons))
colnames(thresholds_t_us_0.4) <- paste0("i.seasons_", range_i.seasons)
thresholds_t_us_0.9 <- thresholds_t_us_0.975 <-
  exceedance_t_us_0.4 <- exceedance_t_us_0.9 <- exceedance_t_us_0.975 <-
  thresholds_t_us_0.4

# run simulation:
for(k in 1:n_sim){
  set.seed(k)
  
  # generate bootstrapped time series:
  inds_train <- sample(1:ncol(dat_for_mem), size = max(range_i.seasons), replace = TRUE)
  dat_train <- dat_for_mem[, inds_train]
  
  # run through different numbers of past seasons used:
  for(a in seq_along(range_i.seasons)){
    
    # restrict to first a seasons:
    dat_train_temp <- dat_train[, 1:range_i.seasons[a]]
    
    # run memmodel:
    # memmodel_temp <- memmodel(dat_train_temp, i.seasons = range_i.seasons[a],
    #                          i.level.intensity = c(0.4, 0.9, 0.975))
    memmodel_temp <- memmodel_lm(dat_train_temp, i.seasons = range_i.seasons[a],
                                 i.level.intensity = c(0.4, 0.9, 0.975),
                                 i.trafo = "log", 
                                 i.n.max = round(30/range_i.seasons[a]))
    
    # store thresholds:
    thresholds_t_us_0.4[k, a] <- memmodel_temp$intensity.thresholds[1]
    thresholds_t_us_0.9[k, a] <- memmodel_temp$intensity.thresholds[2]
    thresholds_t_us_0.975[k, a] <- memmodel_temp$intensity.thresholds[3]
    
    # compute and store exceedance proportions:
    exceedance_t_us_0.4[k, a] <- mean(peaks_test > thresholds_t_us_0.4[k, a])
    exceedance_t_us_0.9[k, a] <- mean(peaks_test > thresholds_t_us_0.9[k, a])
    exceedance_t_us_0.975[k, a] <- mean(peaks_test > thresholds_t_us_0.975[k, a])
  }
  print(k)
}

# save results
save(thresholds_t_us_0.4, thresholds_t_us_0.9, thresholds_t_us_0.975,
     exceedance_t_us_0.4, exceedance_t_us_0.9, exceedance_t_us_0.975,
     file = "Results/results_t_us.rda")
load("Results/results_t_us.rda")

# compute summaries using custom function
summary_thresholds <- sim_summary(thresholds_t_us_0.4, thresholds_t_us_0.9, thresholds_t_us_0.975, range_i.seasons, interval = c(0.1, 0.9))
summary_exceedance <- sim_summary(exceedance_t_us_0.4, exceedance_t_us_0.9, exceedance_t_us_0.975, range_i.seasons, interval = c(0.1, 0.9))
summary_sens_spec <- summarize_sens_spec(thresholds_t_us_0.4, thresholds_t_us_0.9, thresholds_t_us_0.975, peaks_test, range_i.seasons)
summary_ppv_npv <- summarize_ppv_npv(thresholds_t_us_0.4, thresholds_t_us_0.9, thresholds_t_us_0.975, peaks_test, range_i.seasons)


# plot
par(mfrow = c(1, 6))
plot_sim_summary(summary_thresholds, xlab = "# included seasons", ylab = "threshold")


plot_exceedance_summary(summary_exceedance,
                        xlab = "# included seasons", ylab = "av. proportion per category")

plot_sim_summary(summary_sens_spec, xlab = "# included seasons", ylab = "sensitivity",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "sens")
plot_sim_summary(summary_sens_spec, xlab = "# included seasons", ylab = "specificity",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "spec")


plot_sim_summary(summary_ppv_npv, xlab = "# included seasons", ylab = "PPV",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "ppv")
plot_sim_summary(summary_ppv_npv, xlab = "# included seasons", ylab = "NPV",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "npv")

##############################################
### Using one observation per season and t-distribution
##############################################

# initialize matrices to store results (thresholds and exceedance proportions)
thresholds1_t_us_0.4 <- matrix(NA, nrow = n_sim, ncol = length(range_i.seasons))
colnames(thresholds1_t_us_0.4) <- paste0("i.seasons_", range_i.seasons)
thresholds1_t_us_0.9 <- thresholds1_t_us_0.975 <- thresholds1_t_us_0.4

thresholds1_t_us_0.4 <- matrix(NA, nrow = n_sim, ncol = length(range_i.seasons))
colnames(thresholds1_t_us_0.4) <- paste0("i.seasons_", range_i.seasons)
thresholds1_t_us_0.9 <- thresholds1_t_us_0.975 <-
  exceedance1_t_us_0.4 <- exceedance1_t_us_0.9 <- exceedance1_t_us_0.975 <-
  thresholds1_t_us_0.4

for(k in 1:n_sim){
  
  set.seed(k)
  
  # generate bootstrapped time series:
  inds_train <- sample(1:ncol(dat_for_mem), size = max(range_i.seasons), replace = TRUE)
  dat_train <- dat_for_mem[, inds_train]
  
  # run through different numbers of past seasons used:
  for(a in seq_along(range_i.seasons)){
    
    # restrict to first a seasons (could also restrict to last n seasons, does not make a difference):
    dat_train_temp <- dat_train[, 1:range_i.seasons[a]]
    
    # run memmodel
    # memmodel_temp <- memmodel(dat_train, i.seasons = range_i.seasons[a],
    #                           i.level.intensity = c(0.4, 0.9, 0.975),
    #                           i.n.max = 1)
    
    memmodel_temp <- memmodel_lm(dat_train_temp, i.seasons = range_i.seasons[a],
                                 i.level.intensity = c(0.4, 0.9, 0.975),
                                 i.trafo = "log", 
                                 i.n.max = 1)
    
    # store thresholds:
    thresholds1_t_us_0.4[k, a] <- memmodel_temp$intensity.thresholds[1]
    thresholds1_t_us_0.9[k, a] <- memmodel_temp$intensity.thresholds[2]
    thresholds1_t_us_0.975[k, a] <- memmodel_temp$intensity.thresholds[3]
    
    # compute and store exceedance proportions:
    exceedance1_t_us_0.4[k, a] <- mean(peaks_test > thresholds1_t_us_0.4[k, a])
    exceedance1_t_us_0.9[k, a] <- mean(peaks_test > thresholds1_t_us_0.9[k, a])
    exceedance1_t_us_0.975[k, a] <- mean(peaks_test > thresholds1_t_us_0.975[k, a])
  }
  print(k)
}

# save results:
save(thresholds1_t_us_0.4, thresholds1_t_us_0.9, thresholds1_t_us_0.975,
     exceedance1_t_us_0.4, exceedance1_t_us_0.9, exceedance1_t_us_0.975,
     file = "Results/results1_t_us.rda")
load("Results/results1_t_us.rda")

# compute summaries using custom function
summary_thresholds1 <- sim_summary(thresholds1_t_us_0.4, thresholds1_t_us_0.9, thresholds1_t_us_0.975, range_i.seasons, interval = c(0.1, 0.9))
summary_exceedance1 <- sim_summary(exceedance1_t_us_0.4, exceedance1_t_us_0.9, exceedance1_t_us_0.975, range_i.seasons, interval = c(0.1, 0.9))
summary_sens_spec1 <- summarize_sens_spec(thresholds1_t_us_0.4, thresholds1_t_us_0.9, thresholds1_t_us_0.975, peaks_test, range_i.seasons)
summary_ppv_npv1 <- summarize_ppv_npv(thresholds1_t_us_0.4, thresholds1_t_us_0.9, thresholds1_t_us_0.975, peaks_test, range_i.seasons)


# plot:
par(mfrow = c(1, 6))
plot_sim_summary(summary_thresholds1, xlab = "# included seasons", ylab = "threshold")
plot_exceedance_summary(summary_exceedance1,
                        xlab = "# included seasons", ylab = "av. proportion per category")


plot_sim_summary(summary_sens_spec1, xlab = "# included seasons", ylab = "sensitivity",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "sens")
plot_sim_summary(summary_sens_spec1, xlab = "# included seasons", ylab = "specificity",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "spec")


plot_sim_summary(summary_ppv_npv1, xlab = "# included seasons", ylab = "PPV",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "ppv")
plot_sim_summary(summary_ppv_npv1, xlab = "# included seasons", ylab = "NPV",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "npv")

#####################################################
### Using default settings of mem, but without log trafo and t-distribution
#####################################################


# initialize matrices to store results (thresholds and exceedance proportions)
thresholds_nolog_t_us_0.4 <- matrix(NA, nrow = n_sim, ncol = length(range_i.seasons))
colnames(thresholds_nolog_t_us_0.4) <- paste0("i.seasons_", range_i.seasons)
thresholds_nolog_t_us_0.9 <- thresholds_nolog_t_us_0.975 <-
  exceedance_nolog_t_us_0.4 <- exceedance_nolog_t_us_0.9 <- exceedance_nolog_t_us_0.975 <-
  thresholds_nolog_t_us_0.4

# run simulation:
for(k in 1:n_sim){
  set.seed(k)
  
  # generate bootstrapped time series:
  inds_train <- sample(1:ncol(dat_for_mem), size = max(range_i.seasons), replace = TRUE)
  dat_train <- dat_for_mem[, inds_train]
  
  # run through different numbers of past seasons used:
  for(a in seq_along(range_i.seasons)){
    
    # restrict to first a seasons:
    dat_train_temp <- dat_train[, 1:range_i.seasons[a]]
    
    # run memmodel:
    memmodel_temp <- memmodel_lm(dat_train_temp, i.seasons = range_i.seasons[a],
                                 i.level.intensity = c(0.4, 0.9, 0.975),
                                 i.trafo = "identity", 
                                 i.n.max = round(30/range_i.seasons[a]))
    
    # store thresholds:
    thresholds_nolog_t_us_0.4[k, a] <- memmodel_temp$intensity.thresholds[1]
    thresholds_nolog_t_us_0.9[k, a] <- memmodel_temp$intensity.thresholds[2]
    thresholds_nolog_t_us_0.975[k, a] <- memmodel_temp$intensity.thresholds[3]
    
    # compute and store exceedance proportions:
    exceedance_nolog_t_us_0.4[k, a] <- mean(peaks_test > thresholds_nolog_t_us_0.4[k, a])
    exceedance_nolog_t_us_0.9[k, a] <- mean(peaks_test > thresholds_nolog_t_us_0.9[k, a])
    exceedance_nolog_t_us_0.975[k, a] <- mean(peaks_test > thresholds_nolog_t_us_0.975[k, a])
  }
  print(k)
}

# save results
save(thresholds_nolog_t_us_0.4, thresholds_nolog_t_us_0.9, thresholds_nolog_t_us_0.975,
     exceedance_nolog_t_us_0.4, exceedance_nolog_t_us_0.9, exceedance_nolog_t_us_0.975,
     file = "Results/results_nolog_t_us.rda")
load("Results/results_nolog_t_us.rda")

# compute summaries using custom function
summary_thresholds_nolog <- sim_summary(thresholds_nolog_t_us_0.4, thresholds_nolog_t_us_0.9, thresholds_nolog_t_us_0.975, range_i.seasons, interval = c(0.1, 0.9))
summary_exceedance_nolog <- sim_summary(exceedance_nolog_t_us_0.4, exceedance_nolog_t_us_0.9, exceedance_nolog_t_us_0.975, range_i.seasons, interval = c(0.1, 0.9))
summary_sens_spec_nolog <- summarize_sens_spec(thresholds_nolog_t_us_0.4, thresholds_nolog_t_us_0.9, thresholds_nolog_t_us_0.975, peaks_test, range_i.seasons)
summary_ppv_npv_nolog <- summarize_ppv_npv(thresholds_nolog_t_us_0.4, thresholds_nolog_t_us_0.9, thresholds_nolog_t_us_0.975, peaks_test, range_i.seasons)


# plot
par(mfrow = c(1, 6))
plot_sim_summary(summary_thresholds_nolog, xlab = "# included seasons", ylab = "threshold")


plot_exceedance_summary(summary_exceedance_nolog,
                        xlab = "# included seasons", ylab = "av. proportion per category")

plot_sim_summary(summary_sens_spec_nolog, xlab = "# included seasons", ylab = "sensitivity",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "sens")
plot_sim_summary(summary_sens_spec_nolog, xlab = "# included seasons", ylab = "specificity",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "spec")


plot_sim_summary(summary_ppv_npv_nolog, xlab = "# included seasons", ylab = "PPV",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "ppv")
plot_sim_summary(summary_ppv_npv_nolog, xlab = "# included seasons", ylab = "NPV",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "npv")


##############################################
### Using one observation per season, without log trafo, with t-distribution
##############################################

# initialize matrices to store results (thresholds and exceedance proportions)
thresholds1_nolog_t_us_0.4 <- matrix(NA, nrow = n_sim, ncol = length(range_i.seasons))
colnames(thresholds1_nolog_t_us_0.4) <- paste0("i.seasons_", range_i.seasons)
thresholds1_nolog_t_us_0.9 <- thresholds1_nolog_t_us_0.975 <- thresholds1_nolog_t_us_0.4

thresholds1_nolog_t_us_0.4 <- matrix(NA, nrow = n_sim, ncol = length(range_i.seasons))
colnames(thresholds1_nolog_t_us_0.4) <- paste0("i.seasons_", range_i.seasons)
thresholds1_nolog_t_us_0.9 <- thresholds1_nolog_t_us_0.975 <-
  exceedance1_nolog_t_us_0.4 <- exceedance1_nolog_t_us_0.9 <- exceedance1_nolog_t_us_0.975 <-
  thresholds1_nolog_t_us_0.4

for(k in 1:n_sim){
  
  set.seed(k)
  
  # generate bootstrapped time series:
  inds_train <- sample(1:ncol(dat_for_mem), size = max(range_i.seasons), replace = TRUE)
  dat_train <- dat_for_mem[, inds_train]
  
  # run through different numbers of past seasons used:
  for(a in seq_along(range_i.seasons)){
    
    # restrict to first a seasons:
    dat_train_temp <- dat_train[, 1:range_i.seasons[a]]
    
    # run memmodel
    memmodel_temp <- memmodel_lm(dat_train_temp, i.seasons = range_i.seasons[a],
                                 i.level.intensity = c(0.4, 0.9, 0.975),
                                 i.trafo = "identity", i.n.max = 1)
    
    # store thresholds:
    thresholds1_nolog_t_us_0.4[k, a] <- memmodel_temp$intensity.thresholds[1]
    thresholds1_nolog_t_us_0.9[k, a] <- memmodel_temp$intensity.thresholds[2]
    thresholds1_nolog_t_us_0.975[k, a] <- memmodel_temp$intensity.thresholds[3]
    
    # compute and store exceedance proportions:
    exceedance1_nolog_t_us_0.4[k, a] <- mean(peaks_test > thresholds1_nolog_t_us_0.4[k, a])
    exceedance1_nolog_t_us_0.9[k, a] <- mean(peaks_test > thresholds1_nolog_t_us_0.9[k, a])
    exceedance1_nolog_t_us_0.975[k, a] <- mean(peaks_test > thresholds1_nolog_t_us_0.975[k, a])
  }
  print(k)
}

# save results:
save(thresholds1_nolog_t_us_0.4, thresholds1_nolog_t_us_0.9, thresholds1_nolog_t_us_0.975,
     exceedance1_nolog_t_us_0.4, exceedance1_nolog_t_us_0.9, exceedance1_nolog_t_us_0.975,
     file = "Results/results1_nolog_t_us.rda")
load("Results/results1_nolog_t_us.rda")

# compute summaries using custom function
summary_thresholds1_nolog_us <- sim_summary(thresholds1_nolog_t_us_0.4, thresholds1_nolog_t_us_0.9, thresholds1_nolog_t_us_0.975, range_i.seasons, interval = c(0.1, 0.9))
summary_exceedance1_nolog_us <- sim_summary(exceedance1_nolog_t_us_0.4, exceedance1_nolog_t_us_0.9, exceedance1_nolog_t_us_0.975, range_i.seasons, interval = c(0.1, 0.9))
summary_sens_spec1_nolog <- summarize_sens_spec(thresholds1_nolog_t_us_0.4, thresholds1_nolog_t_us_0.9, thresholds1_nolog_t_us_0.975, peaks_test, range_i.seasons)
summary_ppv_npv1_nolog <- summarize_ppv_npv(thresholds1_nolog_t_us_0.4, thresholds1_nolog_t_us_0.9, thresholds1_nolog_t_us_0.975, peaks_test, range_i.seasons)


# plot
par(mfrow = c(1, 6))
plot_sim_summary(summary_thresholds1_nolog_us, xlab = "# included seasons", ylab = "threshold")

plot_exceedance_summary(summary_exceedance1_nolog_us,
                        xlab = "# included seasons", ylab = "av. proportion per category")

plot_sim_summary(summary_sens_spec1_nolog, xlab = "# included seasons", ylab = "sensitivity",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "sens")
plot_sim_summary(summary_sens_spec1_nolog, xlab = "# included seasons", ylab = "specificity",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "spec")


plot_sim_summary(summary_ppv_npv1_nolog, xlab = "# included seasons", ylab = "PPV",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "ppv")
plot_sim_summary(summary_ppv_npv1_nolog, xlab = "# included seasons", ylab = "NPV",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "npv")
