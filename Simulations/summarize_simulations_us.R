############################################
# Summarize simulation results for US data

# set working directory
current_path = rstudioapi::getActiveDocumentContext()$path # get path of this file
setwd(dirname(current_path))

library(zoo)
# load functions:
source("functions_simulations_mem.R")

# the range of number of seasons:
range_i.seasons <- 5:15
interv <- c(0.05, 0.95)

for(country in c("us")){
  load(paste0("Results/results_", country, ".rda"))
  load(paste0("Results/results1_", country, ".rda"))
  load(paste0("Results/results_nolog_", country, ".rda"))
  load(paste0("Results/results1_nolog_", country, ".rda"))
  
  dat_for_mem_temp <- read.csv(paste0("../Data/for_mem/ili_mem_", country, ".csv"))
  # compute peaks of seasons:
  assign(paste0("peaks_test_", country), apply(dat_for_mem_temp, 2, max))
}

# load analytical approximations:
load("Results/approx_expectations.rda")
load("Results/approx_expectations_smoothed.rda")
load("Results/analytical_sens_spec.rda")
load("Results/analytical_sens_spec_t.rda")


lin <- 0.35
ce <- 0.8
ylim_threshold <- c(0, 300)

##############################
# Plots US:
pdf("../Draft/figure/plot_us.pdf", width = 9, height = 11)

par(mar = c(4, 4, 3, 1), las = 1)
layout(matrix(c(1:20, rep(21, 5)), ncol = 5, byrow = TRUE), heights = c(rep(1, 4), 0.7))

five_plots(thresholds1_nolog_us_0.4, thresholds1_nolog_us_0.9, thresholds1_nolog_us_0.975,
           exceedance1_nolog_us_0.4, exceedance1_nolog_us_0.9, exceedance1_nolog_us_0.975,
           peaks_test = peaks_test_us, range_i.seasons = range_i.seasons,
           ylim_threshold = ylim_threshold,  approx_expectations = approx_expectations1_nolog_us,
           approx_sens = approx_sens, approx_spec = approx_spec, approx_ppv = approx_ppv,
           main = "natural scale, n = 1")

five_plots(thresholds1_us_0.4, thresholds1_us_0.9, thresholds1_us_0.975,
           exceedance1_us_0.4, exceedance1_us_0.9, exceedance1_us_0.975,
           peaks_test = peaks_test_us, range_i.seasons = range_i.seasons,
           ylim_threshold = ylim_threshold, approx_expectations = approx_expectations1_us,
           approx_sens = approx_sens, approx_spec = approx_spec, approx_ppv = approx_ppv,
           main = "log-transformed, n = 1")

five_plots(thresholds_nolog_us_0.4, thresholds_nolog_us_0.9, thresholds_nolog_us_0.975,
           exceedance_nolog_us_0.4, exceedance_nolog_us_0.9, exceedance_nolog_us_0.975,
           peaks_test = peaks_test_us, range_i.seasons = range_i.seasons,
           ylim_threshold = ylim_threshold, approx_expectations = approx_expectations_nolog_us,
           main = "natural scale, n = 30/m")

five_plots(thresholds_us_0.4, thresholds_us_0.9, thresholds_us_0.975,
           exceedance_us_0.4, exceedance_us_0.9, exceedance_us_0.975,
           peaks_test = peaks_test_us, range_i.seasons = range_i.seasons,
           ylim_threshold = ylim_threshold, approx_expectations = approx_expectations_us,
           main = "log-transformed, n = 30 / m")

par(mar = c(0, 0, 1, 0))
legend_summary() # legend defined in separate function

dev.off()


# Plots on smoothing:
for(country in c("us")){
  load(paste0("Results/results1_smoothed", 3, "_", country, ".rda"))
  load(paste0("Results/results1_nolog_smoothed", 3, "_", country, ".rda"))
  load(paste0("Results/results1_smoothed", 7, "_", country, ".rda"))
  load(paste0("Results/results1_nolog_smoothed", 7, "_", country, ".rda"))
  
  dat_for_mem_temp <- read.csv(paste0("../Data/for_mem/ili_mem_", country, ".csv"))
  # compute peaks of seasons:
  assign(paste0("peaks_test_", country), apply(dat_for_mem_temp, 2, max))
  
  dat_for_mem_smoothed3_temp <- data.frame(apply(dat_for_mem_temp, MARGIN = 2, FUN = rollmean, k = 3))
  assign(paste0("peaks_test_smoothed3_", country), apply(dat_for_mem_smoothed3_temp, 2, max))
  
  dat_for_mem_smoothed7_temp <- data.frame(apply(dat_for_mem_temp, MARGIN = 2, FUN = rollmean, k = 7))
  assign(paste0("peaks_test_smoothed7_", country), apply(dat_for_mem_smoothed7_temp, 2, max))
}



# 3-day moving average:
pdf("../Draft/figure/plot_smoothing3_us_small.pdf", width = 9, height = 5)

par(mfrow = c(2, 5), mar = c(4, 4, 3, 1), las = 1)

# five_plots(thresholds1_nolog_smoothed3_us_0.4, thresholds1_nolog_smoothed3_us_0.9, thresholds1_nolog_smoothed3_us_0.975,
#            exceedance1_nolog_raw3_us_0.4, exceedance1_nolog_raw3_us_0.9, exceedance1_nolog_raw3_us_0.975,
#            peaks_test = peaks_test_us, range_i.seasons = range_i.seasons,
#            ylim_threshold = ylim_threshold,  approx_expectations = approx_expectations1_nolog_smoothed3_us,
#            approx_sens = NULL, approx_spec = NULL, approx_ppv = NULL,
#            main = "natural scale, n = 1, l = 3, applied to unsmoothed peaks")

five_plots(thresholds1_smoothed3_us_0.4, thresholds1_smoothed3_us_0.9, thresholds1_smoothed3_us_0.975,
           exceedance1_raw3_us_0.4, exceedance1_raw3_us_0.9, exceedance1_raw3_us_0.975,
           peaks_test = peaks_test_us, range_i.seasons = range_i.seasons,
           ylim_threshold = ylim_threshold, approx_expectations = NULL,
           approx_sens = NULL, approx_spec = NULL, approx_ppv = NULL,
           main = "log-transformed, n = 1, l = 3, applied to unsmoothed peaks")

# five_plots(thresholds1_nolog_smoothed3_us_0.4, thresholds1_nolog_smoothed3_us_0.9, thresholds1_nolog_smoothed3_us_0.975,
#            exceedance1_nolog_smoothed3_us_0.4, exceedance1_nolog_smoothed3_us_0.9, exceedance1_nolog_smoothed3_us_0.975,
#            peaks_test = peaks_test_smoothed3_us, range_i.seasons = range_i.seasons,
#            ylim_threshold = ylim_threshold,  approx_expectations = approx_expectations1_nolog_smoothed3_us,
#            approx_sens = approx_sens, approx_spec = approx_spec, approx_ppv = approx_ppv,
#            main = "natural scale, n = 1, l = 3, applied to smoothed peaks")

five_plots(thresholds1_smoothed3_us_0.4, thresholds1_smoothed3_us_0.9, thresholds1_smoothed3_us_0.975,
           exceedance1_smoothed3_us_0.4, exceedance1_smoothed3_us_0.9, exceedance1_smoothed3_us_0.975,
           peaks_test = peaks_test_smoothed3_us, range_i.seasons = range_i.seasons,
           ylim_threshold = ylim_threshold, approx_expectations = NULL,
           approx_sens = approx_sens, approx_spec = approx_spec, approx_ppv = approx_ppv,
           main = "log-transformed, n = 1, l = 3, applied to smoothed peaks")


dev.off()



# 7-day moving average:
pdf("../Draft/figure/plot_smoothing7_us_small.pdf", width = 9, height = 5)

par(mfrow = c(2, 5), mar = c(4, 4, 3, 1), las = 1)

# five_plots(thresholds1_nolog_smoothed7_us_0.4, thresholds1_nolog_smoothed7_us_0.9, thresholds1_nolog_smoothed7_us_0.975,
#            exceedance1_nolog_raw7_us_0.4, exceedance1_nolog_raw7_us_0.9, exceedance1_nolog_raw7_us_0.975,
#            peaks_test = peaks_test_us, range_i.seasons = range_i.seasons,
#            ylim_threshold = ylim_threshold,  approx_expectations = approx_expectations1_nolog_smoothed7_us,
#            approx_sens = NULL, approx_spec = NULL, approx_ppv = NULL,
#            main = "natural scale, n = 1, l = 7, applied to unsmoothed peaks")

five_plots(thresholds1_smoothed7_us_0.4, thresholds1_smoothed7_us_0.9, thresholds1_smoothed7_us_0.975,
           exceedance1_raw7_us_0.4, exceedance1_raw7_us_0.9, exceedance1_raw7_us_0.975,
           peaks_test = peaks_test_us, range_i.seasons = range_i.seasons,
           ylim_threshold = ylim_threshold, approx_expectations = NULL,
           approx_sens = NULL, approx_spec = NULL, approx_ppv = NULL,
           main = "log-transformed, n = 1, l = 7, applied to unsmoothed peaks")

# five_plots(thresholds1_nolog_smoothed7_us_0.4, thresholds1_nolog_smoothed7_us_0.9, thresholds1_nolog_smoothed7_us_0.975,
#            exceedance1_nolog_smoothed7_us_0.4, exceedance1_nolog_smoothed7_us_0.9, exceedance1_nolog_smoothed7_us_0.975,
#            peaks_test = peaks_test_smoothed7_us, range_i.seasons = range_i.seasons,
#            ylim_threshold = ylim_threshold,  approx_expectations = approx_expectations1_nolog_smoothed7_us,
#            approx_sens = approx_sens, approx_spec = approx_spec, approx_ppv = approx_ppv,
#            main = "natural scale, n = 1, l = 7, applied to smoothed peaks")

five_plots(thresholds1_smoothed7_us_0.4, thresholds1_smoothed7_us_0.9, thresholds1_smoothed7_us_0.975,
           exceedance1_smoothed7_us_0.4, exceedance1_smoothed7_us_0.9, exceedance1_smoothed7_us_0.975,
           peaks_test = peaks_test_smoothed7_us, range_i.seasons = range_i.seasons,
           ylim_threshold = ylim_threshold, approx_expectations = NULL,
           approx_sens = approx_sens, approx_spec = approx_spec, approx_ppv = approx_ppv,
           main = "log-transformed, n = 1, l = 7, applied to smoothed peaks")

dev.off()

# Plot for thresholds using CIs:

# load additional US results using confidence intervals:
load(paste0("Results/results_ci_us.rda"))
load(paste0("Results/results1_ci_us.rda"))

summary_thresholds_ci_us <- sim_summary(thresholds_ci_us_0.4, thresholds_ci_us_0.9, thresholds_ci_us_0.975,
                                        range_i.seasons, interval = interv)
summary_exceedance_ci_us <- sim_summary(exceedance_ci_us_0.4, exceedance_ci_us_0.9, exceedance_ci_us_0.975,
                                        range_i.seasons, interval = interv)
summary_thresholds1_ci_us <- sim_summary(thresholds1_ci_us_0.4, thresholds1_ci_us_0.9, thresholds1_ci_us_0.975,
                                         range_i.seasons, interval = interv)
summary_exceedance1_ci_us <- sim_summary(exceedance1_ci_us_0.4, exceedance1_ci_us_0.9, exceedance1_ci_us_0.975,
                                         range_i.seasons, interval = interv)

pdf("../Draft/figure/plot_ci_us.pdf", width = 8, height = 2.5)
par(mfrow = c(1, 4), mar = c(4, 4, 3, 1), las = 1)

plot_sim_summary(summary_thresholds1_ci_us, xlab = "# included seasons m", ylab = "threshold", ylim = c(0, 300))
plot_exceedance_summary(summary_exceedance1_ci_us, xlab = "# included seasons m", ylab = "average share")
mtext("based on CIs, log-transformed, n = 1", 3, at = 0, line = lin, cex = ce)

plot_sim_summary(summary_thresholds_ci_us, xlab = "# included seasons m", ylab = "threshold", ylim = c(0, 300))
plot_exceedance_summary(summary_exceedance_ci_us, xlab = "# included seasons m", ylab = "average share")
mtext("based on CIs, log-transformed, n = 30 / m", 3, at = 0, line = lin, cex = ce)

dev.off()



# Plots on t-distribution:
for(country in c("us")){
  load(paste0("Results/results_t_", country, ".rda"))
  load(paste0("Results/results1_t_", country, ".rda"))
  load(paste0("Results/results_nolog_t_", country, ".rda"))
  load(paste0("Results/results1_nolog_t_", country, ".rda"))
  
  dat_for_mem_temp <- read.csv(paste0("../Data/for_mem/ili_mem_", country, ".csv"))
  # compute peaks of seasons:
  assign(paste0("peaks_test_", country), apply(dat_for_mem_temp, 2, max))
}



# plot:
pdf("../Draft/figure/plot_t_us_small.pdf", width = 9, height = 5)

par(mfrow = c(2, 5), mar = c(4, 4, 3, 1), las = 1)

five_plots(thresholds1_nolog_t_us_0.4, thresholds1_nolog_t_us_0.9, thresholds1_nolog_t_us_0.975,
           exceedance1_nolog_t_us_0.4, exceedance1_nolog_t_us_0.9, exceedance1_nolog_t_us_0.975,
           peaks_test = peaks_test_us, range_i.seasons = range_i.seasons,
           ylim_threshold = ylim_threshold, approx_expectations = approx_expectations1_nolog_t_us,
           approx_sens = approx_sens_t, approx_spec = approx_spec_t, approx_ppv = approx_ppv_t,
           main = "natural scale, n = 1, using t-distribution")

five_plots(thresholds1_t_us_0.4, thresholds1_t_us_0.9, thresholds1_t_us_0.975,
           exceedance1_t_us_0.4, exceedance1_t_us_0.9, exceedance1_t_us_0.975,
           peaks_test = peaks_test_us, range_i.seasons = range_i.seasons,
           ylim_threshold = ylim_threshold, approx_expectations = approx_expectations1_t_us,
           approx_sens = approx_sens_t, approx_spec = approx_spec_t, approx_ppv = approx_ppv_t,
           main = "log-transformed, n = 1, using t-distribution")


# legend_summary() # legend defined in separate function

dev.off()






# Plots on t-distribution plus trends:
for(country in c("us")){
  load(paste0("Results/results1_trend3in_", country, ".rda"))
  load(paste0("Results/results1_trend3out_", country, ".rda"))
  load(paste0("Results/results1_trendm3in_", country, ".rda"))
  load(paste0("Results/results1_trendm3out_", country, ".rda"))
  load(paste0("Results/results1_trend0in_", country, ".rda"))
  load(paste0("Results/results1_trend0out_", country, ".rda"))
  load(paste0("Results/results1_trend7in_", country, ".rda"))
  load(paste0("Results/results1_trend7out_", country, ".rda"))
  load(paste0("Results/results1_trendm7in_", country, ".rda"))
  load(paste0("Results/results1_trendm7out_", country, ".rda"))
  
  dat_for_mem_temp <- read.csv(paste0("../Data/for_mem/ili_mem_", country, ".csv"))
  # compute peaks of seasons:
  assign(paste0("peaks_test_", country), apply(dat_for_mem_temp, 2, max))
}

# compute summaries manually:
summary_thresholds1_trend3in_us <- sim_summary(thresholds1_trend3in_us_0.4, thresholds1_trend3in_us_0.9, thresholds1_trend3in_us_0.975,
                                               range_i.seasons, interval = interv)
summary_exceedance1_trend3in_us <- sim_summary(exceedance1_trend3in_us_0.4, exceedance1_trend3in_us_0.9, exceedance1_trend3in_us_0.975,
                                               range_i.seasons, interval = interv)

summary_thresholds1_trend3out_us <- sim_summary(thresholds1_trend3out_us_0.4, thresholds1_trend3out_us_0.9, thresholds1_trend3out_us_0.975,
                                                range_i.seasons, interval = interv)
summary_exceedance1_trend3out_us <- sim_summary(exceedance1_trend3out_us_0.4, exceedance1_trend3out_us_0.9, exceedance1_trend3out_us_0.975,
                                                range_i.seasons, interval = interv)

summary_thresholds1_trendm3in_us <- sim_summary(thresholds1_trendm3in_us_0.4, thresholds1_trendm3in_us_0.9, thresholds1_trendm3in_us_0.975,
                                                range_i.seasons, interval = interv)
summary_exceedance1_trendm3in_us <- sim_summary(exceedance1_trendm3in_us_0.4, exceedance1_trendm3in_us_0.9, exceedance1_trendm3in_us_0.975,
                                                range_i.seasons, interval = interv)

summary_thresholds1_trendm3out_us <- sim_summary(thresholds1_trendm3out_us_0.4, thresholds1_trendm3out_us_0.9, thresholds1_trendm3out_us_0.975,
                                                 range_i.seasons, interval = interv)
summary_exceedance1_trendm3out_us <- sim_summary(exceedance1_trendm3out_us_0.4, exceedance1_trendm3out_us_0.9, exceedance1_trendm3out_us_0.975,
                                                 range_i.seasons, interval = interv)


summary_thresholds1_trend0in_us <- sim_summary(thresholds1_trend0in_us_0.4, thresholds1_trend0in_us_0.9, thresholds1_trend0in_us_0.975,
                                               range_i.seasons, interval = interv)
summary_exceedance1_trend0in_us <- sim_summary(exceedance1_trend0in_us_0.4, exceedance1_trend0in_us_0.9, exceedance1_trend0in_us_0.975,
                                               range_i.seasons, interval = interv)

summary_thresholds1_trend0out_us <- sim_summary(thresholds1_trend0out_us_0.4, thresholds1_trend0out_us_0.9, thresholds1_trend0out_us_0.975,
                                                range_i.seasons, interval = interv)
summary_exceedance1_trend0out_us <- sim_summary(exceedance1_trend0out_us_0.4, exceedance1_trend0out_us_0.9, exceedance1_trend0out_us_0.975,
                                                range_i.seasons, interval = interv)


summary_thresholds1_trend7in_us <- sim_summary(thresholds1_trend7in_us_0.4, thresholds1_trend7in_us_0.9, thresholds1_trend7in_us_0.975,
                                               range_i.seasons, interval = interv)
summary_exceedance1_trend7in_us <- sim_summary(exceedance1_trend7in_us_0.4, exceedance1_trend7in_us_0.9, exceedance1_trend7in_us_0.975,
                                               range_i.seasons, interval = interv)

summary_thresholds1_trend7out_us <- sim_summary(thresholds1_trend7out_us_0.4, thresholds1_trend7out_us_0.9, thresholds1_trend7out_us_0.975,
                                                range_i.seasons, interval = interv)
summary_exceedance1_trend7out_us <- sim_summary(exceedance1_trend7out_us_0.4, exceedance1_trend7out_us_0.9, exceedance1_trend7out_us_0.975,
                                                range_i.seasons, interval = interv)

summary_thresholds1_trendm7in_us <- sim_summary(thresholds1_trendm7in_us_0.4, thresholds1_trendm7in_us_0.9, thresholds1_trendm7in_us_0.975,
                                                range_i.seasons, interval = interv)
summary_exceedance1_trendm7in_us <- sim_summary(exceedance1_trendm7in_us_0.4, exceedance1_trendm7in_us_0.9, exceedance1_trendm7in_us_0.975,
                                                range_i.seasons, interval = interv)

summary_thresholds1_trendm7out_us <- sim_summary(thresholds1_trendm7out_us_0.4, thresholds1_trendm7out_us_0.9, thresholds1_trendm7out_us_0.975,
                                                 range_i.seasons, interval = interv)
summary_exceedance1_trendm7out_us <- sim_summary(exceedance1_trendm7out_us_0.4, exceedance1_trendm7out_us_0.9, exceedance1_trendm7out_us_0.975,
                                                 range_i.seasons, interval = interv)



# plot:
pdf("../Draft/figure/plot_trend3_us_small.pdf", width = 9, height = 5)

par(mfrow = c(2, 4), mar = c(4, 4, 3, 1), las = 1)

plot_sim_summary(summary_thresholds1_trend3out_us, xlab = "# included seasons m", ylab = "threshold", ylim = ylim_threshold)
plot_exceedance_summary(summary_exceedance1_trend3out_us, xlab = "# included seasons m", ylab = "average share")
mtext("3% yearly growth, n = 1, trend ignored", 3, at = 0, line = lin, cex = ce)

plot_sim_summary(summary_thresholds1_trend3in_us, xlab = "# included seasons m", ylab = "threshold", ylim = ylim_threshold)
plot_exceedance_summary(summary_exceedance1_trend3in_us, xlab = "# included seasons m", ylab = "average share")
mtext("3% yearly growth, n = 1, trend accounted", 3, at = 0, line = lin, cex = ce)


plot_sim_summary(summary_thresholds1_trendm3out_us, xlab = "# included seasons m", ylab = "threshold", ylim = ylim_threshold)
plot_exceedance_summary(summary_exceedance1_trendm3out_us, xlab = "# included seasons m", ylab = "average share")
mtext("-3% yearly growth, n = 1, trend ignored", 3, at = 0, line = lin, cex = ce)

plot_sim_summary(summary_thresholds1_trendm3in_us, xlab = "# included seasons m", ylab = "threshold", ylim = ylim_threshold)
plot_exceedance_summary(summary_exceedance1_trendm3in_us, xlab = "# included seasons m", ylab = "average share")
mtext("-3% yearly growth, n = 1, trend accounted", 3, at = 0, line = lin, cex = ce)

dev.off()


# plot:
pdf("../Draft/figure/plot_trend7_us_small.pdf", width = 9, height = 5)

par(mfrow = c(2, 4), mar = c(4, 4, 3, 1), las = 1)

plot_sim_summary(summary_thresholds1_trend7out_us, xlab = "# included seasons m", ylab = "threshold", ylim = ylim_threshold)
plot_exceedance_summary(summary_exceedance1_trend7out_us, xlab = "# included seasons m", ylab = "average share")
mtext("7% yearly growth, n = 1, trend ignored", 3, at = 0, line = lin, cex = ce)

plot_sim_summary(summary_thresholds1_trend7in_us, xlab = "# included seasons m", ylab = "threshold", ylim = ylim_threshold)
plot_exceedance_summary(summary_exceedance1_trend7in_us, xlab = "# included seasons m", ylab = "average share")
mtext("7% yearly growth, n = 1, trend accounted", 3, at = 0, line = lin, cex = ce)


plot_sim_summary(summary_thresholds1_trendm7out_us, xlab = "# included seasons m", ylab = "threshold", ylim = ylim_threshold)
plot_exceedance_summary(summary_exceedance1_trendm7out_us, xlab = "# included seasons m", ylab = "average share")
mtext("-7% yearly growth, n = 1, trend ignored", 3, at = 0, line = lin, cex = ce)

plot_sim_summary(summary_thresholds1_trendm7in_us, xlab = "# included seasons m", ylab = "threshold", ylim = ylim_threshold)
plot_exceedance_summary(summary_exceedance1_trendm7in_us, xlab = "# included seasons m", ylab = "average share")
mtext("-7% yearly growth, n = 1, trend accounted", 3, at = 0, line = lin, cex = ce)

dev.off()


# compare sensitivity, specificity, PPV when accounting for trend vs when not:

# compute relevant summaries:
summary_sens_spec_trend0in_us <- summarize_sens_spec(thresholds1_trend0in_us_0.4, 
                                                     thresholds1_trend0in_us_0.9,
                                                     thresholds1_trend0in_us_0.975,
                                                     peaks_test_us, range_i.seasons)
summary_ppv_npv_trend0in_us <- summarize_ppv_npv(thresholds1_trend0in_us_0.4, 
                                                 thresholds1_trend0in_us_0.9,
                                                 thresholds1_trend0in_us_0.975,
                                                 peaks_test_us, range_i.seasons)

summary_sens_spec_trend0out_us <- summarize_sens_spec(thresholds1_trend0out_us_0.4, 
                                                      thresholds1_trend0out_us_0.9,
                                                      thresholds1_trend0out_us_0.975,
                                                      peaks_test_us, range_i.seasons)
summary_ppv_npv_trend0out_us <- summarize_ppv_npv(thresholds1_trend0out_us_0.4, 
                                                  thresholds1_trend0out_us_0.9,
                                                  thresholds1_trend0out_us_0.975,
                                                  peaks_test_us, range_i.seasons)

summary_sens_spec_trend3in_us <- summarize_sens_spec(thresholds1_trend3in_us_0.4, 
                                                     thresholds1_trend3in_us_0.9,
                                                     thresholds1_trend3in_us_0.975,
                                                     peaks_test_us, range_i.seasons)
summary_ppv_npv_trend3in_us <- summarize_ppv_npv(thresholds1_trend3in_us_0.4, 
                                                 thresholds1_trend3in_us_0.9,
                                                 thresholds1_trend3in_us_0.975,
                                                 peaks_test_us, range_i.seasons)

summary_sens_spec_trend3out_us <- summarize_sens_spec(thresholds1_trend3out_us_0.4, 
                                                      thresholds1_trend3out_us_0.9,
                                                      thresholds1_trend3out_us_0.975,
                                                      peaks_test_us, range_i.seasons)
summary_ppv_npv_trend3out_us <- summarize_ppv_npv(thresholds1_trend3out_us_0.4, 
                                                  thresholds1_trend3out_us_0.9,
                                                  thresholds1_trend3out_us_0.975,
                                                  peaks_test_us, range_i.seasons)


summary_sens_spec_trend7in_us <- summarize_sens_spec(thresholds1_trend7in_us_0.4, 
                                                     thresholds1_trend7in_us_0.9,
                                                     thresholds1_trend7in_us_0.975,
                                                     peaks_test_us, range_i.seasons)
summary_ppv_npv_trend7in_us <- summarize_ppv_npv(thresholds1_trend7in_us_0.4, 
                                                 thresholds1_trend7in_us_0.9,
                                                 thresholds1_trend7in_us_0.975,
                                                 peaks_test_us, range_i.seasons)

summary_sens_spec_trend7out_us <- summarize_sens_spec(thresholds1_trend7out_us_0.4, 
                                                      thresholds1_trend7out_us_0.9,
                                                      thresholds1_trend7out_us_0.975,
                                                      peaks_test_us, range_i.seasons)
summary_ppv_npv_trend7out_us <- summarize_ppv_npv(thresholds1_trend7out_us_0.4, 
                                                  thresholds1_trend7out_us_0.9,
                                                  thresholds1_trend7out_us_0.975,
                                                  peaks_test_us, range_i.seasons)


pdf("../Draft/figure/plot_cost_trend_us.pdf", width = 9, height = 10)
layout(matrix(c(1:9, rep(10, 3)), ncol = 3, byrow = TRUE), heights = c(rep(1, 3), 0.7))

par(mar = c(4, 4, 3, 1), las = 1)
plot_sim_summary(summary_sens_spec_trend0out_us, xlab = "# included seasons m", ylab = "sensitivity",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "sens", hlines = 0:5/5)
plot_sim_summary(summary_sens_spec_trend0in_us, xlab = "# included seasons m", ylab = "sensitivity",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "sens", hlines = 0:5/5,
                 add = TRUE, pch = 0)


plot_sim_summary(summary_sens_spec_trend0out_us, xlab = "# included seasons m", ylab = "specificity",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "spec", hlines = 0:5/5)
plot_sim_summary(summary_sens_spec_trend0in_us, xlab = "# included seasons m", ylab = "specificity",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "spec", hlines = 0:5/5,
                 add = TRUE, pch = 0)
mtext("no yearly growth, trend accounted vs not accounted, log-transformed, n = 1", line = lin, cex = ce)

plot_sim_summary(summary_ppv_npv_trend0out_us, xlab = "# included seasons m", ylab = "PPV",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "ppv", hlines = 0:5/5)
plot_sim_summary(summary_ppv_npv_trend0in_us, xlab = "# included seasons m", ylab = "PPV",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "ppv", hlines = 0:5/5,
                 add = TRUE, pch = 0)



plot_sim_summary(summary_sens_spec_trend3out_us, xlab = "# included seasons m", ylab = "sensitivity",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "sens", hlines = 0:5/5)
plot_sim_summary(summary_sens_spec_trend3in_us, xlab = "# included seasons m", ylab = "sensitivity",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "sens", hlines = 0:5/5,
                 add = TRUE, pch = 0)


plot_sim_summary(summary_sens_spec_trend3out_us, xlab = "# included seasons m", ylab = "specificity",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "spec", hlines = 0:5/5)
plot_sim_summary(summary_sens_spec_trend3in_us, xlab = "# included seasons m", ylab = "specificity",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "spec", hlines = 0:5/5,
                 add = TRUE, pch = 0)
mtext("3% growth, trend accounted vs not accounted, log-transformed, n = 1", line = lin, cex = ce)

plot_sim_summary(summary_ppv_npv_trend3out_us, xlab = "# included seasons m", ylab = "PPV",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "ppv", hlines = 0:5/5)
plot_sim_summary(summary_ppv_npv_trend3in_us, xlab = "# included seasons m", ylab = "PPV",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "ppv", hlines = 0:5/5,
                 add = TRUE, pch = 0)



plot_sim_summary(summary_sens_spec_trend7out_us, xlab = "# included seasons m", ylab = "sensitivity",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "sens", hlines = 0:5/5)
plot_sim_summary(summary_sens_spec_trend7in_us, xlab = "# included seasons m", ylab = "sensitivity",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "sens", hlines = 0:5/5,
                 add = TRUE, pch = 0)


plot_sim_summary(summary_sens_spec_trend7out_us, xlab = "# included seasons m", ylab = "specificity",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "spec", hlines = 0:5/5)
plot_sim_summary(summary_sens_spec_trend7in_us, xlab = "# included seasons m", ylab = "specificity",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "spec", hlines = 0:5/5,
                 add = TRUE, pch = 0)
mtext("7% yearly growth, trend accounted vs not accounted, log-transformed, n = 1", line = lin, cex = ce)

plot_sim_summary(summary_ppv_npv_trend7out_us, xlab = "# included seasons m", ylab = "PPV",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "ppv", hlines = 0:5/5)
plot_sim_summary(summary_ppv_npv_trend7in_us, xlab = "# included seasons m", ylab = "PPV",
                 ylim = 0:1, show_bands = FALSE, variable_to_plot = "ppv", hlines = 0:5/5,
                 add = TRUE, pch = 0)

par(mar = c(0, 0, 2, 0))
plot(NULL, xlim = 0:1, ylim = 0:1, xlab = "", ylab = "", axes = FALSE)
legend("top", c("Intensity levels:", "very high", "high", "medium",
                "Trend accounted:", "yes", "no", ""),
       col = c(NA, col_very_high, col_high, col_medium, NA, "black", "black", NA),
       lty = c(NA, 1, 1, 1, NA, NA, NA, NA),
       pch = c(NA, NA, NA, NA, NA, 0, 15, NA), ncol = 2, bty = "n")

dev.off()




# Plots of misclassification matrices:

# compute matices for m = 5, 10, 15:
cm1_us_5 <- classification_matrix(thresholds1_us_0.4, thresholds1_us_0.9, thresholds1_us_0.975, peaks_test = peaks_test_us, i.seasons = 5)
cm1_us_10 <- classification_matrix(thresholds1_us_0.4, thresholds1_us_0.9, thresholds1_us_0.975, peaks_test = peaks_test_us, i.seasons = 10)
cm1_us_15 <- classification_matrix(thresholds1_us_0.4, thresholds1_us_0.9, thresholds1_us_0.975, peaks_test = peaks_test_us, i.seasons = 15)

cm_us_5 <- classification_matrix(thresholds_us_0.4, thresholds_us_0.9, thresholds_us_0.975, peaks_test = peaks_test_us, i.seasons = 5)
cm_us_10 <- classification_matrix(thresholds_us_0.4, thresholds_us_0.9, thresholds_us_0.975, peaks_test = peaks_test_us, i.seasons = 10)
cm_us_15 <- classification_matrix(thresholds_us_0.4, thresholds_us_0.9, thresholds_us_0.975, peaks_test = peaks_test_us, i.seasons = 15)

cm1_nolog_us_5 <- classification_matrix(thresholds1_nolog_us_0.4, thresholds1_nolog_us_0.9, thresholds1_nolog_us_0.975, peaks_test = peaks_test_us, i.seasons = 5)
cm1_nolog_us_10 <- classification_matrix(thresholds1_nolog_us_0.4, thresholds1_nolog_us_0.9, thresholds1_nolog_us_0.975, peaks_test = peaks_test_us, i.seasons = 10)
cm1_nolog_us_15 <- classification_matrix(thresholds1_nolog_us_0.4, thresholds1_nolog_us_0.9, thresholds1_nolog_us_0.975, peaks_test = peaks_test_us, i.seasons = 15)

cm_nolog_us_5 <- classification_matrix(thresholds_nolog_us_0.4, thresholds_nolog_us_0.9, thresholds_nolog_us_0.975, peaks_test = peaks_test_us, i.seasons = 5)
cm_nolog_us_10 <- classification_matrix(thresholds_nolog_us_0.4, thresholds_nolog_us_0.9, thresholds_nolog_us_0.975, peaks_test = peaks_test_us, i.seasons = 10)
cm_nolog_us_15 <- classification_matrix(thresholds_nolog_us_0.4, thresholds_nolog_us_0.9, thresholds_nolog_us_0.975, peaks_test = peaks_test_us, i.seasons = 15)


# Plot for thresholds without smoothing:

pdf("../Draft/figure/mosaic_us.pdf", width = 9, height = 11)
# structure plot area:
layout(matrix(c(1, 2, 3, 4,
                5, 6, 7, 8,
                9, 10, 11, 12,
                13, 14, 15, 16,
                17, 17, 17, 17), nrow = 5, byrow = TRUE), widths = c(0.5, 1, 1, 1), heights = c(1, 1, 1, 1, 0.2))
# top row:
# labelling:
par(mar = c(4, 0.5, 3.5, 1), las = 1)
plot(NULL, xlim = 0:1, ylim = 0:1, xlab = "", ylab = "", axes = FALSE)
text(c("very high           ", "high           ", "medium           ", "low           "),
     y = c(0.99, 0.93, 0.65, 0.2), x = 0.5, cex = 1)
mtext(text = c("true category           "), side = 3, at = c(0.6), cex = 0.9)

# mosaic plots:
par(mar = c(4, 0.5, 3.5, 3.5), las = 1)
mosaic(cm1_nolog_us_5)
mtext("m = 5", cex = 0.9)
mosaic(cm1_nolog_us_10)
mtext("m = 10", cex = 0.9)
mtext("natural scale, n = 1", cex = 1, line = 2)
mosaic(cm1_nolog_us_15)
mtext("m = 15", cex = 0.9)

# second row:

# labelling:
par(mar = c(4, 0.5, 3.5, 1), las = 1)
plot(NULL, xlim = 0:1, ylim = 0:1, xlab = "", ylab = "", axes = FALSE)
text(c("very high           ", "high           ", "medium           ", "low           "),
     y = c(0.99, 0.93, 0.65, 0.2), x = 0.5, cex = 1)
mtext(text = c("true category           "), side = 3, at = c(0.6), cex = 0.9)

# mosaic plots:
par(mar = c(4, 0.5, 3.5, 3.5), las = 1)
mosaic(cm1_us_5)
mtext("m = 5", cex = 0.9)
mosaic(cm1_us_10)
mtext("m = 10", cex = 0.9)
mtext("log-transformed, n = 1", cex = 1, line = 2)
mosaic(cm1_us_15)
mtext("m = 15", cex = 0.9)

# third row:

# labelling:
par(mar = c(4, 1, 3.5, 1), las = 1)
plot(NULL, xlim = 0:1, ylim = 0:1, xlab = "", ylab = "", axes = FALSE)
text(c("very high           ", "high           ", "medium           ", "low           "),
     y = c(0.99, 0.93, 0.65, 0.2), x = 0.5, cex = 1)
mtext(text = c("true category           "), side = 3, at = c(0.6), cex = 0.9)

# mosaic plots:
par(mar = c(4, 0.5, 3.5, 3.5), las = 1)
mosaic(cm_nolog_us_5)
mtext("m = 5", cex = 0.9)
mosaic(cm_nolog_us_10)
mtext("m = 10", cex = 0.9)
mtext("natural scale, n = 30/m", cex = 1, line = 2)
mosaic(cm_nolog_us_15)
mtext("m = 15", cex = 0.9)


# bottom row:

# labelling:
par(mar = c(4, 1, 3.5, 1), las = 1)
plot(NULL, xlim = 0:1, ylim = 0:1, xlab = "", ylab = "", axes = FALSE)
text(c("very high           ", "high           ", "medium           ", "low           "),
     y = c(0.99, 0.93, 0.65, 0.2), x = 0.5, cex = 1)
mtext(text = c("true category           "), side = 3, at = c(0.6), cex = 0.9)

# mosaic plots:
par(mar = c(4, 0.5, 3.5, 3.5), las = 1)
mosaic(cm_us_5)
mtext("m = 5", cex = 0.9)
mosaic(cm_us_10)
mtext("m = 10", cex = 0.9)
mtext("log-transformed, n = 30/m", cex = 1, line = 2)
mosaic(cm_us_15)
mtext("m = 15", cex = 0.9)

par(mar = c(0, 0, 0, 0), las = 1)
plot(NULL, xlim = 0:1, ylim = 0:1, xlab = "", ylab = "", axes = FALSE)
legend("top", legend = c("classification:   ", "low", "medium", "high", "very high"),
       pt.bg = c(NA, col_low, col_medium, col_high, col_very_high),
       col = c(NA, rep("black", 4)), pch = 22, ncol = 5, bty = "n", pt.cex = 2)

dev.off()




# Plot for thresholds with smoothing and log trafo:

# compute matices for m = 5, 10, 15:
cm1_smoothed3_raw_us_5 <- classification_matrix(thresholds1_smoothed3_us_0.4, thresholds1_smoothed3_us_0.9, thresholds1_smoothed3_us_0.975, peaks_test = peaks_test_us, i.seasons = 5)
cm1_smoothed3_raw_us_10 <- classification_matrix(thresholds1_smoothed3_us_0.4, thresholds1_smoothed3_us_0.9, thresholds1_smoothed3_us_0.975, peaks_test = peaks_test_us, i.seasons = 10)
cm1_smoothed3_raw_us_15 <- classification_matrix(thresholds1_smoothed3_us_0.4, thresholds1_smoothed3_us_0.9, thresholds1_smoothed3_us_0.975, peaks_test = peaks_test_us, i.seasons = 15)


cm1_smoothed3_smoothed3_us_5 <- classification_matrix(thresholds1_smoothed3_us_0.4, thresholds1_smoothed3_us_0.9, thresholds1_smoothed3_us_0.975, peaks_test = peaks_test_smoothed3_us, i.seasons = 5)
cm1_smoothed3_smoothed3_us_10 <- classification_matrix(thresholds1_smoothed3_us_0.4, thresholds1_smoothed3_us_0.9, thresholds1_smoothed3_us_0.975, peaks_test = peaks_test_smoothed3_us, i.seasons = 10)
cm1_smoothed3_smoothed3_us_15 <- classification_matrix(thresholds1_smoothed3_us_0.4, thresholds1_smoothed3_us_0.9, thresholds1_smoothed3_us_0.975, peaks_test = peaks_test_smoothed3_us, i.seasons = 15)



pdf("../Draft/figure/mosaic_log_smoothed_us.pdf", width = 9, height = 6)
# structure plot area:
layout(matrix(c(1, 2, 3, 4,
                5, 6, 7, 8,
                9, 9, 9, 9), nrow = 3, byrow = TRUE), widths = c(0.5, 1, 1, 1), heights = c(1, 1, 0.2))
# top row:

# labelling:
par(mar = c(4, 0.5, 3.5, 1), las = 1)
plot(NULL, xlim = 0:1, ylim = 0:1, xlab = "", ylab = "", axes = FALSE)
text(c("very high           ", "high           ", "medium           ", "low           "),
     y = c(0.99, 0.93, 0.65, 0.2), x = 0.5, cex = 1)
mtext(text = c("true category           "), side = 3, at = c(0.6), cex = 0.9)

# mosaic plots:
par(mar = c(4, 0.5, 3.5, 3.5), las = 1)
mosaic(cm1_smoothed3_raw_us_5)
mtext("m = 5", cex = 0.9)
mosaic(cm1_smoothed3_raw_us_10)
mtext("m = 10", cex = 0.9)
mtext("log-transformed, n = 1, l = 3, applied to unsmoothed peaks", cex = 1, line = 2)
mosaic(cm1_smoothed3_raw_us_15)
mtext("m = 15", cex = 0.9)

# bottom row:

# labelling:
par(mar = c(4, 0.5, 3.5, 1), las = 1)
plot(NULL, xlim = 0:1, ylim = 0:1, xlab = "", ylab = "", axes = FALSE)
text(c("very high           ", "high           ", "medium           ", "low           "),
     y = c(0.99, 0.93, 0.65, 0.2), x = 0.5, cex = 1)
mtext(text = c("true category           "), side = 3, at = c(0.6), cex = 0.9)

# mosaic plots:
par(mar = c(4, 0.5, 3.5, 3.5), las = 1)
mosaic(cm1_smoothed3_smoothed3_us_5)
mtext("m = 5", cex = 0.9)
mosaic(cm1_smoothed3_smoothed3_us_10)
mtext("m = 10", cex = 0.9)
mtext("log-transformed, n = 1, l = 3, applied to smoothed peaks", cex = 1, line = 2)
mosaic(cm1_smoothed3_smoothed3_us_15)
mtext("m = 15", cex = 0.9)


par(mar = c(0, 0, 0, 0), las = 1)
plot(NULL, xlim = 0:1, ylim = 0:1, xlab = "", ylab = "", axes = FALSE)
legend("top", legend = c("classification:   ", "low", "medium", "high", "very high"),
       pt.bg = c(NA, col_low, col_medium, col_high, col_very_high),
       col = c(NA, rep("black", 4)), pch = 22, ncol = 5, bty = "n", pt.cex = 2)

dev.off()


# Plot for thresholds with t-distribution:

# compute matices for m = 5, 10, 15:

# with log trafo
cm1_t_us_5 <- classification_matrix(thresholds1_t_us_0.4, thresholds1_t_us_0.9, thresholds1_t_us_0.975, peaks_test = peaks_test_us, i.seasons = 5)
cm1_t_us_10 <- classification_matrix(thresholds1_t_us_0.4, thresholds1_t_us_0.9, thresholds1_t_us_0.975, peaks_test = peaks_test_us, i.seasons = 10)
cm1_t_us_15 <- classification_matrix(thresholds1_t_us_0.4, thresholds1_t_us_0.9, thresholds1_t_us_0.975, peaks_test = peaks_test_us, i.seasons = 15)

# without:
cm1_nolog_t_us_5 <- classification_matrix(thresholds1_nolog_t_us_0.4, thresholds1_nolog_t_us_0.9, thresholds1_nolog_t_us_0.975, peaks_test = peaks_test_us, i.seasons = 5)
cm1_nolog_t_us_10 <- classification_matrix(thresholds1_nolog_t_us_0.4, thresholds1_nolog_t_us_0.9, thresholds1_nolog_t_us_0.975, peaks_test = peaks_test_us, i.seasons = 10)
cm1_nolog_t_us_15 <- classification_matrix(thresholds1_nolog_t_us_0.4, thresholds1_nolog_t_us_0.9, thresholds1_nolog_t_us_0.975, peaks_test = peaks_test_us, i.seasons = 15)



pdf("../Draft/figure/mosaic_t_us.pdf", width = 9, height = 6)
# structure plot area:
layout(matrix(c(1, 2, 3, 4,
                5, 6, 7, 8,
                9, 9, 9, 9), nrow = 3, byrow = TRUE), widths = c(0.5, 1, 1, 1), heights = c(1, 1, 0.2))
# top row:

# labelling:
par(mar = c(4, 0.5, 3.5, 1), las = 1)
plot(NULL, xlim = 0:1, ylim = 0:1, xlab = "", ylab = "", axes = FALSE)
text(c("very high           ", "high           ", "medium           ", "low           "),
     y = c(0.99, 0.93, 0.65, 0.2), x = 0.5, cex = 1)
mtext(text = c("true category           "), side = 3, at = c(0.6), cex = 0.9)

# mosaic plots:
par(mar = c(4, 0.5, 3.5, 3.5), las = 1)
mosaic(cm1_nolog_t_us_5)
mtext("m = 5", cex = 0.9)
mosaic(cm1_nolog_t_us_10)
mtext("m = 10", cex = 0.9)
mtext("natural scale, n = 1, using t-distribution", cex = 1, line = 2)
mosaic(cm1_nolog_t_us_15)
mtext("m = 15", cex = 0.9)

# bottom row:

# labelling:
par(mar = c(4, 0.5, 3.5, 1), las = 1)
plot(NULL, xlim = 0:1, ylim = 0:1, xlab = "", ylab = "", axes = FALSE)
text(c("very high           ", "high           ", "medium           ", "low           "),
     y = c(0.99, 0.93, 0.65, 0.2), x = 0.5, cex = 1)
mtext(text = c("true category           "), side = 3, at = c(0.6), cex = 0.9)

# mosaic plots:
par(mar = c(4, 0.5, 3.5, 3.5), las = 1)
mosaic(cm1_t_us_5)
mtext("m = 5", cex = 0.9)
mosaic(cm1_t_us_10)
mtext("m = 10", cex = 0.9)
mtext("log-transformed, n = 1, using t-distribution", cex = 1, line = 2)
mosaic(cm1_t_us_15)
mtext("m = 15", cex = 0.9)


par(mar = c(0, 0, 0, 0), las = 1)
plot(NULL, xlim = 0:1, ylim = 0:1, xlab = "", ylab = "", axes = FALSE)
legend("top", legend = c("classification:   ", "low", "medium", "high", "very high"),
       pt.bg = c(NA, col_low, col_medium, col_high, col_very_high),
       col = c(NA, rep("black", 4)), pch = 22, ncol = 5, bty = "n", pt.cex = 2)

dev.off()