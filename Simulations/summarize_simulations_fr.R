# summarize simulation results

# set working directory
current_path = rstudioapi::getActiveDocumentContext()$path # get path of this file
setwd(dirname(current_path))

library(zoo)
# load functions:
source("functions_simulations_mem.R")

# the range of number of seasons:
range_i.seasons <- 5:15
interv <- c(0.05, 0.95)

for(country in c("fr")){
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


lin <- 0.35
ce <- 0.8

# maybe add code bit to compute all summaries.

# save(summary_thresholds_us, summary_thresholds1_us, summary_thresholds_nolog_us, summary_thresholds1_nolog_us,
#      summary_exceedance_us, summary_exceedance1_us, summary_exceedance_nolog_us, summary_exceedance1_nolog_us,
#      summary_thresholds_ch, summary_thresholds1_ch, summary_thresholds_nolog_ch, summary_thresholds1_nolog_ch,
#      summary_exceedance_ch, summary_exceedance1_ch, summary_exceedance_nolog_ch, summary_exceedance1_nolog_ch,
#      summary_thresholds_es, summary_thresholds1_es, summary_thresholds_nolog_es, summary_thresholds1_nolog_es,
#      summary_exceedance_es, summary_exceedance1_es, summary_exceedance_nolog_es, summary_exceedance1_nolog_es,
#      summary_thresholds_fr, summary_thresholds1_fr, summary_thresholds_nolog_fr, summary_thresholds1_nolog_fr,
#      summary_exceedance_fr, summary_exceedance1_fr, summary_exceedance_nolog_fr, summary_exceedance1_nolog_fr,
#      file = "Results/summaries_results.rda")

##############################
# Plots France:
pdf("../Draft/figure/plot_fr.pdf", width = 9, height = 11)

par(mar = c(4, 4, 3, 1), las = 1)
layout(matrix(c(1:20, rep(21, 5)), ncol = 5, byrow = TRUE), heights = c(rep(1, 4), 0.7))

five_plots(thresholds1_nolog_fr_0.4, thresholds1_nolog_fr_0.9, thresholds1_nolog_fr_0.975,
           exceedance1_nolog_fr_0.4, exceedance1_nolog_fr_0.9, exceedance1_nolog_fr_0.975,
           peaks_test = peaks_test_fr, range_i.seasons = range_i.seasons,
           ylim_threshold = c(0, 500),  approx_expectations = approx_expectations1_nolog_fr,
           approx_sens = approx_sens, approx_spec = approx_spec, approx_ppv = approx_ppv,
           main = "natural scale, n = 1")

five_plots(thresholds1_fr_0.4, thresholds1_fr_0.9, thresholds1_fr_0.975,
           exceedance1_fr_0.4, exceedance1_fr_0.9, exceedance1_fr_0.975,
           peaks_test = peaks_test_fr, range_i.seasons = range_i.seasons,
           ylim_threshold = c(0, 500), approx_expectations = approx_expectations1_fr,
           approx_sens = approx_sens, approx_spec = approx_spec, approx_ppv = approx_ppv,
           main = "log-transformed, n = 1")

five_plots(thresholds_nolog_fr_0.4, thresholds_nolog_fr_0.9, thresholds_nolog_fr_0.975,
           exceedance_nolog_fr_0.4, exceedance_nolog_fr_0.9, exceedance_nolog_fr_0.975,
           peaks_test = peaks_test_fr, range_i.seasons = range_i.seasons,
           ylim_threshold = c(0, 500), approx_expectations = approx_expectations_nolog_fr,
           main = "natural scale, n = 30/m")

five_plots(thresholds_fr_0.4, thresholds_fr_0.9, thresholds_fr_0.975,
           exceedance_fr_0.4, exceedance_fr_0.9, exceedance_fr_0.975,
           peaks_test = peaks_test_fr, range_i.seasons = range_i.seasons,
           ylim_threshold = c(0, 500), approx_expectations = approx_expectations_fr,
           main = "log-transformed, n = 30 / m")

par(mar = c(0, 0, 1, 0))
legend_summary() # legend defined in separate function

dev.off()


# Plots on smoothing:
for(country in c("fr")){
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
pdf("../Draft/figure/plot_smoothing3_fr_small.pdf", width = 9, height = 8)

par(mar = c(4, 4, 3, 1), las = 1)
# layout(matrix(c(1:20, rep(21, 5)), ncol = 5, byrow = TRUE), heights = c(rep(1, 4), 0.7))
layout(matrix(c(1:10, rep(11, 5)), ncol = 5, byrow = TRUE), heights = c(rep(1, 4), 0.7))


# five_plots(thresholds1_nolog_smoothed3_fr_0.4, thresholds1_nolog_smoothed3_fr_0.9, thresholds1_nolog_smoothed3_fr_0.975,
#            exceedance1_nolog_raw3_fr_0.4, exceedance1_nolog_raw3_fr_0.9, exceedance1_nolog_raw3_fr_0.975,
#            peaks_test = peaks_test_fr, range_i.seasons = range_i.seasons,
#            ylim_threshold = c(0, 500),  approx_expectations = approx_expectations1_nolog_smoothed3_fr,
#            approx_sens = NULL, approx_spec = NULL, approx_ppv = NULL,
#            main = "natural scale, n = 1, l = 3, applied to unsmoothed peaks")

five_plots(thresholds1_smoothed3_fr_0.4, thresholds1_smoothed3_fr_0.9, thresholds1_smoothed3_fr_0.975,
           exceedance1_raw3_fr_0.4, exceedance1_raw3_fr_0.9, exceedance1_raw3_fr_0.975,
           peaks_test = peaks_test_fr, range_i.seasons = range_i.seasons,
           ylim_threshold = c(0, 500), approx_expectations = NULL,
           approx_sens = NULL, approx_spec = NULL, approx_ppv = NULL,
           main = "log-transformed, n = 1, l = 3, applied to unsmoothed peaks")

# five_plots(thresholds1_nolog_smoothed3_fr_0.4, thresholds1_nolog_smoothed3_fr_0.9, thresholds1_nolog_smoothed3_fr_0.975,
#            exceedance1_nolog_smoothed3_fr_0.4, exceedance1_nolog_smoothed3_fr_0.9, exceedance1_nolog_smoothed3_fr_0.975,
#            peaks_test = peaks_test_smoothed3_fr, range_i.seasons = range_i.seasons,
#            ylim_threshold = c(0, 500),  approx_expectations = approx_expectations1_nolog_smoothed3_fr,
#            approx_sens = approx_sens, approx_spec = approx_spec, approx_ppv = approx_ppv,
#            main = "natural scale, n = 1, l = 3, applied to smoothed peaks")

five_plots(thresholds1_smoothed3_fr_0.4, thresholds1_smoothed3_fr_0.9, thresholds1_smoothed3_fr_0.975,
           exceedance1_smoothed3_fr_0.4, exceedance1_smoothed3_fr_0.9, exceedance1_smoothed3_fr_0.975,
           peaks_test = peaks_test_smoothed3_fr, range_i.seasons = range_i.seasons,
           ylim_threshold = c(0, 500), approx_expectations = NULL,
           approx_sens = approx_sens, approx_spec = approx_spec, approx_ppv = approx_ppv,
           main = "log-transformed, n = 1, l = 3, applied to smoothed peaks")



legend_summary() # legend defined in separate function

dev.off()



# 7-day moving average:
pdf("../Draft/figure/plot_smoothing7_fr_small.pdf", width = 9, height = 8)

par(mar = c(4, 4, 3, 1), las = 1)
# layout(matrix(c(1:20, rep(21, 5)), ncol = 5, byrow = TRUE), heights = c(rep(1, 4), 0.7))
layout(matrix(c(1:10, rep(11, 5)), ncol = 5, byrow = TRUE), heights = c(rep(1, 4), 0.7))


# five_plots(thresholds1_nolog_smoothed7_fr_0.4, thresholds1_nolog_smoothed7_fr_0.9, thresholds1_nolog_smoothed7_fr_0.975,
#            exceedance1_nolog_raw7_fr_0.4, exceedance1_nolog_raw7_fr_0.9, exceedance1_nolog_raw7_fr_0.975,
#            peaks_test = peaks_test_fr, range_i.seasons = range_i.seasons,
#            ylim_threshold = c(0, 500),  approx_expectations = approx_expectations1_nolog_smoothed7_fr,
#            approx_sens = NULL, approx_spec = NULL, approx_ppv = NULL,
#            main = "natural scale, n = 1, l = 7, applied to unsmoothed peaks")

five_plots(thresholds1_smoothed7_fr_0.4, thresholds1_smoothed7_fr_0.9, thresholds1_smoothed7_fr_0.975,
           exceedance1_raw7_fr_0.4, exceedance1_raw7_fr_0.9, exceedance1_raw7_fr_0.975,
           peaks_test = peaks_test_fr, range_i.seasons = range_i.seasons,
           ylim_threshold = c(0, 500), approx_expectations = NULL,
           approx_sens = NULL, approx_spec = NULL, approx_ppv = NULL,
           main = "log-transformed, n = 1, l = 7, applied to unsmoothed peaks")

# five_plots(thresholds1_nolog_smoothed7_fr_0.4, thresholds1_nolog_smoothed7_fr_0.9, thresholds1_nolog_smoothed7_fr_0.975,
#            exceedance1_nolog_smoothed7_fr_0.4, exceedance1_nolog_smoothed7_fr_0.9, exceedance1_nolog_smoothed7_fr_0.975,
#            peaks_test = peaks_test_smoothed7_fr, range_i.seasons = range_i.seasons,
#            ylim_threshold = c(0, 500),  approx_expectations = approx_expectations1_nolog_smoothed7_fr,
#            approx_sens = approx_sens, approx_spec = approx_spec, approx_ppv = approx_ppv,
#            main = "natural scale, n = 1, l = 7, applied to smoothed peaks")

five_plots(thresholds1_smoothed7_fr_0.4, thresholds1_smoothed7_fr_0.9, thresholds1_smoothed7_fr_0.975,
           exceedance1_smoothed7_fr_0.4, exceedance1_smoothed7_fr_0.9, exceedance1_smoothed7_fr_0.975,
           peaks_test = peaks_test_smoothed7_fr, range_i.seasons = range_i.seasons,
           ylim_threshold = c(0, 500), approx_expectations = NULL,
           approx_sens = approx_sens, approx_spec = approx_spec, approx_ppv = approx_ppv,
           main = "log-transformed, n = 1, l = 7, applied to smoothed peaks")



legend_summary() # legend defined in separate function

dev.off()

# Plot for thresholds using CIs:

# load additional French results using confidence intervals:
load(paste0("Results/results_ci_fr.rda"))
load(paste0("Results/results1_ci_fr.rda"))
# load(paste0("Results/results_ci_nolog_fr.rda"))
# load(paste0("Results/results1_ci_nolog_fr.rda"))

summary_thresholds_ci_fr <- sim_summary(thresholds_ci_fr_0.4, thresholds_ci_fr_0.9, thresholds_ci_fr_0.975,
                                        range_i.seasons, interval = interv)
summary_exceedance_ci_fr <- sim_summary(exceedance_ci_fr_0.4, exceedance_ci_fr_0.9, exceedance_ci_fr_0.975,
                                        range_i.seasons, interval = interv)
summary_thresholds1_ci_fr <- sim_summary(thresholds1_ci_fr_0.4, thresholds1_ci_fr_0.9, thresholds1_ci_fr_0.975,
                                         range_i.seasons, interval = interv)
summary_exceedance1_ci_fr <- sim_summary(exceedance1_ci_fr_0.4, exceedance1_ci_fr_0.9, exceedance1_ci_fr_0.975,
                                         range_i.seasons, interval = interv)
# summary_thresholds_ci_nolog_fr <- sim_summary(thresholds_ci_nolog_fr_0.4, thresholds_ci_nolog_fr_0.9, thresholds_ci_nolog_fr_0.975,
#                                               range_i.seasons, interval = interv)
# summary_exceedance_ci_nolog_fr <- sim_summary(exceedance_ci_nolog_fr_0.4, exceedance_ci_nolog_fr_0.9, exceedance_ci_nolog_fr_0.975,
#                                               range_i.seasons, interval = interv)
# summary_thresholds1_ci_nolog_fr <- sim_summary(thresholds1_ci_nolog_fr_0.4, thresholds1_ci_nolog_fr_0.9, thresholds1_ci_nolog_fr_0.975,
#                                                range_i.seasons, interval = interv)
# summary_exceedance1_ci_nolog_fr <- sim_summary(exceedance1_ci_nolog_fr_0.4, exceedance1_ci_nolog_fr_0.9, exceedance1_ci_nolog_fr_0.975,
#                                                range_i.seasons, interval = interv)


pdf("../Draft/figure/plot_ci_fr.pdf", width = 8, height = 2.5)
par(mfrow = c(1, 4), mar = c(4, 4, 3, 1), las = 1)

plot_sim_summary(summary_thresholds1_ci_fr, xlab = "# included seasons m", ylab = "threshold", ylim = c(0, 300))
plot_exceedance_summary(summary_exceedance1_ci_fr, xlab = "# included seasons m", ylab = "average share")
mtext("based on CIs, log-transformed, n = 1", 3, at = 0, line = lin, cex = ce)

plot_sim_summary(summary_thresholds_ci_fr, xlab = "# included seasons m", ylab = "threshold", ylim = c(0, 300))
plot_exceedance_summary(summary_exceedance_ci_fr, xlab = "# included seasons m", ylab = "average share")
mtext("based on CIs, log-transformed, n = 30 / m", 3, at = 0, line = lin, cex = ce)

dev.off()


# Plots of misclassification matrices:

# compute matices for m = 5, 10, 15:
cm1_fr_5 <- classification_matrix(thresholds1_fr_0.4, thresholds1_fr_0.9, thresholds1_fr_0.975, peaks_test = peaks_test_fr, i.seasons = 5)
cm1_fr_10 <- classification_matrix(thresholds1_fr_0.4, thresholds1_fr_0.9, thresholds1_fr_0.975, peaks_test = peaks_test_fr, i.seasons = 10)
cm1_fr_15 <- classification_matrix(thresholds1_fr_0.4, thresholds1_fr_0.9, thresholds1_fr_0.975, peaks_test = peaks_test_fr, i.seasons = 15)

cm_fr_5 <- classification_matrix(thresholds_fr_0.4, thresholds_fr_0.9, thresholds_fr_0.975, peaks_test = peaks_test_fr, i.seasons = 5)
cm_fr_10 <- classification_matrix(thresholds_fr_0.4, thresholds_fr_0.9, thresholds_fr_0.975, peaks_test = peaks_test_fr, i.seasons = 10)
cm_fr_15 <- classification_matrix(thresholds_fr_0.4, thresholds_fr_0.9, thresholds_fr_0.975, peaks_test = peaks_test_fr, i.seasons = 15)

cm1_nolog_fr_5 <- classification_matrix(thresholds1_nolog_fr_0.4, thresholds1_nolog_fr_0.9, thresholds1_nolog_fr_0.975, peaks_test = peaks_test_fr, i.seasons = 5)
cm1_nolog_fr_10 <- classification_matrix(thresholds1_nolog_fr_0.4, thresholds1_nolog_fr_0.9, thresholds1_nolog_fr_0.975, peaks_test = peaks_test_fr, i.seasons = 10)
cm1_nolog_fr_15 <- classification_matrix(thresholds1_nolog_fr_0.4, thresholds1_nolog_fr_0.9, thresholds1_nolog_fr_0.975, peaks_test = peaks_test_fr, i.seasons = 15)

cm_nolog_fr_5 <- classification_matrix(thresholds_nolog_fr_0.4, thresholds_nolog_fr_0.9, thresholds_nolog_fr_0.975, peaks_test = peaks_test_fr, i.seasons = 5)
cm_nolog_fr_10 <- classification_matrix(thresholds_nolog_fr_0.4, thresholds_nolog_fr_0.9, thresholds_nolog_fr_0.975, peaks_test = peaks_test_fr, i.seasons = 10)
cm_nolog_fr_15 <- classification_matrix(thresholds_nolog_fr_0.4, thresholds_nolog_fr_0.9, thresholds_nolog_fr_0.975, peaks_test = peaks_test_fr, i.seasons = 15)


# Plot for thresholds without smoothing:

pdf("../Draft/figure/mosaic_fr.pdf", width = 9, height = 11)
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
mosaic(cm1_nolog_fr_5)
mtext("m = 5", cex = 0.9)
mosaic(cm1_nolog_fr_10)
mtext("m = 10", cex = 0.9)
mtext("natural scale, n = 1", cex = 1, line = 2)
mosaic(cm1_nolog_fr_15)
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
mosaic(cm1_fr_5)
mtext("m = 5", cex = 0.9)
mosaic(cm1_fr_10)
mtext("m = 10", cex = 0.9)
mtext("log-transformed, n = 1", cex = 1, line = 2)
mosaic(cm1_fr_15)
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
mosaic(cm_nolog_fr_5)
mtext("m = 5", cex = 0.9)
mosaic(cm_nolog_fr_10)
mtext("m = 10", cex = 0.9)
mtext("natural scale, n = 30/m", cex = 1, line = 2)
mosaic(cm_nolog_fr_15)
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
mosaic(cm_fr_5)
mtext("m = 5", cex = 0.9)
mosaic(cm_fr_10)
mtext("m = 10", cex = 0.9)
mtext("log-transformed, n = 30/m", cex = 1, line = 2)
mosaic(cm_fr_15)
mtext("m = 15", cex = 0.9)

par(mar = c(0, 0, 0, 0), las = 1)
plot(NULL, xlim = 0:1, ylim = 0:1, xlab = "", ylab = "", axes = FALSE)
legend("top", legend = c("classification:   ", "low", "medium", "high", "very high"),
       pt.bg = c(NA, col_low, col_medium, col_high, col_very_high),
       col = c(NA, rep("black", 4)), pch = 22, ncol = 5, bty = "n", pt.cex = 2)

dev.off()




# Plot for thresholds with smoothing and log trafo:

# compute matices for m = 5, 10, 15:
cm1_smoothed3_raw_fr_5 <- classification_matrix(thresholds1_smoothed3_fr_0.4, thresholds1_smoothed3_fr_0.9, thresholds1_smoothed3_fr_0.975, peaks_test = peaks_test_fr, i.seasons = 5)
cm1_smoothed3_raw_fr_10 <- classification_matrix(thresholds1_smoothed3_fr_0.4, thresholds1_smoothed3_fr_0.9, thresholds1_smoothed3_fr_0.975, peaks_test = peaks_test_fr, i.seasons = 10)
cm1_smoothed3_raw_fr_15 <- classification_matrix(thresholds1_smoothed3_fr_0.4, thresholds1_smoothed3_fr_0.9, thresholds1_smoothed3_fr_0.975, peaks_test = peaks_test_fr, i.seasons = 15)


cm1_smoothed3_smoothed3_fr_5 <- classification_matrix(thresholds1_smoothed3_fr_0.4, thresholds1_smoothed3_fr_0.9, thresholds1_smoothed3_fr_0.975, peaks_test = peaks_test_smoothed3_fr, i.seasons = 5)
cm1_smoothed3_smoothed3_fr_10 <- classification_matrix(thresholds1_smoothed3_fr_0.4, thresholds1_smoothed3_fr_0.9, thresholds1_smoothed3_fr_0.975, peaks_test = peaks_test_smoothed3_fr, i.seasons = 10)
cm1_smoothed3_smoothed3_fr_15 <- classification_matrix(thresholds1_smoothed3_fr_0.4, thresholds1_smoothed3_fr_0.9, thresholds1_smoothed3_fr_0.975, peaks_test = peaks_test_smoothed3_fr, i.seasons = 15)



pdf("../Draft/figure/mosaic_log_smoothed_fr.pdf", width = 9, height = 6)
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
mosaic(cm1_smoothed3_raw_fr_5)
mtext("m = 5", cex = 0.9)
mosaic(cm1_smoothed3_raw_fr_10)
mtext("m = 10", cex = 0.9)
mtext("log-transformed, n = 1, l = 3, applied to unsmoothed peaks", cex = 1, line = 2)
mosaic(cm1_smoothed3_raw_fr_15)
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
mosaic(cm1_smoothed3_smoothed3_fr_5)
mtext("m = 5", cex = 0.9)
mosaic(cm1_smoothed3_smoothed3_fr_10)
mtext("m = 10", cex = 0.9)
mtext("log-transformed, n = 1, l = 3, applied to smoothed peaks", cex = 1, line = 2)
mosaic(cm1_smoothed3_smoothed3_fr_15)
mtext("m = 15", cex = 0.9)


par(mar = c(0, 0, 0, 0), las = 1)
plot(NULL, xlim = 0:1, ylim = 0:1, xlab = "", ylab = "", axes = FALSE)
legend("top", legend = c("classification:   ", "low", "medium", "high", "very high"),
       pt.bg = c(NA, col_low, col_medium, col_high, col_very_high),
       col = c(NA, rep("black", 4)), pch = 22, ncol = 5, bty = "n", pt.cex = 2)

dev.off()
