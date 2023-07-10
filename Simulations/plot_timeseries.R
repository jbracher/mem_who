# Descriptive plots of data

# set working directory
current_path = rstudioapi::getActiveDocumentContext()$path # get path of this file
setwd(dirname(current_path))

# library for smoothing:
library(zoo)

# a helper function for plotting:
plot_mem_ts <- function(matr, xlim, ylab = "inc / 100,000", ylim = c(0, 400), ...){
  # get years as numeric:
  years <- as.numeric(substr(colnames(matr), start = 1, stop = 4))
  # labels for axes:
  labels_years <- gsub(".", "/", substr(colnames(matr), start = 8, stop = 16), fixed = TRUE)
  
  plot(NULL, xlim = xlim, ylim = ylim, xlab = "", ylab = ylab, axes = FALSE,...)
  # nicer axes:
  axis(1, at = years + 0.5, labels = labels_years)
  axis(2)
  box()
  
  # vertical eparator lines:
  abline(v = seq(from = floor(xlim[1]), to = ceiling(xlim[2])), col = "grey")

  # draw lines for seasons:
  weeks_per_seas <- nrow(matr)
  for(i in 1:ncol(matr)){
    x <- years[i] + (1:weeks_per_seas)/weeks_per_seas
    lines(x, matr[, i])
  }
}

# function for boxplot of observations by rank within season
bplot_ranks <- function(matr, ylab = "", main = ""){
  plot(NULL, xlim = c(0.5, 6.5), ylim = c(0, 1.1*max(matr)), xlab = "", ylab = ylab)
  mtext("value", 2, line = 2.6, cex = 0.6, las = 0)
  mtext("rank in season", 1, line = 1.8, cex = 0.6)
  mtext(main, 3, line = 0.7, cex = 0.7)
  get_nth <- function(vect, n) sort(vect, decreasing = TRUE)[n]
  for(i in 1:6) boxplot(apply(matr, 2, get_nth, n = i), at = i, add = TRUE, axes = FALSE)
}

# function for QQ plot of peak values:
qqnorm_peaks <- function(matr, ylab = "", lim, main = "", ...){
  peaks <- apply(matr, 2, max)
  qqnorm((peaks - mean(peaks))/sd(peaks), xlim = lim, ylim = lim, pch = 15, cex = 0.5,
         xlab = "", ylab = "", main = "", ...)
  mtext("theoretical quantiles", 1, line = 1.8, cex = 0.6)
  mtext("sample quantiles", 2, line = 2.6, cex = 0.6, las = 0)
  mtext(main, 3, line = 0.7, cex = 0.7)
  abline(0:1, col = "lightgrey")
}

# function for boxplot of unsmoothed and smoothed peak values:
bplot_smoothed <- function(matr, ylab = "", xlab = "", main = ""){
  matr_smoothed3 <- data.frame(apply(matr, MARGIN = 2, FUN = rollmean, k = 3))
  matr_smoothed7 <- data.frame(apply(matr, MARGIN = 2, FUN = rollmean, k = 7))
  plot(NULL, xlim = c(0.5, 3.5), ylim = c(0, 1.1*max(matr)), xlab = "", ylab = ylab,
       axes = FALSE)
  axis(1, at = 1:3, labels = c(1, 3, 7))
  axis(2)
  box()
  mtext("smoothing window l", 1, line = 1.8, cex = 0.6)
  mtext("season peak value", 2, line = 2.6, cex = 0.6, las = 0)
  mtext(main, 3, line = 0.7, cex = 0.7)
  get_nth <- function(vect, n) sort(vect, decreasing = TRUE)[n]
  boxplot(apply(matr, 2, get_nth, n = 1), at = 1, add = TRUE, axes = FALSE)
  boxplot(apply(matr_smoothed3, 2, get_nth, n = 1), at = 2, add = TRUE, axes = FALSE)
  boxplot(apply(matr_smoothed7, 2, get_nth, n = 1), at = 3, add = TRUE, axes = FALSE)
}



#########################
# France

# read in data:
dat_fr <- read.csv("../Data/for_mem/ili_mem_fr.csv")

# extract data from two regions (Normandy, Nouvelle Aquitaine)
dat_nouvelle_aquitaine <- dat_fr[, grepl("NOUVELLE.AQUITAINE", colnames(dat_fr))]
colnames(dat_nouvelle_aquitaine) <- gsub("NOUVELLE.AQUITAINE_", "", colnames(dat_nouvelle_aquitaine))
dat_grand_est <- dat_fr[, grepl("GRAND.EST", colnames(dat_fr))]
colnames(dat_grand_est) <- gsub("GRAND.EST_", "", colnames(dat_grand_est))

# generate plot:
pdf("../Draft/figure/plot_data_fr.pdf", width = 9, height = 5.5)

# structure plot area:
par(las = 1, mar = c(3, 4, 1, 1))
layout(matrix(c(1, 1, 1, 1,
                2, 2, 2, 2,
                3, 4, 5, 6), ncol = 4, byrow = TRUE))

# time series Grand Est:
plot_mem_ts(dat_grand_est, xlim = c(1985.5, 2018.5), ylim = c(0, 300))
legend("top", legend = " Grand Est  ", bty = "n", cex = 1.2, bg = "white")

# time series Nouvelle Aquitaine:
plot_mem_ts(dat_nouvelle_aquitaine, xlim = c(1985.5, 2018.5), ylim = c(0, 300))
legend("top", legend = " Nouvelle Aquitaine  ", bty = "n", cex = 1.2, bg = "white")

# boxplots:
bplot_ranks(dat_fr, main = "Distribution by rank in season")
bplot_smoothed(dat_fr, main = "Distribution of smoothed peaks")

# QQ plots:
qqnorm_peaks(dat_fr, lim = c(-3.5, 3.5), main = "Untransformed peaks")
qqnorm_peaks(log(dat_fr), lim = c(-3.5, 3.5), main = "Log-transformed peaks")

dev.off()





#########################
# US

# read in data:
dat_us <- read.csv("../Data/for_mem/ili_mem_us.csv")

# extract data from two regions (Normandy, Nouvelle Aquitaine)
dat_region1 <- dat_us[, grepl("Region1_", colnames(dat_us))]
colnames(dat_region1) <- gsub("Region1_", "", colnames(dat_region1))
dat_region7 <- dat_us[, grepl("Region7_", colnames(dat_us))]
colnames(dat_region7) <- gsub("Region7_", "", colnames(dat_region7))

# generate plot:
pdf("../Draft/figure/plot_data_us.pdf", width = 9, height = 5.5)

# structure plot area:
par(las = 1, mar = c(3, 4, 1, 1))
layout(matrix(c(1, 1, 1, 1,
                2, 2, 2, 2,
                3, 4, 5, 6), ncol = 4, byrow = TRUE))

# time series Region 1:
plot_mem_ts(dat_region1, xlim = c(1998.5, 2017.5), ylim = c(0, 300))
legend("top", legend = " HHS Region 1  ", bty = "n", cex = 1.2, bg = "white")

# time series Region 7:
plot_mem_ts(dat_region7, xlim = c(1998.5, 2017.5), ylim = c(0, 300))
legend("top", legend = " HHS Region 7  ", bty = "n", cex = 1.2, bg = "white")


# boxplots:
bplot_ranks(dat_us, main = "Distribution by rank in season")
bplot_smoothed(dat_us, main = "Distribution of smoothed peaks")

# QQ plots:
qqnorm_peaks(dat_us, lim = c(-3.5, 3.5), main = "Untransformed peaks")
qqnorm_peaks(log(dat_us), lim = c(-3.5, 3.5), main = "Log-transformed peaks")

dev.off()

