# Descriptive plots of data

# set working directory
current_path = rstudioapi::getActiveDocumentContext()$path # get path of this file
setwd(dirname(current_path))

# library for smoothing:
library(zoo)

# define colors for plotting
col_very_high <- "purple4"
col_high <- "tomato2"
col_medium <- rgb(255/255, 178/255, 104/255) # light orange
col_low <- rgb(255/255, 255/255, 180/255) # light red

# read in data:
dat_fr <- read.csv("../Data/for_mem/ili_mem_fr.csv")

# extract data from two regions (Normandy, Nouvelle Aquitaine)
dat_nouvelle_aquitaine <- dat_fr[, grepl("NOUVELLE.AQUITAINE", colnames(dat_fr))]
colnames(dat_nouvelle_aquitaine) <- gsub("NOUVELLE.AQUITAINE_", "", colnames(dat_nouvelle_aquitaine))
dat_grand_est <- dat_fr[, grepl("GRAND.EST", colnames(dat_fr))]
colnames(dat_grand_est) <- gsub("GRAND.EST_", "", colnames(dat_grand_est))

# the portion of the data we will show: last seven years in Grand Est
da <- dat_grand_est[, ncol(dat_grand_est) - (6:0)]
years <- as.numeric(substr(colnames(da), start = 1, stop = 4))
weeks_per_seas <- nrow(da)

# helper function to plot thresholds:
# arguments:
# thresh: a vector with three elements for the medium, high and very high thresholds
# xleft, xright: x-boundaries of rectangle
# ytop: upper end of top rectangle
threshold_rectangles <- function(thresh, xleft, xright, ytop){
  rect(xleft = xleft, xright = xright, ybottom = 0, ytop = ytop, 
       col = col_very_high, border = NA)
  rect(xleft = xleft, xright = xright, ybottom = 0, ytop = thresh[3], 
       col = col_high, border = NA)
  rect(xleft = xleft, xright = xright, ybottom = 0, ytop = thresh[2], 
       col = col_medium, border = NA)
  rect(xleft = xleft, xright = xright, ybottom = 0, ytop = thresh[1], 
       col = col_low, border = NA)
}

# compute MEM and WHO thresholds for example data:
library(mem)

thresh_mem0 <- memmodel(da[, 1:6], i.seasons = 6,
                       i.level.intensity = c(0.4, 0.9, 0.975))
thresh_mem <- thresh_mem0$intensity.thresholds

# WHO needs smoothed data
da_smoothed <- data.frame(apply(da, MARGIN = 2, FUN = rollmean, k = 3))
thresh_who0 <- memmodel(da_smoothed[, 1:6], i.seasons = 6,
                       i.level.intensity = c(0.4, 0.9, 0.975),
                       i.n.max = 1, i.type.intensity = 5)
thresh_who <- thresh_who0$intensity.thresholds


###############
# Figure 1: Comparison of MEM and WHO thresholds
pdf("../Draft/figure/illustration_mem_who.pdf", width = 7, height = 3.5)

# structure plot area
layout(matrix(c(1:4, 5, 5, 5, 5), ncol = 4, byrow = TRUE), widths = c(4, 1, 4, 1), heights = c(4, 1))

# MEM thresholds:

# plot of time series and thresholds:
par(las = 1, mar = c(2.5, 4.1, 3, 0))
yl <- c(0, 1.3) * range(da)
plot(NULL, xlim = range(years) + 0:1, ylim = yl, xlab = "", ylab = "re-scaled incidence", axes = FALSE)
axis(2)
axis(1, at = 2012:2018 + 0.5, labels = 2013:2019)
box()
# add thresholds:
threshold_rectangles(thresh = thresh_mem, xleft = 2012, xright = 2019, ytop = yl[2])
# grey out past seasons:
rect(xleft = 2012, xright = 2018, ybottom = 0, ytop = yl[2], col = rgb(1, 1, 1, 0.5),
     border = NA)
abline(v = c(years, max(years) + 1), col = "lightgrey")

# add time series:
reference_set <- numeric() # vector to collect values in reference set
new_peak <- NA # to store new peak value
for(i in 1:ncol(da)){
  # lines
  y <- da[, i]
  x <- years[i] + (1:weeks_per_seas)/weeks_per_seas
  lty <- ifelse(i == ncol(da), "dotted", "solid")
  lines(x, y, lty = lty)
  # points for reference set
  if(i != ncol(da)){
    ord <- order(y, decreasing = TRUE)
    top_weeks <- x[ord][1:5]
    top_values <- y[ord][1:5]
    points(top_weeks, top_values, pch = 16, cex = 0.9)
    reference_set <- c(reference_set, top_values)
  }else{
    ord <- order(y, decreasing = TRUE)
    top_weeks <- x[ord][1]
    top_values <- y[ord][1]
    points(top_weeks, top_values, pch = 4, cex = 1.3)
    new_peak <- top_values
  }
}
# title:
mtext("MEM default (log-transformed, n = 30/m, l = 1)",
      cex = 0.8, line = 0.7)

# add vertical density plot manually:
par(mar = c(2.5, 0, 3, 0))
plot(NULL, ylim = yl, xlim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
# compute parameters
log_mu <- mean(log(reference_set))
log_sd <- sd(log(reference_set))
y <- seq(from = 0, to = yl[2], by = 1)
# compute density
dens <- 0.5*dnorm(log(y), mean = log_mu, sd = log_sd)
# add polygons: very high
polygon(c(0, dens, 0), c(y[1], y, y[length(y)]), col = col_very_high, border = NA)
# high
inds_high <- which(y < thresh_mem[3])
polygon(c(0, dens[inds_high], 0), 
        c(y[inds_high[1]], y[inds_high], y[inds_high[length(inds_high)]]), col = col_high, border = NA)
# medium
inds_medium <- which(y < thresh_mem[2])
polygon(c(0, dens[inds_medium], 0), 
        c(y[inds_medium[1]], y[inds_medium], y[inds_medium[length(inds_medium)]]), col = col_medium, border = NA)
# low
inds_low <- which(y < thresh_mem[1])
polygon(c(0, dens[inds_low], 0), 
        c(y[inds_low[1]], y[inds_low], y[inds_low[length(inds_low)]]), col = col_low, border = NA)
# add points for reference set and new peak
points(rep(0.02, length(reference_set)), reference_set, pch = 16, cex = 0.9)
points(0.02, new_peak, pch = 4, cex = 1.3)

# black contour:
lines(dens, y, col = "darkgrey")

# labelling
text(x = rep(0.5, 4), y = (c(0, thresh_mem) + c(thresh_mem, yl[2]))/2, 
     labels = c("low", "medium", "high", "very high"), cex = 0.85)



# WHO thresholds:

par(las = 1, mar = c(2.5, 4.1, 3, 0))
yl <- c(0, 1.3) * range(da)
plot(NULL, xlim = range(years) + 0:1, ylim = yl, xlab = "", ylab = "re-scaled incidence", axes = FALSE)
axis(2)
axis(1, at = 2012:2018 + 0.5, labels = 2013:2019)
box()
# add thresholds:
threshold_rectangles(thresh = thresh_who, xleft = 2012, xright = 2019, ytop = yl[2])
# grey out past seasons:
rect(xleft = 2012, xright = 2018, ybottom = 0, ytop = yl[2], col = rgb(1, 1, 1, 0.5),
     border = NA)

abline(v = c(years, max(years) + 1), col = "lightgrey")

# add time series:
reference_set <- numeric()
for(i in 1:ncol(da)){
  # lines, incl smoothed version
  y <- da[, i]
  y_smo <- c(y[1], rollmean(y, k = 3), tail(y, 1))
  x <- years[i] + (1:weeks_per_seas)/weeks_per_seas
  lty <- ifelse(i == ncol(da), "dotted", "dashed")
  # lines(x, y, lty = lty)
  
  # points
  if(i != ncol(da)){
    lines(x, y_smo, lty = lty, col = "black")
    ord <- order(y_smo, decreasing = TRUE)
    top_weeks <- x[ord][1]
    top_values <- y_smo[ord][1]
    points(top_weeks, top_values, pch = 16, cex = 0.9, col = "black")
    reference_set <- c(reference_set, top_values)
  }else{
    lines(x, y, lty = lty, col = "black")
    ord <- order(y, decreasing = TRUE)
    top_weeks <- x[ord][1]
    top_values <- y[ord][1]
    points(top_weeks, top_values, pch = 4, cex = 0.9)
    new_peak <- top_values
  }
}

# title:
mtext("WHO default (natural scale, n = 1, l = 3)", 
      cex = 0.8, line = 0.7)

# add vertical density plot manually:
par(mar = c(2.5, 0, 3, 0))
plot(NULL, ylim = yl, xlim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
# compute parameters
mu <- mean(reference_set)
sd <- sd(reference_set)
y <- seq(from = 0, to = yl[2], by = 1)
# compute density
dens <- 30*dnorm(y, mean = mu, sd = sd)
# add polygons: very high
polygon(c(0, dens, 0), c(y[1], y, y[length(y)]), col = col_very_high, border = NA)
# high
inds_high <- which(y < thresh_who[3])
polygon(c(0, dens[inds_high], 0), 
        c(y[inds_high[1]], y[inds_high], y[inds_high[length(inds_high)]]), col = col_high, border = NA)
# medium
inds_medium <- which(y < thresh_who[2])
polygon(c(0, dens[inds_medium], 0), 
        c(y[inds_medium[1]], y[inds_medium], y[inds_medium[length(inds_medium)]]), col = col_medium, border = NA)
# low
inds_low <- which(y < thresh_who[1])
polygon(c(0, dens[inds_low], 0), 
        c(y[inds_low[1]], y[inds_low], y[inds_low[length(inds_low)]]), col = col_low, border = NA)
# add points for reference set and new peak
points(rep(0.02, length(reference_set)), reference_set, pch = 16, cex = 0.9, col = "black")
points(0.02, new_peak, pch = 4, cex = 1.3)

# black contour:
lines(dens, y, col = "darkgrey")

# labelling
text(x = rep(0.5, 4), y = (c(0, thresh_who) + c(thresh_who, yl[2]))/2, 
     labels = c("low", "medium", "high", "very high"), cex = 0.85)

# legend:
par(mar = c(0, 0, 0, 0))
plot(NULL, axes = FALSE, xlab = "", ylab = "", xlim = 0:1, ylim = 0:1)
legend("center",
       legend = c("unsmoothed historic incidence",
                  "smoothed historic incidence",
                  "unsmoothed new season",
                  "values in reference set",
                  "new season peak"),
       ncol = 2,
       col = c("black", "black", "black", "black", "black"),
       lty = c(1, 2, 3, NA, NA),
       pch = c(NA, NA, NA, 16, 4),
       bty = "n")
# points(0.53, 0.7, pch = 19, col = "darkgrey", cex = 0.9)
# lines(c(0.28, 0.31), c(0.2, 0.2), col = "darkgrey", lty = 3)
dev.off()


###############
### Figure 2 on t-distribution

# get functions for t-distribution an trend
source("mem_lm.R")
library("ciTools")

# compute thresholds under t and normal
thresh_mem_t0 <- memmodel_lm(da[, 1:6], i.seasons = 6,
                             i.level.intensity = c(0.4, 0.9, 0.975), 
                             i.n.max = 1, include_trend = FALSE,
                             include_rank = FALSE, i.trafo = "log")
thresh_mem_t <- thresh_mem_t0$intensity.thresholds

thresh_mem_normal0 <- memmodel(da[, 1:6], i.seasons = 6,
                        i.level.intensity = c(0.4, 0.9, 0.975), 
                        i.n.max = 1)
thresh_mem_normal <- thresh_mem_normal0$intensity.thresholds

# Plot:
pdf("../Draft/figure/illustration_t_normal.pdf", width = 7, height = 3.5)

# structure plot area
layout(matrix(c(1:4, 5, 5, 5, 5), ncol = 4, byrow = TRUE), widths = c(4, 1, 4, 1), heights = c(4, 1))

### normal distribution

# plot of time series and thresholds:
par(las = 1, mar = c(2.5, 4.1, 3, 0))
yl <- c(0, 2) * range(da)
plot(NULL, xlim = range(years) + 0:1, ylim = yl, xlab = "", ylab = "re-scaled incidence", axes = FALSE)
axis(2)
axis(1, at = 2012:2018 + 0.5, labels = 2013:2019)
box()
# add thresholds:
threshold_rectangles(thresh = thresh_mem_normal, xleft = 2012, xright = 2019, ytop = yl[2])
# grey out past seasons:
rect(xleft = 2012, xright = 2018, ybottom = 0, ytop = yl[2], col = rgb(1, 1, 1, 0.5),
     border = NA)
abline(v = c(years, max(years) + 1), col = "lightgrey")

# add time series:
reference_set <- numeric() # to collect values in reference set:
for(i in 1:ncol(da)){
  # lines
  y <- da[, i]
  x <- years[i] + (1:weeks_per_seas)/weeks_per_seas
  lty <- ifelse(i == ncol(da), "dotted", "solid")
  lines(x, y, lty = lty)
  # points for reference set:
  if(i != ncol(da)){
    ord <- order(y, decreasing = TRUE)
    top_weeks <- x[ord][1]
    top_values <- y[ord][1]
    points(top_weeks, top_values, pch = 16, cex = 0.9, col = "black")
    reference_set <- c(reference_set, top_values)
  }else{
    ord <- order(y, decreasing = TRUE)
    top_weeks <- x[ord][1]
    top_values <- y[ord][1]
    points(top_weeks, top_values, pch = 4, cex = 0.9)
    new_peak <- top_values
  }
}

# title:
mtext("log-transformed, n = 1, normal distribution", 
      cex = 0.8, line = 0.7)

# add vertical density plot manually:
par(mar = c(2.5, 0, 3, 0))
plot(NULL, ylim = yl, xlim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
y <- seq(from = 0, to = yl[2], by = 1)
# compute parameters:
log_mu <- mean(log(reference_set))
log_sd <- sd(log(reference_set))
# compute density:
dens <- 0.5*dnorm(log(y), mean = log_mu, sd = log_sd)
# add polygons: very high
polygon(c(0, dens, 0), c(y[1], y, y[length(y)]), col = col_very_high, border = NA)
# high
inds_high <- which(y < thresh_mem_normal[3])
polygon(c(0, dens[inds_high], 0), 
        c(y[inds_high[1]], y[inds_high], y[inds_high[length(inds_high)]]), col = col_high, border = NA)
# medium
inds_medium <- which(y < thresh_mem_normal[2])
polygon(c(0, dens[inds_medium], 0), 
        c(y[inds_medium[1]], y[inds_medium], y[inds_medium[length(inds_medium)]]), col = col_medium, border = NA)
# low
inds_low <- which(y < thresh_mem_normal[1])
polygon(c(0, dens[inds_low], 0), 
        c(y[inds_low[1]], y[inds_low], y[inds_low[length(inds_low)]]), col = col_low, border = NA)
# points for reference set:
points(rep(0.02, length(reference_set)), reference_set, pch = 16, cex = 0.9, col = "black")
points(0.02, new_peak, pch = 4, cex = 1.3)

# black contour:
lines(dens, y, col = "darkgrey")

# labels
text(x = rep(0.5, 4), y = (c(0, thresh_mem_normal) + c(thresh_mem_normal, yl[2]))/2, 
     labels = c("low", "medium", "high", "very high"), cex = 0.85)


### t-distribution

# plot of time series and thresholds:
par(las = 1, mar = c(2.5, 4.1, 3, 0))
yl <- c(0, 2) * range(da)
plot(NULL, xlim = range(years) + 0:1, ylim = yl, xlab = "", ylab = "re-scaled incidence", axes = FALSE)
axis(2)
axis(1, at = 2012:2018 + 0.5, labels = 2013:2019)
box()
# add thresholds:
threshold_rectangles(thresh = thresh_mem_t, xleft = 2012, xright = 2019, ytop = yl[2])
# grey out past seasons:
rect(xleft = 2012, xright = 2018, ybottom = 0, ytop = yl[2], col = rgb(1, 1, 1, 0.5),
     border = NA)
abline(v = c(years, max(years) + 1), col = "lightgrey")

reference_set <- numeric() # to collect reference set
# add time series:
for(i in 1:ncol(da)){
  # lines:
  y <- da[, i]
  x <- years[i] + (1:weeks_per_seas)/weeks_per_seas
  lty <- ifelse(i == ncol(da), "dotted", "solid")
  lines(x, y, lty = lty)
  # points for reference set:
  if(i != ncol(da)){
    ord <- order(y, decreasing = TRUE)
    top_weeks <- x[ord][1]
    top_values <- y[ord][1]
    points(top_weeks, top_values, pch = 16, cex = 0.9, col = "black")
    reference_set <- c(reference_set, top_values)
  }else{
    ord <- order(y, decreasing = TRUE)
    top_weeks <- x[ord][1]
    top_values <- y[ord][1]
    points(top_weeks, top_values, pch = 4, cex = 0.9)
    new_peak <- top_values
  }
}

# title
mtext("log-transformed, n = 1, t-distribution", 
      cex = 0.8, line = 0.7)

# add vertical density plot manually:
par(mar = c(2.5, 0, 3, 0))
plot(NULL, ylim = yl, xlim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
# compute parameters:
log_mu <- mean(log(reference_set))
log_sd <- sd(log(reference_set))
y <- seq(from = 0, to = yl[2], by = 1)
# compute density
dens <- dt((log(y) - log_mu)/log_sd/sqrt(1 + 1/6), df = 5)
# add polygons: very high
polygon(c(0, dens, 0), c(y[1], y, y[length(y)]), col = col_very_high, border = NA)
# high
inds_high <- which(y < thresh_mem_t[3])
polygon(c(0, dens[inds_high], 0), 
        c(y[inds_high[1]], y[inds_high], y[inds_high[length(inds_high)]]), col = col_high, border = NA)
# medium
inds_medium <- which(y < thresh_mem_t[2])
polygon(c(0, dens[inds_medium], 0), 
        c(y[inds_medium[1]], y[inds_medium], y[inds_medium[length(inds_medium)]]), col = col_medium, border = NA)
# low
inds_low <- which(y < thresh_mem_t[1])
polygon(c(0, dens[inds_low], 0), 
        c(y[inds_low[1]], y[inds_low], y[inds_low[length(inds_low)]]), col = col_low, border = NA)
# points for reference set:
points(rep(0.02, length(reference_set)), reference_set, pch = 16, cex = 0.9, col = "black")
points(0.02, new_peak, pch = 4, cex = 1.3)

# black contour:
lines(dens, y, col = "darkgrey")

# labelling
text(x = rep(0.5, 4), y = (c(0, thresh_mem_t) + c(thresh_mem_t, yl[2]))/2, 
     labels = c("low", "medium", "high", "very high"), cex = 0.85)

# legend
par(mar = c(0, 0, 0, 0))
plot(NULL, axes = FALSE, xlab = "", ylab = "", xlim = 0:1, ylim = 0:1)
legend("center",
       legend = c("historic incidence",
                  "new season",
                  "values in reference set",
                  "new season peak"),
       ncol = 2,
       col = c("black", "black", "black", "black"),
       lty = c(1, 3, NA, NA),
       pch = c(NA, NA, 16, 4),
       bty = "n")

dev.off()


#############
# Figure 3 on trends:

# use full data set for trend example
da_trend <- dat_grand_est
years_trend <- as.numeric(substr(colnames(da_trend), start = 1, stop = 4))

# compute thresholds:
thresh_mem_trend0 <- memmodel_lm(da_trend, i.seasons = ncol(da_trend) - 1,
                             i.level.intensity = c(0.4, 0.9, 0.975), 
                             i.n.max = 1, include_trend = TRUE,
                             include_rank = FALSE, i.trafo = "log")
thresh_mem_trend <- thresh_mem_t0$intensity.thresholds

# thresholds over time need to be computed manually:
mod <- thresh_mem_trend0$fit
newdata <- data.frame(t = 1:ncol(da_trend))
q0.4 <- add_quantile(newdata, fit = mod, p = 0.4)
q0.9 <- add_quantile(newdata, fit = mod, p = 0.9)
q0.975 <- add_quantile(newdata, fit = mod, p = 0.975)
q_matr <- cbind(q0.4$quantile0.4, q0.9$quantile0.9, q0.975$quantile0.975)


# Plot:
pdf("../Draft/figure/illustration_trend.pdf", width = 9, height = 4.5)

# structure plot area:
layout(matrix(1:2, ncol = 1), heights = c(4, 1))
par(las = 1, mar = c(2.5, 4.1, 3, 3))
yl <- c(0, 3) * range(da)
plot(NULL, xlim = range(years_trend), ylim = yl, xlab = "", ylab = "re-scaled incidence", axes = FALSE)
axis(2)
labels_trend <- substr(as.character(years_trend + 1), start = 3, stop = 4)
axis(1, at = years_trend + 0.5, labels = labels_trend)
box()
# add thresholds:
for(i in 1:ncol(da_trend)){
  threshold_rectangles(thresh = exp(q_matr[i, ]), xleft = years_trend[i], xright = years_trend[i] + 1, ytop = yl[2])
}
# grey out past seasons:
rect(xleft = 1985, xright = 2018, ybottom = 0, ytop = yl[2], col = rgb(1, 1, 1, 0.5),
     border = NA)
abline(v = c(years_trend, max(years_trend) + 1), col = "lightgrey")

# add time series:
for(i in 1:ncol(da_trend)){
  # lines:
  y <- da_trend[, i]
  x <- years_trend[i] + (1:weeks_per_seas)/weeks_per_seas
  lty <- ifelse(i == ncol(da_trend), "dotted", "solid")
  lines(x, y, lty = lty)
  # points:
  if(i != ncol(da_trend)){
    ord <- order(y, decreasing = TRUE)
    top_weeks <- x[ord][1]
    top_values <- y[ord][1]
    points(top_weeks, top_values, pch = 16, cex = 0.9, col = "black")
    reference_set <- c(reference_set, top_values)
  }else{
    ord <- order(y, decreasing = TRUE)
    top_weeks <- x[ord][1]
    top_values <- y[ord][1]
    points(top_weeks, top_values, pch = 4, cex = 0.9)
    new_peak <- top_values
  }
}
# title
mtext("log-transformed, n = 1, t-distribution, accounting for trend")

# legend:
par(mar = c(0, 0, 0, 0))
plot(NULL, axes = FALSE, xlab = "", ylab = "", xlim = 0:1, ylim = 0:1)
legend("center",
       legend = c("historic incidence",
                  "new season",
                  "values in reference set",
                  "new season peak"),
       ncol = 2,
       col = c("black", "black", "black", "black"),
       lty = c(1, 3, NA, NA),
       pch = c(NA, NA, 16, 4),
       bty = "n")

dev.off()
