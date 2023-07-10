# define colours:
col_very_high <- "darkorchid4"
col_high <- "red"
col_medium <- "orange"
col_low <- "yellow"


# function to summarize simulation results
# results_0.4: simulation results for medium thereshold
# results_0.9: simulation results for high thereshold
# results_0.975: simulation results for very high thereshold
# range_i.seasons: vector containing the number of seasons included (refers to different entries in results_0.4 etc)
# interval: which uncertainty intervals are to be computed?
sim_summary <- function(results_0.4, results_0.9, results_0.975, range_i.seasons,
                        interval = c(0.025, 0.975)){
  # check input:
  if(length(interval) != 2) stop("interval needs to be of length 2.")
  interval <- sort(interval)

  # initialize return list:
  ret <- list()

  # to store average thresholds:
  ret$averages <- matrix(ncol = 3, nrow = length(range_i.seasons))
  colnames(ret$averages) <- c("medium_0.4", "high_0.9", "very_high_0.975")
  rownames(ret$averages) <- paste0("i.seasons_", range_i.seasons)

  # fill in averages thresholds:
  ret$averages[, "medium_0.4"] <- colMeans(results_0.4, na.rm = TRUE)
  ret$averages[, "high_0.9"] <- colMeans(results_0.9, na.rm = TRUE)
  ret$averages[, "very_high_0.975"] <- colMeans(results_0.975, na.rm = TRUE)

  # add quantiles forming uncertainty intervals:
  for(p in interval){
    ret[[paste0("q_", p)]] <- NA*ret$averages
    # compute quantiles of thresholds at different levels:
    ret[[paste0("q_", p)]][, "medium_0.4"] <- apply(results_0.4, MARGIN = 2, FUN = quantile, na.rm = TRUE, p = p)
    ret[[paste0("q_", p)]][, "high_0.9"] <- apply(results_0.9, MARGIN = 2, FUN = quantile, na.rm = TRUE, p = p)
    ret[[paste0("q_", p)]][, "very_high_0.975"] <- apply(results_0.975, MARGIN = 2, FUN = quantile, na.rm = TRUE, p = p)
  }
  # also store range_i.seasons with results
  ret$range_i.seasons <- range_i.seasons
  # return
  return(ret)
}

# function to summarize sensitivity and specificity
# results_0.4: simulation results for medium thereshold
# results_0.9: simulation results for high thereshold
# results_0.975: simulation results for very high thereshold
# peaks_test: the actual peak values
# range_i.seasons: vector containing the number of seasons included (refers to different entries in results_0.4 etc)
summarize_sens_spec <- function(results_0.4, results_0.9, results_0.975, peaks_test, range_i.seasons){
  # function for sensitivity
  sens <- function(threshold, peaks_test, level){
    true_quantile <- quantile(peaks_test, probs = level)
    peaks_test_above <- peaks_test[peaks_test >= true_quantile]
    sens <- mean(peaks_test_above >= threshold)
    return(sens)
  }
  # function for specificity
  spec <- function(threshold, peaks_test, level){
    true_quantile <- quantile(peaks_test, probs = level)
    peaks_test_below <- peaks_test[peaks_test < true_quantile]
    spec <- mean(peaks_test_below < threshold)
    return(spec)
  }
  
  # set up data.frames to store results:
  result_sens <- result_spec <- data.frame(row.names = paste0("i.seasons_", range_i.seasons))
  
  # run through number of seasons used:
  for(i in seq_along(range_i.seasons)){
    # compute and store average sensitivities:
    result_sens$medium_0.4[i] <- mean(sapply(results_0.4[, i], sens, peaks_test = peaks_test, level = 0.4))
    result_sens$high_0.9[i] <- mean(sapply(results_0.9[, i], sens, peaks_test = peaks_test, level = 0.9))
    result_sens$very_high_0.975[i] <- mean(sapply(results_0.975[, i], sens, peaks_test = peaks_test, level = 0.975))
    
    # compute and store average specificities:
    result_spec$medium_0.4[i] <- mean(sapply(results_0.4[, i], spec, peaks_test = peaks_test, level = 0.4))
    result_spec$high_0.9[i] <- mean(sapply(results_0.9[, i], spec, peaks_test = peaks_test, level = 0.9))
    result_spec$very_high_0.975[i] <- mean(sapply(results_0.975[, i], spec, peaks_test = peaks_test, level = 0.975))
  }
  
  # return:
  return(list(sens = result_sens,
              spec = result_spec,
              range_i.seasons = range_i.seasons))
}

# function to summarize positive and negative predictive values
# results_0.4: simulation results for medium thereshold
# results_0.9: simulation results for high thereshold
# results_0.975: simulation results for very high thereshold
# peaks_test: the actual peak values
# range_i.seasons: vector containing the number of seasons included (refers to different entries in results_0.4 etc)
summarize_ppv_npv <- function(results_0.4, results_0.9, results_0.975, peaks_test, range_i.seasons){
  # some helper functions to summarize positives and negatives:
  total_positives  <- function(threshold, peaks_test){
    sum(peaks_test >= threshold)
  }
  
  true_positives <-  function(threshold, peaks_test, level){
    true_quantile <- quantile(peaks_test, probs = level)
    sum(peaks_test >= threshold & peaks_test >= true_quantile)
  }
  
  total_negatives  <- function(threshold, peaks_test){
    sum(peaks_test < threshold)
  }
  
  true_negatives <-  function(threshold, peaks_test, level){
    true_quantile <- quantile(peaks_test, probs = level)
    sum(peaks_test < threshold & peaks_test < true_quantile)
  }
  
  # data.frame to store results:
  result_true_positives <- result_total_positives <-
    result_true_negatives <- result_total_negatives <-
    data.frame(row.names = paste0("i.seasons_", range_i.seasons))
  
  # run through number of seasons used, compute true and false positives and negatives:
  for(i in seq_along(range_i.seasons)){
    result_true_positives$medium_0.4[i] <- sum(sapply(results_0.4[, i], true_positives, peaks_test = peaks_test, level = 0.4))
    result_true_positives$high_0.9[i] <- sum(sapply(results_0.9[, i], true_positives, peaks_test = peaks_test, level = 0.9))
    result_true_positives$very_high_0.975[i] <- sum(sapply(results_0.975[, i], true_positives, peaks_test = peaks_test, level = 0.975))
    
    result_total_positives$medium_0.4[i] <- sum(sapply(results_0.4[, i], total_positives, peaks_test = peaks_test))
    result_total_positives$high_0.9[i] <- sum(sapply(results_0.9[, i], total_positives, peaks_test = peaks_test))
    result_total_positives$very_high_0.975[i] <- sum(sapply(results_0.975[, i], total_positives, peaks_test = peaks_test))
    
    result_true_negatives$medium_0.4[i] <- sum(sapply(results_0.4[, i], true_negatives, peaks_test = peaks_test, level = 0.4))
    result_true_negatives$high_0.9[i] <- sum(sapply(results_0.9[, i], true_negatives, peaks_test = peaks_test, level = 0.9))
    result_true_negatives$very_high_0.975[i] <- sum(sapply(results_0.975[, i], true_negatives, peaks_test = peaks_test, level = 0.975))
    
    result_total_negatives$medium_0.4[i] <- sum(sapply(results_0.4[, i], total_negatives, peaks_test = peaks_test))
    result_total_negatives$high_0.9[i] <- sum(sapply(results_0.9[, i], total_negatives, peaks_test = peaks_test))
    result_total_negatives$very_high_0.975[i] <- sum(sapply(results_0.975[, i], total_negatives, peaks_test = peaks_test))
  }
  
  # compute PPV and NPV from this:
  result_ppv <- result_true_positives / result_total_positives
  result_npv <- result_true_negatives / result_total_negatives
  
  return(list(ppv = result_ppv,
              npv = result_npv,
              range_i.seasons = range_i.seasons))
}

# modify the alpha value of a given color (to generate transparent versions for prediction bands)
modify_alpha <- function(col, alpha){
  x <- col2rgb(col)/255
  rgb(x[1], x[2], x[3], alpha = alpha)
}

# plotting function for thresholds
plot_sim_summary <- function(sim_summary, variable_to_plot = "averages", ylim = NULL, show_bands = TRUE, hlines = NULL, ...){
  # get indices of lower and upper interval ends:
  ind_lower <- names(sim_summary)[grepl("q_", names(sim_summary))][1]
  ind_upper <- names(sim_summary)[grepl("q_", names(sim_summary))][2]

  # determine ylim if not provided
  if(is.null(ylim)){
    yl <- c(0, 1.5*max(sim_summary[[variable_to_plot]][, "very_high_0.975"], na.rm = TRUE))
  }else{
    yl <- ylim
  }

  # plot for medium threshold::
  plot(sim_summary$range_i.seasons,
       sim_summary[[variable_to_plot]][, "medium_0.4"], ylim = yl,
       col = col_medium, pch = 15, ...)
  abline(h = hlines, col = "lightgrey")
  # add intervals if desired:
  if(show_bands){
    polygon(c(sim_summary$range_i.seasons, rev(sim_summary$range_i.seasons)),
            c(sim_summary[[ind_lower]][, "medium_0.4"], rev(sim_summary[[ind_upper]][, "medium_0.4"])), border = NA,
            col = modify_alpha(col_medium, 0.2))
  }

  # add high thresholds:
  points(sim_summary$range_i.seasons, sim_summary[[variable_to_plot]][, "high_0.9"], col = col_high, pch = 15)
  # intervals if desired
  if(show_bands){
    polygon(c(sim_summary$range_i.seasons, rev(sim_summary$range_i.seasons)),
            c(sim_summary[[ind_lower]][, "high_0.9"], rev(sim_summary[[ind_upper]][, "high_0.9"])), border = NA,
            col = modify_alpha(col_high, 0.2))
  }

  # very high thresholds
  points(sim_summary$range_i.seasons, sim_summary[[variable_to_plot]][, "very_high_0.975"], col = col_very_high, pch = 15)
  # intervals if desired
  if(show_bands){
    polygon(c(sim_summary$range_i.seasons, rev(sim_summary$range_i.seasons)),
            c(sim_summary[[ind_lower]][, "very_high_0.975"], rev(sim_summary[[ind_upper]][, "very_high_0.975"])), border = NA,
            col = modify_alpha(col_very_high, 0.2))
  }
}

# plotting function for exceedance proportions
plot_exceedance_summary <- function(summary_exceedance, ...){
  # extract numbers of seasons used:
  range_i.seasons <- as.numeric(gsub("i.seasons_", "", rownames(summary_exceedance$averages), 1))
  # helper vector for plotting:
  x <- rep(range_i.seasons, each = 2) + c(-0.5, 0.5)
  
  # initialize plot:
  plot(NULL, xlim = range(range_i.seasons) + c(-0.5, 0.5), ylim = 0:1, ...)
  # very high:
  polygon(c(x, rev(x)),
          rep(0:1, each = length(x)), col = col_very_high, border = NA)
  # high:
  y_orange <- rep(1 - c(summary_exceedance$averages[, "very_high_0.975"], rep(1, nrow(summary_exceedance$averages))), each = 2)
  polygon(c(x, rev(x)),
          y_orange,
          col = col_high, border = NA)
  # medium:
  y_yellow <- rep(1 - c(summary_exceedance$averages[, "high_0.9"], rep(1, nrow(summary_exceedance$averages))), each = 2)
  polygon(c(x, rev(x)),
          y_yellow,
          col = col_medium, border = NA)
  # low:
  y_green <- rep(1 - c(summary_exceedance$averages[, "medium_0.4"], rep(1, nrow(summary_exceedance$averages))), each = 2)
  polygon(c(x, rev(x)),
          y_green,
          col = col_low, border = NA)
  
  abline(v = 0.5 + (0:16), col = "white", lwd = 4)
  abline(h = c(0.4, 0.9, 0.975), lty = 3)
  box()
}

# custon legend for plots
legend_summary <- function(){
  plot(NULL, xlim = 0:1, ylim = 0:1, axes = FALSE, xlab = "", ylab = "")
  legend("top", legend = c("Intensity levels:", "very high", "high", "medium", "low", 
                           "First column:", "mean threshold", "empirical 5% and 95% quantiles", "approximation of mean", "threshold under normality",
                           "Second column:",  "intended shares of intensity levels", "empirical shares", "", "",
                           "Third to fifth columns:", "empirical sensitivity /", "specificity / PPV", "approximation under normality", "(where applicable)"), 
         col = c(NA, col_very_high, col_high, col_medium, col_low,
                 NA, "black", "grey", "black", NA,
                 NA, "black", col_very_high, NA, NA,
                 NA, "black", NA, "black", NA),
         pch = c(NA, 15, 15, 15, 15, 
                 NA, 15, 15, NA, NA,
                 NA, NA, 15, NA, NA,
                 NA, 15, NA, NA, NA), 
         lty = c(NA, NA, NA, NA, NA,
                 NA, NA, NA, 1, NA,
                 NA, 3, NA, NA, NA,
                 NA, NA, NA, 1, NA),
         pt.cex = c(1, 1, 1, 1, 1,
                 1, 1, 2, 0.6, 1,
                 1, 1, 2, 1, 1,
                 1, 1, 1, 1, 1),
         ncol = 4, bty = "n")
}


# function to add analytical approximations as lines
lines_approx_expectations <- function(approx_expectations, lty = 1, pch = NA, cex = 0.6, black = FALSE){
  lines(approx_expectations[, "m"], approx_expectations[, "medium_0.4"], type = "l", lty = lty, pch = pch, cex = cex,
        col = ifelse(black, "black", col_medium))
  lines(approx_expectations[, "m"], approx_expectations[, "high_0.9"], type = "l", lty = lty, pch = pch, cex = cex,
        col = ifelse(black, "black", col_high))
  lines(approx_expectations[, "m"], approx_expectations[, "very_high_0.975"], type = "l", lty = lty, pch = pch, cex = cex,
        col = ifelse(black, "black", col_very_high))
}

# wrapper function to make five plots:
five_plots <- function(thresholds_0.4, thresholds_0.9, thresholds_0.975,
                       exceedance_0.4, exceedance_0.9, exceedance_0.975,
                       peaks_test, range_i.seasons, ylim_threshold,
                       approx_expectations, 
                       approx_sens = NULL, approx_spec = NULL, approx_ppv = NULL,
                       main = "", line = lin, cex = ce,
                       main2 = "", line2 = lin + 1.5, cex2 = ce){
  
  # compute summaries:
  summary_thresholds <- sim_summary(thresholds_0.4, thresholds_0.9, thresholds_0.975,
                                    range_i.seasons, interval = interv)
  summary_exceedance <- sim_summary(exceedance_0.4, exceedance_0.9, exceedance_0.975,
                                    range_i.seasons, interval = interv)
  summary_sens_spec <- summarize_sens_spec(thresholds_0.4, thresholds_0.9, thresholds_0.975,
                                           peaks_test, range_i.seasons)
  summary_ppv_npv <- summarize_ppv_npv(thresholds_0.4, thresholds_0.9, thresholds_0.975,
                                       peaks_test, range_i.seasons)
  
  # make plots:
  
  # summary of thresolds:
  plot_sim_summary(summary_thresholds, xlab = "# included seasons m", ylab = "threshold", ylim = ylim_threshold)
  if(!is.null(approx_expectations)){
    lines_approx_expectations(approx_expectations)
  }
  
  # summary of exceedance probabilities:
  plot_exceedance_summary(summary_exceedance, xlab = "# included seasons m", ylab = "average share")
  
  # summary of sensitivity:
  plot_sim_summary(summary_sens_spec, xlab = "# included seasons", ylab = "sensitivity",
                   ylim = 0:1, show_bands = FALSE, variable_to_plot = "sens", hlines = 0:5/5)
  if(!is.null(approx_sens)){
    lines_approx_expectations(approx_sens)
  }
  mtext(side = 3, text = main, cex = cex, line = line)
  mtext(side = 3, text = main2, cex = cex2, line = line2, font = 2)
  
  # summary of specificity:
  plot_sim_summary(summary_sens_spec, xlab = "# included seasons", ylab = "specificity",
                   ylim = 0:1, show_bands = FALSE, variable_to_plot = "spec", hlines = 0:5/5)
  if(!is.null(approx_spec)){
    lines_approx_expectations(approx_spec)
  }
  
  # summary of PPVs:
  plot_sim_summary(summary_ppv_npv, xlab = "# included seasons", ylab = "PPV",
                   ylim = 0:1, show_bands = FALSE, variable_to_plot = "ppv", hlines = 0:5/5)
  if(!is.null(approx_ppv)){
    lines_approx_expectations(approx_ppv)
  }
}


# function to get (mis-)classificatin matrix of the four categories
# results_0.4, results_0.9, results_0.975: the thresholds
# peaks_test: the true peak values
# i.seasons: the number of historical seasons (i.e., column in results_0.4 etc) to use
classification_matrix <- function(results_0.4, results_0.9, results_0.975, peaks_test, i.seasons){
  # true classification of peaks:
  true_quantiles <- quantile(peaks_test, probs = c(0.4, 0.9, 0.975))
  list_peaks_test <- list()
  list_peaks_test$low <- peaks_test[peaks_test < true_quantiles[1]]
  list_peaks_test$medium <- peaks_test[peaks_test >= true_quantiles[1] & peaks_test < true_quantiles[2]]
  list_peaks_test$high <- peaks_test[peaks_test >= true_quantiles[2] & peaks_test < true_quantiles[3]]
  list_peaks_test$very_high <- peaks_test[peaks_test >= true_quantiles[3]]
  
  # label to handle indices later
  label_i.seasons <- paste0("i.seasons_", i.seasons)
  
  # matrices to store stuff:
  classifications <- matrix(ncol = 4, nrow = 4)
  rownames(classifications) <- colnames(classifications) <- c("low", "medium", "high", "very_high")
  exceedance <- classifications
  
  # function to compare a vector of peak values and a single threshold (to be put into sapply)
  at_least <- function(threshold, peaks_test) sum(peaks_test >= threshold)
  
  # compute exceedance probabilities for different true classes:
  for(true_level in names(list_peaks_test)){
    # low (handle using multiplication by 0):
    exceedance[true_level, "low"] <-  sum(sapply(0*results_0.4[, label_i.seasons], at_least, peaks_test = list_peaks_test[[true_level]]))
    # medium:
    exceedance[true_level, "medium"] <- sum(sapply(results_0.4[, label_i.seasons], at_least, peaks_test = list_peaks_test[[true_level]]))
    # high
    exceedance[true_level, "high"] <- sum(sapply(results_0.9[, label_i.seasons], at_least, peaks_test = list_peaks_test[[true_level]]))
    # very high
    exceedance[true_level, "very_high"] <- sum(sapply(results_0.975[, label_i.seasons], at_least, peaks_test = list_peaks_test[[true_level]]))
  }
  
  # get classification probabilities via differences
  classifications[, "very_high"] <- exceedance[, "very_high"]
  classifications[, "high"] <- (exceedance[, "high"] - exceedance[, "very_high"])
  classifications[, "medium"] <- (exceedance[, "medium"] - exceedance[, "high"])
  classifications[, "low"] <- (exceedance[, "low"] - exceedance[, "medium"])
  
  # re-scale
  classifications <- classifications/sum(classifications)
  
  # return
  return(classifications)
}


# function to generate mosaic plot of classifications:
mosaic <- function(classification_matrix){
  # initialize plot:
  plot(NULL, xlim = 0:1, ylim = 0:1, xlab = "fractions of assigned intensity levels", ylab = "", axes = FALSE)
  # put axis on the right:
  axis(1)
  axis(4)
  box()

  # figure out coordinates and colors:
  y_coords <- c(0, cumsum(rowSums(classification_matrix)))
  cols <- c(col_low, col_medium, col_high, col_very_high, las = 1)
  
  # run through rows and columns of classification matrix and add rectangles:
  for(row in 1:4){
    for(col in 1:4){
      x_coords_temp <- c(0, cumsum(classification_matrix[row, ]))/sum(classification_matrix[row, ])
      rect(x_coords_temp[col], y_coords[row], x_coords_temp[col + 1], y_coords[row + 1], col = cols[col], border = "white")
    }
  }
}



# classification_matrix <- function(results_0.4, results_0.9, results_0.975, peaks_test, i.seasons){
#   # true classification of peaks:
#   peaks_test <- sort(peaks_test)
#   n_peaks_test <- length(peaks_test)
#   
#   # label to handle indices later
#   label_i.seasons <- paste0("i.seasons_", i.seasons)
#   
#   # matrices to store stuff:
#   classifications <- matrix(ncol = 4, nrow = n_peaks_test)
#   colnames(classifications) <- c("low", "medium", "high", "very_high")
#   exceedance <- classifications
#   
#   # function to compare a vector of peak values and a single threshold (to be put into sapply)
#   at_least <- function(new, thresholds) mean(new >= thresholds)
#   
#   # compute exceedance probabilities for different true classes:
#   for (i in 1:n_peaks_test) {
#       # low (handle using multiplication by 0):
#       exceedance[i, "low"] <-  1
#       # medium:
#       exceedance[i, "medium"] <-  at_least(new = peaks_test[i], thresholds = results_0.4[, label_i.seasons])
#       # high
#       exceedance[i, "high"] <-  at_least(new = peaks_test[i], thresholds = results_0.9[, label_i.seasons])
#       # very high
#       exceedance[i, "very_high"] <-  at_least(new = peaks_test[i], thresholds = results_0.975[, label_i.seasons])
#   }
#   
#   # get classification probabilities via differences
#   classifications[, "very_high"] <- exceedance[, "very_high"]
#   classifications[, "high"] <- (exceedance[, "high"] - exceedance[, "very_high"])
#   classifications[, "medium"] <- (exceedance[, "medium"] - exceedance[, "high"])
#   classifications[, "low"] <- (exceedance[, "low"] - exceedance[, "medium"])
#   
#   # return
#   return(classifications)
# }
# 
# 
# mosaic <- function(classification_matrix){
#   # initialize plot:
#   par(las = 1)
#   plot(NULL, xlim = 0:1, ylim = 0:1, xlab = "fractions of assigned intensity levels", ylab = "", axes = FALSE)
#   # put axis on the right:
#   axis(1)
#   axis(4)
#   par(las = 0)
#   mtext("true quantile level of peak", side = 4, cex = 0.65, line = 2.5)
#   box()
#   
#   # figure out coordinates and colors:
#   y_coords <- rep(seq(from = 0, to = 1, length.out = nrow(classification_matrix) + 1), each = 2)
#   x_coords <- matrix(NA, ncol = 4, nrow = 2*nrow(classification_matrix) + 2)
#   colnames(x_coords) <- c("low", "medium", "high", "very_high")
#   x_coords[, "very_high"] <- c(0, rep(1, 2*nrow(classification_matrix)), 0)
#   x_coords[, "high"] <- c(0, rep(1 - classification_matrix[, "very_high"], each = 2), 0)
#   x_coords[, "medium"] <- c(0, rep(classification_matrix[, "low"] + classification_matrix[, "medium"], each = 2), 0)
#   x_coords[, "low"] <- c(0, rep(classification_matrix[, "low"], each = 2), 0)
#   
#   cols <- c(col_low, col_medium, col_high, col_very_high, las = 1)
#   
#   # run through rows and columns of classification matrix and add rectangles:
#   polygon(x = x_coords[, "very_high"], y = y_coords, col = col_very_high, border = NA)
#   polygon(x = x_coords[, "high"], y = y_coords, col = col_high, border = NA)
#   polygon(x = x_coords[, "medium"], y = y_coords, col = col_medium, border = NA)
#   polygon(x = x_coords[, "low"], y = y_coords, col = col_low, border = NA)
#   
#   abline(h = c(0.4, 0.9, 0.975), col = "black", lty = 3)
# }
