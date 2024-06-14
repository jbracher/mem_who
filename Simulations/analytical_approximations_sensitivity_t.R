# Approximating the sensitivity, specificity, PPV and NPV analytically

# set working directory
current_path = rstudioapi::getActiveDocumentContext()$path # get path of this file
setwd(dirname(current_path))

# analytical approximation of sensitivity:
# n: number of available observations
# alpha: the nominal quantile level at the threshold
sens_t <- function(n, alpha){
  # approximate expectation and variance of threshold
  mean_q <- qt(alpha, df = n - 1)*sqrt(1 + 1/n)
  var_q <- 1/n + qt(alpha, df = n - 1)^2*(1 + 1/n)/(2*(n - 1))
  
  # function to integrate
  f <- function(x){
    dnorm(x)*(pnorm(x, mean = mean_q, sd = sqrt(var_q)))
  }
  
  # integrate:
  integrate(f, lower = qnorm(alpha), upper = 10)$value/(1 - alpha)
}

# simulation-based computation of sensitivity:
# n: number of available observations
# alpha: the nominal quantile level at the threshold
# n_sim: number of simulation samples
sens_t_sim <- function(n, alpha, n_sim = 50000){
  # generate samples
  Y0 <- rnorm(n = n_sim*n)
  # arrange into sequences of length n
  Y <- matrix(Y0, ncol = n, nrow = n_sim)
  # sample new extreme observations, which are above the
  # alpha quantile of the true distribution:
  Y_new_high <- qnorm(runif(n_sim, alpha, 1))
  
  # compute alpha thresholds:
  means <- rowMeans(Y)
  sds <- apply(Y, MARGIN = 1, FUN = sd)
  vars <- apply(Y, MARGIN = 1, FUN = var)
  q <- means + qt(alpha, df = n - 1)*sqrt(1 + 1/n)*sds
  
  # compute sensitivity
  mean(Y_new_high > q)
}

# analytical approximation of specificity:
# n: number of available observations
# alpha: the nominal quantile level at the threshold
spec_t <- function(n, alpha){
  # approximate expectation and variance of threshold
  mean_q <- qt(alpha, df = n - 1)*sqrt(1 + 1/n)
  var_q <- 1/n + qt(alpha, df = n - 1)^2*(1 + 1/n)/(2*(n - 1))
  
  # function to integrate
  f <- function(x){
    dnorm(x)*(1 - pnorm(x, mean = mean_q, sd = sqrt(var_q)))
  }
  
  # integrate:
  integrate(f, lower = -10, upper = qnorm(alpha))$value/(alpha)
}

# simulation-based computation of specificity:
# n: number of available observations
# alpha: the nominal quantile level at the threshold
# n_sim: number of simulation samples
spec_t_sim <- function(n, alpha, n_sim = 50000){
  # generate samples
  Y0 <- rnorm(n = n_sim*n)
  # arrange into sequences of length n
  Y <- matrix(Y0, ncol = n, nrow = n_sim)
  # sample new extreme observations, which are above the
  # alpha quantile of the true distribution:
  Y_new_low <- qnorm(runif(n_sim, 0, alpha))
  
  # compute alpha thresholds:
  means <- rowMeans(Y)
  sds <- apply(Y, MARGIN = 1, FUN = sd)
  vars <- apply(Y, MARGIN = 1, FUN = var)
  q <- means + qt(alpha, df = n - 1)*sqrt(1 + 1/n)*sds
  
  # compute specificity:
  mean(Y_new_low < q)
}

# apply and plot:
values_m <- 5:30

# load("Results/analytical_sens_spec_t.rda")

# medium threshold:
values_sens_t_0.4 <- sapply(values_m, sens_t, alpha = 0.4)
values_sens_t_sim_0.4 <- sapply(values_m, sens_t_sim, alpha = 0.4)
values_spec_t_0.4 <- sapply(values_m, spec_t, alpha = 0.4)
values_spec_t_sim_0.4 <- sapply(values_m, spec_t_sim, alpha = 0.4)
# compute PPVs via usual formula:
values_ppv_t_0.4 <- values_sens_t_0.4*0.6/(values_sens_t_0.4*0.6 + (1 - values_spec_t_0.4)*0.4)
values_ppv_t_sim_0.4 <- values_sens_t_sim_0.4*0.6/(values_sens_t_sim_0.4*0.6 + (1 - values_spec_t_sim_0.4)*0.4)

# high threshold:
values_sens_t_0.9 <- sapply(values_m, sens_t, alpha = 0.9)
values_sens_t_sim_0.9 <- sapply(values_m, sens_t_sim, alpha = 0.9)
values_spec_t_0.9 <- sapply(values_m, spec_t, alpha = 0.9)
values_spec_t_sim_0.9 <- sapply(values_m, spec_t_sim, alpha = 0.9)
# compute PPVs via usual formula:
values_ppv_t_0.9 <- values_sens_t_0.9*0.1/(values_sens_t_0.9*0.1 + (1 - values_spec_t_0.9)*0.9)
values_ppv_t_sim_0.9 <- values_sens_t_sim_0.9*0.1/(values_sens_t_sim_0.9*0.1 + (1 - values_spec_t_sim_0.9)*0.9)

# very high threshold:
values_sens_t_0.975 <- sapply(values_m, sens_t, alpha = 0.975)
values_sens_t_sim_0.975 <- sapply(values_m, sens_t_sim, alpha = 0.975)
values_spec_t_0.975 <- sapply(values_m, spec_t, alpha = 0.975)
values_spec_t_sim_0.975 <- sapply(values_m, spec_t_sim, alpha = 0.975)
# compute PPVs via usual formula:
values_ppv_t_0.975 <- values_sens_t_0.975*0.025/(values_sens_t_0.975*0.025 + (1 - values_spec_t_0.975)*0.975)
values_ppv_t_sim_0.975 <- values_sens_t_sim_0.975*0.025/(values_sens_t_sim_0.975*0.025 + (1 - values_spec_t_sim_0.975)*0.975)


# summarize results in data.frames:
approx_sens_t <- data.frame(m = values_m,
                          medium_0.4 = values_sens_t_0.4,
                          high_0.9 = values_sens_t_0.9, 
                          very_high_0.975 = values_sens_t_0.975)

approx_spec_t <- data.frame(m = values_m,
                          medium_0.4 = values_spec_t_0.4,
                          high_0.9 = values_spec_t_0.9, 
                          very_high_0.975 = values_spec_t_0.975)

approx_ppv_t <- data.frame(m = values_m,
                          medium_0.4 = values_ppv_t_0.4,
                          high_0.9 = values_ppv_t_0.9, 
                          very_high_0.975 = values_ppv_t_0.975)

# write out:
save(approx_sens_t, approx_spec_t, approx_ppv_t,
     values_sens_t_0.4, values_sens_t_0.9, values_sens_t_0.975,
     values_spec_t_0.4, values_spec_t_0.9, values_spec_t_0.975,
     values_ppv_t_0.4, values_ppv_t_0.9, values_ppv_t_0.975,
     values_sens_t_sim_0.4, values_sens_t_sim_0.9, values_sens_t_sim_0.975,
     values_spec_t_sim_0.4, values_spec_t_sim_0.9, values_spec_t_sim_0.975,
     values_ppv_t_sim_0.4, values_ppv_t_sim_0.9, values_ppv_t_sim_0.975,
     file = "Results/analytical_sens_spec_t.rda")

# load("Results/analytical_sens_spec_t.rda")

# Plot:

pdf("../Draft/figure/analytical_sens_spec_t.pdf", width = 8, height = 3)
par(mfrow = c(1, 3))
# sensitivity:
plot(values_m, values_sens_t_0.975, type = "l", ylim = c(0, 1), col = "darkorchid4", lty = 2,
     xlab = "n", ylab = "Sensitivity")
abline(h = 0:5/5, col = "lightgrey")

lines(values_m, values_sens_t_sim_0.975, col = "darkorchid4")

lines(values_m, values_sens_t_0.9, col = "red", lty = 2)
lines(values_m, values_sens_t_sim_0.9, col = "red")

lines(values_m, values_sens_t_0.4, col = "orange", lty = 2)
lines(values_m, values_sens_t_sim_0.4, col = "orange")

legend("bottom", col = "black", lty = c(2, 1),
       legend = c("theoretical approximation", "simulation"),
       bty = "n")

# specificity:
plot(values_m, values_spec_t_0.975, type = "l", ylim = c(0, 1), col = "darkorchid4", lty = 2,
     xlab = "m", ylab = "Specificity")
abline(h = 0:5/5, col = "lightgrey")

lines(values_m, values_spec_t_sim_0.975, col = "darkorchid4")

lines(values_m, values_spec_t_0.9, col = "red", lty = 2)
lines(values_m, values_spec_t_sim_0.9, col = "red")

lines(values_m, values_spec_t_0.4, col = "orange", lty = 2)
lines(values_m, values_spec_t_sim_0.4, col = "orange")

mtext("Using quantiles of the t-distribution", line = 1.5)

# PPV:
plot(values_m, values_ppv_t_0.975, type = "l", ylim = c(0, 1), col = "darkorchid4", lty = 2,
     ylab = "PPV", xlab = "m")
abline(h = 0:5/5, col = "lightgrey")
lines(values_m, values_ppv_t_sim_0.975, col = "darkorchid4")


lines(values_m, values_ppv_t_0.9, col = "red", lty = 2)
lines(values_m, values_ppv_t_sim_0.9, col = "red")

lines(values_m, values_ppv_t_0.4, col = "orange", lty = 2)
lines(values_m, values_ppv_t_sim_0.4, col = "orange")

legend("bottom", col = c("darkorchid4", "red", "orange"),
       legend = c(expression(alpha == 0.975~(very~ high)), 
                  expression(alpha==0.9~(high)),
                  expression(alpha==0.4~(medium))),
       lty = 1, bty = "n")

dev.off()
