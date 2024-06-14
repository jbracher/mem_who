# approximating the exceedance probability when using normal rather than 
# t distributions

n_sim <- 100000

# generate samples
n <- 5
alpha <- 0.9

Y0 <- rnorm(n = n_sim*n)
# arrange into sequences of length n
Y <- matrix(Y0, ncol = n, nrow = n_sim)

# compute alpha thresholds:
means <- rowMeans(Y)
sds <- apply(Y, MARGIN = 1, FUN = sd)

q_t <- means + qt(alpha, df = n - 1)*sqrt(1 + 1/n)*sds
mean(1 - pnorm(q_t))

q_norm <- means + qnorm(alpha)*sds
mean(1 - pnorm(q_norm))

1 - pt(qnorm(alpha)/sqrt(1 + 1/n), df = n - 1)

values_n <- 5:30

exceedance_0.4 <- 1 - pt(qnorm(0.4)/sqrt(1 + 1/values_n), df = values_n - 1)
exceedance_0.9 <- 1 - pt(qnorm(0.9)/sqrt(1 + 1/values_n), df = values_n - 1)
exceedance_0.975 <- 1 - pt(qnorm(0.975)/sqrt(1 + 1/values_n), df = values_n - 1)

pdf("../Draft/figure/exceedance_prob_normal.pdf", width = 6, height = 3.5)
par(las = 1, mar = c(4, 4, 1, 1))
plot(NULL, NULL, type = "l", xlim = c(5, 30), ylim = c(0, 0.6),
     xlab = expression(n %*% m), ylab = "exceedance probability")
abline(h = c(0.025, 0.1, 0.6), lty = "dotted")
lines(values_n, exceedance_0.4, col = "orange")
lines(values_n, exceedance_0.9, col = "red")
lines(values_n, exceedance_0.975, col = "darkorchid4")
legend("center", legend = c("Threshold levels:", "very high", "high", "medium",
                         "", "intended excedance", "probabilities", ""), ncol = 2,
       col = c(NA, "darkorchid4", "red", "orange", NA, "black", NA, NA),
       lty = c(NA, 1, 1, 1, NA, 3, NA, NA), bty = "n", cex = 0.75)
dev.off()