memmodel_lm <- function(i.data, 
                        i.seasons = ncol(i.data),
                        i.n.max = 1,
                        include_trend = FALSE,
                        include_rank = FALSE,
                        i.level.intensity = c(0.4, 0.9, 0.975),
                        i.trafo = "log"){
  
  # check arguments:
  if(!i.trafo %in% c("identity", "log")) stop("i.trafo needs to be 'identity' or 'log'.")
  
  # restrict data:
  i.data <- i.data[, seq(to = ncol(i.data), length.out = i.seasons, by = 1)]
  
  # sort within each season:
  for(i in 1:ncol(i.data)){
    i.data[, i] <- sort(i.data[, i], decreasing = TRUE)
  }
  
  # restrict to last i.n.max highest values:
  i.data <- i.data[1:i.n.max, ]
  i.rank <- matrix(1:i.n.max, nrow = i.n.max, ncol = i.seasons)
  i.t <- matrix(1:i.seasons, nrow = i.n.max, ncol = i.seasons, byrow = TRUE)
  
  # pool into vector:
  x <- as.vector(unlist(i.data))
  rank <- as.factor(as.vector(i.rank))
  t <- as.vector(i.t)
  
  # transform:
  if(i.trafo == "log"){
    y <- log(x)
  }
  if(i.trafo == "identity"){
    y <- x
  }
  n <- length(y)
  
  # compute thresholds on tansformed scale:
  
  # choose formula:
  formula <- y ~ 1
  if(include_trend){
    formula <- y ~ 1 + t
  }
  if(include_rank){
    formula <- y ~ 1 + rank
  }
  if(include_rank & include_trend){
    formula <- y ~ 1 + t + rank
  }
  
  # fit linear model:
  m <- lm(formula = formula)
  # obtain predictive quantiles:
  newdata <- data.frame(t = i.seasons + 1, rank = rank[1])
  q_y <- numeric(0)
  for (i in seq_along(i.level.intensity)){
    pred <- add_quantile(newdata, fit = m, p = i.level.intensity[i])
    q_y[i] <- pred[1, ncol(pred)]
  }
  
  # transform back:
  if(i.trafo == "log"){
    q_x <- exp(q_y)
  }
  if(i.trafo == "identity"){
    q_x <- q_y
  }
  
  result <- list(intensity.thresholds = q_x, fit = m)
}
