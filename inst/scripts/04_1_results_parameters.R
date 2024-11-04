###################################################
### Section 4 - RESULTS FULL MODEL (PARAMETERS) ###
###################################################

# input  data: "data/model.rds"

library("ggplot2")
library("jtools")

# FUNCTIONS
# This function plots the coefficients
plot_coefs <- function (params, 
                        ci_level = 0.9, 
                        intercept = 0,
                        colors = "black", 
                        point.size = 5, 
                        legend.title = "Model") {
  
  coef.names <- colnames(params)
  
  tidies <- data.frame(estimate = colMeans(params))
  tidies$conf.low  <- apply(params, 2, quantile, prob = (1 - ci_level) / 2)
  tidies$conf.high <- apply(params, 2, quantile, prob = (1 + ci_level) / 2)
  tidies$model     <- as.factor("M1")
  tidies$term      <- factor(coef.names, levels = rev(coef.names))
  
  dh <- as.numeric(!TRUE) * 0.5
  
  p <- ggplot(data = tidies, 
              aes(y = term, x = estimate, xmin = conf.low, xmax = conf.high, colour = model))
  
  p <- p + ggplot2::geom_pointrange(
    aes(y = term, x = estimate, xmin = conf.low, xmax = conf.high, colour = model, shape = model), 
    position = ggplot2::position_dodge(width = dh), 
    fill = "white", fatten = 1)
  
  p <- p + geom_vline(xintercept = intercept, linetype = 2) + 
    scale_colour_manual(values = colors, limits = rev(levels(tidies$model)), 
                        breaks = rev(levels(tidies$model)), 
                        labels = rev(levels(tidies$model))) + 
    drop_y_gridlines() + 
    theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 10)) + 
    xlab(ifelse(FALSE, no = "Estimate", yes = "exp(Estimate)"))
  
  p <- p + 
    scale_y_discrete(limits = levels(tidies$term), name = "Parameters") +
    scale_x_continuous(breaks=seq(-2,2,.2))
  
  yrange <- ggplot_build(p)$layout$panel_params[[1]]$y.range
  xrange <- ggplot_build(p)$layout$panel_params[[1]]$x.range
  
  if (yrange[2] <= (length(unique(tidies$term)) + 0.8)) {
    upper_y <- length(unique(tidies$term)) + 0.8
  } else {
    upper_y <- yrange[2]
  }
  
  lower_y <- 0.8
  p <- p + coord_cartesian(ylim = c(lower_y, upper_y), xlim = xrange, expand = FALSE) +
    theme_minimal() + theme(legend.position = "none")  
  
  return(p)
}

# This function returns the mean and quantiles by columns in X
my.summary <- function(X, prob = c(0.05, 0.95), round = 4) {
  
  if (is.matrix(X)) {
    x <- colMeans(X)
    y <- apply(X, 2, quantile, prob = prob)
    return(t(round(rbind("mean" = x, y), round)))
  } else {
    x <- mean(X)
    y <- quantile(X, prob = prob)
    return(round(c("mean" = x, y), round))
  }
}

# This function edits the output from my.summary to latex notation for a table
in.tex <- function(x) {
  internt.fun <- function(x) {
    paste0("\texttt{", x[4], "} & $", format(x[1], nsmall=4), "$ & $(", format(x[2], nsmall=4), ",", format(x[3], nsmall=2), ")$ \\")
  }
  if (is.matrix(x)) {
    x <- data.frame(x, rownames(x))
    rownames(x) <- NULL
    return(apply(x, 1, internt.fun))
  } else {
    x[4] <- "x"
    return(internt.fun(x))
  }
}


# MODEL PARAMETERS
## beta
beta <- rbind(model[,1]$params$beta, model[,2]$params$beta)
my_names <- c(
  "poly(trend, 2)1"          , "poly(trend, 2)2"          ,
  "sine"                     , "cosi"                     ,
  "log(dist)"                , "lag1"                     ,
  "lag2"                     , "lag1:lag2"                ,
  "sine:poly(trend, 2)1"     , "cosi:poly(trend, 2)1"     ,
  "sine:poly(trend, 2)2"     , "cosi:poly(trend, 2)2"     ,
  "log(dist):poly(trend, 2)1", "log(dist):poly(trend, 2)2",
  "trend:lag1"               , "trend:lag2"               ,
  "trend:lag1:lag2"          , "log(dist):lag1"           , 
  "log(dist):lag2"           , "log(dist):lag1:lag2"      )
my_order <- match(my_names, colnames(beta))
colnames(beta)[my_order] <- c(
  "trend1"            , "trend2"             ,
  "sin"               , "cos"                ,
  "log(dist)"         , "lag1"               ,
  "lag2"              , "lag1:lag2"          ,
  "sin:trend1"        , "cos:trend1"         ,
  "sin:trend2"        , "cos:trend2"         ,
  "log(dist):trend1"  , "log(dist):trend2"   ,
  "log(t-1):lag1"     , "log(t-1):lag2"      ,
  "log(t-1):lag1:lag2", "log(dist):lag1"     , 
  "log(dist):lag2"    , "log(dist):lag1:lag2")

## plot coefficients
p <- plot_coefs(beta[, my_order])
ggsave("inst/img/MAIN_coeffs.pdf", p, width = 5.5, height = 3.9)

## table coefficients
beta <- rbind(model[,1]$params$beta, model[,2]$params$beta)
beta <- sweep(beta, 2, attr(model[,1]$x, "scaled:scale")[colnames(beta)], FUN = '/')
coef <- attr(poly(model[,1]$model$trend, 2), "coefs")

plot(poly(model[,1]$model$trend, 2)[1:61,1])
lines((model[,1]$model$trend[1:61] - coef$alpha[1]) / sqrt(coef$norm2)[3], col = "red")

plot(poly(model[,1]$model$trend, 2)[1:61,2])
lines(((model[,1]$model$trend[1:61] - coef$alpha[2]) * (model[,1]$model$trend[1:61] - coef$alpha[1]) - coef$norm2[3] / length(model[,1]$model$trend)) / sqrt(coef$norm2)[4], col = "red")

c0 <- 1 / sqrt(coef$norm2)[3]
c1 <- 1 / sqrt(coef$norm2)[4]
c2 <- -c1 * sum(coef$alpha)

beta[,"poly(trend, 2)1"] <- c0 * beta[,"poly(trend, 2)1"] + c2 * beta[,"poly(trend, 2)2"]
beta[,"poly(trend, 2)2"] <- c1 * beta[,"poly(trend, 2)2"]
beta[,"sine:poly(trend, 2)1"] <- c0 * beta[,"sine:poly(trend, 2)1"] + c2 * beta[,"sine:poly(trend, 2)2"]
beta[,"sine:poly(trend, 2)2"] <- c1 * beta[,"sine:poly(trend, 2)2"]
beta[,"cosi:poly(trend, 2)1"] <- c0 * beta[,"cosi:poly(trend, 2)1"] + c2 * beta[,"cosi:poly(trend, 2)2"]
beta[,"cosi:poly(trend, 2)2"] <- c1 * beta[,"cosi:poly(trend, 2)2"]
beta[,"log(dist):poly(trend, 2)1"] <- c0 * beta[,"log(dist):poly(trend, 2)1"] + c2 * beta[,"log(dist):poly(trend, 2)2"]
beta[,"log(dist):poly(trend, 2)2"] <- c1 * beta[,"log(dist):poly(trend, 2)2"]

in.tex(my.summary(beta[, my_order]))

## hyperparameters
hp <- rbind(model[,1]$params$hp, model[,2]$params$hp)
in.tex(my.summary(hp[,c(1,3,5,8)]))
quantile(3 / hp[,"decay0"], prob = c(0.05, 0.95)) # effective range
in.tex(my.summary(exp(- hp[,"decay0"] * max(dist(coords))), round = 2)) # correlation
in.tex(my.summary(1 / hp[,c(2,4,6,7,9,10)], round = 2)) # sigma2
in.tex(my.summary(1 / sqrt(hp[,c(2,4,6,7,9,10)]))) # sigma

## beta l = 1,2
betal <- cbind(
  rbind(model[,1]$params$betal1, model[,2]$params$betal1),
  rbind(model[,1]$params$betal2, model[,2]$params$betal2))
betal <- sweep(betal, 2, attr(model[,1]$x, "scaled:scale")[colnames(betal)], FUN = '/')
betal[,c(1,3)] <- c0 * betal[,c(1,3)] 

in.tex(my.summary(betal))

# TREND
l <- c(15, 196)
t <- 2:62
logdist <- log(c(stations$dist[2], stations$dist[25]))

beta <- colMeans(beta)

# Madrid / January 15 / I_{t,l-1} = 0, I_{t,l-2} = 0
X <- cbind(sin(2 * pi * l[1] / 365), cos(2 * pi * l[1] / 365), 
           logdist[1], log(t - 1), log(t - 1)^2, 
           0, 0, 
           sin(2 * pi * l[1] / 365) * log(t - 1), sin(2 * pi * l[1] / 365) * log(t - 1)^2, 
           cos(2 * pi * l[1] / 365) * log(t - 1), cos(2 * pi * l[1] / 365) * log(t - 1)^2, 
           logdist[1] * log(t - 1), logdist[1] * log(t - 1)^2, 0, 0, 0, 0, 0, 0, 0)
plot(X %*% beta, ylim = c(-9-7, 2-6.5), type = "l")
# Madrid / July 15 / I_{t,l-1} = 0, I_{t,l-2} = 0
X <- cbind(sin(2 * pi * l[2] / 365), cos(2 * pi * l[2] / 365), 
           logdist[1], log(t - 1), log(t - 1)^2, 
           0, 0, 
           sin(2 * pi * l[2] / 365) * log(t - 1), sin(2 * pi * l[2] / 365) * log(t - 1)^2, 
           cos(2 * pi * l[2] / 365) * log(t - 1), cos(2 * pi * l[2] / 365) * log(t - 1)^2, 
           logdist[1] * log(t - 1), logdist[1] * log(t - 1)^2, 0, 0, 0, 0, 0, 0, 0)
lines(X %*% beta, col = "gray")
# Madrid / January 15 / I_{t,l-1} = 1, I_{t,l-2} = 0
X <- cbind(sin(2 * pi * l[1] / 365), cos(2 * pi * l[1] / 365), 
           logdist[1], log(t - 1), log(t - 1)^2, 
           1, 0, 
           sin(2 * pi * l[1] / 365) * log(t - 1), sin(2 * pi * l[1] / 365) * log(t - 1)^2, 
           cos(2 * pi * l[1] / 365) * log(t - 1), cos(2 * pi * l[1] / 365) * log(t - 1)^2, 
           logdist[1] * log(t - 1), logdist[1] * log(t - 1)^2, 
           log(t - 1), logdist[1], 0, 0, 0, 0, 0)
lines(X %*% beta, lty = 2)
# Madrid / July 15 / I_{t,l-1} = 1, I_{t,l-2} = 0
X <- cbind(sin(2 * pi * l[2] / 365), cos(2 * pi * l[2] / 365), 
           logdist[1], log(t - 1), log(t - 1)^2, 
           1, 0, 
           sin(2 * pi * l[2] / 365) * log(t - 1), sin(2 * pi * l[2] / 365) * log(t - 1)^2, 
           cos(2 * pi * l[2] / 365) * log(t - 1), cos(2 * pi * l[2] / 365) * log(t - 1)^2, 
           logdist[1] * log(t - 1), logdist[1] * log(t - 1)^2, 
           log(t - 1), logdist[1], 0, 0, 0, 0, 0)
lines(X %*% beta, col = "gray", lty = 2)
