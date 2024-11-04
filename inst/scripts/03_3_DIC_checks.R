#######################################################
### Section 4 - DIC and POSTERIOR PREDICTIVE CHECKS ###
#######################################################

# input  data: data, "data/model.rds", "data/model_others.RData"

### DIC
# M0
metrics(model[,1]$y, 
        matrix(rep(1 / 2:62, 365 * 40), 1000, 61*365*40, byrow = TRUE), 
        metric = "DIC")

# M1
p <- matrix(nrow = 1000, ncol = 890600)
set.seed(12345)
p[1:500, ] <- 
  predict.rom(
    object  = M1[, 1], 
    newdata = data,
    type = "sample")
p[501:1000, ] <- 
  predict.rom(
    object  = M1[, 2], 
    newdata = data,
    type = "sample")
metrics(M1[,1]$y, p, metric = "DIC")

# M2
p <- matrix(nrow = 1000, ncol = 890600)
set.seed(12345)
p[1:500, ] <- 
  predict.rom(
    object  = M2[, 1], 
    newdata = data,
    type = "sample")
p[501:1000, ] <- 
  predict.rom(
    object  = M2[, 2], 
    newdata = data,
    type = "sample")
metrics(M2[,1]$y, p, metric = "DIC")

# M3
p <- matrix(nrow = 1000, ncol = 890600)
set.seed(12345)
p[1:500, ] <- 
  predict.rom(
    object  = M3[, 1], 
    newdata = data,
    type = "sample")
p[501:1000, ] <- 
  predict.rom(
    object  = M3[, 2], 
    newdata = data,
    type = "sample")
metrics(M3[,1]$y, p, metric = "DIC")

# M4
p <- matrix(nrow = 1000, ncol = 890600)
set.seed(12345)
p[1:500, ] <- 
  predict.rom(
    object  = M4[, 1], 
    newdata = data,
    type = "sample")
p[501:1000, ] <- 
  predict.rom(
    object  = M4[, 2], 
    newdata = data,
    type = "sample")
metrics(M4[,1]$y, p, metric = "DIC")

# M5
p <- matrix(nrow = 1000, ncol = 890600)
set.seed(12345)
p[1:500, ] <- 
  predict.rom(
    object  = model[, 1], 
    newdata = data,
    type = "sample")
p[501:1000, ] <- 
  predict.rom(
    object  = model[, 2], 
    newdata = data,
    type = "sample")
metrics(model[,1]$y, p, metric = "DIC")



### Posterior predictive distribution
I <- matrix(rbinom(1000 * 890600, 1, p), 1000, 890600)
Ipre <- array(I, dim = c(1000, 61, 365, 40))
Iobs <- array(model[,1]$y, dim = c(61, 365, 40))

# Preparation
# Seasons
season <- list(
  DJF = c(1:59, 335:365), 
  MAM = 60:151, 
  JJA = 152:243, 
  SON = 244:334)
# Months
month <- list(
  jan = 1:31,
  feb = 32:59,
  mar = 60:90,
  apr = 91:120,
  may = 121:151,
  jun = 152:181,
  jul = 182:212,
  aug = 213:243,
  sep = 244:273,
  oct = 274:304,
  nov = 305:334,
  dec = 335:365)

# Features
# Nall
Nallobs <- rep(NA, 12 * 40)
Nallpre <- matrix(ncol = 1000, nrow = 12 * 40)
j <- 0; t <- 1:61
for (i in 1:40) {
  for (m in 1:12) {
    j <- j + 1
    l <- month[[m]]
    L <- length(l)
    Nallobs[j]   <- sum(Iobs[t, l, i]) / L
    Nallpre[j, ] <- apply(Ipre[, t, l, i], 1, sum) / L
  }
}

# N21st
N21stobs <- rep(NA, 12 * 40)
N21stpre <- matrix(ncol = 1000, nrow = 12 * 40)
j <- 0; t <- 40:61
for (i in 1:40) {
  for (m in 1:12) {
    j <- j + 1
    l <- month[[m]]
    L <- length(l)
    N21stobs[j]   <- sum(Iobs[t, l, i]) / L
    N21stpre[j, ] <- apply(Ipre[, t, l, i], 1, sum) / L
  }
}

# R
E <- sum(1 / 53:62)
Robs <- rep(NA, 12 * 40)
Rpre <- matrix(ncol = 1000, nrow = 12 * 40)
j <- 0; t <- 52:61
for (i in 1:40) {
  for (m in 1:12) {
    j <- j + 1
    l <- month[[m]]
    L <- length(l)
    Robs[j]   <- sum(Iobs[t, l, i]) / (L * E)
    Rpre[j, ] <- apply(Ipre[, t, l, i], 1, sum) / (L * E)
  }
}

# ERS
ERSobs <- rep(NA, 4 * 61)
ERSpre <- matrix(ncol = 1000, nrow = 4 * 61)
j <- 0
for (t in 1:61) {
  for (m in 1:4) {
    j <- j + 1
    l <- season[[m]]
    L <- length(l)
    ERSobs[j]   <- (t + 1) * sum(Iobs[t, l, ]) / (40 * L)
    ERSpre[j, ] <- (t + 1) * apply(Ipre[, t, l, ], 1, sum) / (40 * L)
  }
}

# PIT
PIThist <- function(I, p, J = 10,
                    title = NULL,
                    picture.name = "photo.pdf",
                    save = TRUE) {
  
  func <- function(u, I, p) {
    
    p0 <- rowMeans(I <  p)
    p1 <- rowMeans(I <= p)
    
    ifelse(u <= p0, 0, ifelse(u > p1, 1, (u - p0) / (p1 - p0)))
  }
  
  f <- rep(0, J); f[J] <- 1
  for (j in 1:(J - 1)) {
    f[j] <- mean(func(j / J, I, p))
  }
  
  f <- diff(c(0, f))
  
  df <- data.frame("PIT" = 1:J / J, "Relative Frequency" = J * f)
  
  gg <- ggplot2::ggplot(df, ggplot2::aes(PIT, Relative.Frequency)) + 
    ggplot2::geom_bar(stat = 'identity', just = 1, width = 1 / J, col = "black", fill = "white") +
    ggplot2::theme_classic() +
    ggplot2::labs(x="PIT", y="Relative Frequency", title=title) +
    ggplot2::ylim(c(0, 1.75))
  
  if (save) {
    ggplot2::ggsave(picture.name, gg, width = 4, height = 4)
  }
  
  gg
}

PIThist(Nallobs, Nallpre, 
        title = bquote(bar(N)["62,month"](s)), 
        picture.name = "inst/img/SUPP_PIT_Nall.pdf")
PIThist(N21stobs, N21stpre, 
        title = bquote(bar(N)["41:62,month"](s)), 
        picture.name = "inst/img/SUPP_PIT_N21st.pdf")
PIThist(Robs, Rpre, 
        title = bquote(R["53:62,month"](s)), 
        picture.name = "inst/img/SUPP_PIT_R.pdf")
PIThist(ERSobs, ERSpre, 
        title = expression(t %*% widehat(bar(ERS))["t,season"](D)), 
        picture.name = "inst/img/SUPP_PIT_ERS.pdf")



# Observed vs predicted 
OPplot <- function(I, p,
                   lower.lim = 0,
                   title = NULL,
                   picture.name = "photo.pdf",
                   save = TRUE) {
  
  df <- data.frame(
    cbind(I, rowMeans(p),
          t(apply(p, 1, quantile, prob = c(0.05, 0.95), type = 2))
    )
  )
  colnames(df) <- c("Xobs", "Xpre", "q05", "q95")
  
  df$col <- apply(df[,c(1,3,4)], 1, function(x) (x[1] - x[2]) * (x[1] - x[3]) < 0)
  print(mean(df$col))
  
  gg <- ggplot2::ggplot(df) +
    ggplot2::geom_abline(intercept=0, slope=1) +
    ggplot2::geom_linerange(ggplot2::aes(y=Xobs, xmin=q05, xmax=q95, color=col), alpha=0.1) +
    ggplot2::geom_point(ggplot2::aes(x=Xpre, y=Xobs, color=col), alpha=0.8) +
    ggplot2::scale_color_manual(values=c("red", "black")) +
    ggplot2::labs(x='Predicted', y='Observed', title=title) +
    ggplot2::xlim(c(lower.lim, ceiling(max(df[,1:4])))) + 
    ggplot2::ylim(c(lower.lim, ceiling(max(df[,1:4])))) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position="none") +
    ggplot2::coord_fixed()
  
  if (save) {
    ggplot2::ggsave(picture.name, gg, width = 4, height = 4)
  }
  
  gg
}

OPplot(Nallobs + 1, Nallpre + 1, 
       lower.lim = 3,
       title = bquote(bar(N)["62,month"](s)), 
       picture.name = "inst/img/SUPP_OP_Nall.pdf")
OPplot(N21stobs, N21stpre, 
       title = bquote(bar(N)["41:62,month"](s)), 
       picture.name = "inst/img/SUPP_OP_N21st.pdf")
OPplot(Robs, Rpre, 
        title = bquote(R["53:62,month"](s)), 
        picture.name = "inst/img/SUPP_OP_R.pdf")
OPplot(ERSobs, ERSpre, 
        title = expression(t %*% widehat(bar(ERS))["t,season"](D)), 
        picture.name = "inst/img/SUPP_OP_ERS.pdf")
