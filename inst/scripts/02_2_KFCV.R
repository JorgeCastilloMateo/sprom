############################################
### Section 4 - 10-FOLD CROSS-VALIDATION ###
############################################

# input  data: coords, data, I365, (updated)stations, Y365
# output data: "data/KFCV.RData", "data/KFCVfullmodel.rds"

# Clear workspace
rm(list = setdiff(ls(), c("background", "coords", "data", "grid", 
                          "I365", "limits", "stations", "Y365")))

# Indexes
TT <- nrow(I365)
LL <- ncol(I365)
SS <- dim(I365)[3]
K <- 10

# Choose folds at random
set.seed(23)
folds <- matrix(nrow = SS / K, ncol = K)
remain <- 1:SS
for (k in 1:K) {
  folds[,k] <- sort(sample(remain, SS / K))
  remain    <- remain[!(remain %in% folds[,k])]
}

### Parameters: [iter, chain, params, n]
KFCVmodels    <- list()
KFCVfullmodel <- list()
KFCVmetric <- matrix(nrow = 6, ncol = 67,
  dimnames = list(
    model  = c(paste0("M", 0:5)), 
    metric = c("BS", "BS1", "BS2", "AUC", "AUC1", "AUC2", paste0("ADt", 1:61)))
)

form <- Y ~ (sine + cosi + log(dist)) * poly(trend, 2) + (trend + log(dist)) * lag1 * lag2

METRICS <- matrix(nrow = K, ncol = 67)



### M0
# Metrics
set.seed(12345)
for (k in 1:K) {
  METRICS[k,1:6] <- metrics(
    data$Y[data$site %in% folds[,k]], 
    rep(1 / 2:TT, LL * SS / K), 
    data$year[data$site %in% folds[,k]],
    data$day[data$site %in% folds[,k]],
    data$site[data$site %in% folds[,k]],
    metric = "point")
  METRICS[k,7:67] <- metrics(
    data$Y[data$site %in% folds[,k]], 
    rep(1 / 2:TT, LL * SS / K), 
    data$year[data$site %in% folds[,k]],
    data$day[data$site %in% folds[,k]],
    data$site[data$site %in% folds[,k]],
    metric = "N")
}
(KFCVmetric[1,] <- colMeans(METRICS))

save(KFCVmodels, KFCVmetric, file = "data/KFCV.RData")



### M1 FIXED EFFECTS
# Chain 1 (10 folds in parallel): < 4 hours
time <- Sys.time()
cl <- parallel::makeCluster(K)
parallel::clusterExport(cl, c("form", "data", "folds"))
KFCVmodels[[1]] <- parallel::parSapply(cl = cl, X = 1:K, 
  FUN = function(iter) {set.seed(iter); sprom::rom(
    formula = form, data = data[!(data$site %in% folds[, iter]), ],
    model = "1",
    n.report = 10001, n.burnin = 10000, n.sims = 10000, n.thin = 20)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Chain 2 (10 folds in parallel): < 4 hours
time <- Sys.time()
cl <- parallel::makeCluster(K)
parallel::clusterExport(cl, c("form", "data", "folds"))
KFCVmodels[[2]] <- parallel::parSapply(cl = cl, X = 1:K, 
  FUN = function(iter) {set.seed(iter + 1000); sprom::rom(
    formula = form, data = data[!(data$site %in% folds[, iter]), ],
    model = "1",
    n.report = 10001, n.burnin = 10000, n.sims = 10000, n.thin = 20)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Checking convergence
v <- matrix(nrow = K, ncol = 2)
for (k in 1:K) {
  v[k,1] <- coda::gelman.diag(coda::mcmc.list(
    coda::mcmc(Reduce(cbind, KFCVmodels[[1]][,k]$params)),
    coda::mcmc(Reduce(cbind, KFCVmodels[[2]][,k]$params))))$mpsrf
  v[k,2] <- min((
    coda::effectiveSize(Reduce(cbind, KFCVmodels[[1]][,k]$params)) +
    coda::effectiveSize(Reduce(cbind, KFCVmodels[[2]][,k]$params))))
}
v

# Metrics
set.seed(12345)
for (k in 1:K) {
  print(k)
  p1 <- sprom::predict.rom(
    object  = KFCVmodels[[1]][,k], 
    newdata = data[data$site %in% folds[,k],])
  p2 <- sprom::predict.rom(
    object  = KFCVmodels[[2]][,k], 
    newdata = data[data$site %in% folds[,k],])
  METRICS[k,1:6] <- sprom::metrics(
    data$Y[data$site %in% folds[,k]], 
    rbind(p1, p2),
    data$year[data$site %in% folds[,k]],
    data$day[data$site %in% folds[,k]],
    data$site[data$site %in% folds[,k]],
    metric = "point")
  METRICS[k,7:67] <- sprom::metrics(
    data$Y[data$site %in% folds[,k]], 
    rbind(p1, p2),
    data$year[data$site %in% folds[,k]],
    data$day[data$site %in% folds[,k]],
    data$site[data$site %in% folds[,k]],
    metric = "N")
}
(KFCVmetric[2,] <- colMeans(METRICS))

save(KFCVmodels, KFCVmetric, file = "data/KFCV.RData")



### M2
# Chain 1 (10 folds in parallel): < 24 hours
time <- Sys.time()
cl <- parallel::makeCluster(K)
parallel::clusterExport(cl, c("form", "data", "folds", "coords", "SS"))
KFCVmodels[[3]] <- parallel::parSapply(cl = cl, X = 1:K, 
  FUN = function(iter) {set.seed(iter); sprom::rom(
    formula = form, data = data[!(data$site %in% folds[, iter]), ],
    coords = coords[!(1:SS %in% folds[, iter]), ],
    site = data$site[!(data$site %in% folds[, iter])],
    year = data$year[!(data$site %in% folds[, iter])],
    day  = data$day[!(data$site %in% folds[, iter])],
    model = "2",
    decay.prior = "unif",
    prior = list("beta" = 1e-04, "sigma" = c(2, 1), "decay" = c(3 / 906.8, 3 / 9.1)),
    n.report = 100001, n.burnin = 50000, n.sims = 50000, n.thin = 100)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Chain 2 (10 folds in parallel): < 24 hours
time <- Sys.time()
cl <- parallel::makeCluster(K)
parallel::clusterExport(cl, c("form", "data", "folds", "coords", "SS"))
KFCVmodels[[4]] <- parallel::parSapply(cl = cl, X = 1:K, 
  FUN = function(iter) {set.seed(iter + 100); sprom::rom(
    formula = form, data = data[!(data$site %in% folds[, iter]), ],
    coords = coords[!(1:SS %in% folds[, iter]), ],
    site = data$site[!(data$site %in% folds[, iter])],
    year = data$year[!(data$site %in% folds[, iter])],
    day  = data$day[!(data$site %in% folds[, iter])],
    model = "2",
    decay.prior = "unif",
    prior = list("beta" = 1e-04, "sigma" = c(2, 1), "decay" = c(3 / 906.8, 3 / 9.1)),
    n.report = 100001, n.burnin = 50000, n.sims = 50000, n.thin = 100)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Checking convergence (beta's for illustration)
v <- matrix(nrow = K, ncol = 2)
for (k in 1:K) {
  v[k,1] <- coda::gelman.diag(coda::mcmc.list(
    coda::mcmc(KFCVmodels[[3]][,k]$params$beta),
    coda::mcmc(KFCVmodels[[4]][,k]$params$beta)))$mpsrf
  v[k,2] <- min(
    coda::effectiveSize(KFCVmodels[[3]][,k]$params$beta) +
    coda::effectiveSize(KFCVmodels[[4]][,k]$params$beta))
}
v

# Metrics
set.seed(12345)
for (k in 1:K) {
  print(k)
  p1 <- sprom::predict.rom(
    object  = KFCVmodels[[3]][,k], 
    newdata = data[data$site %in% folds[,k],],
    newcoords = coords[folds[,k],],
    newsite = data[data$site %in% folds[,k],]$site,
    newyear = data[data$site %in% folds[,k],]$year,
    newday  = data[data$site %in% folds[,k],]$day)
  p2 <- sprom::predict.rom(
    object  = KFCVmodels[[4]][,k], 
    newdata = data[data$site %in% folds[,k],],
    newcoords = coords[folds[,k],],
    newsite = data[data$site %in% folds[,k],]$site,
    newyear = data[data$site %in% folds[,k],]$year,
    newday  = data[data$site %in% folds[,k],]$day)
  METRICS[k,1:6] <- sprom::metrics(
    data$Y[data$site %in% folds[,k]], 
    rbind(p1, p2),
    data$year[data$site %in% folds[,k]],
    data$day[data$site %in% folds[,k]],
    data$site[data$site %in% folds[,k]],
    metric = "point")
  METRICS[k,7:67] <- sprom::metrics(
    data$Y[data$site %in% folds[,k]], 
    rbind(p1, p2),
    data$year[data$site %in% folds[,k]],
    data$day[data$site %in% folds[,k]],
    data$site[data$site %in% folds[,k]],
    metric = "N")
}
(KFCVmetric[3,] <- colMeans(METRICS))

save(KFCVmodels, KFCVmetric, file = "data/KFCV.RData")



### M3
# Chain 1 (10 folds in parallel): < 24 hours
time <- Sys.time()
cl <- parallel::makeCluster(K)
parallel::clusterExport(cl, c("form", "data", "folds", "coords", "SS"))
KFCVmodels[[5]] <- parallel::parSapply(cl = cl, X = 1:K, 
  FUN = function(iter) {set.seed(iter); sprom::rom(
    formula = form, data = data[!(data$site %in% folds[, iter]), ],
    coords = coords[!(1:SS %in% folds[, iter]), ],
    site = data$site[!(data$site %in% folds[, iter])],
    year = data$year[!(data$site %in% folds[, iter])],
    day  = data$day[!(data$site %in% folds[, iter])],
    model = "3",
    decay.prior = "unif",
    prior = list("beta" = 1e-04, "sigma" = c(2, 1), "decay" = c(3 / 906.8, 3 / 9.1)),
    n.report = 100001, n.burnin = 50000, n.sims = 50000, n.thin = 100)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Chain 2 (10 folds in parallel): < 24 hours
time <- Sys.time()
cl <- parallel::makeCluster(K)
parallel::clusterExport(cl, c("form", "data", "folds", "coords", "SS"))
KFCVmodels[[6]] <- parallel::parSapply(cl = cl, X = 1:K, 
  FUN = function(iter) {set.seed(iter + 100); sprom::rom(
    formula = form, data = data[!(data$site %in% folds[, iter]), ],
    coords = coords[!(1:SS %in% folds[, iter]), ],
    site = data$site[!(data$site %in% folds[, iter])],
    year = data$year[!(data$site %in% folds[, iter])],
    day  = data$day[!(data$site %in% folds[, iter])],
    model = "3",
    decay.prior = "unif",
    prior = list("beta" = 1e-04, "sigma" = c(2, 1), "decay" = c(3 / 906.8, 3 / 9.1)),
    n.report = 100001, n.burnin = 50000, n.sims = 50000, n.thin = 100)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Checking convergence (beta's for illustration)
v <- matrix(nrow = K, ncol = 2)
for (k in 1:K) {
  v[k,1] <- coda::gelman.diag(coda::mcmc.list(
    coda::mcmc(KFCVmodels[[5]][,k]$params$beta),
    coda::mcmc(KFCVmodels[[6]][,k]$params$beta)))$mpsrf
  v[k,2] <- min(
    coda::effectiveSize(KFCVmodels[[5]][,k]$params$beta) +
    coda::effectiveSize(KFCVmodels[[6]][,k]$params$beta))
}
v

# Metrics
set.seed(12345)
for (k in 1:K) {
  print(k)
  p1 <- sprom::predict.rom(
    object  = KFCVmodels[[5]][,k], 
    newdata = data[data$site %in% folds[,k],],
    newcoords = coords[folds[,k],],
    newsite = data[data$site %in% folds[,k],]$site,
    newyear = data[data$site %in% folds[,k],]$year,
    newday  = data[data$site %in% folds[,k],]$day)
  p2 <- sprom::predict.rom(
    object  = KFCVmodels[[6]][,k], 
    newdata = data[data$site %in% folds[,k],],
    newcoords = coords[folds[,k],],
    newsite = data[data$site %in% folds[,k],]$site,
    newyear = data[data$site %in% folds[,k],]$year,
    newday  = data[data$site %in% folds[,k],]$day)
  METRICS[k,1:6] <- sprom::metrics(
    data$Y[data$site %in% folds[,k]], 
    rbind(p1, p2),
    data$year[data$site %in% folds[,k]],
    data$day[data$site %in% folds[,k]],
    data$site[data$site %in% folds[,k]],
    metric = "point")
  METRICS[k,7:67] <- sprom::metrics(
    data$Y[data$site %in% folds[,k]], 
    rbind(p1, p2),
    data$year[data$site %in% folds[,k]],
    data$day[data$site %in% folds[,k]],
    data$site[data$site %in% folds[,k]],
    metric = "N")
}
(KFCVmetric[4,] <- colMeans(METRICS))

save(KFCVmodels, KFCVmetric, file = "data/KFCV.RData")



### M4 
# Chain 1 (10 folds in parallel): < 24 hours
time <- Sys.time()
cl <- parallel::makeCluster(K)
parallel::clusterExport(cl, c("form", "data", "folds", "coords", "SS"))
KFCVmodels[[7]] <- parallel::parSapply(cl = cl, X = 1:K, 
  FUN = function(iter) {set.seed(iter); sprom::rom(
    formula = form, data = data[!(data$site %in% folds[, iter]), ],
    coords = coords[!(1:SS %in% folds[, iter]), ],
    site = data$site[!(data$site %in% folds[, iter])],
    year = data$year[!(data$site %in% folds[, iter])],
    day  = data$day[!(data$site %in% folds[, iter])],
    model = "4",
    decay.prior = "unif",
    prior = list("beta" = 1e-04, "sigma" = c(2, 1), "decay" = c(3 / 906.8, 3 / 9.1)),
    n.report = 100001, n.burnin = 50000, n.sims = 50000, n.thin = 100)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Chain 2 (10 folds in parallel): < 24 hours
time <- Sys.time()
cl <- parallel::makeCluster(K)
parallel::clusterExport(cl, c("form", "data", "folds", "coords", "SS"))
KFCVmodels[[8]] <- parallel::parSapply(cl = cl, X = 1:K, 
  FUN = function(iter) {set.seed(iter + 100); sprom::rom(
    formula = form, data = data[!(data$site %in% folds[, iter]), ],
    coords = coords[!(1:SS %in% folds[, iter]), ],
    site = data$site[!(data$site %in% folds[, iter])],
    year = data$year[!(data$site %in% folds[, iter])],
    day  = data$day[!(data$site %in% folds[, iter])],
    model = "4",
    decay.prior = "unif",
    prior = list("beta" = 1e-04, "sigma" = c(2, 1), "decay" = c(3 / 906.8, 3 / 9.1)),
    n.report = 100001, n.burnin = 50000, n.sims = 50000, n.thin = 100)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Checking convergence (beta's for illustration)
v <- matrix(nrow = K, ncol = 2)
for (k in 1:K) {
  v[k,1] <- coda::gelman.diag(coda::mcmc.list(
    coda::mcmc(KFCVmodels[[7]][,k]$params$beta),
    coda::mcmc(KFCVmodels[[8]][,k]$params$beta)))$mpsrf
  v[k,2] <- min(
    coda::effectiveSize(KFCVmodels[[7]][,k]$params$beta) +
    coda::effectiveSize(KFCVmodels[[8]][,k]$params$beta))
}
v

# Metrics
set.seed(12345)
for (k in 1:K) {
  print(k)
  p1 <- sprom::predict.rom(
    object  = KFCVmodels[[7]][,k], 
    newdata = data[data$site %in% folds[,k],],
    newcoords = coords[folds[,k],],
    newsite = data[data$site %in% folds[,k],]$site,
    newyear = data[data$site %in% folds[,k],]$year,
    newday  = data[data$site %in% folds[,k],]$day)
  p2 <- sprom::predict.rom(
    object  = KFCVmodels[[8]][,k], 
    newdata = data[data$site %in% folds[,k],],
    newcoords = coords[folds[,k],],
    newsite = data[data$site %in% folds[,k],]$site,
    newyear = data[data$site %in% folds[,k],]$year,
    newday  = data[data$site %in% folds[,k],]$day)
  METRICS[k,1:6] <- sprom::metrics(
    data$Y[data$site %in% folds[,k]], 
    rbind(p1, p2),
    data$year[data$site %in% folds[,k]],
    data$day[data$site %in% folds[,k]],
    data$site[data$site %in% folds[,k]],
    metric = "point")
  METRICS[k,7:67] <- sprom::metrics(
    data$Y[data$site %in% folds[,k]], 
    rbind(p1, p2),
    data$year[data$site %in% folds[,k]],
    data$day[data$site %in% folds[,k]],
    data$site[data$site %in% folds[,k]],
    metric = "N")
}
(KFCVmetric[5,] <- colMeans(METRICS))

save(KFCVmodels, KFCVmetric, file = "data/KFCV.RData")



### M5 FULL MODEL
# Chain 1 (10 folds in parallel): ~ 78 hours
time <- Sys.time()
cl <- parallel::makeCluster(K)
parallel::clusterExport(cl, c("form", "data", "folds", "coords", "SS"))
KFCVfullmodel[[1]] <- parallel::parSapply(cl = cl, X = 1:K, 
  FUN = function(iter) {set.seed(iter); sprom::rom(
    formula = form, data = data[!(data$site %in% folds[, iter]), ],
    coords = coords[!(1:SS %in% folds[, iter]), ],
    site = data$site[!(data$site %in% folds[, iter])],
    year = data$year[!(data$site %in% folds[, iter])],
    day  = data$day[!(data$site %in% folds[, iter])],
    model = "5",
    n.report = 100001, n.burnin = 100000, n.sims = 100000, n.thin = 200)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Chain 2 (10 folds in parallel): ~ 78 hours
time <- Sys.time()
cl <- parallel::makeCluster(K)
parallel::clusterExport(cl, c("form", "data", "folds", "coords", "SS"))
KFCVfullmodel[[2]] <- parallel::parSapply(cl = cl, X = 1:K, 
  FUN = function(iter) {set.seed(iter + 100); sprom::rom(
    formula = form, data = data[!(data$site %in% folds[, iter]), ],
    coords = coords[!(1:SS %in% folds[, iter]), ],
    site = data$site[!(data$site %in% folds[, iter])],
    year = data$year[!(data$site %in% folds[, iter])],
    day  = data$day[!(data$site %in% folds[, iter])],
    model = "5",
    n.report = 100001, n.burnin = 100000, n.sims = 100000, n.thin = 200)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Checking convergence (beta's for illustration)
v <- matrix(nrow = K, ncol = 2)
for (k in 1:K) {
  v[k,1] <- coda::gelman.diag(coda::mcmc.list(
    coda::mcmc(KFCVfullmodel[[1]][,k]$params$beta),
    coda::mcmc(KFCVfullmodel[[2]][,k]$params$beta)))$mpsrf
  v[k,2] <- min(
    coda::effectiveSize(KFCVfullmodel[[1]][,k]$params$beta) +
    coda::effectiveSize(KFCVfullmodel[[2]][,k]$params$beta))
}
v

saveRDS(KFCVfullmodel, file = "data/KFCVfullmodel.rds")

# Metrics
set.seed(12345)
p <- matrix(nrow = 1000, ncol = SS * (TT - 1) * LL)
for (k in 1:K) {
  print(k)
  p[1:500, 1:(SS/K * (TT - 1) * LL) + (k - 1) * SS/K * (TT - 1) * LL] <- 
    p1 <- sprom::predict.rom(
      object  = KFCVfullmodel[[1]][,k], 
      newdata = data[data$site %in% folds[,k],],
      newcoords = coords[folds[,k],],
      newsite = data[data$site %in% folds[,k],]$site,
      newyear = data[data$site %in% folds[,k],]$year,
      newday  = data[data$site %in% folds[,k],]$day)
  p[501:1000, 1:(SS/K * (TT - 1) * LL) + (k - 1) * SS/K * (TT - 1) * LL] <- 
    p2 <- sprom::predict.rom(
      object  = KFCVfullmodel[[2]][,k], 
      newdata = data[data$site %in% folds[,k],],
      newcoords = coords[folds[,k],],
      newsite = data[data$site %in% folds[,k],]$site,
      newyear = data[data$site %in% folds[,k],]$year,
      newday  = data[data$site %in% folds[,k],]$day)
  METRICS[k,1:6] <- sprom::metrics(
    data$Y[data$site %in% folds[,k]], 
    rbind(p1, p2),
    data$year[data$site %in% folds[,k]],
    data$day[data$site %in% folds[,k]],
    data$site[data$site %in% folds[,k]],
    metric = "point")
  METRICS[k,7:67] <- sprom::metrics(
    data$Y[data$site %in% folds[,k]], 
    rbind(p1, p2),
    data$year[data$site %in% folds[,k]],
    data$day[data$site %in% folds[,k]],
    data$site[data$site %in% folds[,k]],
    metric = "N")
}
(KFCVmetric[6,] <- colMeans(METRICS))

save(KFCVmodels, KFCVmetric, file = "data/KFCV.RData")



### AD metric plot
pdf("inst/img/SUPP_ADmetric.pdf", width = 5.6, height = 4)
plot(x = 1:62,
     y = c(0, KFCVmetric[1, 7:67]), 
     xlab = "t (year)", ylab = "",
     xaxt='n', ylim = c(0, 2),
     type = "l", lty = 6)
axis(side = 1, at = c(1, 21, 41, 61), 
     labels = c("1 (1960)", "21 (1980)", "41 (2000)", "61 (2020)"))
title(ylab = expression(AD[t]), 
      line = 2, cex.lab = 1.2)
lines(c(0, KFCVmetric[2, 7:67]), lty = 5)
lines(c(0, KFCVmetric[3, 7:67]), lty = 4)
lines(c(0, KFCVmetric[4, 7:67]), lty = 3)
lines(c(0, KFCVmetric[5, 7:67]), lty = 2)
lines(c(0, KFCVmetric[6, 7:67]), lty = 1)
legend("bottomright", 
       legend=c(expression(M[0]), expression(M[1]), expression(M[2]), 
                expression(M[3]), expression(M[4]), expression(M[5])),
       lty=6:1, cex=0.8)
dev.off()

### AUC by decade and location in M5
paux <- array(colMeans(p), dim = c(61, 365, 40))
Iaux <- array(data$Y,      dim = c(61, 365, 40))[,,c(folds)]

paux[!(Iaux %in% 0:1)] <- NA
Iaux[!(Iaux %in% 0:1)] <- NA

# Preparation
# Seasons
season <- list(
  DJF = c(1:59, 335:365), 
  MAM = 60:151, 
  JJA = 152:243, 
  SON = 244:334)
# Decades
decade <- list(
  D1 = 1:10,
  D2 = 11:20,
  D3 = 21:30,
  D4 = 31:40,
  D5 = 41:50,
  D6 = 51:61)

# Compute AUC
AUC <- array(dim = c(6, 4, 40))
for (t in 1:6) {
  for (l in 1:4) {
    for (i in 1:40) {
      AUC[t, l, i] <- pROC::auc(
        c(Iaux[decade[[t]], season[[l]], i]), 
        c(paux[decade[[t]], season[[l]], i]))
    }
  }
}

stations$NAME2 <-  
  c("BADAJOZ", "MADRID (RETIRO)", "MALAGA", "NAVACERRADA", "SALAMANCA", 
    "SAN SEBASTIAN", "TORTOSA", "VALENCIA", "ZARAGOZA", "BARCELONA (FABRA)",
    "ALBACETE", "BURGOS", "CIUDAD REAL", "CORUNA", "MURCIA", "SEVILLA", 
    "SORIA", "BILBAO", "SANTIAGO", "PONFERRADA", "LEON", "LOGRONO", "ZAMORA", 
    "REUS", "BARCELONA (AIR)", "MADRID (TORREJON)", "VITORIA", "ALMERIA", 
    "GIJON", "CACERES", "SANTANDER", "CASTELLON", "HUELVA", "LLEIDA", 
    "MADRID (BARAJAS)", "MADRID (4VIENTOS)", "MADRID (GETAFE)", "MORON", 
    "VALLADOLID", "DAROCA")

name_seasons <- c("Winter", "Spring", "Summer", "Autumn")

# Plot/Table with the AUC
pdf("inst/img/SUPP_AUCmatrix.pdf", width = 1.2 * 8.27, height = 0.6 * 11.69)
par(mfrow = c(4, 1))
par(mar = c(0,0,0,0), oma = c(0,0,6,0), las = 1)
for (l in 1:4) {
  aux <- AUC[,l,]
  rownames(aux) <- c("1961-1970", "1971-1980", "1981-1990",
                     "1991-2000", "2001-2010", "2011-2021")
  colnames(aux) <- stations$NAME2[c(folds)]
  aux <- aux[,order(stations$HGHT[c(folds)])]
  
  corrplot::corrplot(
    aux, method = 'number', 
    col = colorRampPalette(RColorBrewer::brewer.pal(8, "RdBu"))(8), 
    col.lim = c(0.5, 1),
    is.corr = FALSE,
    number.cex = 0.8,
    tl.col = "black",
    tl.srt = 45,
    tl.pos = "l",
    tl.offset = 1)
  mtext(name_seasons[l], at = -1.5, line = -1)
  
  if (l == 1) 
    text(x = 1:40, y = 7, labels = colnames(aux), srt = 45, adj = 0, col = "black", xpd = NA)
}
dev.off()

mean(AUC >= 0.8)
mean(AUC >= 0.9)
