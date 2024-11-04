##################################################
### Section 4 - MODEL FITTING (M1, M2, M3, M4) ###
##################################################

# input  data: coords, data
# output data: "data/model_others.RData"

# M1 - 2 chains in parallel < 4 hours
time <- Sys.time()
cl <- parallel::makeCluster(2)
parallel::clusterExport(cl, c("data", "coords"))
M1 <- parallel::parSapply(cl = cl, X = 1:2, 
  FUN = function(iter) {set.seed(123 * iter + 1234); sprom::rom(
    formula = Y ~ (sine + cosi + log(dist)) * poly(trend, 2) + (trend + log(dist)) * lag1 * lag2, 
    data = data,
    model = "1",
    n.report = 10001, n.burnin = 10000, n.sims = 10000, n.thin = 20)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Checking convergence (PSRF and ESS)
v <- rep(NA, 2)
# beta
v[1] <- coda::gelman.diag(coda::mcmc.list(
  coda::mcmc(M1[,1]$params$beta),
  coda::mcmc(M1[,2]$params$beta)))$mpsrf
v[2] <- min(
  coda::effectiveSize(M1[,1]$params$beta) +
  coda::effectiveSize(M1[,2]$params$beta))
v
# betal1
v[1] <- coda::gelman.diag(coda::mcmc.list(
  coda::mcmc(M1[,1]$params$betal1),
  coda::mcmc(M1[,2]$params$betal1)))$mpsrf
v[2] <- min(
  coda::effectiveSize(M1[,1]$params$betal1) +
  coda::effectiveSize(M1[,2]$params$betal1))
v
# betal2
v[1] <- coda::gelman.diag(coda::mcmc.list(
  coda::mcmc(M1[,1]$params$betal2),
  coda::mcmc(M1[,2]$params$betal2)))$mpsrf
v[2] <- min(
  coda::effectiveSize(M1[,1]$params$betal2) +
  coda::effectiveSize(M1[,2]$params$betal2))
v

# Save data
save(M1, file = "data/model_others.RData")



# M2 - 2 chains in parallel < 18 hours 
time <- Sys.time()
cl <- parallel::makeCluster(2)
parallel::clusterExport(cl, c("data", "coords"))
M2 <- parallel::parSapply(cl = cl, X = 1:2, 
  FUN = function(iter) {set.seed(123 * iter + 1234); sprom::rom(
    formula = Y ~ (sine + cosi + log(dist)) * poly(trend, 2) + (trend + log(dist)) * lag1 * lag2, 
    data = data,
    coords = coords,
    site = data$site,
    year = data$year,
    day  = data$day,
    model = "2",
    decay.prior = "unif",
    prior = list("beta" = 1e-04, "sigma" = c(2, 1), "decay" = c(3 / 906.8, 3 / 9.1)),
    n.report = 100001, n.burnin = 50000, n.sims = 50000, n.thin = 100)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Checking convergence (PSRF and ESS)
v <- rep(NA, 2)
# beta
v[1] <- coda::gelman.diag(coda::mcmc.list(
  coda::mcmc(M2[,1]$params$beta),
  coda::mcmc(M2[,2]$params$beta)))$mpsrf
v[2] <- min(
  coda::effectiveSize(M2[,1]$params$beta) +
    coda::effectiveSize(M2[,2]$params$beta))
v
# betal1
v[1] <- coda::gelman.diag(coda::mcmc.list(
  coda::mcmc(M2[,1]$params$betal1),
  coda::mcmc(M2[,2]$params$betal1)))$mpsrf
v[2] <- min(
  coda::effectiveSize(M2[,1]$params$betal1) +
    coda::effectiveSize(M2[,2]$params$betal1))
v
# betal2
v[1] <- coda::gelman.diag(coda::mcmc.list(
  coda::mcmc(M2[,1]$params$betal2),
  coda::mcmc(M2[,2]$params$betal2)))$mpsrf
v[2] <- min(
  coda::effectiveSize(M2[,1]$params$betal2) +
    coda::effectiveSize(M2[,2]$params$betal2))
v
# ws
v[1] <- coda::gelman.diag(coda::mcmc.list(
  coda::mcmc(M2[,1]$params$ws),
  coda::mcmc(M2[,2]$params$ws)))$mpsrf
v[2] <- min(
  coda::effectiveSize(M2[,1]$params$ws) +
    coda::effectiveSize(M2[,2]$params$ws))
v
# hp
v[1] <- coda::gelman.diag(coda::mcmc.list(
  coda::mcmc(M2[,1]$params$hp),
  coda::mcmc(M2[,2]$params$hp)))$mpsrf
v[2] <- min(
  coda::effectiveSize(M2[,1]$params$hp) +
    coda::effectiveSize(M2[,2]$params$hp))
v

# Save data
save(M1, M2, file = "data/model_others.RData")



# M3 - 2 chains in parallel < 18 hours 
time <- Sys.time()
cl <- parallel::makeCluster(2)
parallel::clusterExport(cl, c("data", "coords"))
M3 <- parallel::parSapply(cl = cl, X = 1:2, 
  FUN = function(iter) {set.seed(123 * iter + 1234); sprom::rom(
    formula = Y ~ (sine + cosi + log(dist)) * poly(trend, 2) + (trend + log(dist)) * lag1 * lag2, 
    data = data,
    coords = coords,
    site = data$site,
    year = data$year,
    day  = data$day,
    model = "3",
    decay.prior = "unif",
    prior = list("beta" = 1e-04, "sigma" = c(2, 1), "decay" = c(3 / 906.8, 3 / 9.1)),
    n.report = 100001, n.burnin = 50000, n.sims = 50000, n.thin = 100)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Checking convergence (PSRF and ESS)
v <- rep(NA, 2)
# beta
v[1] <- coda::gelman.diag(coda::mcmc.list(
  coda::mcmc(M3[,1]$params$beta),
  coda::mcmc(M3[,2]$params$beta)))$mpsrf
v[2] <- min(
  coda::effectiveSize(M3[,1]$params$beta) +
    coda::effectiveSize(M3[,2]$params$beta))
v
# betal1
v[1] <- coda::gelman.diag(coda::mcmc.list(
  coda::mcmc(M3[,1]$params$betal1),
  coda::mcmc(M3[,2]$params$betal1)))$mpsrf
v[2] <- min(
  coda::effectiveSize(M3[,1]$params$betal1) +
    coda::effectiveSize(M3[,2]$params$betal1))
v
# betal2
v[1] <- coda::gelman.diag(coda::mcmc.list(
  coda::mcmc(M3[,1]$params$betal2),
  coda::mcmc(M3[,2]$params$betal2)))$mpsrf
v[2] <- min(
  coda::effectiveSize(M3[,1]$params$betal2) +
    coda::effectiveSize(M3[,2]$params$betal2))
v
# ws
v[1] <- coda::gelman.diag(coda::mcmc.list(
  coda::mcmc(M3[,1]$params$ws),
  coda::mcmc(M3[,2]$params$ws)))$mpsrf
v[2] <- min(
  coda::effectiveSize(M3[,1]$params$ws) +
    coda::effectiveSize(M3[,2]$params$ws))
v
# wt
v[1] <- coda::gelman.diag(coda::mcmc.list(
  coda::mcmc(M3[,1]$params$wt),
  coda::mcmc(M3[,2]$params$wt)))$mpsrf
v[2] <- min(
  coda::effectiveSize(M3[,1]$params$wt) +
    coda::effectiveSize(M3[,2]$params$wt))
v
# hp
v[1] <- coda::gelman.diag(coda::mcmc.list(
  coda::mcmc(M3[,1]$params$hp),
  coda::mcmc(M3[,2]$params$hp)))$mpsrf
v[2] <- min(
  coda::effectiveSize(M3[,1]$params$hp) +
    coda::effectiveSize(M3[,2]$params$hp))
v

# Save data
save(M1, M2, M3, file = "data/model_others.RData")



# M4 - 2 chains in parallel < 18 hours 
time <- Sys.time()
cl <- parallel::makeCluster(2)
parallel::clusterExport(cl, c("data", "coords"))
M4 <- parallel::parSapply(cl = cl, X = 1:2, 
  FUN = function(iter) {set.seed(123 * iter + 1234); sprom::rom(
    formula = Y ~ (sine + cosi + log(dist)) * poly(trend, 2) + (trend + log(dist)) * lag1 * lag2, 
    data = data,
    coords = coords,
    site = data$site,
    year = data$year,
    day  = data$day,
    model = "4",
    decay.prior = "unif",
    prior = list("beta" = 1e-04, "sigma" = c(2, 1), "decay" = c(3 / 906.8, 3 / 9.1)),
    n.report = 100001, n.burnin = 50000, n.sims = 50000, n.thin = 100)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

# Checking convergence (PSRF and ESS)
v <- rep(NA, 2)
# beta
v[1] <- coda::gelman.diag(coda::mcmc.list(
  coda::mcmc(M4[,1]$params$beta),
  coda::mcmc(M4[,2]$params$beta)))$mpsrf
v[2] <- min(
  coda::effectiveSize(M4[,1]$params$beta) +
    coda::effectiveSize(M4[,2]$params$beta))
v
# betal1
v[1] <- coda::gelman.diag(coda::mcmc.list(
  coda::mcmc(M4[,1]$params$betal1),
  coda::mcmc(M4[,2]$params$betal1)))$mpsrf
v[2] <- min(
  coda::effectiveSize(M4[,1]$params$betal1) +
    coda::effectiveSize(M4[,2]$params$betal1))
v
# betal2
v[1] <- coda::gelman.diag(coda::mcmc.list(
  coda::mcmc(M4[,1]$params$betal2),
  coda::mcmc(M4[,2]$params$betal2)))$mpsrf
v[2] <- min(
  coda::effectiveSize(M4[,1]$params$betal2) +
    coda::effectiveSize(M4[,2]$params$betal2))
v
# ws
v[1] <- coda::gelman.diag(coda::mcmc.list(
  coda::mcmc(M4[,1]$params$ws),
  coda::mcmc(M4[,2]$params$ws)))$mpsrf
v[2] <- min(
  coda::effectiveSize(M4[,1]$params$ws) +
    coda::effectiveSize(M4[,2]$params$ws))
v
# wtl
v[1] <- coda::gelman.diag(coda::mcmc.list(
  coda::mcmc(M4[,1]$params$wtl),
  coda::mcmc(M4[,2]$params$wtl)))$mpsrf
v[2] <- min(
  coda::effectiveSize(M4[,1]$params$wtl) +
    coda::effectiveSize(M4[,2]$params$wtl))
v
# hp
v[1] <- coda::gelman.diag(coda::mcmc.list(
  coda::mcmc(M4[,1]$params$hp),
  coda::mcmc(M4[,2]$params$hp)))$mpsrf
v[2] <- min(
  coda::effectiveSize(M4[,1]$params$hp) +
    coda::effectiveSize(M4[,2]$params$hp))
v

# Save data
save(M1, M2, M3, M4, file = "data/model_others.RData")
