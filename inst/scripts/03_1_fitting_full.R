##############################################
### Section 4 - MODEL FITTING (FULL MODEL) ###
##############################################

# input  data: coords, data
# output data: "data/model.rds"

# 2 Chains in parallel: ~ 75 hours 
time <- Sys.time()
cl <- parallel::makeCluster(2)
parallel::clusterExport(cl, c("data", "coords"))
model <- parallel::parSapply(cl = cl, X = 1:2, 
  FUN = function(iter) {set.seed(123 * iter + 1234); sprom::rom(
    formula = Y ~ (sine + cosi + log(dist)) * poly(trend, 2) + (trend + log(dist)) * lag1 * lag2, 
    data = data,
    coords = coords,
    site = data$site,
    year = data$year,
    day  = data$day,
    model = "5",
    n.report = 100001, n.burnin = 100000, n.sims = 100000, n.thin = 200)}
)
parallel::stopCluster(cl = cl)
Sys.time() - time

saveRDS(model, file = "data/model.rds")



# CHECKING CONVERGENCE (PSRF and ESS)
v <- rep(NA, 2)
## beta
v[1] <- coda::gelman.diag(coda::mcmc.list(
  coda::mcmc(model[,1]$params$beta),
  coda::mcmc(model[,2]$params$beta)))$mpsrf
v[2] <- min(
  coda::effectiveSize(model[,1]$params$beta) +
  coda::effectiveSize(model[,2]$params$beta))
v
## betal1
v[1] <- coda::gelman.diag(coda::mcmc.list(
  coda::mcmc(model[,1]$params$betal1),
  coda::mcmc(model[,2]$params$betal1)))$mpsrf
v[2] <- min(
  coda::effectiveSize(model[,1]$params$betal1) +
  coda::effectiveSize(model[,2]$params$betal1))
v
## betal2
v[1] <- coda::gelman.diag(coda::mcmc.list(
  coda::mcmc(model[,1]$params$betal2),
  coda::mcmc(model[,2]$params$betal2)))$mpsrf
v[2] <- min(
  coda::effectiveSize(model[,1]$params$betal2) +
  coda::effectiveSize(model[,2]$params$betal2))
v
## hp
v[1] <- coda::gelman.diag(coda::mcmc.list(
  coda::mcmc(model[,1]$params$hp),
  coda::mcmc(model[,2]$params$hp)))$mpsrf
v[2] <- min(
  coda::effectiveSize(model[,1]$params$hp) +
  coda::effectiveSize(model[,2]$params$hp))
v
## wtl
v <- c(1, 10000)
for (i in 1:22265) {
  v[1] <- max(v[1], coda::gelman.diag(coda::mcmc.list(
    coda::mcmc(model[,1]$params$wtl[,i]),
    coda::mcmc(model[,2]$params$wtl[,i])))$psrf[1])
  v[2] <- min(v[2],
    coda::effectiveSize(model[,1]$params$wtl[,i]) +
    coda::effectiveSize(model[,2]$params$wtl[,i]))
}
v
## wtls
v <- c(1, 10000)
for (i in 1:890600) {
  v[1] <- max(v[1], coda::gelman.diag(coda::mcmc.list(
    coda::mcmc(model[,1]$params$wtls[,i]),
    coda::mcmc(model[,2]$params$wtls[,i])))$psrf[1])
  v[2] <- min(v[2],
    coda::effectiveSize(model[,1]$params$wtls[,i]) +
    coda::effectiveSize(model[,2]$params$wtls[,i]))
}
v

# CHECKING CONVERGENCE (trace-plots)
## beta
name <- c("sin", "cos", "log(dist)", "trend1", "trend2", "lag1", "lag2",
          "sin:trend1", "sin:trend2", "cos:trend1", "cos:trend2", 
          "log(dist):trend1", "log(dist):trend2", "log(t-1):lag1", 
          "log(dist):lag1", "log(t-1):lag2", "log(dist):lag2", "lag1:lag2",
          "log(t-1):lag1:lag2", "log(dist):lag1:lag2")   
for (i in 1:20) {
  pdf(paste0("inst/img/SUPP_traceplot", gsub(':', '', colnames(model[,1]$params$beta)[i]), ".pdf"),
      width = 6, height = 3)
  y_min <- min(min(model[,1]$params$beta[,i]), min(model[,2]$params$beta[,i]))
  y_max <- max(max(model[,1]$params$beta[,i]), max(model[,2]$params$beta[,i]))
  plot(model[,1]$params$beta[,i], type = "l", 
       xlab = "Iteration", ylab = "Parameter", ylim = c(y_min, y_max),
       main = name[i])
  lines(model[,2]$params$beta[,i], col = "gray")
  dev.off()
}
## betal1
name <- c("trend1 (\u2113 = 1)", "lag1 (\u2113 = 1)")
for (i in 1:2) {
  Cairo::CairoPDF(paste0("inst/img/SUPP_traceplot", gsub(':', '', colnames(model[,1]$params$betal1)[i]), "l1.pdf"),
      width = 6, height = 3)
  y_min <- min(min(model[,1]$params$betal1[,i]), min(model[,2]$params$betal1[,i]))
  y_max <- max(max(model[,1]$params$betal1[,i]), max(model[,2]$params$betal1[,i]))
  plot(model[,1]$params$betal1[,i], type = "l", 
       xlab = "Iteration", ylab = "Parameter", ylim = c(y_min, y_max),
       main = name[i])
  lines(model[,2]$params$betal1[,i], col = "gray")
  dev.off()
}
## betal2
name <- c("trend1 (\u2113 = 2)", "lag1 (\u2113 = 2)")
for (i in 1:2) {
  Cairo::CairoPDF(paste0("inst/img/SUPP_traceplot", gsub(':', '', colnames(model[,1]$params$betal2)[i]), "l2.pdf"),
      width = 6, height = 3)
  y_min <- min(min(model[,1]$params$betal2[,i]), min(model[,2]$params$betal2[,i]))
  y_max <- max(max(model[,1]$params$betal2[,i]), max(model[,2]$params$betal2[,i]))
  plot(model[,1]$params$betal2[,i], type = "l", 
       xlab = "Iteration", ylab = "Parameter", ylim = c(y_min, y_max),
       main = name[i])
  lines(model[,2]$params$betal2[,i], col = "gray")
  dev.off()
}
## hp (hyper parameters)
name <- c("(Intercept)", expression(1 / sigma[0]^2), expression(phi[0]), expression(1 / sigma[1]^2),
          "(Intercept) (\u2113 = 1)", expression(1 / sigma["0,1"]^2), expression(1 / sigma["1,1"]^2),
          "(Intercept) (\u2113 = 2)", expression(1 / sigma["0,2"]^2), expression(1 / sigma["1,2"]^2))
for (i in 1:10) {
  Cairo::CairoPDF(paste0("inst/img/SUPP_traceplot", gsub(':', '', colnames(model[,1]$params$hp)[i]), ".pdf"),
      width = 6, height = 3)
  y_min <- min(min(model[,1]$params$hp[,i]), min(model[,2]$params$hp[,i]))
  y_max <- max(max(model[,1]$params$hp[,i]), max(model[,2]$params$hp[,i]))
  plot(model[,1]$params$hp[,i], type = "l", 
       xlab = "Iteration", ylab = "Parameter", ylim = c(y_min, y_max),
       main = name[i])
  lines(model[,2]$params$hp[,i], col = "gray")
  dev.off()
}
