##########################################
### Section 6 (SUPPLEMENTARY MATERIAL) ###
### SIMULATION STUDY                   ###
##########################################

# Additional R packages: "RecordTest", "spTReg", "truncnorm"
# spTReg (GitHub: https://github.com/JorgeCastilloMateo/spTReg)

# Clear workspace
rm(list = ls())


# Setting
TT <- 62
LL <- 365
cc <- 0.035
tt <- rep(1:TT, LL)
sigma <- 3.56


### Normal (A)
set.seed(543108)
p.true <- RecordTest::p.record(matrix(cc * (1:TT) + sigma * rnorm(TT * 1e+06), nrow = TT, ncol = 1e+06))[TT]
p.hat1 <- p.hat2 <- rep(NA, 100)
for (i in 1:100) {

  # simulation
  epsilon <- rnorm(TT * LL)
  Y <- cc * tt + sigma * epsilon

  # approach 1
  I <- RecordTest::I.record(matrix(Y, nrow = TT, ncol = LL))

  glm.fit <- glm(c(I[-1,]) ~ poly(log(rep(2:TT, LL) - 1), 2), family = binomial)

  p.hat1[i] <- glm.fit$fitted.values[TT - 1]

  # approach 2
  lm.fit <- spTReg::iidm(Y ~ I(tt) + I(tt^2), verbose = FALSE)

  b0.hat <- lm.fit$p.params.samples[,1]
  b1.hat <- lm.fit$p.params.samples[,2]
  b2.hat <- lm.fit$p.params.samples[,3]
  sd.hat <- lm.fit$p.params.samples[,4]

  p.hat <- rep(NA, 1000)
  for (b in 1:1000) {
    Y.hat <- matrix(rnorm(TT * LL, mean = b0.hat[b] + b1.hat[b] * tt + b2.hat[b] * tt^2, sd = sd.hat[b]))
    p.hat[b] <- RecordTest::p.record(matrix(Y.hat, nrow = TT, ncol = LL))[TT]
  }

  p.hat2[i] <- mean(p.hat)
}

metrics_sims(p.true, p.hat1, p.hat2)



### Student's t (B)
set.seed(534259)
for (nu in c(1, 2, 3, 5, 10, 20, 50)) {
  print(paste("Value of nu:", nu))
  p.true <- RecordTest::p.record(matrix(cc * (1:TT) + sigma * rt(TT * 1e+06, df = nu), nrow = TT, ncol = 1e+06))[TT]
  p.hat1 <- p.hat2 <- rep(NA, 100)
  for (i in 1:100) {

    # simulation
    epsilon <- rt(TT * LL, df = nu)
    Y <- cc * tt + sigma * epsilon

    # approach 1
    I <- RecordTest::I.record(matrix(Y, nrow = TT, ncol = LL))

    glm.fit <- glm(c(I[-1,]) ~ poly(log(rep(2:TT, LL) - 1), 2), family = binomial)

    p.hat1[i] <- glm.fit$fitted.values[TT - 1]

    # approach 2
    lm.fit <- spTReg::iidm(Y ~ I(tt) + I(tt^2), verbose = FALSE)

    b0.hat <- lm.fit$p.params.samples[,1]
    b1.hat <- lm.fit$p.params.samples[,2]
    b2.hat <- lm.fit$p.params.samples[,3]
    sd.hat <- lm.fit$p.params.samples[,4]

    p.hat <- rep(NA, 1000)
    for (b in 1:1000) {
      Y.hat <- matrix(rnorm(TT * LL, mean = b0.hat[b] + b1.hat[b] * tt + b2.hat[b] * tt^2, sd = sd.hat[b]))
      p.hat[b] <- RecordTest::p.record(matrix(Y.hat, nrow = TT, ncol = LL))[TT]
    }

    p.hat2[i] <- mean(p.hat)
  }

  print(metrics_sims(p.true, p.hat1, p.hat2))
}



### Skew normal (C)
set.seed(938371)
for (alpha in c(-4, -2, -1, 1, 2, 4)) {
  print(paste("Value of alpha:", alpha))
  p.true <- RecordTest::p.record(
    matrix(cc * (1:TT) + sigma * (alpha * abs(rnorm(TT * 1e+06)) + rnorm(TT * 1e+06)) / sqrt(1 + alpha^2),
           nrow = TT, ncol = 1e+06))[TT]
  p.hat1 <- p.hat2 <- rep(NA, 100)
  for (i in 1:100) {

    # simulation
    epsilon <- (alpha * abs(rnorm(TT * LL)) + rnorm(TT * LL)) / sqrt(1 + alpha^2)
    Y <- cc * tt + sigma * epsilon

    # approach 1
    I <- RecordTest::I.record(matrix(Y, nrow = TT, ncol = LL))

    glm.fit <- glm(c(I[-1,]) ~ poly(log(rep(2:TT, LL) - 1), 2), family = binomial)

    p.hat1[i] <- glm.fit$fitted.values[TT - 1]

    # approach 2
    lm.fit <- spTReg::iidm(Y ~ I(tt) + I(tt^2), verbose = FALSE)

    b0.hat <- lm.fit$p.params.samples[,1]
    b1.hat <- lm.fit$p.params.samples[,2]
    b2.hat <- lm.fit$p.params.samples[,3]
    sd.hat <- lm.fit$p.params.samples[,4]

    p.hat <- rep(NA, 1000)
    for (b in 1:1000) {
      Y.hat <- matrix(rnorm(TT * LL, mean = b0.hat[b] + b1.hat[b] * tt + b2.hat[b] * tt^2, sd = sd.hat[b]))
      p.hat[b] <- RecordTest::p.record(matrix(Y.hat, nrow = TT, ncol = LL))[TT]
    }

    p.hat2[i] <- mean(p.hat)
  }

  print(metrics_sims(p.true, p.hat1, p.hat2))
}



### Normal bulk and truncated normal tail (D)
rmixture <- function(n, tau0 = 0.95, mean1 = 0, mean2 = 0, sd1 = 1, sd2 = 1) {

  q_tau0 <- sd1 * qnorm(tau)

  fx <- truncnorm::dtruncnorm(q_tau0, b = q_tau0, mean = mean1, sd = sd1)
  gx <- truncnorm::dtruncnorm(q_tau0, a = q_tau0, mean = mean2, sd = sd2)

  U <- runif(n)

  tau <- gx / (fx + gx)

  sel <- (U < tau)

  U[sel]  <- truncnorm::qtruncnorm((U / tau)[sel],
    b = q_tau0, mean = mean1[sel], sd = sd1)
  U[!sel] <- truncnorm::qtruncnorm(((U - tau) / (1 - tau))[!sel],
    a = q_tau0, mean = mean2[!sel], sd = sd2)

  return(U)
}

set.seed(938371)
for (tau in c(0.90, 0.95, 0.99)) {
  print(paste("Value of tau:", tau))
  for (cc_ in c(0, 0.5 * cc, 0.75 * cc, cc, 1.5 * cc, 2 * cc)) {
    print(paste("Value of cc_:", cc_))

    p.true <- RecordTest::p.record(
      matrix(rmixture(TT * 1e+06, tau, rep(cc * (1:TT), 1e+06), rep(cc_ * (1:TT), 1e+06), sigma, sigma),
             nrow = TT, ncol = 1e+06))[TT]
    p.hat1 <- p.hat2 <- rep(NA, 100)
    for (i in 1:100) {

      # simulation
      Y <- rmixture(TT * LL, tau, cc * tt, cc_ * tt, sigma, sigma)

      # approach 1
      I <- RecordTest::I.record(matrix(Y, nrow = TT, ncol = LL))

      glm.fit <- glm(c(I[-1,]) ~ poly(log(rep(2:TT, LL) - 1), 2), family = binomial)

      p.hat1[i] <- glm.fit$fitted.values[TT - 1]

      # approach 2
      lm.fit <- spTReg::iidm(Y ~ I(tt) + I(tt^2), verbose = FALSE)

      b0.hat <- lm.fit$p.params.samples[,1]
      b1.hat <- lm.fit$p.params.samples[,2]
      b2.hat <- lm.fit$p.params.samples[,3]
      sd.hat <- lm.fit$p.params.samples[,4]

      p.hat <- rep(NA, 1000)
      for (b in 1:1000) {
        Y.hat <- matrix(rnorm(TT * LL, mean = b0.hat[b] + b1.hat[b] * tt + b2.hat[b] * tt^2, sd = sd.hat[b]))
        p.hat[b] <- RecordTest::p.record(matrix(Y.hat, nrow = TT, ncol = LL))[TT]
      }

      p.hat2[i] <- mean(p.hat)
    }

    print(metrics_sims(p.true, p.hat1, p.hat2))

  }
}

