#################################################
### Section 2 - MISSING DATA and WEAK RECORDS ###
#################################################

# input  data: "data/Y365.rda"

# Clear workspace
rm(list = setdiff(ls(), c("background", "grid", "limits", "stations", "Y365")))

# Read Tmax data. First 01_1_download_Tmax_data
#load("data/Y365.rda")



# MISSING DATA
table(apply(Y365, 3, function(x) sum(is.na(x))))

# Simulation missing data
TT <- nrow(Y365)
tt <- 1:TT
hatsd <- mean(apply(Y365, 2:3, function(x) summary(lm(x ~ tt))$sigma)) # average sd for every day within year and site
hatb1 <- mean(apply(Y365, 2:3, function(x) lm(x ~ tt)$coef[2])) # average linear trend across years
myNA <- is.na(Y365) # position of NA's in data

set.seed(23)
difference <- matrix(nrow = 10000, ncol = 4)
for (b in 1:10000) {
  cat(paste0("..", b))
  X <- matrix(rnorm(905200, mean = tt * hatb1, sd = hatsd), 
              nrow = TT, ncol = 14600)
  I1 <- apply(X, 2, sprom::I.weak.record)
  X[myNA] <- NA
  I2 <- apply(X, 2, sprom::I.weak.record)
  difference[b, 1] <- sum(I1 - I2)
  difference[b, 2] <- sum(I1 != I2)
  difference[b, 3] <- sum(I1[32:62,] - I2[32:62,])
  difference[b, 4] <- sum(I1[32:62,] != I2[32:62,])
}

mean(difference[,1])
quantile(difference[,1], c(0.05, 0.95))
mean(difference[,2])
quantile(difference[,2], c(0.05, 0.95))
mean(difference[,3])
quantile(difference[,3], c(0.05, 0.95))
mean(difference[,4])
quantile(difference[,4], c(0.05, 0.95))

par(mfrow = c(1, 2))
hist(difference[, 1], freq = FALSE)
hist(difference[, 3], freq = FALSE)



# WEAK RECORDS
# Compute indicators
TT <- nrow(Y365)
LL <- ncol(Y365)
SS <- dim(Y365)[3]
I365 <- array(dim = c(TT, LL, SS))
for (i in 1:SS) {
  I365[,,i] <- apply(Y365[,,i], 2, sprom::I.weak.record)
}

X <- apply(I365[-1,,], 3, function(x) 100 * sum(!(x %in% 0:1)) / sum(x != 0))
mean(X)
range(X)

sum(I365[-1,,] == 1)
sum(I365[-1,,] == 1 / 2)
sum(I365[-1,,] == 1 / 3)
sum(I365[-1,,] == 1 / 4)
sum(I365[-1,,] == 1 / 5)
sum(I365[-1,,] == 1 / 6)

sum(!(I365[-1,,] %in% 0:1)) / sum(I365[-1,,] != 0)
sum(!(I365[3:12,,] %in% 0:1)) / sum(I365[3:12,,] != 0)
sum(!(I365[13:22,,] %in% 0:1)) / sum(I365[13:22,,] != 0)
sum(!(I365[23:32,,] %in% 0:1)) / sum(I365[23:32,,] != 0)
sum(!(I365[33:42,,] %in% 0:1)) / sum(I365[33:42,,] != 0)
sum(!(I365[43:52,,] %in% 0:1)) / sum(I365[43:52,,] != 0)
sum(!(I365[53:62,,] %in% 0:1)) / sum(I365[53:62,,] != 0)

# Simulating weak records
set.seed(23)
Y1 <- rnorm(1000000, mean = 30 * hatb1, sd = hatsd)
Y2 <- rnorm(1000000, mean = 40 * hatb1, sd = hatsd)
e1 <- Y1 - round(Y1, 1)
e2 <- Y2 - round(Y2, 1)

ks.test(e1, e2)

par(mfrow = c(1, 2))
hist(e1[round(Y1, 1) == 6], freq = FALSE)
hist(e2[round(Y2, 1) == 6], freq = FALSE)
