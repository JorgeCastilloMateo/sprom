###################################
### Section 4 - DATA PROCESSING ###
###################################

# input  data: (updated)stations, "data/Y365.rda"
# output data: I365, data, coords

# Clear workspace
rm(list = setdiff(ls(), c("background", "grid", "limits", "stations", "Y365")))

# Load the R-package
library("sprom")

# Read Tmax data. First 01_1_download_Tmax_data
#load("data/Y365.rda")

# Indexes
TT <- nrow(Y365)
LL <- ncol(Y365)
SS <- dim(Y365)[3]
tt <- 2:TT
ll <- 1:LL

# Compute indicators
I365 <- array(dim = c(TT, LL, SS))
for (i in 1:SS) {
  I365[,,i] <- apply(Y365[,,i], 2, I.weak.record)
}

# Main effects / indexes
data <- data.frame(
  Y       = c(I365[-1, , ]),
  year    = c(array(tt, dim = c(TT - 1, LL, SS))),
  day     = c(aperm(array(ll, dim = c(LL, TT - 1, SS)), perm = c(2, 1, 3))),
  yearday = NA,
  site    = c(aperm(array(1:SS, dim = c(SS, TT - 1, LL)), perm = c(2, 3, 1))),
  trend   = c(array(log(tt - 1), dim = c(TT - 1, LL, SS))),
  sine    = c(aperm(array(sin(2 * pi * ll / 365), dim = c(LL, TT - 1, SS)), perm = c(2, 1, 3))),
  cosi    = c(aperm(array(cos(2 * pi * ll / 365), dim = c(LL, TT - 1, SS)), perm = c(2, 1, 3))),
  lon     = c(aperm(array(stations$LON, dim = c(SS, TT - 1, LL)), perm = c(2, 3, 1))),
  lat     = c(aperm(array(stations$LAT, dim = c(SS, TT - 1, LL)), perm = c(2, 3, 1))),
  dist    = c(aperm(array(stations$DIST, dim = c(SS, TT - 1, LL)), perm = c(2, 3, 1))),
  elev    = c(aperm(array(stations$HGHT, dim = c(SS, TT - 1, LL)), perm = c(2, 3, 1))),
  lag1    = c(abind::abind(I365[-TT, 365, , drop = FALSE], I365[-1, -365, ], along = 2)),
  lag2    = c(abind::abind(I365[-TT, 364:365, ], I365[-1, -c(364:365), ], along = 2))
)

data$yearday <- interaction(data$year, data$day)

coords <- sf::st_coordinates(stations$geometry) / 1000
