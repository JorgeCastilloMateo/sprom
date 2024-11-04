##################################################
### Section 4 - RESULTS FULL MODEL (INFERENCE) ###
##################################################

# input  data: grid, "data/model.rds"

core <- function(b, T, n, newn,
                 lList, tList, DJF, MAM, JJA,
                 sigma1, mu1, sigma2, mu2, sigma, mu,
                 correct1, correct2, correct,
                 beta1, beta2, beta, decay,
                 Wtl, wtls, prec0l1, prec0l2, prec0,
                 d, trendOrtho, lag1, lag2, trend, sine, cosi,
                 newdata, Terms, namesBeta = colnames(object$x)) {
  # Statistics of interest:
  # prob of interest
  P <- matrix(0, nrow = newn, ncol = length(lList))
  # N total (1), N 41:62 - 1 by season (4), N 53:62 - 1 by season (4)
  statisticsS <- matrix(0, nrow = newn, ncol = 9)
  # ERS total and by season (5)
  statisticsT <- matrix(0, nrow = T, ncol = 5)
  
  ### matrices
  # l = 1
  beta1 <- beta1 / sigma1
  correct1 <- - sum(beta1 * mu1)
  # l = 2
  beta2 <- beta2 / sigma2
  correct2 <- - sum(beta2 * mu2)
  # l = 3,...,L
  beta <- beta / sigma
  correct <- - sum(beta * mu)
  R22 <- exp(- decay * d[newn + 1:n, newn + 1:n])
  R22inv <- solve(R22)
  R11 <- exp(- decay * d[1:newn, 1:newn])
  R12 <- exp(- decay * d[1:newn, newn + 1:n])
  R12R22inv <- R12 %*% R22inv
  R <- t(chol(R11 - R12R22inv %*% t(R12))) #lower triangular
  lListCount <- 0
  
  for (t in 1:T){
    # l = 1
    wtl <- Wtl[t]
    w <- 
      wtl + R12R22inv %*% (wtls[t + 0:(n - 1) * (365 * T)] - wtl) + 
      R %*% rnorm(newn) / sqrt(prec0l1)
    
    p <- 1 / (1 + exp(- (cbind(trendOrtho[t,1], lag1) %*% beta1 + w + correct1)))
    lag2 <- rbinom(newn, size = 1, prob = p)
    
    # N total (1), N 41:62 - 1 by season (4), N 53:62 - 1 by season (4)
    statisticsS[,1] <- statisticsS[,1] + lag2 # total
    if (t %in% (41:62 - 1)) {
      statisticsS[,2] <- statisticsS[,2] + lag2 # Winter
    }
    if (t %in% (53:62 - 1)) {
      statisticsS[,6] <- statisticsS[,6] + lag2 # Winter
    }
    # ERS total and by season (5)
    statisticsT[t,1] <- statisticsT[t,1] + mean(lag2) # total
    statisticsT[t,2] <- statisticsT[t,2] + mean(lag2) # Winter

    # l = 2
    wtl <- Wtl[t + T]
    w <- 
      wtl + R12R22inv %*% (wtls[t + T + 0:(n - 1) * (365 * T)] - wtl) + 
      R %*% rnorm(newn) / sqrt(prec0l2)
    
    p <- 1 / (1 + exp(- (cbind(trendOrtho[t,1], lag2) %*% beta2 + w + correct2)))
    lag1 <- rbinom(newn, size = 1, prob = p)
    
    # N total (1), N 41:62 - 1 by season (4), N 53:62 - 1 by season (4)
    statisticsS[,1] <- statisticsS[,1] + lag1 # total
    if (t %in% (41:62 - 1)) {
      statisticsS[,2] <- statisticsS[,2] + lag1 # Winter
    }
    if (t %in% (53:62 - 1)) {
      statisticsS[,6] <- statisticsS[,6] + lag1 # Winter
    }
    # ERS total and by season (5)
    statisticsT[t,1] <- statisticsT[t,1] + mean(lag1) # total
    statisticsT[t,2] <- statisticsT[t,2] + mean(lag1) # Winter
    
    newdata$x$trend <- rep(trend[t], newn)
    
    # l = 3,...,L
    for (l in 3:365) {
      newdata$x$sine <- rep(sine[l], newn)
      newdata$x$cosi <- rep(cosi[l], newn)
      newdata$x$lag1 <- lag1
      newdata$x$lag2 <- lag2
      
      m <- model.frame(Terms, newdata$x)
      X <- model.matrix(Terms, m)[,namesBeta]
      
      wtl <- Wtl[t + (l - 1) * T]
      w <- 
        wtl + R12R22inv %*% (wtls[t + (l - 1) * T + 0:(n - 1) * (365 * T)] - wtl) + 
        R %*% rnorm(newn) / sqrt(prec0)
      
      p <- 1 / (1 + exp(- (X %*% beta + w + correct)))
      
      if ((t %in% tList) & (l %in% lList)) {
        lListCount <- lListCount + 1
        P[,lListCount] <- p
      }
      
      lag2 <- lag1
      lag1 <- rbinom(newn, size = 1, prob = p)
      
      # N total (1), N 41:62 - 1 by season (4), N 53:62 - 1 by season (4)
      statisticsS[,1] <- statisticsS[,1] + lag1 # total
      if (t %in% (41:62 - 1)) {
        if (l %in% DJF) {
          statisticsS[,2] <- statisticsS[,2] + lag1
        } else if (l %in% MAM) {
          statisticsS[,3] <- statisticsS[,3] + lag1
        } else if (l %in% JJA) {
          statisticsS[,4] <- statisticsS[,4] + lag1
        } else {
          statisticsS[,5] <- statisticsS[,5] + lag1
        }
        if (t %in% (53:62 - 1)) {
          if (l %in% DJF) {
            statisticsS[,6] <- statisticsS[,6] + lag1
          } else if (l %in% MAM) {
            statisticsS[,7] <- statisticsS[,7] + lag1
          } else if (l %in% JJA) {
            statisticsS[,8] <- statisticsS[,8] + lag1
          } else {
            statisticsS[,9] <- statisticsS[,9] + lag1
          }
        }
      }
      
      # ERS total and by season (5)
      statisticsT[t,1] <- statisticsT[t,1] + mean(lag1) # total
      if (l %in% DJF) {
        statisticsT[t,2] <- statisticsT[t,2] + mean(lag1)
      } else if (l %in% MAM) {
        statisticsT[t,3] <- statisticsT[t,3] + mean(lag1)
      } else if (l %in% JJA) {
        statisticsT[t,4] <- statisticsT[t,4] + mean(lag1)
      } else {
        statisticsT[t,5] <- statisticsT[t,5] + mean(lag1)
      }
    }
  }
  
  return(list(P, statisticsS, statisticsT))
}

# object matrix of two chains
# newdata coords and logdist of the grid
# tList is the position of the year (without counting the trivial year)
# lList is the vector of days
predictModelPaper <- function(object,
                              newdata, 
                              tList, lList) {
  
  mu     <- attr(object$x, "scaled:center")
  sigma  <- attr(object$x, "scaled:scale")
  mu1    <- attr(object$x, "scaled:center")[colnames(object$params$betal1)]
  sigma1 <- attr(object$x, "scaled:scale")[colnames(object$params$betal1)]
  mu2    <- attr(object$x, "scaled:center")[colnames(object$params$betal2)]
  sigma2 <- attr(object$x, "scaled:scale")[colnames(object$params$betal2)]
  
  tt <- terms(object)
  Terms <- delete.response(tt)
  
  # seasons
  DJF <- c(1:59, 335:365)
  MAM <- 60:151
  JJA <- 152:243
  SON <- 244:334
  
  B <- nrow(object$params$beta)
  n <- nrow(object$coords)
  T <- length(unique(object$date$year))
  newn <- nrow(newdata$coords)
  
  poly_coefs <- attr(poly(object$model$trend, degree = 2), "coefs")
  trend <- log(1:T)
  trendOrtho <- poly(trend, degree = 2, coefs = poly_coefs)
  sine <- sin(2 * pi * 1:365 / 365)
  cosi <- cos(2 * pi * 1:365 / 365)
  
  lag1 <- rep(1, newn)
  lag2 <- rep(NA, newn)
  
  distance <- dist(rbind(newdata$coords, object$coords))
  d <- matrix(0, nrow = newn + n, ncol = newn + n)
  d[lower.tri(d)] <- distance
  d <- d + t(d)
  
  newdata$x <- data.frame("trend" = rep(NA, newn),
                          "sine"  = rep(NA, newn),
                          "cosi"  = rep(NA, newn),
                          "dist"  = newdata$dist,
                          "lag1"  = lag1,
                          "lag2"  = lag2)
 
  # para dias concretos
  cl <- parallel::makeCluster(10)
  parallel::clusterExport(cl, c("core"))
  statistics <- parallel::parSapply(cl = cl, X = 1:B, 
    FUN = function(iter) {set.seed(123 * iter + 1234); core(
      b = iter, T, n, newn,
      lList, tList, DJF, MAM, JJA,
      sigma1, mu1, sigma2, mu2, sigma, mu,
      correct1, correct2, correct,
      beta1 = object$params$betal1[iter,], beta2 = object$params$betal2[iter,], 
      beta = object$params$beta[iter,], decay = object$params$hp[iter,"decay0"],
      Wtl = object$params$wtl[iter,], wtls = object$params$wtls[iter,], 
      prec0l1 = object$params$hp[iter,"prec0l1"], 
      prec0l2 = object$params$hp[iter,"prec0l2"], 
      prec0 = object$params$hp[iter,"prec0"],
      d, trendOrtho, lag1, lag2, trend, sine, cosi,
      newdata, Terms, namesBeta = colnames(object$x))}
  )
  parallel::stopCluster(cl = cl)
  
  # Statistics of interest:
  # prob of interest
  P <- array(dim = c(newn, length(lList), B))
  # N total (1), N 41:62 - 1 by season (4), N 53:62 - 1 by season (4)
  statisticsS <- array(0, dim = c(newn, 9, B))
  # ERS total and by season (5)
  statisticsT <- array(0, dim = c(T, 5, B))
  
  for (b in 1:B) {
    P[,,b] <- statistics[,b][[1]]
    statisticsS[,,b] <- statistics[,b][[2]]
    statisticsT[,,b] <- statistics[,b][[3]]
  }
  
  return(list(P, statisticsS, statisticsT))
}

newdata <- list()
newdata$dist <- grid$dist
newdata$coords <- sf::st_coordinates(grid) / 1000

models <- model[,1]
models$params$beta   <- rbind(model[,1]$params$beta,   model[,2]$params$beta)
models$params$wtls   <- rbind(model[,1]$params$wtls,   model[,2]$params$wtls)
models$params$wtl    <- rbind(model[,1]$params$wtl,    model[,2]$params$wtl)
models$params$wt     <- rbind(model[,1]$params$wt,     model[,2]$params$wt)
models$params$hp     <- rbind(model[,1]$params$hp,     model[,2]$params$hp)
models$params$betal1 <- rbind(model[,1]$params$betal1, model[,2]$params$betal1)
models$params$betal2 <- rbind(model[,1]$params$betal2, model[,2]$params$betal2)
rm(model)

# 10 clusters in parallel: < 75 hours 
system.time(
  pred <- predictModelPaper(models,
                            newdata, 
                            61, 222:231)
)

saveRDS(pred, file = "data/pred.rds")




library("sf")
library("sp")
library("stars")
library("smoothr")
library("tidyverse")
library("geomtextpath")
library("rnaturalearth")
library("rnaturalearthdata")

#' @param Z The values
#' @param coords_limit Boundaries of the map (typical xlim and ylim, better to
#'   leave them by default)
#' @param ref number: reference white color in the blue to red continuous color 
#'   scale (not factor)
#' @param zlim limits of the Z if continuous (not factor)
#' @param picture.name name of the .png to be saved
#' @param title map title
#' @param legend.name legend title
#' @param save saves the map in the directory
#' @param contour one of c("none", "auto", "smooth")
#' @param breaks the contour levels will be drawn at these values * ref
#' @param label label for each break
#' @param smoothness level of smoothness if contour = "smooth"
#' @param threshold minimum length of the contour line to be drawn in meters
#' @param threshold2 if the length of the line in meters is smaller, the label
#'   is not drawn
#' @param dist vector of distance to the coast (if < 2 km the pixel is not 
#'   considered in the contour line)
#' @return A map

mapSpain <- function(Z, 
                     coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
                     ref = .5,  
                     zlim = NULL, 
                     picture.name = "photo.png",
                     title = expression(p["tl"](s)), 
                     legend.name = "",
                     save = TRUE,
                     contour = c("none", "auto", "smooth"),
                     breaks = 8:12 / 10,
                     label = "",
                     n.contour = 1,
                     smoothness = 5,
                     threshold = 1e+05,
                     threshold2 = 4 * threshold,
                     dist) {
  
  contour <- match.arg(contour)
  
  coords_limit <- st_transform(
    as(
      SpatialPointsDataFrame(coords = coords_limit, 
                             data = coords_limit,
                             proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")),
      'sf'), 
    2062)
  
  #define grid of Spain
  spain <- ne_countries(scale = "large", country = "Spain", returnclass = "sf")
  spain <- st_transform(spain, 2062)
  
  spain_coords <- Polygons(
    list(Polygon(st_coordinates(spain)[st_coordinates(spain)[,"L2"] == 3,1:2])),
    ID = "spain")
  spain_coords <- SpatialPolygons(list(spain_coords))
  spain_coords <- as(spain_coords, "sf")
  st_crs(spain_coords) <- st_crs(spain)
  
  grid <- st_make_grid(spain, cellsize = 10 * 1000, what = "centers")
  grid <- st_intersection(grid, spain_coords)
  
  # background
  world_map <- ne_countries(scale = "large", returnclass = 'sf')
  european_union <- c("Andorra","France","Morocco","Portugal","Spain","Gibraltar","Algeria")
  background <- 
    world_map %>% 
    filter(name %in% european_union) %>%
    st_transform(2062)
  
  map <- ggplot(data = background) + 
    geom_sf(fill="antiquewhite") + 
    xlab("Longitude") + ylab("Latitude") + ggtitle(title) +
    theme(panel.background = element_rect(fill = "aliceblue"),
          axis.text.x=element_text(size=6),
          axis.text.y=element_text(size=6,angle=90),
          axis.title=element_text(size=10,face="bold")) + 
    geom_tile(data = grid, ggplot2::aes(x = st_coordinates(grid)[,1], y = st_coordinates(grid)[,2], fill = Z)) + 
    scale_fill_gradient2(midpoint = ref, low = scales::muted("blue"), mid = "white", high = scales::muted("red"), space = "Lab", limits = zlim, name = legend.name)
  
  if (contour == "auto") {
    map <- map + 
      geom_textcontour(data = grid, 
                       mapping = ggplot2::aes(x = sf::st_coordinates(grid)[,1], y = sf::st_coordinates(grid)[,2], z = Z),
                       breaks = ref * breaks)
  } else if (contour == "smooth") {
    raster <- st_rasterize(st_sf(geometry = grid[dist > 2], value = Z[dist > 2]))
    raster$value[raster$value == 0] <- NA
    contour1 <- stars::st_contour(stars::st_as_stars(raster), contour_lines = TRUE, breaks = ref * breaks)
    contour2 <- smoothr::smooth(contour1, method = "ksmooth", smoothness = smoothness)
    units(threshold) <- units(threshold2) <- "m"
    contour3 <- contour2[st_length(contour2) > threshold,]
    label <- label[match(round(contour3$value, 3), round(ref * breaks, 3))]
    label[st_length(contour3) < threshold2] <- ""
    
    map <- map + 
      geom_textsf(data = contour3, label = label, 
                  color = "black", linecolor = "black",
                  size = 2, linewidth = 0.275)
  }
  
  map <- map + coord_sf(xlim = st_coordinates(coords_limit)[,1], ylim = st_coordinates(coords_limit)[,2])# +
  #ggplot2::geom_point(ggplot2::aes(x = X, y = Y), data = data.frame(coords * 1000), color = "black")
  
  if (save) {
    ggplot2::ggsave(picture.name, map, width = 8.27 / 2, height = 11.69 / 4)
  }
  
  map
}



# map number of records
mapSpain(rowMeans(pred[[2]][,1,]) / 365 + 1, 
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
         ref = sum(1 / 1:62),  
         zlim = c(4.55, 6), 
         picture.name = "inst/img/MAIN_Nt62s.png",
         title = expression(paste(bar(N)["62"](s), " (mean)")), 
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 1:4 / 20,
         label = paste0(10 * 1:4 / 2, "%"),
         smoothness = 5,
         threshold = 1.5 * 1e+05,
         threshold2 = 1.5 * 1e+05,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,1,], 1, quantile, prob = 0.05) / 365 + 1, 
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
         ref = sum(1 / 1:62),  
         zlim = c(4.55, 6), 
         picture.name = "inst/img/MAIN_Nt62s-q05.png",
         title = expression(paste(bar(N)["62"](s), " (0.05 quantile)")), 
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 0:3 / 20,
         label =  c("CRM", paste0(10 * 1:3 / 2, "%")),
         smoothness = 5,
         threshold = 1.5 * 1e+05,
         threshold2 = 1.5 * 1e+05,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,1,], 1, quantile, prob = 0.95) / 365 + 1, 
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
         ref = sum(1 / 1:62),  
         zlim = c(4.55, 6), 
         picture.name = "inst/img/MAIN_Nt62s-q95.png",
         title = expression(paste(bar(N)["62"](s), " (0.95 quantile)")), 
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 1:5 / 20,
         label = paste0(10 * 1:5 / 2, "%"),
         smoothness = 5,
         threshold = 1.5 * 1e+05,
         threshold2 = 1.5 * 1e+05,
         dist = grid$dist)

summary(rowMeans(pred[[2]][,1,]) / 365 + 1, digits = 3)

interest <- colMeans((pred[[2]][,1,] / 365 + 1) > sum(1 / 1:62))
mean(interest)
quantile(interest, prob = c(0.05, 0.95))

mean(apply(pred[[2]][,1,] / 365 + 1, 1, quantile, prob = 0.05) > sum(1 / 1:62))
mean(apply(pred[[2]][,1,] / 365 + 1, 1, quantile, prob = 0.95) < sum(1 / 1:62))



# map number of records 21st century
mapSpain(rowMeans(pred[[2]][,2,]) / 90, 
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
         ref = sum(1 / 41:62),  
         zlim = c(0, 1.3), 
         picture.name = "inst/img/MAIN_Nt41-62lDJFs.png",
         title = expression(paste(bar(N)["41:62,DJF"](s), " (mean)")), 
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 3:7 / 10,
         label = paste0(10 * 3:7, "%"),
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,2,], 1, quantile, prob = 0.05) / 90, 
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
         ref = sum(1 / 41:62),  
         zlim = c(0, 1.3), 
         picture.name = "inst/img/MAIN_Nt41-62lDJFs-q05.png",
         title = expression(paste(bar(N)["41:62,DJF"](s), " (0.05 quantile)")), 
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 0:4 / 10,
         label = c("CRM", paste0(10 * 1:4, "%")),
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,2,], 1, quantile, prob = 0.95) / 90, 
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
         ref = sum(1 / 41:62),  
         zlim = c(0, 1.3), 
         picture.name = "inst/img/MAIN_Nt41-62lDJFs-q95.png",
         title = expression(paste(bar(N)["41:62,DJF"](s), " (0.95 quantile)")), 
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 6:10 / 10,
         label = paste0(10 * 6:10, "%"),
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)

mapSpain(rowMeans(pred[[2]][,3,]) / 92, 
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
         ref = sum(1 / 41:62),  
         zlim = c(0, 1.3), 
         picture.name = "inst/img/MAIN_Nt41-62lMAMs.png",
         title = expression(paste(bar(N)["41:62,MAM"](s), " (mean)")), 
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 2:9 / 10,
         label = paste0(10 * 2:9, "%"),
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,3,], 1, quantile, prob = 0.05) / 92, 
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
         ref = sum(1 / 41:62),  
         zlim = c(0, 1.3), 
         picture.name = "inst/img/MAIN_Nt41-62lMAMs-q05.png",
         title = expression(paste(bar(N)["41:62,MAM"](s), " (0.05 quantile)")), 
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 0:6 / 10,
         label = c("CRM", paste0(10 * 1:6, "%")),
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,3,], 1, quantile, prob = 0.95) / 92, 
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
         ref = sum(1 / 41:62),  
         zlim = c(0, 1.3), 
         picture.name = "inst/img/MAIN_Nt41-62lMAMs-q95.png",
         title = expression(paste(bar(N)["41:62,MAM"](s), " (0.95 quantile)")), 
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 5:12 / 10,
         label = paste0(10 * 5:12, "%"),
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)

mapSpain(rowMeans(pred[[2]][,4,]) / 92, 
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
         ref = sum(1 / 41:62),  
         zlim = c(0, 1.3), 
         picture.name = "inst/img/MAIN_Nt41-62lJJAs.png",
         title = expression(paste(bar(N)["41:62,JJA"](s), " (mean)")), 
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 5:14 / 10,
         label = paste0(10 * 5:14, "%"),
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,4,], 1, quantile, prob = 0.05) / 92, 
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
         ref = sum(1 / 41:62),  
         zlim = c(0, 1.3), 
         picture.name = "inst/img/MAIN_Nt41-62lJJAs-q05.png",
         title = expression(paste(bar(N)["41:62,JJA"](s), " (0.05 quantile)")), 
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 1:11 / 10,
         label = paste0(10 * 1:11, "%"),
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,4,], 1, quantile, prob = 0.95) / 92, 
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
         ref = sum(1 / 41:62),  
         zlim = c(0, 1.3), 
         picture.name = "inst/img/MAIN_Nt41-62lJJAs-q95.png",
         title = expression(paste(bar(N)["41:62,JJA"](s), " (0.95 quantile)")), 
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 8:18 / 10,
         label = paste0(10 * 8:18, "%"),
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)

mapSpain(rowMeans(pred[[2]][,5,]) / 91, 
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
         ref = sum(1 / 41:62),  
         zlim = c(0, 1.3), 
         picture.name = "inst/img/MAIN_Nt41-62lSONs.png",
         title = expression(paste(bar(N)["41:62,SON"](s), " (mean)")), 
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 3:8 / 10,
         label = paste0(10 * 3:8, "%"),
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,5,], 1, quantile, prob = 0.05) / 91, 
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
         ref = sum(1 / 41:62),  
         zlim = c(0, 1.3), 
         picture.name = "inst/img/MAIN_Nt41-62lSONs-q05.png",
         title = expression(paste(bar(N)["41:62,SON"](s), " (0.05 quantile)")), 
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 0:5 / 10,
         label = c("CRM", paste0(10 * 1:5, "%")),
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,5,], 1, quantile, prob = 0.95) / 91, 
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
         ref = sum(1 / 41:62),  
         zlim = c(0, 1.3), 
         picture.name = "inst/img/MAIN_Nt41-62lSONs-q95.png",
         title = expression(paste(bar(N)["41:62,SON"](s), " (0.95 quantile)")), 
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 6:11 / 10,
         label = paste0(10 * 6:11, "%"),
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)

## all year
interest <- colMeans((pred[[2]][,2,] + pred[[2]][,3,] + pred[[2]][,4,] + pred[[2]][,5,]) / 365)
round(mean(interest), 2)
round(quantile(interest, prob = c(0.05, 0.95)), 2)
## winter
interest <- colMeans(pred[[2]][,2,] / 90)
round(mean(interest), 2)
round(quantile(interest, prob = c(0.05, 0.95)), 2)
## spring
interest <- colMeans(pred[[2]][,3,] / 92)
round(mean(interest), 2)
round(quantile(interest, prob = c(0.05, 0.95)), 2)
## summer
interest <- colMeans(pred[[2]][,4,] / 92)
round(mean(interest), 2)
round(quantile(interest, prob = c(0.05, 0.95)), 2)
## autumn
interest <- colMeans(pred[[2]][,5,] / 91)
round(mean(interest), 2)
round(quantile(interest, prob = c(0.05, 0.95)), 2)

interest <- colMeans(((pred[[2]][,2,] + pred[[2]][,3,] + pred[[2]][,4,] + pred[[2]][,5,]) / 365) > sum(1 / 41:62))
mean(apply(((pred[[2]][,2,] + pred[[2]][,3,] + pred[[2]][,4,] + pred[[2]][,5,]) / 365), 1, quantile, prob = 0.025) > sum(1 / 41:62))



# maps of ratios in the last decade
mapSpain(rowMeans(pred[[2]][,6,]) / 90 / sum(1 / 53:62), 
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
         ref = 1,  
         zlim = c(0, 3.8), 
         picture.name = "inst/img/MAIN_Rt53-62lDJF.png",
         title = expression(paste(R["53:62,DJF"](s), " (mean)")), 
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 3:5 / 5,
         label = 1 + 3:5 / 5,
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,6,], 1, quantile, prob = 0.05) / 90 / sum(1 / 53:62), 
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
         ref = 1,  
         zlim = c(0, 3.8), 
         picture.name = "inst/img/MAIN_Rt53-62lDJF-q05.png",
         title = expression(paste(R["53:62,DJF"](s), " (0.05 quantile)")), 
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 1:2 / 5,
         label = 1 + 1:2 / 5,
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,6,], 1, quantile, prob = 0.95) / 90 / sum(1 / 53:62), 
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
         ref = 1,  
         zlim = c(0, 3.8), 
         picture.name = "inst/img/MAIN_Rt53-62lDJF-q95.png",
         title = expression(paste(R["53:62,DJF"](s), " (0.95 quantile)")), 
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 6:9 / 5,
         label = 1 + 6:9 / 5,
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)

mapSpain(rowMeans(pred[[2]][,7,]) / 92 / sum(1 / 53:62), 
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
         ref = 1,  
         zlim = c(0, 3.8), 
         picture.name = "inst/img/MAIN_Rt53-62lMAM.png",
         title = expression(paste(R["53:62,MAM"](s), " (mean)")), 
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 2:6 / 5,
         label = 1 + 2:6 / 5,
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,7,], 1, quantile, prob = 0.05) / 92 / sum(1 / 53:62), 
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
         ref = 1,  
         zlim = c(0, 3.8), 
         picture.name = "inst/img/MAIN_Rt53-62lMAM-q05.png",
         title = expression(paste(R["53:62,MAM"](s), " (0.05 quantile)")), 
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + -1:4 / 5,
         label = 1 + -1:4 / 5,
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,7,], 1, quantile, prob = 0.95) / 92 / sum(1 / 53:62), 
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
         ref = 1,  
         zlim = c(0, 3.8), 
         picture.name = "inst/img/MAIN_Rt53-62lMAM-q95.png",
         title = expression(paste(R["53:62,MAM"](s), " (0.95 quantile)")), 
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 4:9 / 5,
         label = 1 + 4:9 / 5,
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)

mapSpain(rowMeans(pred[[2]][,8,]) / 92 / sum(1 / 53:62), 
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
         ref = 1,  
         zlim = c(0, 3.8), 
         picture.name = "inst/img/MAIN_Rt53-62lJJA.png",
         title = expression(paste(R["53:62,JJA"](s), " (mean)")),
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 3:9 / 5,
         label = 1 + 3:9 / 5,
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,8,], 1, quantile, prob = 0.05) / 92 / sum(1 / 53:62), 
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
         ref = 1,  
         zlim = c(0, 3.8), 
         picture.name = "inst/img/MAIN_Rt53-62lJJA-q05.png",
         title = expression(paste(R["53:62,JJA"](s), " (0.05 quantile)")), 
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 0:6 / 5,
         label = 1 + 0:6 / 5,
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,8,], 1, quantile, prob = 0.95) / 92 / sum(1 / 53:62), 
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
         ref = 1,  
         zlim = c(0, 3.8), 
         picture.name = "inst/img/MAIN_Rt53-62lJJA-q95.png",
         title = expression(paste(R["53:62,JJA"](s), " (0.95 quantile)")), 
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 5:13 / 5,
         label = 1 + 5:13 / 5,
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)

mapSpain(rowMeans(pred[[2]][,9,]) / 91 / sum(1 / 53:62), 
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
         ref = 1,  
         zlim = c(0, 3.8), 
         picture.name = "inst/img/MAIN_Rt53-62lSON.png",
         title = expression(paste(R["53:62,SON"](s), " (mean)")), 
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 3:6 / 5,
         label = 1 + 3:6 / 5,
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,9,], 1, quantile, prob = 0.05) / 91 / sum(1 / 53:62), 
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
         ref = 1,  
         zlim = c(0, 3.8), 
         picture.name = "inst/img/MAIN_Rt53-62lSON-q05.png",
         title = expression(paste(R["53:62,SON"](s), " (0.05 quantile)")), 
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 0:8 / 5,
         label = 1 + 0:8 / 5,
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)
mapSpain(apply(pred[[2]][,9,], 1, quantile, prob = 0.95) / 91 / sum(1 / 53:62), 
         coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
         ref = 1,  
         zlim = c(0, 3.8), 
         picture.name = "inst/img/MAIN_Rt53-62lSON-q95.png",
         title = expression(paste(R["53:62,SON"](s), " (0.95 quantile)")), 
         legend.name = "",
         save = TRUE,
         contour = "smooth",
         breaks = 1 + 6:14 / 5,
         label = 1 + 6:14 / 5,
         smoothness = 15,
         threshold = 1e+05,
         threshold2 = 1e+05 * 2,
         dist = grid$dist)

## all year
interest <- colMeans((pred[[2]][,6,] + pred[[2]][,7,] + pred[[2]][,8,] + pred[[2]][,9,]) / 365 / sum(1 / 53:62))
round(mean(interest), 2)
round(quantile(interest, prob = c(0.05, 0.95)), 2)
## winter
interest <- colMeans(pred[[2]][,6,] / 90 / sum(1 / 53:62))
round(mean(interest), 2)
round(quantile(interest, prob = c(0.05, 0.95)), 2)
## spring
interest <- colMeans(pred[[2]][,7,] / 92 / sum(1 / 53:62))
round(mean(interest), 2)
round(quantile(interest, prob = c(0.05, 0.95)), 2)
## summer
interest <- colMeans(pred[[2]][,8,] / 92 / sum(1 / 53:62))
round(mean(interest), 2)
round(quantile(interest, prob = c(0.05, 0.95)), 2)
## autumn
interest <- colMeans(pred[[2]][,9,] / 91 / sum(1 / 53:62))
round(mean(interest), 2)
round(quantile(interest, prob = c(0.05, 0.95)), 2)

interest <- colMeans(((pred[[2]][,6,] + pred[[2]][,7,] + pred[[2]][,8,] + pred[[2]][,9,]) / 365) > sum(1 / 53:62))
round(mean(interest), 2)
round(quantile(interest, prob = c(0.05, 0.95)), 2)
mean(apply(((pred[[2]][,6,] + pred[[2]][,7,] + pred[[2]][,8,] + pred[[2]][,9,]) / 365), 1, quantile, prob = 0.05) > sum(1 / 53:62))
mean(apply((pred[[2]][,6,] / 90), 1, quantile, prob = 0.05) > sum(1 / 53:62))
mean(apply((pred[[2]][,7,] / 92), 1, quantile, prob = 0.05) > sum(1 / 53:62))
mean(apply((pred[[2]][,8,] / 92), 1, quantile, prob = 0.05) > sum(1 / 53:62))
mean(apply((pred[[2]][,9,] / 91), 1, quantile, prob = 0.05) > sum(1 / 53:62))



# ERS
ggERS <- function(y, xchar, ychar, ylim, xnumbreaks, xlabbreaks,
                  title = NULL,
                  picture.name = "photo.png",
                  save = TRUE) {
  
  n <- nrow(y) + 1
  df <- data.frame(
    y = c(1, 2:n * rowMeans(y)),
    t = 1:n,
    CI1 = c(1, 2:n * apply(y, 1, quantile, prob = 0.05)),
    CI2 = c(1, 2:n * apply(y, 1, quantile, prob = 0.95))
  )
  
  gg <- ggplot(data = df,
         mapping =  aes(x = t, y = y)) +
    geom_hline(yintercept=1,
               color = "gray") +
    theme_bw() +
    theme(legend.position="none") +
    ylab(ychar) +
    xlab(xchar) +
    ylim(ylim) +
    scale_x_continuous(breaks=xnumbreaks,
                       labels=xlabbreaks) +
    geom_ribbon(ggplot2::aes(
      ymin = CI1, ymax = CI2), 
      alpha = 0.2) +
    geom_line(size = 0.2)
  
  if (!is.null(title)) {
    gg <- gg + ggtitle(title)
  }
  
  if (save) {
    ggplot2::ggsave(picture.name, gg, width = 8.27 / 2, height = 11.69 / 4)
  }
  
  gg
}

ggERS(pred[[3]][,1,] / 365,
      "t (year)",
      expression(t %*% widehat(bar(ERS))[t](D)),
      c(0, 3.3),
      c(1, 21, 41, 61),
      c("1 (1960)", "21 (1980)", "41 (2000)", "61 (2020)"),
      picture.name = "inst/img/MAIN_ERS.pdf")

ggERS(pred[[3]][,2,] / 90,
      "t (year)",
      expression(t %*% widehat(bar(ERS))["t,DJF"](D)),
      c(0, 5.5),
      c(1, 21, 41, 61),
      c("1 (1960)", "21 (1980)", "41 (2000)", "61 (2020)"),
      title = "Winter",
      picture.name = "inst/img/SUPP_ERS-DJF.pdf")

ggERS(pred[[3]][,3,] / 92,
      "t (year)",
      expression(t %*% widehat(bar(ERS))["t,MAM"](D)),
      c(0, 5.5),
      c(1, 21, 41, 61),
      c("1 (1960)", "21 (1980)", "41 (2000)", "61 (2020)"),
      title = "Spring",
      picture.name = "inst/img/SUPP_ERS-MAM.pdf")

ggERS(pred[[3]][,4,] / 92,
      "t (year)",
      expression(t %*% widehat(bar(ERS))["t,JJA"](D)),
      c(0, 5.5),
      c(1, 21, 41, 61),
      c("1 (1960)", "21 (1980)", "41 (2000)", "61 (2020)"),
      title = "Summer",
      picture.name = "inst/img/SUPP_ERS-JJA.pdf")

ggERS(pred[[3]][,5,] / 91,
      "t (year)",
      expression(t %*% widehat(bar(ERS))["t,SON"](D)),
      c(0, 5.5),
      c(1, 21, 41, 61),
      c("1 (1960)", "21 (1980)", "41 (2000)", "61 (2020)"),
      title = "Autumn",
      picture.name = "inst/img/SUPP_ERS-SON.pdf")



# map daily
for (l in 1:8) {
  p <- mapSpain(rowMeans(pred[[1]][,l,]), 
                coords_limit = data.frame(X = c(-10, 4), Y = c(35.5, 44)), 
                ref = .5,  
                zlim = c(0, 1), 
                title = bquote(p[.(paste0("62,",l+221))](s)), 
                legend.name = "",
                save = FALSE,
                contour = c("none", rep("smooth", 7))[l],
                breaks = c(.1, .2, .5, .8, .9) / .5,
                label = c(.1, .2, .5, .8, .9),
                smoothness = 5,
                threshold = 1e+05 / 2,
                threshold2 = 1e+05,
                dist = grid$dist) + 
    geom_point(aes(x = sf::st_coordinates(stations$geometry)[,1],
               y = sf::st_coordinates(stations$geometry)[,2]),
               data = data.frame(X = 1:40),
               color = c("blue", "gray", "red")[round(2*data$Y[data$year == 62 & data$day == l+221]+1, 0)],
               size = 0.75, alpha = 0.75)
  
  ggplot2::ggsave(paste0("inst/img/MAIN_pt62l", l+221, "sOBSERVED.png"), 
                  p, width = 8.27 / 2, height = 11.69 / 4)
}
