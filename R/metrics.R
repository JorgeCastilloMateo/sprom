#' Metrics for Validation
#' 
#' @description 
#'   This function is used to validate the predicted probabilities and number
#'   of records with the actual values.
#' 
#' @details 
#'   Values different from 0 or 1 (tied records) are removed from the metrics.
#'   
#'   The DIC uses mean values instead of removing or sampling values different 
#'   from 0 or 1.
#' 
#' @param Y The observed values (vector)
#' @param p The estimated probabilities (vector or matrix with rows being 
#'   replicates)
#' @param year,day,site The year, day, and site associated with the values 
#'   (vectors) (\code{years} must be in order from lower to higher)
#' @param metric A character string \code{"point"}, \code{"N"}, or \code{"DIC"}
#' @param D1,D2 Two shorter year periods on which to calculate the metrics, in 
#'   addition to the period of all years
#' @param B If \code{p} is a vector, the number of replicates to calculate AD's
#' @return If \code{metric = "point"}; the metrics (vector):
#'   \item{BS}{Brier Score}
#'   \item{BS1}{BS in \code{D1}}
#'   \item{BS2}{BS in \code{D2}}
#'   \item{AUC}{Area under the ROC curve}
#'   \item{AUC1}{AUC in \code{D1}}
#'   \item{AUC2}{AUC in \code{D2}}
#'   
#'   If \code{metric = "N"}; a vector of AD's averaged grouped by year
#'   
#'   If \code{metric = "DIC"}; a vector with DIC, D (deviance), and penalty (pD)
#'   
#' @author Jorge Castillo-Mateo
#'   
#' @importFrom stats rbinom
#' @importFrom pROC auc
#'   
#' @export 
metrics <- function(Y, p, year, day, site, metric = c("point", "N", "DIC"),
                    D1 = 2:(62 / 2), D2 = (62 / 2 + 1):62, B = 2000) {
  
  metric <- match.arg(metric)
  
  if (metric == "DIC") {
    p <- t(p)
    D <- -2 * mean(colSums(Y * log(p) + (1 - Y) * log(1 - p)))
    
    p <- rowMeans(p)
    pD <- D + 2 * sum(Y * log(p) + (1 - Y) * log(1 - p))
    
    return(c("DIC" = D + pD, "D" = D, "pD" = pD))
  }
  
  weak <- !(Y %in% 0:1)
  Y    <- Y[!weak]
  year <- year[!weak]
  day  <- day[!weak]
  site <- site[!weak]
  
  if (metric == "N") { # AD
    
    daysite <- interaction(day, site)
    lengthY <- length(Y)
  
    if (is.matrix(p)) {
      p <- p[,!weak]
      B <- nrow(p)
      N <- matrix(nrow = lengthY, ncol = B)
      for (b in 1:B) {
        N[,b] <- unlist(tapply(rbinom(lengthY, 1, p[b,]), daysite, cumsum))
      }
    } else {
      p <- p[!weak]
      N <- matrix(nrow = lengthY, ncol = B)
      for (b in 1:B) {
        N[,b] <- unlist(tapply(rbinom(lengthY, 1, p), daysite, cumsum))
      }
    }
  
    Nobs <- unlist(tapply(Y, daysite, cumsum))

    metric <- tapply(c(abs(N - Nobs)), rep(year, B), mean)
    names(metric) <- paste0("ADt", unique(year))
    
  } else {
    
    if (is.matrix(p)) p <- colMeans(p)
    p <- p[!weak]
    
    metric <- rep(NA, 6)
    names(metric) <- c("BS", "BS1", "BS2", "AUC", "AUC1", "AUC2")
    
    # Brier Score
    r2 <- (Y - p)^2
    metric[1] <- 100 * mean(r2)
    metric[2] <- 100 * mean(r2[year %in% D1])
    metric[3] <- 100 * mean(r2[year %in% D2])
  
    # AUC
    metric[4] <- pROC::auc(Y, p)
    metric[5] <- pROC::auc(Y[year %in% D1], p[year %in% D1])
    metric[6] <- pROC::auc(Y[year %in% D2], p[year %in% D2])
    
  }
  
  return(metric)
}
