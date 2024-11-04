#' Predict Method for Record Occurrence Model Fits
#' 
#' @description 
#'   Predicted values based on rom (record occurrence model) object.
#' 
#' @details 
#'   Predicted values are on the scale of the response (probabilities).
#' 
#' @param object Object of class inheriting from "rom"
#' @param newdata An optional list with \code{x} in which to look for variables 
#'   with which to predict, \code{coords} and \code{site} in which to predict, 
#'   and \code{day} or \code{year}. 
#'   If omitted, the fitted values are used.
#' @param newcoords New coordinates for the sites in \code{newdata$site} or 
#'   \code{newsite}.
#' @param newsite,newyear,newday A vector where for each observation includes 
#'   the site, year or day.
#' @param type the type of prediction required. For \eqn{K}-fold 
#'   cross-validation \code{"KFCV"} uses observed indicators from \code{data}
#'   and samples weak records from their prior probability.
#' @return A \code{"rom"} list with elements:
#'   \item{params}{Matrix where rows are MCMC simulations and cols are 
#'     parameters:
#'   \deqn{\beta_1,\ldots,\beta_k}}
#'   \item{\code{y}}{Data fitted}
#'   \item{\code{x}}{Covariates}
#'   
#' @author Jorge Castillo-Mateo
#'   
#' @importFrom stats delete.response model.frame model.matrix model.response terms
#'   
#' @export 
predict.rom <- function(object, 
                        newdata, 
                        newcoords,
                        newsite,
                        newyear,
                        newday,
                        type = c("KFCV", "sample")) {
  
  type <- match.arg(type)
  
  if (missing(newsite)) newsite <- newdata$site
  if (missing(newyear)) newyear <- newdata$year
  if (missing(newday))  newday  <- newdata$day
  
  newsite <- match(newsite, sort(unique(newsite)))
  newyear <- match(newyear, sort(unique(newyear)))
  newday  <- match(newday,  sort(unique(newday)))
  
  n <- length(unique(object$date$site))
  T <- length(unique(object$date$year))
  L <- length(unique(object$date$day))
  newn <- length(unique(newsite))
  newN <- nrow(newdata)
  
  tt <- terms(object)
  Y <- model.response(model.frame(tt, newdata))
  Terms <- delete.response(tt)
  m <- model.frame(Terms, newdata)
  m0 <- m1 <- m
  m0[!(m0[,"lag1"] %in% 0:1),"lag1"] <- 0
  m0[!(m0[,"lag2"] %in% 0:1),"lag2"] <- 0
  m1[!(m1[,"lag1"] %in% 0:1),"lag1"] <- 1
  m1[!(m1[,"lag2"] %in% 0:1),"lag2"] <- 1
  X <- model.matrix(Terms, m)[,colnames(object$x)]
  X0 <- model.matrix(Terms, m0)[,colnames(object$x)]
  X1 <- model.matrix(Terms, m1)[,colnames(object$x)]
  
  if (type == "KFCV" & object$number %in% c("2", "3", "4", "5")) {
    distance <- stats::dist(rbind(newcoords, object$coords))
    dAux <- matrix(0, nrow = newn + n, ncol = newn + n)
    dAux[lower.tri(dAux)] <- distance
    dAux <- dAux + t(dAux)
  }
  
  indLag1  <- grep("lag1", colnames(X))
  indLag2  <- grep("lag2", colnames(X))
  indLag12 <- intersect(indLag1, indLag2)
  indLag1  <- indLag1[!(indLag1 %in% indLag12)]
  indLag2  <- indLag2[!(indLag2 %in% indLag12)]
  
  if (object$scale.cov) {
    X  <- sweep(X,  2, attr(object$x, "scaled:center"), FUN = '-')
    X0 <- sweep(X0, 2, attr(object$x, "scaled:center"), FUN = '-')
    X1 <- sweep(X1, 2, attr(object$x, "scaled:center"), FUN = '-')
    X  <- sweep(X,  2, attr(object$x, "scaled:scale"), FUN = '/')
    X0 <- sweep(X0, 2, attr(object$x, "scaled:scale"), FUN = '/')
    X1 <- sweep(X1, 2, attr(object$x, "scaled:scale"), FUN = '/')
    
    if ("(Intercept)" %in% colnames(X)) {
      X[,"(Intercept)"] <- 1
      X0[,"(Intercept)"] <- 1
      X1[,"(Intercept)"] <- 1
    } 
  }
  
  B <- nrow(object$params$beta)
  
  indl  <- which(newday %in% 3:365)
  indl1 <- which(newday == 1)
  indl2 <- which(newday == 2)
  
  Xb <- matrix(nrow = B, ncol = newN)
  
  weak <- which(!(Y %in% 0:1))
  Nweak <- length(weak)
  
  day1 <- newday[weak] + 1
  day2 <- newday[weak] + 2
  year1 <- newyear[weak]
  year2 <- newyear[weak]
  year1[day1 > 365] <- year1[day1 > 365] + 1
  day1[day1 > 365]  <- day1[day1 > 365] - 365
  year2[day2 > 365] <- year2[day2 > 365] + 1
  day2[day2 > 365]  <- day2[day2 > 365] - 365
  
  weak1 <- match(apply(apply(cbind(day1, year1, newsite[weak]), 2, format, width = 3), 1, paste, collapse = ""),
                 apply(apply(cbind(newday, newyear, newsite), 2, format, width = 3), 1, paste, collapse = ""))
  weak2 <- match(apply(apply(cbind(day2, year2, newsite[weak]), 2, format, width = 3), 1, paste, collapse = ""),
                 apply(apply(cbind(newday, newyear, newsite), 2, format, width = 3), 1, paste, collapse = ""))
  
  indSubmodel <- which(colnames(X) %in% colnames(object$params$betal1))
  
  if (object$number == "1") {
    
    Xb <- predGlmBerKFCV(
      B, Nweak,
      X, X0, X1, Xb,
      object$params$beta, object$params$betal1, object$params$betal2,
      weak - 1, weak1 - 1, weak2 - 1,
      !is.na(weak1), !is.na(weak2),
      weak1[!is.na(weak1)] - 1, weak2[!is.na(weak2)] - 1, 
      indLag1 - 1, indLag2 - 1, indLag12 - 1,
      indl - 1, indl1 - 1, indl2 - 1, indSubmodel - 1,
      Y)
    
  } else if (object$number == "2") {
    
    site2 <- newsite
    site2[newday == 2] <- newsite[newday == 2] + newn
    site2[newday > 2] <- newsite[newday > 2] + 2 * newn
    
    if (type == "KFCV") {
      
      Ws <- matrix(nrow = B, ncol = 3 * newn)
      
      Xb <- predGlmBerGPtl(
        B, Nweak,
        n, newn, T, L,
        X, X0, X1, 
        site2 - 1,
        dAux, Xb, Ws,
        object$params$beta, object$params$betal1, object$params$betal2,
        object$params$ws, 
        object$params$hp[, c("(Intercept)", "prec0", "decay0", 
                             "(Intercept)l1", "prec0l1", 
                             "(Intercept)l2", "prec0l2")],
        weak - 1, weak1 - 1, weak2 - 1,
        !is.na(weak1), !is.na(weak2),
        weak1[!is.na(weak1)] - 1, weak2[!is.na(weak2)] - 1, 
        indLag1 - 1, indLag2 - 1, indLag12 - 1,
        indl - 1, indl1 - 1, indl2 - 1, indSubmodel - 1,
        Y)
      
    } else {
      
      Xb <- predGlmBerKFCV(
        B, Nweak,
        X, X0, X1, Xb,
        object$params$beta, object$params$betal1, object$params$betal2,
        weak - 1, weak1 - 1, weak2 - 1,
        !is.na(weak1), !is.na(weak2),
        weak1[!is.na(weak1)] - 1, weak2[!is.na(weak2)] - 1, 
        indLag1 - 1, indLag2 - 1, indLag12 - 1,
        indl - 1, indl1 - 1, indl2 - 1, indSubmodel - 1,
        Y) + object$params$ws[,site2]
      
    }
    
  } else if (object$number == "3") {
    
    site2 <- newsite
    site2[newday == 2] <- newsite[newday == 2] + newn
    site2[newday > 2] <- newsite[newday > 2] + 2 * newn
    time2 <- newyear
    time2[newday == 2] <- newyear[newday == 2] + T
    time2[newday > 2] <- newyear[newday > 2] + 2 * T
    
    if (type == "KFCV") {
      
      Ws <- matrix(nrow = B, ncol = 3 * newn)
      
      Xb <- predGlmBerGPtl(
        B, Nweak,
        n, newn, T, L,
        X, X0, X1, 
        site2 - 1,
        dAux, Xb, Ws,
        object$params$beta, object$params$betal1, object$params$betal2,
        object$params$ws, 
        object$params$hp[, c("(Intercept)", "prec0", "decay0", 
                             "(Intercept)l1", "prec0l1", 
                             "(Intercept)l2", "prec0l2")],
        weak - 1, weak1 - 1, weak2 - 1,
        !is.na(weak1), !is.na(weak2),
        weak1[!is.na(weak1)] - 1, weak2[!is.na(weak2)] - 1, 
        indLag1 - 1, indLag2 - 1, indLag12 - 1,
        indl - 1, indl1 - 1, indl2 - 1, indSubmodel - 1,
        Y) + object$params$wt[,time2]
      
    } else {

      Xb <- predGlmBerKFCV(
        B, Nweak,
        X, X0, X1, Xb,
        object$params$beta, object$params$betal1, object$params$betal2,
        weak - 1, weak1 - 1, weak2 - 1,
        !is.na(weak1), !is.na(weak2),
        weak1[!is.na(weak1)] - 1, weak2[!is.na(weak2)] - 1, 
        indLag1 - 1, indLag2 - 1, indLag12 - 1,
        indl - 1, indl1 - 1, indl2 - 1, indSubmodel - 1,
        Y) + object$params$ws[,site2] + object$params$wt[,time2]
      
    }
    
  } else if (object$number == "4") {
    
    site2 <- newsite
    site2[newday == 2] <- newsite[newday == 2] + newn
    site2[newday > 2] <- newsite[newday > 2] + 2 * newn
    time2 <- interaction(newyear, newday)
    time2 <- match(time2, sort(unique(time2)))
    
    if (type == "KFCV") {
      
      Ws <- matrix(nrow = B, ncol = newN)
      
      Xb <- predGlmBerGPtl(
        B, Nweak,
        n, newn, T, L,
        X, X0, X1, 
        site2 - 1,
        dAux, Xb, Ws,
        object$params$beta, object$params$betal1, object$params$betal2,
        object$params$ws, 
        object$params$hp[, c("(Intercept)", "prec0", "decay0", 
                             "(Intercept)l1", "prec0l1", 
                             "(Intercept)l2", "prec0l2")],
        weak - 1, weak1 - 1, weak2 - 1,
        !is.na(weak1), !is.na(weak2),
        weak1[!is.na(weak1)] - 1, weak2[!is.na(weak2)] - 1, 
        indLag1 - 1, indLag2 - 1, indLag12 - 1,
        indl - 1, indl1 - 1, indl2 - 1, indSubmodel - 1,
        Y) + object$params$wtl[,time2]
      
    } else {
      
      Xb <- predGlmBerKFCV(
        B, Nweak,
        X, X0, X1, Xb,
        object$params$beta, object$params$betal1, object$params$betal2,
        weak - 1, weak1 - 1, weak2 - 1,
        !is.na(weak1), !is.na(weak2),
        weak1[!is.na(weak1)] - 1, weak2[!is.na(weak2)] - 1, 
        indLag1 - 1, indLag2 - 1, indLag12 - 1,
        indl - 1, indl1 - 1, indl2 - 1, indSubmodel - 1,
        Y) + object$params$ws[,site2] + object$params$wtl[,time2]
      
    }
    
  } else if (object$number == "5") {
    
    if (type == "KFCV") {
      
      Wtls <- matrix(nrow = B, ncol = newN)
      
      Xb <- predGlmBerModelPaperKFCV(
        B, Nweak,
        n, newn, T, L,
        X, X0, X1, 
        object$date$year - 1, object$date$day - 1, newyear - 1, newday - 1,
        dAux, Xb, Wtls,
        object$params$beta, object$params$betal1, object$params$betal2,
        object$params$wtls, object$params$wtl, object$params$hp[, c("prec0", "decay0", "prec0l1", "prec0l2")],
        weak - 1, weak1 - 1, weak2 - 1,
        !is.na(weak1), !is.na(weak2),
        weak1[!is.na(weak1)] - 1, weak2[!is.na(weak2)] - 1, 
        indLag1 - 1, indLag2 - 1, indLag12 - 1,
        indl - 1, indl1 - 1, indl2 - 1, indSubmodel - 1,
        Y)
      
    } else {
      
      Xb <- predGlmBerKFCV(
        B, Nweak,
        X, X0, X1, Xb,
        object$params$beta, object$params$betal1, object$params$betal2,
        weak - 1, weak1 - 1, weak2 - 1,
        !is.na(weak1), !is.na(weak2),
        weak1[!is.na(weak1)] - 1, weak2[!is.na(weak2)] - 1, 
        indLag1 - 1, indLag2 - 1, indLag12 - 1,
        indl - 1, indl1 - 1, indl2 - 1, indSubmodel - 1,
        Y) + object$params$wtls
      
    }
  
  }
  
  Y <- 1 / (1 + exp(-Xb))
  
  return(Y)
}
