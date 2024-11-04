#' Record Occurrence Model Fitting
#' 
#' @description 
#'   This function fits the Record Occurrence Model in a Bayesian framework 
#'   using MCMC.
#'    
#' @details 
#'   The fitting algorithm is a data-augmentation Metropolis-within-Gibbs MCMC 
#'   algorithm. 
#'   
#'   \code{formula} requires two variables, \code{lag1} and \code{lag2}, 
#'   that are the record indicator one and two days before. 
#'   
#'   \code{data} requires column names \code{trend} (annual trend), 
#'   \code{lag1}, and \code{lag2}.
#'   
#'   The response variable must be \code{Y}. It can be 0's and 1's indicating 
#'   the occurrence of a record, but values \eqn{1/r} (\eqn{r=1,,2\ldots}) are 
#'   allowed indicating the probability of a weak record being a strong record.
#'   (If all values are 0's and 1's the function could return an error.)
#'   
#'   The models for the first and second days of the year are 
#'   \code{Y ~ trend + lag1}. 
#' 
#' @note It is necessary to have data ordered as using if c() over a 
#'   3-dimensional array where the first dimension is year, the second is day, 
#'   and the third is site. This is, first all data for the first site, then 
#'   for the first day and it progresses by years. (This will also be necessary 
#'   for prediction.)
#' 
#' @param formula an object of class "formula": a symbolic description of the 
#'   model to be fitted. (Has some requirements with \code{trend}, \code{lag1}, 
#'   and \code{lag2}).
#' @param data a data frame (\eqn{N \times p}) containing the variables in the 
#'   model.
#' @param coords a matrix or data frame (\eqn{n \times 2}) containing the 
#'   coordinates of the sites.
#' @param site,year,day vectors (\eqn{N \times 1}) containing the site, year, 
#'   and day the measurements were taken.
#' @param model a character indicating the model to fit from the paper.
#' @param scale.cov a boolean indicating whether or not the columns of the 
#'   design matrix must be scaled to have zero mean and unit variance.
#' @param decay.prior a character string indicating the prior for the decay 
#'   parameter, gamma or uniform in an interval.
#' @param inits a vector of initial values (if NULL, random values).
#' @param prior a list containing the parameters of the prior distributions.
#' @param n.sims,n.thin,n.burnin,n.report 
#'   (i) Number of iterations not discarded. 
#'   (ii) Thinning rate. 
#'   (iii) Number of iterations discarded at the beginning. 
#'   (iv) Report the number of iterations rate.
#' @return A list with many elements. Samples from the model parameters are in 
#'   \code{params} where rows are MCMC simulations and cols are parameters.
#'  
#' @author Jorge Castillo-Mateo
#'  
#' @importFrom stats model.matrix model.response rgamma rnorm runif
#'   
#' @export 
rom <- function(formula,
                data,
                coords,
                site,
                year,
                day,
                model = c("1", "2", "3", "4", "5"),
                scale.cov = TRUE,
                decay.prior = c("gamma", "unif"),
                inits = NULL, 
                prior = list("beta" = 1e-04, "sigma" = c(2, 1), "decay" = c(2, 1)),
                n.sims = 1000,
                n.thin = 1,
                n.burnin = 1000,
                n.report = 100
) {
  
  model <- match.arg(model)
  
  if (missing(site)) site <- data$site
  if (missing(year)) year <- data$year
  if (missing(day))  day  <- data$day
  site <- match(site, sort(unique(site)))
  year <- match(year, sort(unique(year)))
  day  <- match(day,  sort(unique(day)))
  n <- length(unique(site))
  T <- length(unique(year))
  L <- length(unique(day))
  
  if (L != 365) stop("365 days needed to fit the model")
  
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- al <- mf[c(1L, m)]
  al[[1L]] <- quote(stats::alias)
  names(al)[2] <- "object"
  al <- eval(al, parent.frame())
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "numeric")
  
  mf0 <- mf1 <- mf
  mf0[!(mf0[,"lag1"] %in% 0:1),"lag1"] <- 0
  mf0[!(mf0[,"lag2"] %in% 0:1),"lag2"] <- 0
  mf1[!(mf1[,"lag1"] %in% 0:1),"lag1"] <- 1
  mf1[!(mf1[,"lag2"] %in% 0:1),"lag2"] <- 1
  
  if (is.null(attr(al$Complete, "dimnames"))) {
    X  <- model.matrix(mt, mf )
    X0 <- model.matrix(mt, mf0)
    X1 <- model.matrix(mt, mf1)
  } else {
    X  <- model.matrix(mt, mf )[, attr(al$Complete, "dimnames")[[2]]]
    X0 <- model.matrix(mt, mf0)[, attr(al$Complete, "dimnames")[[2]]]
    X1 <- model.matrix(mt, mf1)[, attr(al$Complete, "dimnames")[[2]]]
  }
  
  if (model %in% c("2", "3", "4", "5")) {
    if (!("(Intercept)" %in% colnames(X))) 
      stop("Intercept must be included in the model")
    X <- X[,-1]
    X0 <- X0[,-1]
    X1 <- X1[,-1]
    
    distance <- stats::dist(coords)
    drange <- range(distance)
    dAux <- matrix(0, nrow = n, ncol = n)
    dAux[lower.tri(dAux)] <- distance
    dAux <- dAux + t(dAux)
    
    decay.prior <- match.arg(decay.prior)
    decay.prior <- (decay.prior == "unif")
  }
  
  indLag1  <- grep("lag1", colnames(X))
  indLag2  <- grep("lag2", colnames(X))
  indLag12 <- intersect(indLag1, indLag2)
  indLag1  <- indLag1[!(indLag1 %in% indLag12)]
  indLag2  <- indLag2[!(indLag2 %in% indLag12)]
  
  if (scale.cov) {
    X <- scale(X)
    X0 <- sweep(X0, 2, attr(X, "scaled:center"), FUN = '-')
    X1 <- sweep(X1, 2, attr(X, "scaled:center"), FUN = '-')
    X0 <- sweep(X0, 2, attr(X, "scaled:scale"), FUN = '/')
    X1 <- sweep(X1, 2, attr(X, "scaled:scale"), FUN = '/')
    
    if ("(Intercept)" %in% colnames(X)) {
      X[,"(Intercept)"] <- 1
      X0[,"(Intercept)"] <- 1
      X1[,"(Intercept)"] <- 1
    } 
  }
  
  k <- ncol(X)
  N <- length(Y)
  abLim <- rep(-1, N)
  abLim[Y == 1] <- 1
  weak <- which(!(Y %in% 0:1))
  Nweak <- length(weak)
  
  day1 <- day[weak] + 1
  day2 <- day[weak] + 2
  year1 <- year[weak]
  year2 <- year[weak]
  year1[day1 > 365] <- year1[day1 > 365] + 1
  day1[day1 > 365]  <- day1[day1 > 365] - 365
  year2[day2 > 365] <- year2[day2 > 365] + 1
  day2[day2 > 365]  <- day2[day2 > 365] - 365
  
  weak1 <- match(apply(apply(cbind(day1, year1, site[weak]), 2, format, width = 3), 1, paste, collapse = ""),
                 apply(apply(cbind(day, year, site), 2, format, width = 3), 1, paste, collapse = ""))
  weak2 <- match(apply(apply(cbind(day2, year2, site[weak]), 2, format, width = 3), 1, paste, collapse = ""),
                 apply(apply(cbind(day, year, site), 2, format, width = 3), 1, paste, collapse = ""))
  
  if ("(Intercept)" %in% colnames(X)) {
    if ("trend" %in% colnames(X)) {
      indSubmodel <- match(c("(Intercept)", "trend", "lag1"), colnames(X))
    } else if ("poly(trend, 2)1" %in% colnames(X)) {
      indSubmodel <- match(c("(Intercept)", "poly(trend, 2)1", "lag1"), colnames(X))
    } else {
      stop("Only trend with names 'trend' or 'poly(trend, 2)1' are supported")
    }
  } else {
    if ("trend" %in% colnames(X)) {
      indSubmodel <- match(c("trend", "lag1"), colnames(X))
    } else if ("poly(trend, 2)1" %in% colnames(X)) {
      indSubmodel <- match(c("poly(trend, 2)1", "lag1"), colnames(X))
    } else {
      stop("Only trend with names 'trend' or 'poly(trend, 2)1' are supported")
    }
  }

  indl  <- which(day %in% 3:365)
  indl1 <- which(day == 1)
  indl2 <- which(day == 2)
  Nl <- length(indl1)
  kl <- length(indSubmodel)
  
  if (model == "1") {
    
    keep <- matrix(nrow = n.sims / n.thin, ncol = k + 2 * kl)
    if (is.null(inits)) {
      inits <- rnorm(k + 2 * kl) #rep(0, k + 2 * kl)
    } else if ((k + 2 * kl) != length(inits)) {
      stop("'inits' does not have the proper length 'k + 2 * kl'")
    }
    
    keep <- glmBer(
      N, k, Nl, kl, Nweak, X, X0, X1, 
      inits[1:k], inits[1:kl + k], inits[1:kl + k + kl], 
      weak - 1, weak1 - 1, weak2 - 1,
      !is.na(weak1), !is.na(weak2),
      weak1[!is.na(weak1)] - 1, weak2[!is.na(weak2)] - 1, 
      indLag1 - 1, indLag2 - 1, indLag12 - 1,
      indl - 1, indl1 - 1, indl2 - 1, indSubmodel - 1,
      Y, abLim,
      prior$beta * diag(k), keep, 
      n.sims, n.thin, n.burnin, n.report
    )
    
    colnames(keep) <- c(colnames(X), 
                        colnames(X)[indSubmodel],
                        colnames(X)[indSubmodel])
    
    keep.list <- list()
    keep.list$params$beta   <- keep[,1:k]
    keep.list$params$betal1 <- keep[,1:kl + k]
    keep.list$params$betal2 <- keep[,1:kl + k + kl]
    keep.list$date$day <- day
    keep.list$date$year <- year
    keep.list$date$site <- site
    keep.list$y <- Y
    keep.list$x <- X
    keep.list$call <- cl
    keep.list$terms <- mt
    keep.list$model <- mf
    keep.list$scale.cov <- scale.cov
    keep.list$number <- model
    
  } else if (model == "2") {
    
    # beta + ws (l1 l2 others)
    d <- k + 3 * n
    # [+ beta0 + prec0 + decay0 + (l1 l2)]
    keep <- matrix(nrow = n.sims / n.thin, ncol = d + 3 + 2 * kl + 4 + 4) 
    if (is.null(inits)) {
      inits <- rnorm(d + 3 + 2 * kl + 4)
      inits[d + c(2,kl+5,2*kl+7)] <- rgamma(3, 2, 1)
      if (decay.prior) {
        inits[d + 3] <- runif(1, prior$decay[1], prior$decay[2])
      } else {
        inits[d + 3] <- runif(1, 3 / drange[2], 3 / drange[1])
      }
    } else if ((d + 3 + 2 * kl + 4) != length(inits)) {
      stop("'inits' does not have the proper length")
    }
    
    site2 <- site
    site2[day == 2] <- site[day == 2] + n
    site2[day > 2] <- site[day > 2] + 2 * n
    time2 <- rep(1, N)
    
    keep <- glmBerGPtl(
      N, k, Nl, kl, Nweak, X, X0, X1, 
      inits[1:k], inits[1:kl + d + 3], inits[1:kl + d + 5 + kl],
      c(inits[kl + 4:5 + d], 0), c(inits[2*kl + 6:7 + d], 0),
      inits[k + 1:(3*n)], c(inits[d + 1:3], 0), 
      dAux, decay.prior, n, 
      0, T, L,
      site2 - 1, time2 - 1,
      weak - 1, weak1 - 1, weak2 - 1,
      !is.na(weak1), !is.na(weak2),
      weak1[!is.na(weak1)] - 1, weak2[!is.na(weak2)] - 1, 
      indLag1 - 1, indLag2 - 1, indLag12 - 1,
      indl - 1, indl1 - 1, indl2 - 1, indSubmodel - 1,
      Y, abLim,
      prior$beta * diag(k), prior$beta, 
      prior$sigma[1], prior$sigma[2], prior$decay[1], prior$decay[2],
      keep,
      n.sims, n.thin, n.burnin, n.report
    )
    
    colnames(keep) <- c(colnames(X),
                        paste0("ws", 1:n, "l1"), paste0("ws", 1:n, "l2"), paste0("ws", 1:n), "removeV",
                        "(Intercept)", "prec0", "decay0", "removeV",
                        colnames(X)[indSubmodel], "(Intercept)l1", "prec0l1", "removeV",
                        colnames(X)[indSubmodel], "(Intercept)l2", "prec0l2", "removeV")
    
    keep <- keep[,!colnames(keep) == "removeV"]
    
    keep.list <- list()
    keep.list$params$beta   <- keep[,1:k]
    keep.list$params$ws     <- keep[,1:(3*n) + k]
    keep.list$params$hp     <- keep[,c(1:3,kl+4:5,2*kl+6:7) + d]
    keep.list$params$betal1 <- keep[,1:kl + d + 3]
    keep.list$params$betal2 <- keep[,1:kl + d + 5 + kl]
    keep.list$date$day  <- day
    keep.list$date$year <- year
    keep.list$date$site <- site
    keep.list$coords <- coords
    keep.list$y <- Y
    keep.list$x <- X
    keep.list$call <- cl
    keep.list$terms <- mt
    keep.list$model <- mf
    keep.list$scale.cov <- scale.cov
    keep.list$number <- model
    
  } else if (model == "3") {
    
    # beta + ws (l1 l2 others) + wt (l1 l2 others)
    d <- k + 3 * n + 3 * T
    # [+ beta0 + prec0 + decay0 + prec1 + (l1 l2)]
    keep <- matrix(nrow = n.sims / n.thin, ncol = d + 4 + 2 * kl + 6) 
    if (is.null(inits)) {
      inits <- rnorm(d + 4 + 2 * kl + 6)
      inits[d + c(2,4,kl+6:7,2*kl+9:10)] <- rgamma(6, 2, 1)
      if (decay.prior) {
        inits[d + 3] <- runif(1, prior$decay[1], prior$decay[2])
      } else {
        inits[d + 3] <- runif(1, 3 / drange[2], 3 / drange[1])
      }
    } else if ((d + 4 + 2 * kl + 6) != length(inits)) {
      stop("'inits' does not have the proper length")
    }
    
    site2 <- site
    site2[day == 2] <- site[day == 2] + n
    site2[day > 2] <- site[day > 2] + 2 * n
    time2 <- year
    time2[day == 2] <- year[day == 2] + T
    time2[day > 2] <- year[day > 2] + 2 * T
    
    keep <- glmBerGPtl(
      N, k, Nl, kl, Nweak, X, X0, X1, 
      inits[1:k], inits[1:kl + d + 4], inits[1:kl + d + 7 + kl],
      inits[kl + 5:7 + d], inits[2*kl + 8:10 + d],
      inits[k + 1:(3*n)], inits[d + 1:4], 
      dAux, decay.prior, n,
      inits[(k+3*n) + 1:(3*T)], T, L,
      site2 - 1, time2 - 1,
      weak - 1, weak1 - 1, weak2 - 1,
      !is.na(weak1), !is.na(weak2),
      weak1[!is.na(weak1)] - 1, weak2[!is.na(weak2)] - 1, 
      indLag1 - 1, indLag2 - 1, indLag12 - 1,
      indl - 1, indl1 - 1, indl2 - 1, indSubmodel - 1,
      Y, abLim,
      prior$beta * diag(k), prior$beta, 
      prior$sigma[1], prior$sigma[2], prior$decay[1], prior$decay[2],
      keep,
      n.sims, n.thin, n.burnin, n.report
    )
    
    colnames(keep) <- c(colnames(X),
                        paste0("ws", 1:n, "l1"), paste0("ws", 1:n, "l2"), paste0("ws", 1:n), 
                        paste0("wt", 1:T, "l1"), paste0("wt", 1:T, "l2"), paste0("wt", 1:T),
                        "(Intercept)", "prec0", "decay0", "prec1",
                        colnames(X)[indSubmodel], "(Intercept)l1", "prec0l1", "prec1l1",
                        colnames(X)[indSubmodel], "(Intercept)l2", "prec0l2", "prec1l2")
    
    keep.list <- list()
    keep.list$params$beta   <- keep[,1:k]
    keep.list$params$ws     <- keep[,1:(3*n) + k]
    keep.list$params$wt     <- keep[,1:(3*T) + k + 3*n]
    keep.list$params$hp     <- keep[,c(1:4,kl+5:7,2*kl+8:10) + d]
    keep.list$params$betal1 <- keep[,1:kl + d + 4]
    keep.list$params$betal2 <- keep[,1:kl + d + 7 + kl]
    keep.list$date$day  <- day
    keep.list$date$year <- year
    keep.list$date$site <- site
    keep.list$coords <- coords
    keep.list$y <- Y
    keep.list$x <- X
    keep.list$call <- cl
    keep.list$terms <- mt
    keep.list$model <- mf
    keep.list$scale.cov <- scale.cov
    keep.list$number <- model
    
  } else if (model == "4") {
    
    # beta + ws (l1 l2 others) + wtl
    d <- k + 3 * n + T * L
    # [+ beta0 + prec0 + decay0 + prec1 + (l1 l2)]
    keep <- matrix(nrow = n.sims / n.thin, ncol = d + 4 + 2 * kl + 6) 
    if (is.null(inits)) {
      inits <- rnorm(d + 4 + 2 * kl + 6)
      inits[d + c(2,4,kl+6:7,2*kl+9:10)] <- rgamma(6, 2, 1)
      if (decay.prior) {
        inits[d + 3] <- runif(1, prior$decay[1], prior$decay[2])
      } else {
        inits[d + 3] <- runif(1, 3 / drange[2], 3 / drange[1])
      }
    } else if ((d + 4 + 2 * kl + 6) != length(inits)) {
      stop("'inits' does not have the proper length")
    }
    
    site2 <- site
    site2[day == 2] <- site[day == 2] + n
    site2[day > 2] <- site[day > 2] + 2 * n
    time2 <- interaction(year, day)
    time2 <- match(time2, sort(unique(time2)))
    
    keep <- glmBerGPtl(
      N, k, Nl, kl, Nweak, X, X0, X1, 
      inits[1:k], inits[1:kl + d + 4], inits[1:kl + d + 7 + kl],
      inits[kl + 5:7 + d], inits[2*kl + 8:10 + d],
      inits[k + 1:(3*n)], inits[d + 1:4], 
      dAux, decay.prior, n,
      inits[(k+3*n) + 1:(T*L)], T, L,
      site2 - 1, time2 - 1,
      weak - 1, weak1 - 1, weak2 - 1,
      !is.na(weak1), !is.na(weak2),
      weak1[!is.na(weak1)] - 1, weak2[!is.na(weak2)] - 1, 
      indLag1 - 1, indLag2 - 1, indLag12 - 1,
      indl - 1, indl1 - 1, indl2 - 1, indSubmodel - 1,
      Y, abLim,
      prior$beta * diag(k), prior$beta, 
      prior$sigma[1], prior$sigma[2], prior$decay[1], prior$decay[2],
      keep,
      n.sims, n.thin, n.burnin, n.report
    )
    
    colnames(keep) <- c(colnames(X),
                        paste0("ws", 1:n, "l1"), paste0("ws", 1:n, "l2"), paste0("ws", 1:n),
                        paste0("wt", 1:T, "l", rep(1:L, each = T)),
                        "(Intercept)", "prec0", "decay0", "prec1",
                        colnames(X)[indSubmodel], "(Intercept)l1", "prec0l1", "prec1l1",
                        colnames(X)[indSubmodel], "(Intercept)l2", "prec0l2", "prec1l2")
    
    keep.list <- list()
    keep.list$params$beta   <- keep[,1:k]
    keep.list$params$ws     <- keep[,1:(3*n) + k]
    keep.list$params$wtl    <- keep[,1:(T*L) + k + 3*n]
    keep.list$params$hp     <- keep[,c(1:4,kl+5:7,2*kl+8:10) + d]
    keep.list$params$betal1 <- keep[,1:kl + d + 4]
    keep.list$params$betal2 <- keep[,1:kl + d + 7 + kl]
    keep.list$date$day  <- day
    keep.list$date$year <- year
    keep.list$date$site <- site
    keep.list$coords <- coords
    keep.list$y <- Y
    keep.list$x <- X
    keep.list$call <- cl
    keep.list$terms <- mt
    keep.list$model <- mf
    keep.list$scale.cov <- scale.cov
    keep.list$number <- model
    
  } else if (model == "5") {
    
    # beta + wtls + wtl
    d <- k + n * T * L + T * L
    # [+ beta0 + prec0 + decay0 + prec1 + (l1 l2)]
    keep <- matrix(nrow = n.sims / n.thin, ncol = d + 4 + 2 * kl + 6) 
    if (is.null(inits)) {
      inits <- rnorm(d + 4 + 2 * kl + 6) # rep(0, d + 4 + 2 * kl + 6)
      inits[d + c(2,4,kl+6:7,2*kl+9:10)] <- rgamma(6, 2, 1) # 1
      if (decay.prior) {
        inits[d + 3] <- runif(1, prior$decay[1], prior$decay[2])
      } else {
        inits[d + 3] <- runif(1, 3 / drange[2], 3 / drange[1]) # 0.5 * 3 / drange[2]
      }
    } else if ((d + 4 + 2 * kl + 6) != length(inits)) {
      stop("'inits' does not have the proper length")
    }
    
    keep <- glmBerModelPaper(
      N, k, Nl, kl, Nweak, X, X0, X1, 
      inits[1:k], inits[1:kl + d + 4], inits[1:kl + d + 7 + kl],
      inits[kl + 5:7 + d], inits[2*kl + 8:10 + d],
      inits[k + 1:(n*T*L)], inits[d + 1:4], 
      dAux, decay.prior, n,
      inits[(k+n*T*L) + 1:(T*L)], T, L,
      site - 1, year - 1, day - 1,
      weak - 1, weak1 - 1, weak2 - 1,
      !is.na(weak1), !is.na(weak2),
      weak1[!is.na(weak1)] - 1, weak2[!is.na(weak2)] - 1, 
      indLag1 - 1, indLag2 - 1, indLag12 - 1,
      indl - 1, indl1 - 1, indl2 - 1, indSubmodel - 1,
      Y, abLim,
      prior$beta * diag(k), prior$beta, 
      prior$sigma[1], prior$sigma[2], prior$decay[1], prior$decay[2],
      keep,
      n.sims, n.thin, n.burnin, n.report
    )
    
    colnames(keep) <- c(colnames(X),
                        paste0("wt", 1:T, "l", rep(1:L, each = T), "s", rep(1:n, each = T * L)),
                        paste0("wt", 1:T, "l", rep(1:L, each = T)),
                        "(Intercept)", "prec0", "decay0", "prec1",
                        colnames(X)[indSubmodel], "(Intercept)l1", "prec0l1", "prec1l1",
                        colnames(X)[indSubmodel], "(Intercept)l2", "prec0l2", "prec1l2")
    
    keep.list <- list()
    keep.list$params$beta   <- keep[,1:k]
    keep.list$params$wtls   <- keep[,1:(T*L*n) + k]
    keep.list$params$wtl    <- keep[,1:(T*L) + k + T*L*n]
    keep.list$params$hp     <- keep[,c(1:4,kl+5:7,2*kl+8:10) + d]
    keep.list$params$betal1 <- keep[,1:kl + d + 4]
    keep.list$params$betal2 <- keep[,1:kl + d + 7 + kl]
    keep.list$date$day  <- day
    keep.list$date$year <- year
    keep.list$date$site <- site
    keep.list$coords <- coords
    keep.list$y <- Y
    keep.list$x <- X
    keep.list$call <- cl
    keep.list$terms <- mt
    keep.list$model <- mf
    keep.list$scale.cov <- scale.cov
    keep.list$number <- model
    
  } else {
    
    stop(paste("'model'", model, "is not implemented"))
    
  }
  
  class(keep.list) <- "rom"
  
  return(keep.list)
}
