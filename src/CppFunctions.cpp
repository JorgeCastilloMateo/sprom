#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// @description Random generation for the truncated normal distribution with
//   mean equal to mu, standard deviation equal to sigma and lower and upper
//   limits (0, Inf).
// @param N number of observations.
// @param mu vector of means.
// @param sigma vector of standard deviations.
// [[Rcpp::export]]
arma::vec rtnorm(const int N, arma::vec mu, arma::vec sigma) {
  
  double a, alpha, z, g;
  arma::vec x(N);
  
  for (int i = 0; i < N; ++i) {
    a = - mu(i) / sigma(i);
    if (a > 0.67448975) { // Robert 1995
      alpha = (a + sqrt(a * a + 4)) / 2;
      do {
        z = R::rexp(1 / alpha) + a;
        g = exp(-pow(z - alpha, 2) / 2);
      } while (R::runif(0, 1) > g);
    } else { // inverse method
      z = -R::qnorm(R::pnorm(-a, 0, 1, 1, 0) * R::runif(0, 1), 0, 1, 1, 0);
    }
    
    x(i) = z * sigma(i) + mu(i);
  }
  
  return x;
}

// @description Not vectorized truncated normal in (a, b).
// [[Rcpp::export]]
double rtnorm1(double mu, double sigma, double a, double b) {

  double pAlpha = R::pnorm(a, mu, sigma, 1, 0);
  double pBeta  = R::pnorm(b, mu, sigma, 1, 0);
  double x = R::qnorm(pAlpha + R::runif(0, 1) * (pBeta - pAlpha), mu, sigma, 1, 0);
  
  return x;
}

// [[Rcpp::export]]
double dtnorm1(double x, double mu, double sigma, double a, double b) { // log
  return R::dnorm(x, mu, sigma, 1) - std::log(R::pnorm(b, mu, sigma, 1, 0) - R::pnorm(a, mu, sigma, 1, 0));
}



// @description Random generation for the KS distribution.
// @details 
//   t2 = t * t
//   tP = pi * pi / (8 * t2)
//   c1 = 4 * t * t / (pi * pi)
//   c2 = 4 * exp(-6 * t * t)
// @param N number of observations.
// [[Rcpp::export]]
arma::vec rkolmogorov(const int N) {
  
  static const double p  = 0.372833;
  static const double pi = 3.141593;
  static const double t2 = 0.5625;
  static const double tP = 2.193245;
  static const double c1 = 0.2279727;
  static const double c2 = 0.1368725;
  
  arma::vec output(N);
  
  bool next, accept;
  arma::vec E(2);
  double G, X, W, Z, P, Q, U;
  int n, n2;
  double auxG;
  
  for (int i = 0; i < N; ++i) {
    
    next = false;
    
    if (R::runif(0, 1) < p) {
      
      while (!next) {
        accept = false;
        while (!accept) {
          E = Rcpp::rexp(2);
          E(0) /=  1 - 1 / (2 * tP);
          E(1) *= 2;
          G = tP + E(0);
          accept = (E(0) * E(0) <= tP * E(1) * (G + tP));
          if (!accept) {
            auxG = G / tP;
            accept = (auxG - 1 - log(auxG) <= E(1));
          }
        }
        U = R::runif(0, 1);
        X = pi / sqrt(8 * G);
        if (U >= c1) {
          output(i) = X;
          break;
        }
        W = 0;
        Z = 1 / (2 * G);
        P = exp(-G);
        n = 1;
        Q = 1;
        while (U >= W) {
          W += Z * Q;
          if (U >= W) {
            output(i) = X;
            next = true;
            break;
          }
          n += 2;
          n2 = n * n;
          Q = pow(P, n2 - 1);
          W -= n2 * Q;
        }
      }
      
    } else {
      
      while (!next) {
        E(0) = R::rexp(1);
        U = R::runif(0, 1);
        X = sqrt(t2 + E(0) / 2);
        if (U >= c2) {
          output(i) = X;
          break;
        }
        W = 0;
        n = 1;
        Z = exp(-2 * X * X);
        while (U > W) {
          ++n;
          n2 = n * n;
          W += n2 * pow(Z, n2 - 1);
          if (U >= W) {
            output(i) = X;
            next = true;
            break;
          }
          ++n;
          n2 = n * n;
          W -= n2 * pow(Z, n2 - 1);
        }
      }
      
    }
    
  }
  
  return output;
}

// @description Random generation for the lambda's in the augmented logistic 
//   regression model.
// @param N number of observations.
// @param old old vector of lambda's.
// @param e vector of eta = (Y - Xb)'s.
// [[Rcpp::export]]
arma::vec MH(const int N, arma::vec old, const arma::vec e) {
  
  static const double p  = 0.372833;
  static const double pi = 3.141593;
  static const double t2 = 0.5625;
  static const double tP = 2.193245;
  static const double c1 = 0.2279727;
  static const double c2 = 0.1368725;
  
  bool next, accept;
  arma::vec E(2);
  double A, G, X, W, Z, P, Q, U;
  int n, n2;
  double auxG;
  
  for (int i = 0; i < N; ++i) {
    
    next = false;
    
    if (R::runif(0, 1) < p) {
      
      while (!next) {
        accept = false;
        while (!accept) {
          E = Rcpp::rexp(2);
          E(0) /=  1 - 1 / (2 * tP);
          E(1) *= 2;
          G = tP + E(0);
          accept = (E(0) * E(0) <= tP * E(1) * (G + tP));
          if (!accept) {
            auxG = G / tP;
            accept = (auxG - 1 - log(auxG) <= E(1));
          }
        }
        U = R::runif(0, 1);
        X = pi / sqrt(8 * G);
        if (U >= c1) {
          break;
        }
        W = 0;
        Z = 1 / (2 * G);
        P = exp(-G);
        n = 1;
        Q = 1;
        while (U >= W) {
          W += Z * Q;
          if (U >= W) {
            next = true;
            break;
          }
          n += 2;
          n2 = n * n;
          Q = pow(P, n2 - 1);
          W -= n2 * Q;
        }
      }
      
    } else {
      
      while (!next) {
        E(0) = R::rexp(1);
        U = R::runif(0, 1);
        X = sqrt(t2 + E(0) / 2);
        if (U >= c2) {
          break;
        }
        W = 0;
        n = 1;
        Z = exp(-2 * X * X);
        while (U > W) {
          ++n;
          n2 = n * n;
          W += n2 * pow(Z, n2 - 1);
          if (U >= W) {
            next = true;
            break;
          }
          ++n;
          n2 = n * n;
          W -= n2 * pow(Z, n2 - 1);
        }
      }
      
    }
    
    X = pow(2 * X, 2);
    A = (log(old(i) / X) + e(i) * e(i) * (1 / old(i) - 1 / X)) / 2;
    if (log(R::runif(0, 1)) <= A) { old(i) = X; }
  }
  
  return old;
}    

// @description Iterations of the Metropolis-within-Gibbs algorithm for the 
//   logistic regression model.
// @param N number of observations.
// @param k number of covariates.
// @param Nl number of observations in l = 1 or l = 2.
// @param kl number of covariates for l = 1 or l = 2 models.
// @param Nweak number of weak records.
// @param X,X0,X1 matrix of covariates (X0 and X1 when weak records are all 0 
//   or all 1, respectively).
// @param beta,betal1,betal2 vector initial value of beta, betal1 and betal2.
// @param weak,weak1,weak2 vector of positions of weak records (weak1 and weak2
//   position in lag1 and lag2).
// @param indLag1,indLag2,indLag12 vector of positions of columns with 
//   lag1, lag2, or lag1:lag2.
// @param prob vector of probabilities that the observed indicator is 1.
// @param abLim vector indicating a record (1) or non-record (-1).
// @param V matrix prior var-cov of beta.
// @param keep matrix of beta's to keep.
// @param ... other arguments.
// [[Rcpp::export]]
arma::mat glmBer(
  const int N, const int k, const int Nl, const int kl, const int Nweak, 
  arma::mat X, const arma::mat X0, const arma::mat X1, 
  arma::vec beta, arma::vec betal1, arma::vec betal2, 
  const arma::uvec weak, const arma::uvec weak1, const arma::uvec weak2,
  const arma::vec noNA1, const arma::vec noNA2, 
  const arma::uvec Weak1, const arma::uvec Weak2, 
  const arma::uvec indLag1, const arma::uvec indLag2, const arma::uvec indLag12, 
  const arma::uvec indl, const arma::uvec indl1, const arma::uvec indl2, 
  const arma::uvec indSubmodel, 
  const arma::vec prob,
  arma::vec abLim,
  const arma::mat V, arma::mat keep,
  const int nSims, const int nThin, const int nBurnin, const int nReport) {

  arma::vec Xb(N), Xl(k), xl(kl);
  arma::vec Z(N); // the first sampled (no need to specify)
  arma::vec lambda = rkolmogorov(N);
  const arma::vec zerok(k, arma::fill::zeros);
  
  arma::mat OmegaBeta(k, k);
  arma::mat OmegaBetal(kl, kl);
  arma::uword ind; 
  arma::uvec ind1, ind2;
  
  for (int b = 1 - nBurnin; b <= nSims; ++b) {
    
    if ((b != 0) && (b % nReport == 0)) {
      Rcpp::Rcout << "The value of b : " << b << "\n";
    }
    
    // weak records
    X.submat(Weak1, indLag12) = X1.submat(Weak1, indLag12);
    X.submat(Weak2, indLag12) = X1.submat(Weak2, indLag12);
    for (int i = 0; i < Nweak; i++) {
      ind = weak(i);
      if (R::runif(0, 1) < prob(ind)) {
        abLim(ind) = 1;
        if (noNA1(i)) {
          ind1  = weak1(i);
          X.submat(ind1, indLag1) = X1.submat(ind1, indLag1);
        } 
        if (noNA2(i)) {
          ind2  = weak2(i);
          X.submat(ind2, indLag2) = X1.submat(ind2, indLag2);
        } 
      } else {
        abLim(ind) = -1;
        if (noNA1(i)) {
          ind1 = weak1(i);
          X.submat(ind1, indLag1)  = X0.submat(ind1, indLag1);
          X.submat(ind1, indLag12) = X0.submat(ind1, indLag12);
        }
        if (noNA2(i)) {
          ind2 = weak2(i);
          X.submat(ind2, indLag2)  = X0.submat(ind2, indLag2);
          X.submat(ind2, indLag12) = X0.submat(ind2, indLag12);
        }
      }
    }
    
    // Z and lambda
    Xb.rows(indl)  = X.rows(indl)  * beta;
    Xb.rows(indl1) = X.submat(indl1, indSubmodel) * betal1;
    Xb.rows(indl2) = X.submat(indl2, indSubmodel) * betal2;
    Z = abLim % rtnorm(N, abLim % Xb, sqrt(lambda));
    lambda = MH(N, lambda, Z - Xb);
    
    // beta
    OmegaBeta = V;
    beta = zerok;
    for (int i = 0; i < (N - 2 * Nl); ++i) {
      ind = indl(i);
      Xl = X.row(ind).t() / lambda(ind);
      OmegaBeta += Xl * X.row(ind);
      beta += Xl * Z(ind);
    }
    OmegaBeta = arma::inv_sympd(OmegaBeta);
    beta = OmegaBeta * beta + arma::chol(OmegaBeta, "lower") * arma::randn(k);
    
    // beta l = 1
    OmegaBetal = V.submat(0, 0, kl - 1, kl - 1);
    betal1 = zerok.head(kl);
    for (int i = 0; i < Nl; ++i) {
      ind  = indl1(i);
      ind1 = indl1(i);
      xl = X.submat(ind1, indSubmodel).t() / lambda(ind);
      OmegaBetal += xl * X.submat(ind1, indSubmodel);
      betal1 += xl * Z(ind);
    }
    OmegaBetal = arma::inv_sympd(OmegaBetal);
    betal1 = OmegaBetal * betal1 + arma::chol(OmegaBetal, "lower") * arma::randn(kl);
    
    // beta l = 2
    OmegaBetal = V.submat(0, 0, kl - 1, kl - 1);
    betal2 = zerok.head(kl);
    for (int i = 0; i < Nl; ++i) {
      ind  = indl2(i);
      ind2 = indl2(i);
      xl = X.submat(ind2, indSubmodel).t() / lambda(ind);
      OmegaBetal += xl * X.submat(ind2, indSubmodel);
      betal2 += xl * Z(ind);
    }
    OmegaBetal = arma::inv_sympd(OmegaBetal);
    betal2 = OmegaBetal * betal2 + arma::chol(OmegaBetal, "lower") * arma::randn(kl);
    
    if ((b > 0) && (b % nThin == 0)) {
      keep(b / nThin - 1, arma::span(0, k - 1)) = beta.t();
      keep(b / nThin - 1, arma::span(k, k + kl - 1)) = betal1.t();
      keep(b / nThin - 1, arma::span(k + kl, k + 2 * kl - 1)) = betal2.t();
    }
  }
  
  return keep;
}

// [[Rcpp::export]]
arma::mat predGlmBerKFCV(
    const int B, const int Nweak,
    arma::mat X, const arma::mat X0, const arma::mat X1, arma::mat Xb,
    const arma::mat beta, const arma::mat betal1, const arma::mat betal2,
    const arma::uvec weak, const arma::uvec weak1, const arma::uvec weak2,
    const arma::vec noNA1, const arma::vec noNA2, 
    const arma::uvec Weak1, const arma::uvec Weak2, 
    const arma::uvec indLag1, const arma::uvec indLag2, const arma::uvec indLag12, 
    const arma::uvec indl, const arma::uvec indl1, const arma::uvec indl2, 
    const arma::uvec indSubmodel,
    const arma::vec prob
  ) {
  
  arma::uword ind; 
  arma::uvec ind1, ind2, indAux;
  
  for (int b = 0; b < B; ++b) {
    
    // weak records
    X.submat(Weak1, indLag12) = X1.submat(Weak1, indLag12);
    X.submat(Weak2, indLag12) = X1.submat(Weak2, indLag12);
    for (int i = 0; i < Nweak; i++) {
      ind = weak(i);
      if (R::runif(0, 1) < prob(ind)) {
        if (noNA1(i)) {
          ind1  = weak1(i);
          X.submat(ind1, indLag1) = X1.submat(ind1, indLag1);
        } 
        if (noNA2(i)) {
          ind2  = weak2(i);
          X.submat(ind2, indLag2) = X1.submat(ind2, indLag2);
        } 
      } else {
        if (noNA1(i)) {
          ind1 = weak1(i);
          X.submat(ind1, indLag1)  = X0.submat(ind1, indLag1);
          X.submat(ind1, indLag12) = X0.submat(ind1, indLag12);
        }
        if (noNA2(i)) {
          ind2 = weak2(i);
          X.submat(ind2, indLag2)  = X0.submat(ind2, indLag2);
          X.submat(ind2, indLag12) = X0.submat(ind2, indLag12);
        }
      }
    }
    
    // Xb
    indAux = b;
    Xb.submat(indAux, indl)  = beta.row(b) * X.rows(indl).t();
    Xb.submat(indAux, indl1) = betal1.row(b) * X.submat(indl1, indSubmodel).t();
    Xb.submat(indAux, indl2) = betal2.row(b) * X.submat(indl2, indSubmodel).t();
  }
  
  return Xb;
}

// @description Iterations of the Metropolis-within-Gibbs algorithm for the 
//   logistic regression model with spatial GPs (old version, not used).
// @param N number of observations.
// @param k number of covariates.
// @param Y vector of responses.
// @param X matrix of covariates.
// @param beta vector initial value of beta
// @param aLim,bLim vectors of the limits for the Z's
// @param V matrix prior var-cov of beta
// @param keep matrix of beta's to keep
// @param ... other arguments
// [[Rcpp::export]]
arma::mat glmBerGP(
    const int N, const int k, const arma::mat X, 
    arma::vec beta, const arma::vec abLim,
    arma::vec alpha, arma::vec hp, arma::mat Rinv, arma::uvec site, int n,
    const arma::mat V, arma::mat keep,
    arma::vec psi, arma::vec hpPsi, arma::uvec year, int T,
    bool randomEffectYear, bool randomEffectDay,
    const int nSims, const int nThin, const int nBurnin, const int nReport) {
  
  arma::vec Xb(N), Xl(k);
  arma::vec Z(N, arma::fill::zeros);
  arma::vec auxZ = Z - (X * beta + alpha(site) + psi(year));
  arma::vec lambda = rkolmogorov(N);
  arma::vec lambdaInv(N);
  const arma::vec zerok(k, arma::fill::zeros);
  
  arma::mat OmegaBeta(k, k);
  arma::mat OmegaAlpha(n, n);
  
  int count0, count1; 
  int Nn = N / n;
  
  // prior
  const double na = 0;
  const double nb = 0.0001;
  const double ga = 2;
  const double gb = 1;
  
  // for the GPs
  arma::vec process(n);
  arma::vec onen(n, arma::fill::ones);
  double oneRone = arma::accu(Rinv);
  double delta, chi;

  //arma::uvec indexPsi(N / T);
  arma::uvec indexPsiuvec(N / T);
  arma::umat indexPsi(N / T, T);
  for (int t = 1; t <= T; ++t) {
    indexPsi.col(t - 1) = arma::find(year == t);
  }
  
  for (int b = 1 - nBurnin; b <= nSims; ++b) {
    
    if ((b != 0) && (b % nReport == 0)) {
      Rcpp::Rcout << "The value of b : " << b << "\n";
    }
    
    // lambda and Z
    lambda = MH(N, lambda, auxZ);
    auxZ -= Z;
    Z = abLim % rtnorm(N, -abLim % auxZ, sqrt(lambda));
    auxZ += Z;
    lambdaInv = 1 / lambda;
    
    // beta
    auxZ += X * beta;
    OmegaBeta = V;
    beta = zerok;
    for (int i = 0; i < N; ++i) {
      Xl = X.row(i).t() * lambdaInv(i);
      OmegaBeta += Xl * X.row(i);
      beta += Xl * auxZ(i);
    }
    OmegaBeta = arma::inv_sympd(OmegaBeta);
    beta = OmegaBeta * beta + arma::chol(OmegaBeta, "lower") * arma::randn(k);
    auxZ -= X * beta;
    
    // alpha s
    auxZ += alpha(site);
    OmegaAlpha = hp(1) * Rinv;
    alpha = OmegaAlpha * onen * hp(0);
    count1 = 0;
    // IMPORTANT: data ordered by site (first all s1, then s2, etc.)
    for (int i = 0; i < n; ++i) {
      count0 = count1;
      count1 += Nn;
      OmegaAlpha(i, i) += arma::accu(lambdaInv(arma::span(count0, count1 - 1)));
      alpha(i) += arma::accu(auxZ(arma::span(count0, count1 - 1)) % lambdaInv(arma::span(count0, count1 - 1)));
    }
    OmegaAlpha = arma::inv_sympd(OmegaAlpha);
    alpha = OmegaAlpha * alpha + arma::chol(OmegaAlpha, "lower") * arma::randn(n);
    auxZ -= alpha(site);
    
    // hiperPriors
    //// alpha
    delta = 1 / (oneRone * hp(1) + nb);
    chi   = (onen.t() * Rinv * alpha).eval()(0,0) * hp(1) + na * nb;
    hp(0) = R::rnorm(delta * chi, sqrt(delta));
    
    //// prec
    process = alpha - hp(0);
    hp(1) = R::rgamma(n / 2 + ga,
      1 / ((process.t() * Rinv * process).eval()(0,0) / 2 + gb));
    
    if (randomEffectYear) {
      // psi t
      auxZ += psi(year);
      for (int t = 1; t < T; ++t) {
        indexPsi = arma::find(year == t);
        delta = 1 / (arma::accu(lambdaInv(indexPsi)) + (1 + hpPsi(0) * hpPsi(0)) * hpPsi(1));
        chi   = arma::accu(auxZ(indexPsi) % lambdaInv(indexPsi)) + hpPsi(0) * (psi(t-1) + psi(t+1)) * hpPsi(1);
        psi(t) = R::rnorm(delta * chi, sqrt(delta));
      }
      indexPsi = arma::find(year == T);
      delta = 1 / (arma::accu(lambdaInv(indexPsi)) + hpPsi(1));
      chi   = arma::accu(auxZ(indexPsi) % lambdaInv(indexPsi)) + hpPsi(0) * psi(T-1) * hpPsi(1);
      psi(T) = R::rnorm(delta * chi, sqrt(delta));
      auxZ -= psi(year);
      
      // hiperPriors Psi
      //// rho
      delta = (psi(arma::span(1, T)).t() * psi(arma::span(0, T - 1))).eval()(0);
      chi = arma::accu(pow(psi(arma::span(1, T - 1)), 2));
      hpPsi(0) = rtnorm1(delta / chi, 1 / sqrt(hpPsi(1) * chi), -1, 1);
      
      //// prec
      hpPsi(1) = R::rgamma(T / 2 + ga,
            1 / (arma::accu(pow(psi(arma::span(1, T)) - hpPsi(0) * psi(arma::span(0, T - 1)), 2)) / 2 + gb));
    }
    
    if (randomEffectDay) {
      // psi t
      auxZ += psi(year);
      for (int t = 1; t <= T; ++t) {
        //indexPsi = arma::find(year == t);
        indexPsiuvec = indexPsi.col(t - 1);
        delta = 1 / (arma::accu(lambdaInv(indexPsiuvec)) + hpPsi(1));
        chi   = arma::accu(auxZ(indexPsiuvec) % lambdaInv(indexPsiuvec));
        psi(t) = R::rnorm(delta * chi, sqrt(delta));
      }
      auxZ -= psi(year);
      
      // hiperPriors Psi
      //// prec
      hpPsi(1) = R::rgamma(T / 2 + ga,
            1 / (arma::accu(pow(psi(arma::span(1, T)), 2)) / 2 + gb));
    }
    
      
    if ((b > 0) && (b % nThin == 0)) {
      keep(b / nThin - 1, arma::span(0, k - 1)) = beta.t();
      keep(b / nThin - 1, arma::span(k, k + n - 1)) = alpha.t();
      keep(b / nThin - 1, arma::span(k + n, k + n + T - 1)) = psi(arma::span(1, T)).t();
      keep(b / nThin - 1, arma::span(k + n + T, k + n + T + 1)) = hp.t();
      keep(b / nThin - 1, arma::span(k + n + T + 2, k + n + T + 3)) = hpPsi.t();
    }
  }
  
  return keep;
}

// @description Iterations of the MCMC algorithm for the 
//   logistic regression model with a spatial GP and optionally 
//   daily random effects.
// [[Rcpp::export]]
arma::mat glmBerGPtl(
    const int N, const int k, const int Nl, const int kl, const int Nweak, 
    arma::mat X, const arma::mat X0, const arma::mat X1, 
    arma::vec beta, arma::vec betal1, arma::vec betal2, 
    arma::vec hpl1, arma::vec hpl2,
    arma::vec Ws, arma::vec hp0, 
    arma::mat dist, bool decayPrior, int n,
    arma::vec wtl, int T, int L,
    arma::uvec site, arma::uvec time,
    const arma::uvec weak, const arma::uvec weak1, const arma::uvec weak2,
    const arma::vec noNA1, const arma::vec noNA2, 
    const arma::uvec Weak1, const arma::uvec Weak2, 
    const arma::uvec indLag1, const arma::uvec indLag2, const arma::uvec indLag12, 
    const arma::uvec indl, const arma::uvec indl1, const arma::uvec indl2, 
    const arma::uvec indSubmodel,
    const arma::vec prob, arma::vec abLim,
    const arma::mat V, const double nb, 
    const double ga, const double gb, const double da, const double db,
    arma::mat keep,
    const int nSims, const int nThin, const int nBurnin, const int nReport) {
  
  const int TL2 = T * (L - 2); // T * 363
  const int TIME = arma::unique(time).eval().n_elem; // wtl = 0 -> 1 | wtl = wt -> T * 3 | wtl = wtl -> T * 365
  const int TIME2 = Rcpp::max(Rcpp::NumericVector::create(1, TIME - 2 * T)); // wtl = 0 -> 1 | wtl = wt -> T | wtl = wtl -> T * 363
  const bool randomEffectsTime = (TIME != 1);
  arma::vec Xb(N), Xl(k), xl(kl);
  Xb.rows(indl)  = X.rows(indl)  * beta;
  Xb.rows(indl1) = X.submat(indl1, indSubmodel) * betal1;
  Xb.rows(indl2) = X.submat(indl2, indSubmodel) * betal2;
  
  arma::vec Z(N, arma::fill::zeros); // the first sampled (no need to specify)
  arma::vec auxZ = Z - (Xb + Ws(site) + wtl(time));
  arma::vec lambda = rkolmogorov(N);
  arma::vec lambdaInv = 1 / lambda;
  const arma::vec zerok(k, arma::fill::zeros);
  arma::vec vn(n);
  arma::vec vTIME2(TIME2);
  arma::vec vT(T);
  arma::vec vN(N);
  
  arma::mat OmegaBeta(k, k);
  arma::mat OmegaBetal(kl, kl);
  arma::mat OmegaW(n, n);
  
  // prior
  const double na = 0;
  
  // for the GPs
  arma::vec process(n);
  arma::vec onen(n, arma::fill::ones);
  double delta, chi;
  
  arma::uvec indWUvec(TL2);
  arma::umat indWUmat(TL2, n);
  arma::uvec indWl1Uvec(T);
  arma::umat indWl1Umat(T, n);
  arma::uvec indWl2Uvec(T);
  arma::umat indWl2Umat(T, n);
  for (int i = 0; i < n; ++i) { 
    indWUmat.col(i) = arma::find(site == i + 2 * n);
    indWl1Umat.col(i) = arma::find(site == i);
    indWl2Umat.col(i) = arma::find(site == i + n);
  }
  
  arma::uvec indwUvec(n * TL2 / TIME2);
  arma::umat indwUmat(n * TL2 / TIME2, TIME2);
  arma::uvec indwl1Uvec(n);
  arma::umat indwl1Umat(n, T);
  arma::uvec indwl2Uvec(n);
  arma::umat indwl2Umat(n, T);
  if (randomEffectsTime) {
    for (int t = 0; t < TIME2; ++t) { 
      indwUmat.col(t) = arma::find(time == t + 2 * T);
    }
    for (int t = 0; t < T; ++t) { 
      indwl1Umat.col(t) = arma::find(time == t);
      indwl2Umat.col(t) = arma::find(time == t + T);
    }
  }
  
  // decay (starter)
  int accept = 0;
  int total = 0;
  double r;
  double sd = 1;
  double lsd = 0; //log(sd);
  
  arma::mat Rinv = arma::inv_sympd(exp(- hp0(2) * dist));
  double Rlogdet = arma::log_det_sympd(Rinv);
  double oneRone = arma::accu(Rinv);
  
  double decay_aux, ldecay_aux;
  double ldecay = log(hp0(2));
  arma::mat Rinv_aux(n, n);
  double Rlogdet_aux;
  
  arma::mat ZtRZ(3, 1), ZtRZ_aux(3, 1); // l = 1, l = 2, others
  double A;
  
  arma::uword ind; 
  arma::uvec ind1, ind2;
  
  // DATA is ordered!!! All data in s1 first, all data in l = 1 first, move from t
  
  for (int b = 1 - nBurnin; b <= nSims; ++b) {
    
    // check
    if (b % nReport == 0) {
      Rcpp::Rcout << "The value of b : " << b << "\n";
      if ((b < 1) && (total != 0)) {
        Rcpp::Rcout << "Acceptance rate : " << r << "\n";
      } else if (b < 1) {
        Rcpp::Rcout << "Acceptance rate : " << 0 << "\n";
      } else {
        Rcpp::Rcout << "Acceptance rate : " <<  (double)accept / (double)b << "\n";
      }
    }
    
    // weak records
    auxZ += Xb;
    X.submat(Weak1, indLag12) = X1.submat(Weak1, indLag12);
    X.submat(Weak2, indLag12) = X1.submat(Weak2, indLag12);
    for (int i = 0; i < Nweak; i++) {
      ind = weak(i);
      if (R::runif(0, 1) < prob(ind)) {
        abLim(ind) = 1;
        if (noNA1(i)) {
          ind1  = weak1(i);
          X.submat(ind1, indLag1) = X1.submat(ind1, indLag1);
        } 
        if (noNA2(i)) {
          ind2  = weak2(i);
          X.submat(ind2, indLag2) = X1.submat(ind2, indLag2);
        } 
      } else {
        abLim(ind) = -1;
        if (noNA1(i)) {
          ind1 = weak1(i);
          X.submat(ind1, indLag1)  = X0.submat(ind1, indLag1);
          X.submat(ind1, indLag12) = X0.submat(ind1, indLag12);
        }
        if (noNA2(i)) {
          ind2 = weak2(i);
          X.submat(ind2, indLag2)  = X0.submat(ind2, indLag2);
          X.submat(ind2, indLag12) = X0.submat(ind2, indLag12);
        }
      }
    }
    Xb.rows(indl)  = X.rows(indl)  * beta;
    Xb.rows(indl1) = X.submat(indl1, indSubmodel) * betal1;
    Xb.rows(indl2) = X.submat(indl2, indSubmodel) * betal2;
    auxZ -= Xb;
    
    // lambda and Z
    auxZ -= Z;
    Z = abLim % rtnorm(N, -abLim % auxZ, sqrt(lambda));
    auxZ += Z;
    lambda = MH(N, lambda, auxZ);
    lambdaInv = 1 / lambda;
    
    // beta
    auxZ += Xb;
    OmegaBeta = V;
    beta = zerok;
    for (int i = 0; i < (N - 2 * Nl); ++i) {
      ind = indl(i);
      Xl = X.row(ind).t() * lambdaInv(ind);
      OmegaBeta += Xl * X.row(ind);
      beta += Xl * auxZ(ind);
    }
    OmegaBeta = arma::inv_sympd(OmegaBeta);
    beta = OmegaBeta * beta + arma::chol(OmegaBeta, "lower") * arma::randn(k);
    
    // beta l = 1
    OmegaBetal = V.submat(0, 0, kl - 1, kl - 1);
    betal1 = zerok.head(kl);
    for (int i = 0; i < Nl; ++i) {
      ind  = indl1(i);
      ind1 = indl1(i);
      xl = X.submat(ind1, indSubmodel).t() * lambdaInv(ind);
      OmegaBetal += xl * X.submat(ind1, indSubmodel);
      betal1 += xl * auxZ(ind);
    }
    OmegaBetal = arma::inv_sympd(OmegaBetal);
    betal1 = OmegaBetal * betal1 + arma::chol(OmegaBetal, "lower") * arma::randn(kl);
    
    // beta l = 2
    OmegaBetal = V.submat(0, 0, kl - 1, kl - 1);
    betal2 = zerok.head(kl);
    for (int i = 0; i < Nl; ++i) {
      ind  = indl2(i);
      ind2 = indl2(i);
      xl = X.submat(ind2, indSubmodel).t() * lambdaInv(ind);
      OmegaBetal += xl * X.submat(ind2, indSubmodel);
      betal2 += xl * auxZ(ind);
    }
    OmegaBetal = arma::inv_sympd(OmegaBetal);
    betal2 = OmegaBetal * betal2 + arma::chol(OmegaBetal, "lower") * arma::randn(kl);
    Xb.rows(indl)  = X.rows(indl) * beta;
    Xb.rows(indl1) = X.submat(indl1, indSubmodel) * betal1;
    Xb.rows(indl2) = X.submat(indl2, indSubmodel) * betal2;
    auxZ -= Xb;
    
    // W(s)
    //// IMPORTANT: data ordered by site, day, year (first all s1, then s2, etc.)
    auxZ += Ws(site);
    vN = auxZ % lambdaInv;
    
    //// l = 3,...,L
    OmegaW = hp0(1) * Rinv;
    Ws(arma::span(2 * n, 3 * n - 1)) = OmegaW * onen * hp0(0);
    for (int i = 0; i < n; ++i) {
      indWUvec = indWUmat.col(i);
      OmegaW(i, i)  += arma::accu(lambdaInv(indWUvec));
      Ws(i + 2 * n) += arma::accu(vN(indWUvec));
    }
    OmegaW = arma::inv_sympd(OmegaW);
    Ws(arma::span(2 * n, 3 * n - 1)) = 
      OmegaW * Ws(arma::span(2 * n, 3 * n - 1)) + arma::chol(OmegaW, "lower") * arma::randn(n);
      
    //// l = 1
    OmegaW = hpl1(1) * Rinv;
    Ws(arma::span(0, n - 1)) = OmegaW * onen * hpl1(0);
    for (int i = 0; i < n; ++i) {
      indWl1Uvec = indWl1Umat.col(i);
      OmegaW(i, i) += arma::accu(lambdaInv(indWl1Uvec));
      Ws(i)        += arma::accu(vN(indWl1Uvec));
    }
    OmegaW = arma::inv_sympd(OmegaW);
    Ws(arma::span(0, n - 1)) = 
      OmegaW * Ws(arma::span(0, n - 1)) + arma::chol(OmegaW, "lower") * arma::randn(n);
      
    //// l = 2
    OmegaW = hpl2(1) * Rinv;
    Ws(arma::span(n, 2 * n - 1)) = OmegaW * onen * hpl2(0);
    for (int i = 0; i < n; ++i) {
      indWl2Uvec = indWl2Umat.col(i);
      OmegaW(i, i) += arma::accu(lambdaInv(indWl2Uvec));
      Ws(i + n)    += arma::accu(vN(indWl2Uvec));
    }
    OmegaW = arma::inv_sympd(OmegaW);
    Ws(arma::span(n, 2 * n - 1)) = 
      OmegaW * Ws(arma::span(n, 2 * n - 1)) + arma::chol(OmegaW, "lower") * arma::randn(n);
    
    auxZ -= Ws(site);
    
    if (randomEffectsTime) {
      // w time
      auxZ += wtl(time);
      vN = auxZ % lambdaInv;
      
      //// l = 3,...,L
      for (int t = 0; t < TIME2; ++t) {
        indwUvec = indwUmat.col(t);
        delta = 1 / (arma::accu(lambdaInv(indwUvec)) + hp0(3));
        chi   = arma::accu(vN(indwUvec)); // + 0 * hp0(3);
        wtl(t + 2 * T) = R::rnorm(delta * chi, sqrt(delta));
      }
      
      //// l = 1
      for (int t = 0; t < T; ++t) {
        indwl1Uvec = indwl1Umat.col(t);
        delta = 1 / (arma::accu(lambdaInv(indwl1Uvec)) + hpl1(2));
        chi   = arma::accu(vN(indwl1Uvec)); // + 0 * hpl1(2);
        wtl(t) = R::rnorm(delta * chi, sqrt(delta));
      }
      
      //// l = 2
      for (int t = 0; t < T; ++t) {
        indwl2Uvec = indwl2Umat.col(t);
        delta = 1 / (arma::accu(lambdaInv(indwl2Uvec)) + hpl2(2));
        chi   = arma::accu(vN(indwl2Uvec)); // + 0 * hpl2(2);
        wtl(t + T) = R::rnorm(delta * chi, sqrt(delta));
      }
      
      auxZ -= wtl(time);
      
      //// prec1
      vTIME2 = wtl(arma::span(2 * T, TIME - 1));
      hp0(3) = R::rgamma(TIME2 / 2 + ga,
          1 / (arma::dot(vTIME2, vTIME2) / 2 + gb));
      
      //// prec1 l = 1
      vT = wtl(arma::span(0, T - 1));
      hpl1(2) = R::rgamma(T / 2 + ga,
           1 / (arma::dot(vT, vT) / 2 + gb));
      
      //// prec1 l = 2
      vT = wtl(arma::span(T, 2 * T - 1));
      hpl2(2) = R::rgamma(T / 2 + ga,
           1 / (arma::dot(vT, vT) / 2 + gb));
    }
    
    // hiperPriors
    //// beta0
    delta = 1 / (oneRone * hp0(1) + nb);
    chi   = arma::as_scalar(onen.t() * Rinv * Ws(arma::span(2 * n, 3 * n - 1))) * hp0(1) + na * nb;
    hp0(0) = R::rnorm(delta * chi, sqrt(delta));
    
    //// beta0 l = 1
    delta = 1 / (oneRone * hpl1(1) + nb);
    chi   = arma::as_scalar(onen.t() * Rinv * Ws(arma::span(0, n - 1))) * hpl1(1) + na * nb;
    hpl1(0) = R::rnorm(delta * chi, sqrt(delta));
    
    //// beta0 l = 2
    delta = 1 / (oneRone * hpl2(1) + nb);
    chi   = arma::as_scalar(onen.t() * Rinv * Ws(arma::span(n, 2 * n - 1))) * hpl2(1) + na * nb;
    hpl2(0) = R::rnorm(delta * chi, sqrt(delta));
    
    //// decay0
    if (decayPrior) {
      decay_aux   = rtnorm1(hp0(2), sd, da, db);
      Rinv_aux    = arma::inv_sympd(exp(- decay_aux * dist));
      Rlogdet_aux = arma::log_det_sympd(Rinv_aux);
      
      vn = Ws(arma::span(2 * n, 3 * n - 1)) - hp0(0);
      ZtRZ_aux.row(2) = vn.t() * Rinv_aux * vn;
      ZtRZ.row(2)     = vn.t() * Rinv * vn;
      
      vn = Ws(arma::span(0, n - 1)) - hpl1(0);
      ZtRZ_aux.row(0) = vn.t() * Rinv_aux * vn;
      ZtRZ.row(0)     = vn.t() * Rinv * vn;
      
      vn = Ws(arma::span(n, 2 * n - 1)) - hpl2(0);
      ZtRZ_aux.row(1) = vn.t() * Rinv_aux * vn;
      ZtRZ.row(1)     = vn.t() * Rinv * vn;
      
      A = 
        (3 * Rlogdet_aux - hp0(1) * ZtRZ_aux(2, 0) - hpl1(1) * ZtRZ_aux(0, 0) - hpl2(1) * ZtRZ_aux(1, 0)) / 2 - 
        (3 * Rlogdet - hp0(1) * ZtRZ(2, 0) - hpl1(1) * ZtRZ(0, 0) - hpl2(1) * ZtRZ(1, 0)) / 2 + 
        dtnorm1(hp0(2), decay_aux, sd, da, db) - 
        dtnorm1(decay_aux, hp0(2), sd, da, db);
    } else {
      ldecay_aux  = R::rnorm(ldecay, sd);
      decay_aux   = exp(ldecay_aux);
      Rinv_aux    = arma::inv_sympd(exp(- decay_aux * dist));
      Rlogdet_aux = arma::log_det_sympd(Rinv_aux);
      
      vn = Ws(arma::span(2 * n, 3 * n - 1)) - hp0(0);
      ZtRZ_aux.row(2) = vn.t() * Rinv_aux * vn;
      ZtRZ.row(2)     = vn.t() * Rinv * vn;
      
      vn = Ws(arma::span(0, n - 1)) - hpl1(0);
      ZtRZ_aux.row(0) = vn.t() * Rinv_aux * vn;
      ZtRZ.row(0)     = vn.t() * Rinv * vn;
      
      vn = Ws(arma::span(n, 2 * n - 1)) - hpl2(0);
      ZtRZ_aux.row(1) = vn.t() * Rinv_aux * vn;
      ZtRZ.row(1)     = vn.t() * Rinv * vn;
      
      A = 
        (3 * Rlogdet_aux - hp0(1) * ZtRZ_aux(2, 0) - hpl1(1) * ZtRZ_aux(0, 0) - hpl2(1) * ZtRZ_aux(1, 0)) / 2 + 
        da * ldecay_aux - db * decay_aux - 
        ((3 * Rlogdet - hp0(1) * ZtRZ(2, 0) - hpl1(1) * ZtRZ(0, 0) - hpl2(1) * ZtRZ(1, 0)) / 2 +
        da * ldecay - db * hp0(2));
    }
    if (log(R::runif(0, 1)) <= A) {
      ++accept;
      hp0(2) = decay_aux;
      ldecay = ldecay_aux;
      Rinv = Rinv_aux;
      Rlogdet = Rlogdet_aux;
      oneRone = arma::accu(Rinv);
      ZtRZ = ZtRZ_aux;
    }
    // tune sd of the proposal for decay
    if ((b < 1) && (++total % 25 == 0)) {
      r = (double)accept / (double)total;
      if (r > 0.33) {
        lsd += 1 / sqrt((nBurnin + b) / 25);
      } else {
        lsd -= 1 / sqrt((nBurnin + b) / 25);
      }
      sd = exp(lsd);
      //accept = 0;
      //total = 0;
      Rcpp::Rcout << "Tuning now sd : " << sd << "\n";
    } else if (b == 0) {
      accept = 0;
      total = 0;
    }
    
    //// prec0
    hp0(1) = R::rgamma(n / 2 + ga, 
        1 / (ZtRZ(2, 0) / 2 + gb));
    
    //// prec0 l = 1
    hpl1(1) = R::rgamma(n / 2 + ga, 
         1 / (ZtRZ(0, 0) / 2 + gb));
    
    //// prec0 l = 2
    hpl2(1) = R::rgamma(n / 2 + ga, 
         1 / (ZtRZ(1, 0) / 2 + gb));
    
    // keep
    if ((b > 0) && (b % nThin == 0)) {
      keep(b / nThin - 1, arma::span(0,                             k - 1)) = beta.t();
      keep(b / nThin - 1, arma::span(k,                             k + n * 3 - 1)) = Ws.t();
      if (randomEffectsTime) {
        keep(b / nThin - 1, arma::span(k + n * 3,                   k + n * 3 + TIME - 1)) = wtl.t();
      }
      keep(b / nThin - 1, arma::span(k + n * 3 + TIME,              k + n * 3 + TIME + 3)) = hp0.t();
      keep(b / nThin - 1, arma::span(k + n * 3 + TIME + 4,          k + n * 3 + TIME + 3 + kl)) = betal1.t();
      keep(b / nThin - 1, arma::span(k + n * 3 + TIME + 4 + kl,     k + n * 3 + TIME + 6 + kl)) = hpl1.t();
      keep(b / nThin - 1, arma::span(k + n * 3 + TIME + 7 + kl,     k + n * 3 + TIME + 6 + 2 * kl)) = betal2.t();
      keep(b / nThin - 1, arma::span(k + n * 3 + TIME + 7 + 2 * kl, k + n * 3 + TIME + 9 + 2 * kl)) = hpl2.t();
    }
  }
  
  return keep;
}

// [[Rcpp::export]]
arma::mat predGlmBerGPtl(
    const int B, const int Nweak, 
    const int n, const int newn, const int T, const int L,
    arma::mat X, const arma::mat X0, const arma::mat X1, 
    arma::uvec newsite,
    arma::mat dist, arma::mat Xb, arma::mat Ws0,
    const arma::mat beta, const arma::mat betal1, const arma::mat betal2,
    const arma::mat Ws, const arma::mat hp,
    const arma::uvec weak, const arma::uvec weak1, const arma::uvec weak2,
    const arma::vec noNA1, const arma::vec noNA2, 
    const arma::uvec Weak1, const arma::uvec Weak2, 
    const arma::uvec indLag1, const arma::uvec indLag2, const arma::uvec indLag12, 
    const arma::uvec indl, const arma::uvec indl1, const arma::uvec indl2, 
    const arma::uvec indSubmodel,
    const arma::vec prob
) {
  
  arma::mat Sigma11(newn, newn);
  arma::mat Sigma12(newn, n);
  arma::mat Sigma22inv(n, n);
  arma::mat SigmaAux(n, newn); // transposed
  arma::mat Sigma(newn, newn);
  
  arma::uword ind; 
  arma::uvec ind1, ind2, indAux;
  
  for (int b = 0; b < B; ++b) {
    
    // weak records
    X.submat(Weak1, indLag12) = X1.submat(Weak1, indLag12);
    X.submat(Weak2, indLag12) = X1.submat(Weak2, indLag12);
    for (int i = 0; i < Nweak; i++) {
      ind = weak(i);
      if (R::runif(0, 1) < prob(ind)) {
        if (noNA1(i)) {
          ind1  = weak1(i);
          X.submat(ind1, indLag1) = X1.submat(ind1, indLag1);
        } 
        if (noNA2(i)) {
          ind2  = weak2(i);
          X.submat(ind2, indLag2) = X1.submat(ind2, indLag2);
        } 
      } else {
        if (noNA1(i)) {
          ind1 = weak1(i);
          X.submat(ind1, indLag1)  = X0.submat(ind1, indLag1);
          X.submat(ind1, indLag12) = X0.submat(ind1, indLag12);
        }
        if (noNA2(i)) {
          ind2 = weak2(i);
          X.submat(ind2, indLag2)  = X0.submat(ind2, indLag2);
          X.submat(ind2, indLag12) = X0.submat(ind2, indLag12);
        }
      }
    }
    
    // Xb
    indAux = b;
    Xb.submat(indAux, indl)  = beta.row(b) * X.rows(indl).t();
    Xb.submat(indAux, indl1) = betal1.row(b) * X.submat(indl1, indSubmodel).t();
    Xb.submat(indAux, indl2) = betal2.row(b) * X.submat(indl2, indSubmodel).t();
    
    // Ws
    Sigma11    = exp(- hp(b, 2) * dist.submat(0, 0, newn - 1, newn - 1));
    Sigma12    = exp(- hp(b, 2) * dist.submat(0, newn, newn - 1, newn + n - 1));
    Sigma22inv = arma::inv_sympd(exp(- hp(b, 2) * dist.submat(newn, newn, newn + n - 1, newn + n - 1)));
    SigmaAux   = (Sigma12 * Sigma22inv).t();
    Sigma      = Sigma11 - SigmaAux.t() * Sigma12.t();
    Sigma      = arma::chol(Sigma, "lower");
    
    // l = 3,...,L
    Ws0(b, arma::span(2 * newn, 3 * newn - 1)) = 
      hp(b, 0) + 
      (Ws(b, arma::span(2 * n, 3 * n - 1)) - hp(b, 0)) * SigmaAux +
      (Sigma * arma::randn(newn)).t() / sqrt(hp(b, 1));

    // l = 1, 2 (0, 1)
    Ws0(b, arma::span(0, newn - 1)) = 
      hp(b, 3) + 
      (Ws(b, arma::span(0, n - 1)) - hp(b, 3)) * SigmaAux +
      (Sigma * arma::randn(newn)).t() / sqrt(hp(b, 4));
    
    Ws0(b, arma::span(newn, 2 * newn - 1)) = 
      hp(b, 5) + 
      (Ws(b, arma::span(n, 2 * n - 1)) - hp(b, 5)) * SigmaAux +
      (Sigma * arma::randn(newn)).t() / sqrt(hp(b, 6));
  }
  
  return Xb + Ws0.cols(newsite);
}

// @description Iterations of the MCMC algorithm for the 
//   specific logistic regression model in the PAPER!!!
// [[Rcpp::export]]
arma::mat glmBerModelPaper(
    const int N, const int k, const int Nl, const int kl, const int Nweak, 
    arma::mat X, const arma::mat X0, const arma::mat X1, 
    arma::vec beta, arma::vec betal1, arma::vec betal2, 
    arma::vec hpl1, arma::vec hpl2,
    arma::vec Wtls, arma::vec hp0, 
    arma::mat dist, bool decayPrior, int n,
    arma::vec wtl, int T, int L,
    arma::uvec site, arma::uvec year, arma::uvec day,
    const arma::uvec weak, const arma::uvec weak1, const arma::uvec weak2,
    const arma::vec noNA1, const arma::vec noNA2, 
    const arma::uvec Weak1, const arma::uvec Weak2, 
    const arma::uvec indLag1, const arma::uvec indLag2, const arma::uvec indLag12, 
    const arma::uvec indl, const arma::uvec indl1, const arma::uvec indl2, 
    const arma::uvec indSubmodel,
    const arma::vec prob, arma::vec abLim,
    const arma::mat V, const double nb, 
    const double ga, const double gb, const double da, const double db,
    arma::mat keep,
    const int nSims, const int nThin, const int nBurnin, const int nReport) {
  
  const int TL = T * L;
  arma::vec Xb(N), Xl(k), xl(kl);
  Xb.rows(indl)  = X.rows(indl)  * beta;
  Xb.rows(indl1) = X.submat(indl1, indSubmodel) * betal1;
  Xb.rows(indl2) = X.submat(indl2, indSubmodel) * betal2;
  
  arma::vec Z(N, arma::fill::zeros); // the first sampled (no need to specify)
  arma::vec auxZ = Z - (Xb + Wtls);
  arma::vec lambda = rkolmogorov(N);
  arma::vec lambdaInv = 1 / lambda;
  const arma::vec zerok(k, arma::fill::zeros);
  arma::vec vn(n);
  arma::vec vTL2(T * (L - 2));
  arma::vec vT(T);
  arma::vec vN(N);
  
  arma::mat OmegaBeta(k, k);
  arma::mat OmegaBetal(kl, kl);
  arma::mat OmegaW(n, n);
  
  // prior
  const double na = 0;
  
  // for the GPs
  arma::vec process(n);
  arma::vec onen(n, arma::fill::ones);
  double delta, chi;
  
  arma::uvec indWUvec(n);
  arma::umat indWUmat(n, TL);
  for (int l = 0; l < L; ++l) { 
    for (int t = 0; t < T; ++t) {
      indWUmat.col(l * T + t) = arma::find((year == t) && (day == l));
    }
  }
  
  // decay (starter)
  int accept = 0;
  int total = 0;
  double r;
  double sd = 1;
  double lsd = 0; //log(sd);
  
  arma::mat Rinv = arma::inv_sympd(exp(- hp0(2) * dist));
  double Rlogdet = arma::log_det_sympd(Rinv);
  double oneRone = arma::accu(Rinv);
  
  double decay_aux, ldecay_aux;
  double ldecay = log(hp0(2));
  arma::mat Rinv_aux(n, n);
  double Rlogdet_aux;
  
  arma::mat ZtRZ(1, 1), ZtRZ_aux(1, 1), ZtRZl(2, 1), ZtRZl_aux(2, 1);
  double A;
  
  arma::uword ind; 
  arma::uvec ind1, ind2;
  
  // DATA is ordered!!! All data in s1 first, all data in l = 1 first, move from t
  
  for (int b = 1 - nBurnin; b <= nSims; ++b) {
    
    // check
    if (b % nReport == 0) {
      Rcpp::Rcout << "The value of b : " << b << "\n";
      if ((b < 1) && (total != 0)) {
        Rcpp::Rcout << "Acceptance rate : " << r << "\n";
      } else if (b < 1) {
        Rcpp::Rcout << "Acceptance rate : " << 0 << "\n";
      } else {
        Rcpp::Rcout << "Acceptance rate : " <<  (double)accept / (double)b << "\n";
      }
    }
    
    // weak records
    auxZ += Xb;
    X.submat(Weak1, indLag12) = X1.submat(Weak1, indLag12);
    X.submat(Weak2, indLag12) = X1.submat(Weak2, indLag12);
    for (int i = 0; i < Nweak; i++) {
      ind = weak(i);
      if (R::runif(0, 1) < prob(ind)) {
        abLim(ind) = 1;
        if (noNA1(i)) {
          ind1  = weak1(i);
          X.submat(ind1, indLag1) = X1.submat(ind1, indLag1);
        } 
        if (noNA2(i)) {
          ind2  = weak2(i);
          X.submat(ind2, indLag2) = X1.submat(ind2, indLag2);
        } 
      } else {
        abLim(ind) = -1;
        if (noNA1(i)) {
          ind1 = weak1(i);
          X.submat(ind1, indLag1)  = X0.submat(ind1, indLag1);
          X.submat(ind1, indLag12) = X0.submat(ind1, indLag12);
        }
        if (noNA2(i)) {
          ind2 = weak2(i);
          X.submat(ind2, indLag2)  = X0.submat(ind2, indLag2);
          X.submat(ind2, indLag12) = X0.submat(ind2, indLag12);
        }
      }
    }
    Xb.rows(indl)  = X.rows(indl)  * beta;
    Xb.rows(indl1) = X.submat(indl1, indSubmodel) * betal1;
    Xb.rows(indl2) = X.submat(indl2, indSubmodel) * betal2;
    auxZ -= Xb;
    
    // lambda and Z
    auxZ -= Z;
    Z = abLim % rtnorm(N, -abLim % auxZ, sqrt(lambda));
    auxZ += Z;
    lambda = MH(N, lambda, auxZ);
    lambdaInv = 1 / lambda;
    
    // beta
    auxZ += Xb;
    OmegaBeta = V;
    beta = zerok;
    for (int i = 0; i < (N - 2 * Nl); ++i) {
      ind = indl(i);
      Xl = X.row(ind).t() * lambdaInv(ind);
      OmegaBeta += Xl * X.row(ind);
      beta += Xl * auxZ(ind);
    }
    OmegaBeta = arma::inv_sympd(OmegaBeta);
    beta = OmegaBeta * beta + arma::chol(OmegaBeta, "lower") * arma::randn(k);

    // beta l = 1
    OmegaBetal = V.submat(0, 0, kl - 1, kl - 1);
    betal1 = zerok.head(kl);
    for (int i = 0; i < Nl; ++i) {
      ind  = indl1(i);
      ind1 = indl1(i);
      xl = X.submat(ind1, indSubmodel).t() * lambdaInv(ind);
      OmegaBetal += xl * X.submat(ind1, indSubmodel);
      betal1 += xl * auxZ(ind);
    }
    OmegaBetal = arma::inv_sympd(OmegaBetal);
    betal1 = OmegaBetal * betal1 + arma::chol(OmegaBetal, "lower") * arma::randn(kl);
    
    // beta l = 2
    OmegaBetal = V.submat(0, 0, kl - 1, kl - 1);
    betal2 = zerok.head(kl);
    for (int i = 0; i < Nl; ++i) {
      ind  = indl2(i);
      ind2 = indl2(i);
      xl = X.submat(ind2, indSubmodel).t() * lambdaInv(ind);
      OmegaBetal += xl * X.submat(ind2, indSubmodel);
      betal2 += xl * auxZ(ind);
    }
    OmegaBetal = arma::inv_sympd(OmegaBetal);
    betal2 = OmegaBetal * betal2 + arma::chol(OmegaBetal, "lower") * arma::randn(kl);
    Xb.rows(indl)  = X.rows(indl)  * beta;
    Xb.rows(indl1) = X.submat(indl1, indSubmodel) * betal1;
    Xb.rows(indl2) = X.submat(indl2, indSubmodel) * betal2;
    auxZ -= Xb;
    
    // Wtl(s)
    // IMPORTANT: data ordered by site, day, year (first all s1, then s2, etc.)
    auxZ += Wtls;
    vN = auxZ % lambdaInv;
    
    // l = 3,...,L
    delta = 1 / (oneRone * hp0(1) + hp0(3));
    for (int l = 2; l < L; ++l) {
      for (int t = 0; t < T; ++t) {
        indWUvec = indWUmat.col(l * T + t);
        // W tl s
        OmegaW = hp0(1) * Rinv;
        Wtls(indWUvec) = OmegaW * onen * wtl(l * T + t) + vN(indWUvec);
        OmegaW.diag() += lambdaInv(indWUvec);
        OmegaW = arma::inv_sympd(OmegaW);
        Wtls(indWUvec) = OmegaW * Wtls(indWUvec) + arma::chol(OmegaW, "lower") * arma::randn(n);
        
        // w tl
        chi   = arma::as_scalar(onen.t() * Rinv * Wtls(indWUvec)) * hp0(1) + hp0(0) * hp0(3);
        wtl(l * T + t) = R::rnorm(delta * chi, sqrt(delta));
      }
    }
    
    // l = 1 (0)
    delta = 1 / (oneRone * hpl1(1) + hpl1(2));
    for (int t = 0; t < T; ++t) {
      indWUvec = indWUmat.col(t);
      // W tl s
      OmegaW = hpl1(1) * Rinv;
      Wtls(indWUvec) = OmegaW * onen * wtl(t) + vN(indWUvec);
      OmegaW.diag() += lambdaInv(indWUvec);
      OmegaW = arma::inv_sympd(OmegaW);
      Wtls(indWUvec) = OmegaW * Wtls(indWUvec) + arma::chol(OmegaW, "lower") * arma::randn(n);
        
      // w tl
      chi   = arma::as_scalar(onen.t() * Rinv * Wtls(indWUvec)) * hpl1(1) + hpl1(0) * hpl1(2);
      wtl(t) = R::rnorm(delta * chi, sqrt(delta));
    }
    
    // l = 2 (1)
    delta = 1 / (oneRone * hpl2(1) + hpl2(2));
    for (int t = 0; t < T; ++t) {
      indWUvec = indWUmat.col(T + t);
      // W tl s
      OmegaW = hpl2(1) * Rinv;
      Wtls(indWUvec) = OmegaW * onen * wtl(T + t) + vN(indWUvec);
      OmegaW.diag() += lambdaInv(indWUvec);
      OmegaW = arma::inv_sympd(OmegaW);
      Wtls(indWUvec) = OmegaW * Wtls(indWUvec) + arma::chol(OmegaW, "lower") * arma::randn(n);
      
      // w tl
      chi   = arma::as_scalar(onen.t() * Rinv * Wtls(indWUvec)) * hpl2(1) + hpl2(0) * hpl2(2);
      wtl(T + t) = R::rnorm(delta * chi, sqrt(delta));
    }
    
    auxZ -= Wtls;
    
    // hiperPriors
    //// beta0
    delta = 1 / (T * (L - 2) * hp0(3) + nb);
    chi   = arma::accu(wtl(arma::span(2 * T, TL - 1))) * hp0(3) + na * nb;
    hp0(0) = R::rnorm(delta * chi, sqrt(delta));
    
    //// beta0 l = 1
    delta = 1 / (T * hpl1(2) + nb);
    chi   = arma::accu(wtl(arma::span(0, T - 1))) * hpl1(2) + na * nb;
    hpl1(0) = R::rnorm(delta * chi, sqrt(delta));
    
    //// beta0 l = 2
    delta = 1 / (T * hpl2(2) + nb);
    chi   = arma::accu(wtl(arma::span(T, 2 * T - 1))) * hpl2(2) + na * nb;
    hpl2(0) = R::rnorm(delta * chi, sqrt(delta));
    
    //// decay0
    if (decayPrior) {
      decay_aux   = rtnorm1(hp0(2), sd, da, db);
      Rinv_aux    = arma::inv_sympd(exp(- decay_aux * dist));
      Rlogdet_aux = arma::log_det_sympd(Rinv_aux);
      ZtRZ_aux.zeros(); 
      ZtRZ.zeros(); 
      ZtRZl_aux.zeros(); 
      ZtRZl.zeros(); 
      for (int l = 2; l < L; ++l) {
        for (int t = 0; t < T; ++t) {
          indWUvec = indWUmat.col(l * T + t);
          vn = Wtls(indWUvec) - wtl(l * T + t);
          ZtRZ_aux += vn.t() * Rinv_aux * vn;
          ZtRZ     += vn.t() * Rinv * vn;
        }
      }
      for (int l = 0; l < 2; ++l) {
        for (int t = 0; t < T; ++t) {
          indWUvec = indWUmat.col(l * T + t);
          vn = Wtls(indWUvec) - wtl(l * T + t);
          ZtRZl_aux.row(l) += vn.t() * Rinv_aux * vn;
          ZtRZl.row(l)     += vn.t() * Rinv * vn;
        }
      }
      A = 
        (TL * Rlogdet_aux - hp0(1) * arma::as_scalar(ZtRZ_aux) - hpl1(1) * ZtRZl_aux(0, 0) - hpl2(1) * ZtRZl_aux(1, 0)) / 2 - 
        (TL * Rlogdet - hp0(1) * arma::as_scalar(ZtRZ) - hpl1(1) * ZtRZl(0, 0) - hpl2(1) * ZtRZl(1, 0)) / 2 + 
        dtnorm1(hp0(2), decay_aux, sd, da, db) - 
        dtnorm1(decay_aux, hp0(2), sd, da, db);
    } else {
      ldecay_aux  = R::rnorm(ldecay, sd);
      decay_aux   = exp(ldecay_aux);
      Rinv_aux    = arma::inv_sympd(exp(- decay_aux * dist));
      Rlogdet_aux = arma::log_det_sympd(Rinv_aux);
      ZtRZ_aux.zeros(); 
      ZtRZ.zeros(); 
      ZtRZl_aux.zeros(); 
      ZtRZl.zeros(); 
      for (int l = 2; l < L; ++l) {
        for (int t = 0; t < T; ++t) {
          indWUvec = indWUmat.col(l * T + t);
          vn = Wtls(indWUvec) - wtl(l * T + t);
          ZtRZ_aux += vn.t() * Rinv_aux * vn;
          ZtRZ     += vn.t() * Rinv * vn;
        }
      }
      for (int l = 0; l < 2; ++l) {
        for (int t = 0; t < T; ++t) {
          indWUvec = indWUmat.col(l * T + t);
          vn = Wtls(indWUvec) - wtl(l * T + t);
          ZtRZl_aux.row(l) += vn.t() * Rinv_aux * vn;
          ZtRZl.row(l)     += vn.t() * Rinv * vn;
        }
      }
      A = 
        (TL * Rlogdet_aux - hp0(1) * arma::as_scalar(ZtRZ_aux) - hpl1(1) * ZtRZl_aux(0, 0) - hpl2(1) * ZtRZl_aux(1, 0)) / 2 + 
        da * ldecay_aux - db * decay_aux - 
        ((TL * Rlogdet - hp0(1) * arma::as_scalar(ZtRZ) - hpl1(1) * ZtRZl(0, 0) - hpl2(1) * ZtRZl(1, 0)) / 2 +
        da * ldecay - db * hp0(2));
    }
    if (log(R::runif(0, 1)) <= A) {
      ++accept;
      hp0(2) = decay_aux;
      ldecay = ldecay_aux;
      Rinv = Rinv_aux;
      Rlogdet = Rlogdet_aux;
      oneRone = arma::accu(Rinv);
      ZtRZ = ZtRZ_aux;
      ZtRZl = ZtRZl_aux;
    }
    // tune sd of the proposal for decay
    if ((b < 1) && (++total % 25 == 0)) {
      r = (double)accept / (double)total;
      if (r > 0.33) {
        lsd += 1 / sqrt((nBurnin + b) / 25);
      } else {
        lsd -= 1 / sqrt((nBurnin + b) / 25);
      }
      sd = exp(lsd);
      //accept = 0;
      //total = 0;
      Rcpp::Rcout << "Tuning now sd : " << sd << "\n";
    } else if (b == 0) {
      accept = 0;
      total = 0;
    }
    
    //// prec0
    hp0(1) = R::rgamma(n * T * (L - 2) / 2 + ga, 
      1 / (arma::as_scalar(ZtRZ) / 2 + gb));
    
    //// prec0 l = 1
    hpl1(1) = R::rgamma(n * T / 2 + ga, 
        1 / (ZtRZl(0, 0) / 2 + gb));
    
    //// prec0 l = 2
    hpl2(1) = R::rgamma(n * T / 2 + ga, 
        1 / (ZtRZl(1, 0) / 2 + gb));
    
    //// prec1
    vTL2 = wtl(arma::span(2 * T, TL - 1)) - hp0(0);
    hp0(3) = R::rgamma(T * (L - 2) / 2 + ga,
         1 / (arma::dot(vTL2, vTL2) / 2 + gb));
    
    //// prec1 l = 1
    vT = wtl(arma::span(0, T - 1)) - hpl1(0);
    hpl1(2) = R::rgamma(T / 2 + ga,
         1 / (arma::dot(vT, vT) / 2 + gb));
    
    //// prec1 l = 2
    vT = wtl(arma::span(T, 2 * T - 1)) - hpl2(0);
    hpl2(2) = R::rgamma(T / 2 + ga,
         1 / (arma::dot(vT, vT) / 2 + gb));
    
    // keep
    if ((b > 0) && (b % nThin == 0)) {
      keep(b / nThin - 1, arma::span(0, k - 1)) = beta.t();
      keep(b / nThin - 1, arma::span(k, k + n * TL - 1)) = Wtls.t();
      keep(b / nThin - 1, arma::span(k + n * TL, k + n * TL + TL - 1)) = wtl.t();
      keep(b / nThin - 1, arma::span(k + n * TL + TL, k + n * TL + TL + 3)) = hp0.t();
      keep(b / nThin - 1, arma::span(k + n * TL + TL + 4, k + n * TL + TL + 3 + kl)) = betal1.t();
      keep(b / nThin - 1, arma::span(k + n * TL + TL + 4 + kl, k + n * TL + TL + 6 + kl)) = hpl1.t();
      keep(b / nThin - 1, arma::span(k + n * TL + TL + 7 + kl, k + n * TL + TL + 6 + 2 * kl)) = betal2.t();
      keep(b / nThin - 1, arma::span(k + n * TL + TL + 7 + 2 * kl, k + n * TL + TL + 9 + 2 * kl)) = hpl2.t();
    }
  }
  
  return keep;
}

// [[Rcpp::export]]
arma::mat predGlmBerModelPaperKFCV(
    const int B, const int Nweak, 
    const int n, const int newn, const int T, const int L,
    arma::mat X, const arma::mat X0, const arma::mat X1, 
    arma::vec year, arma::vec day, arma::vec newyear, arma::vec newday, 
    arma::mat dist, arma::mat Xb, arma::mat Wtls0,
    const arma::mat beta, const arma::mat betal1, const arma::mat betal2,
    const arma::mat Wtls, const arma::mat wtl, const arma::mat hp,
    const arma::uvec weak, const arma::uvec weak1, const arma::uvec weak2,
    const arma::vec noNA1, const arma::vec noNA2, 
    const arma::uvec Weak1, const arma::uvec Weak2, 
    const arma::uvec indLag1, const arma::uvec indLag2, const arma::uvec indLag12, 
    const arma::uvec indl, const arma::uvec indl1, const arma::uvec indl2, 
    const arma::uvec indSubmodel,
    const arma::vec prob
) {
  
  arma::mat Sigma11(newn, newn);
  arma::mat Sigma12(newn, n);
  arma::mat Sigma22inv(n, n);
  arma::mat SigmaAux(n, newn); // transposed
  arma::mat Sigma(newn, newn);
  
  arma::uvec indWUvec(n);
  arma::umat indWUmat(n, T * L);
  arma::uvec indW0Uvec(newn);
  arma::umat indW0Umat(newn, T * L);
  for (int l = 0; l < L; ++l) { 
    for (int t = 0; t < T; ++t) {
      indWUmat.col(l * T + t) = arma::find((year == t) && (day == l));
      indW0Umat.col(l * T + t) = arma::find((newyear == t) && (newday == l));
    }
  }
  
  double sigma0, sigma0l1, sigma0l2;
  
  arma::uword ind; 
  arma::uvec ind1, ind2, indAux;
  
  for (int b = 0; b < B; ++b) {
    
    // weak records
    X.submat(Weak1, indLag12) = X1.submat(Weak1, indLag12);
    X.submat(Weak2, indLag12) = X1.submat(Weak2, indLag12);
    for (int i = 0; i < Nweak; i++) {
      ind = weak(i);
      if (R::runif(0, 1) < prob(ind)) {
        if (noNA1(i)) {
          ind1  = weak1(i);
          X.submat(ind1, indLag1) = X1.submat(ind1, indLag1);
        } 
        if (noNA2(i)) {
          ind2  = weak2(i);
          X.submat(ind2, indLag2) = X1.submat(ind2, indLag2);
        } 
      } else {
        if (noNA1(i)) {
          ind1 = weak1(i);
          X.submat(ind1, indLag1)  = X0.submat(ind1, indLag1);
          X.submat(ind1, indLag12) = X0.submat(ind1, indLag12);
        }
        if (noNA2(i)) {
          ind2 = weak2(i);
          X.submat(ind2, indLag2)  = X0.submat(ind2, indLag2);
          X.submat(ind2, indLag12) = X0.submat(ind2, indLag12);
        }
      }
    }
    
    // Xb
    indAux = b;
    Xb.submat(indAux, indl)  = beta.row(b) * X.rows(indl).t();
    Xb.submat(indAux, indl1) = betal1.row(b) * X.submat(indl1, indSubmodel).t();
    Xb.submat(indAux, indl2) = betal2.row(b) * X.submat(indl2, indSubmodel).t();
    
    // Wtls
    Sigma11    = exp(- hp(b, 1) * dist.submat(0, 0, newn - 1, newn - 1));
    Sigma12    = exp(- hp(b, 1) * dist.submat(0, newn, newn - 1, newn + n - 1));
    Sigma22inv = arma::inv_sympd(exp(- hp(b, 1) * dist.submat(newn, newn, newn + n - 1, newn + n - 1)));
    SigmaAux   = (Sigma12 * Sigma22inv).t();
    Sigma      = Sigma11 - SigmaAux.t() * Sigma12.t();
    Sigma      = arma::chol(Sigma, "lower");
    // l = 3,...,L
    sigma0 = 1 / sqrt(hp(b, 0));
    for (int l = 2; l < L; ++l) { 
      for (int t = 0; t < T; ++t) {
        indWUvec = indWUmat.col(l * T + t);
        indW0Uvec = indW0Umat.col(l * T + t);
        Wtls0.submat(indAux, indW0Uvec) = 
          wtl(b, l * T + t) + 
          (Wtls.submat(indAux, indWUvec) - wtl(b, l * T + t)) * SigmaAux +
          sigma0 * (Sigma * arma::randn(newn)).t();
      }
    }
    // l = 1, 2 (0, 1)
    sigma0l1 = 1 / sqrt(hp(b, 2));
    sigma0l2 = 1 / sqrt(hp(b, 3));
    for (int t = 0; t < T; ++t) {
      
      indWUvec = indWUmat.col(t);
      indW0Uvec = indW0Umat.col(t);
      Wtls0.submat(indAux, indW0Uvec) = 
        wtl(b, t) + 
        (Wtls.submat(indAux, indWUvec) - wtl(b, t)) * SigmaAux +
        sigma0l1 * (Sigma * arma::randn(newn)).t();
      
      indWUvec = indWUmat.col(T + t);
      indW0Uvec = indW0Umat.col(T + t);
      Wtls0.submat(indAux, indW0Uvec) = 
        wtl(b, T + t) + 
        (Wtls.submat(indAux, indWUvec) - wtl(b, T + t)) * SigmaAux +
        sigma0l2 * (Sigma * arma::randn(newn)).t();
    }
    
  }
  
  return Xb + Wtls0;
}