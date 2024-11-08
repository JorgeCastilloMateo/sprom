# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

rtnorm <- function(N, mu, sigma) {
    .Call('_sprom_rtnorm', PACKAGE = 'sprom', N, mu, sigma)
}

rtnorm1 <- function(mu, sigma, a, b) {
    .Call('_sprom_rtnorm1', PACKAGE = 'sprom', mu, sigma, a, b)
}

dtnorm1 <- function(x, mu, sigma, a, b) {
    .Call('_sprom_dtnorm1', PACKAGE = 'sprom', x, mu, sigma, a, b)
}

rkolmogorov <- function(N) {
    .Call('_sprom_rkolmogorov', PACKAGE = 'sprom', N)
}

MH <- function(N, old, e) {
    .Call('_sprom_MH', PACKAGE = 'sprom', N, old, e)
}

glmBer <- function(N, k, Nl, kl, Nweak, X, X0, X1, beta, betal1, betal2, weak, weak1, weak2, noNA1, noNA2, Weak1, Weak2, indLag1, indLag2, indLag12, indl, indl1, indl2, indSubmodel, prob, abLim, V, keep, nSims, nThin, nBurnin, nReport) {
    .Call('_sprom_glmBer', PACKAGE = 'sprom', N, k, Nl, kl, Nweak, X, X0, X1, beta, betal1, betal2, weak, weak1, weak2, noNA1, noNA2, Weak1, Weak2, indLag1, indLag2, indLag12, indl, indl1, indl2, indSubmodel, prob, abLim, V, keep, nSims, nThin, nBurnin, nReport)
}

predGlmBerKFCV <- function(B, Nweak, X, X0, X1, Xb, beta, betal1, betal2, weak, weak1, weak2, noNA1, noNA2, Weak1, Weak2, indLag1, indLag2, indLag12, indl, indl1, indl2, indSubmodel, prob) {
    .Call('_sprom_predGlmBerKFCV', PACKAGE = 'sprom', B, Nweak, X, X0, X1, Xb, beta, betal1, betal2, weak, weak1, weak2, noNA1, noNA2, Weak1, Weak2, indLag1, indLag2, indLag12, indl, indl1, indl2, indSubmodel, prob)
}

glmBerGP <- function(N, k, X, beta, abLim, alpha, hp, Rinv, site, n, V, keep, psi, hpPsi, year, T, randomEffectYear, randomEffectDay, nSims, nThin, nBurnin, nReport) {
    .Call('_sprom_glmBerGP', PACKAGE = 'sprom', N, k, X, beta, abLim, alpha, hp, Rinv, site, n, V, keep, psi, hpPsi, year, T, randomEffectYear, randomEffectDay, nSims, nThin, nBurnin, nReport)
}

glmBerGPtl <- function(N, k, Nl, kl, Nweak, X, X0, X1, beta, betal1, betal2, hpl1, hpl2, Ws, hp0, dist, decayPrior, n, wtl, T, L, site, time, weak, weak1, weak2, noNA1, noNA2, Weak1, Weak2, indLag1, indLag2, indLag12, indl, indl1, indl2, indSubmodel, prob, abLim, V, nb, ga, gb, da, db, keep, nSims, nThin, nBurnin, nReport) {
    .Call('_sprom_glmBerGPtl', PACKAGE = 'sprom', N, k, Nl, kl, Nweak, X, X0, X1, beta, betal1, betal2, hpl1, hpl2, Ws, hp0, dist, decayPrior, n, wtl, T, L, site, time, weak, weak1, weak2, noNA1, noNA2, Weak1, Weak2, indLag1, indLag2, indLag12, indl, indl1, indl2, indSubmodel, prob, abLim, V, nb, ga, gb, da, db, keep, nSims, nThin, nBurnin, nReport)
}

predGlmBerGPtl <- function(B, Nweak, n, newn, T, L, X, X0, X1, newsite, dist, Xb, Ws0, beta, betal1, betal2, Ws, hp, weak, weak1, weak2, noNA1, noNA2, Weak1, Weak2, indLag1, indLag2, indLag12, indl, indl1, indl2, indSubmodel, prob) {
    .Call('_sprom_predGlmBerGPtl', PACKAGE = 'sprom', B, Nweak, n, newn, T, L, X, X0, X1, newsite, dist, Xb, Ws0, beta, betal1, betal2, Ws, hp, weak, weak1, weak2, noNA1, noNA2, Weak1, Weak2, indLag1, indLag2, indLag12, indl, indl1, indl2, indSubmodel, prob)
}

glmBerModelPaper <- function(N, k, Nl, kl, Nweak, X, X0, X1, beta, betal1, betal2, hpl1, hpl2, Wtls, hp0, dist, decayPrior, n, wtl, T, L, site, year, day, weak, weak1, weak2, noNA1, noNA2, Weak1, Weak2, indLag1, indLag2, indLag12, indl, indl1, indl2, indSubmodel, prob, abLim, V, nb, ga, gb, da, db, keep, nSims, nThin, nBurnin, nReport) {
    .Call('_sprom_glmBerModelPaper', PACKAGE = 'sprom', N, k, Nl, kl, Nweak, X, X0, X1, beta, betal1, betal2, hpl1, hpl2, Wtls, hp0, dist, decayPrior, n, wtl, T, L, site, year, day, weak, weak1, weak2, noNA1, noNA2, Weak1, Weak2, indLag1, indLag2, indLag12, indl, indl1, indl2, indSubmodel, prob, abLim, V, nb, ga, gb, da, db, keep, nSims, nThin, nBurnin, nReport)
}

predGlmBerModelPaperKFCV <- function(B, Nweak, n, newn, T, L, X, X0, X1, year, day, newyear, newday, dist, Xb, Wtls0, beta, betal1, betal2, Wtls, wtl, hp, weak, weak1, weak2, noNA1, noNA2, Weak1, Weak2, indLag1, indLag2, indLag12, indl, indl1, indl2, indSubmodel, prob) {
    .Call('_sprom_predGlmBerModelPaperKFCV', PACKAGE = 'sprom', B, Nweak, n, newn, T, L, X, X0, X1, year, day, newyear, newday, dist, Xb, Wtls0, beta, betal1, betal2, Wtls, wtl, hp, weak, weak1, weak2, noNA1, noNA2, Weak1, Weak2, indLag1, indLag2, indLag12, indl, indl1, indl2, indSubmodel, prob)
}

