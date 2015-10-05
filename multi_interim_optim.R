###############################################################################
# Multiple Interim Optimization
# Author
#   Xiaoqi Lu: lx2170@cumc.columbia.edu
#   Qiqi Deng: qiqi.deng@beohringer-ingelheim.com
# Version
#   07/19/2015
#   07/28/2015  set default values for some arguments
#   07/30/2015  fix a fatal error in ProbPass()
#   08/14/2015  heatmap and power contour completed
###############################################################################

###############################################################################
# To-do list:
# 1. Add revenue() and cost() as input, rather than bcRatio. (as actual cost)
# 2. Add constraints for fix t, N0
# 3. Consider d as distribution
# 4. Design() as interactive
###############################################################################

library(optimx)

###############################################################################
# Information Fraction
# Input
#   n0    a vector, number(s) of subjects in group 0 at interim(s)
#   n1    a vector, number(s) of subjects in group 1 at interim(s)
#   N0    sample size of group 0
#   N1    sample size of group 1
#   s0    standard deviation of group 0
#         DEFAULT = 1
#   s1    standard deviation of group 1
#         DEFAULT = 1
# Output
#         a vector, time(s) at interim(s)
###############################################################################
Info <- function(n0, n1, N0, N1, s0 = 1, s1 = 1) {
  v0 <- s0 ^ 2
  v1 <- s1 ^ 2
  t <- (v0 / N0 + v1 / N1) / (v0 / n0 + v1 / n1)
  return(t)
}

###############################################################################
# Sample Size at Interim
# Input
#   t     a vector, time(s) at interim(s)
#   N0    sample size of group 0
#   r     n0 : n1 = N0 : N1 = 1 : r
#         DEFAULT = 1
#   s0    standard deviation of group 0
#         DEFAULT = 1
#   s1    standard deviation of group 1
#         DEFAULT = 1
# Output
#   n0    a vector, number(s) of subjects in group 0 at interim(s)
#   n1    a vector, number(s) of subjects in group 1 at interim(s)
#   N0    sample size of group 0
#   N1    sample size of group 1
###############################################################################
NInterim <- function(t, N0, r = 1, s0 = 1, s1 = 1) {
  v0 <- s0 ^ 2
  v1 <- s1 ^ 2
  N1 <- ceiling(N0 * r)
  n0 <- ceiling(t * (v0 + v1 / r) / (v0 / N0 + v1 / N1))
  n1 <- ceiling(t * (v0 * r + v1) / (v0 / N0 + v1 / N1))
  return(list(n0 = n0, n1 = n1, N0 = N0, N1 = N1))
}

###############################################################################
# Drift of Brownian Motion
# Input
#   d     true difference in means
#   N0    sample size of group 0
#   N1    sample size of group 1
#   s0    standard deviation of group 0
#         DEFAULT = 1
#   s1    standard deviation of group 1
#         DEFAULT = 1
# Output
#         drift of Brownian Motion
###############################################################################
Drift <- function(d, N0, N1, s0 = 1, s1 = 1) {
  mu <- d / sqrt(s0 ^ 2 / N0 + s1 ^ 2 / N1)
  return(mu)
}

###############################################################################
# Probability of Passing Critical Values (Continuing)
# Input
#   c     a vector, critical value(s) for B-value(s)
#   t     a vector, time(s) at interim(s)
#   mu    drift of Brownian Motion
#   pass  a vector, whether B passes critical value(s) at interim(s)
#         default = rep(TRUE, length(c))
# Output
#         probability of (not) passing
# Notes
#   In the following cases, extend c and t with final analysis.
#   overall power: pass = rep(TRUE, length(c))
#   wrong continuing: pass = c(rep(TRUE, length(c) - 1), FALSE)
###############################################################################
ProbPass <- function(c, t, mu, pass) {
  if (length(c) == 1) {
    p <- pnorm(mu * sqrt(t) - c / sqrt(t), lower.tail = pass)
  } else {
    fun <- function(b1, c, t, pass) {
      f <- ProbPass(c[-1] - b1, t[-1] - t[1], mu, pass[-1]) * dnorm(b1, mu * t[1], sqrt(t[1]))
      return(f)
    }
    funApply <- function(b1, c, t, pass) {
      f <- sapply(b1, fun, c = c, t = t, pass = pass)
      return(f)
    }
    if (pass[1]) {
      p <- integrate(funApply, c[1], Inf, c = c, t = t, pass = pass)$value
    } else {
      p <- integrate(funApply, -Inf, c[1], c = c, t = t, pass = pass)$value
    }
  }
  return(p)
}

###############################################################################
# Conditional Power
# Input
#   muA   drift (standardized difference) under alternative hypothesis
#   c     a vector, critical value(s) for B-value(s)
#   t     a vector, time(s) at interim(s)
#   z     critical value for final analysis
# Output
#   trend a vector, conditional power at each interim, using current trend
#   pred  a vector, contitional power at each interim, using Bayesian 
###############################################################################
CondPower <- function(muA, c, t, z) {
  n <- length(t)
  if (n == 0) return(list(alt = numeric(0),
                          trend = numeric(0),
                          pred = numeric(0)))
  mu <- c / t
  sd <- 1 / sqrt(t)
  cExt <- c(c, z)
  tExt <- c(t, 1)
  alt <- trend <- pred <- numeric(n)
  for (i in seq_len(n)) {
    alt[i] <- ProbPass(cExt[(i + 1) : (n + 1)] - cExt[i],
                       tExt[(i + 1) : (n + 1)] - tExt[i],
                       muA,
                       rep(TRUE, n - i + 1))
    trend[i] <- ProbPass(cExt[(i + 1) : (n + 1)] - cExt[i],
                         tExt[(i + 1) : (n + 1)] - tExt[i],
                         mu[i],
                         rep(TRUE, n - i + 1))
    Fun <- function(x) {
      FunSingle <- function(xSingle) {
        p <- ProbPass(cExt[(i + 1) : (n + 1)] - cExt[i],
                      tExt[(i + 1) : (n + 1)] - tExt[i],
                      xSingle,
                      rep(TRUE, n - i + 1))
        return(p)
      }
      return(sapply(x, FunSingle) * dnorm(x, mu[i], sd[i]))
    }
    pred[i] <- integrate(Fun, -Inf, Inf)$value
  }
  return(list(alt = alt, trend = trend, pred = pred))
}

###############################################################################
# Cost of Clinical Trial
# Input
#   n0    number of subjects in group 0
#   n1    number of subjects in group 1
# Output
#         cost of clinical trial
###############################################################################
Cost <- function(n0, n1) {
  u <- n0 + n1
  return(u)
}

###############################################################################
# Utility
# Input
#   ben   benefit of clinical trial, if success
#   d     true difference in means
#   c     a vector, critical value(s) for B-value(s)
#   z     critical value for final analysis
#   n0    a vector, number(s) of subjects in group 0 at interim(s)
#   n1    a vector, number(s) of subjects in group 1 at interim(s)
#   N0    sample size of group 0
#   N1    sample size of group 1
#   s0    standard deviation of group 0
#         DEFAULT = 1
#   s1    standard deviation of group 1
#         DEFAULT = 1
# Output
#   u     expected utility
#   pS    a vector, probability of each scenario
#   nS0   a vector, comsumed sample size (group 0) of each scenario
#   nS1   a vector, comsumed sample size (group 1) of each scenario
#   uS    a vector, utility of each scenario
###############################################################################
Utility <- function(ben, d, c, z, n0, n1, N0, N1, s0 = 1, s1 = 1) {
  mu <- Drift(d, N0, N1, s0, s1)
  t <- Info(n0, n1, N0, N1, s0, s1)
  k <- length(c)
  pStop <- uStop <- numeric(k)
  for (i in seq_len(k)) {
    pStop[i] <- ProbPass(c[seq_len(i)], t[seq_len(i)], mu, c(rep(TRUE, i - 1), FALSE))
    uStop[i] <- - Cost(n0[i], n1[i])
  }
  pFail <- ProbPass(c(c, z), c(t, 1), mu, c(rep(TRUE, k), FALSE))
  uFail <- - Cost(N0, N1)
  pSuccess <- ProbPass(c(c, z), c(t, 1), mu, rep(TRUE, k + 1))
  uSuccess <- - Cost(N0, N1) + ben
  u <- sum(pStop * uStop) + pFail * uFail + pSuccess * uSuccess
  pS <- c(pStop, pFail, pSuccess)
  nS0 <- c(n0, N0, N0)
  nS1 <- c(n1, N1, N1)
  uS <- c(uStop, uFail, uSuccess)
  return(list(u = u, pS = pS, nS0 = nS0, nS1 = nS1, uS = uS))
}

###############################################################################
# Utility Wrapper
#   Instead of n0, n1, N0, N1, the wrapper takes t, r and N0 as input.
# Input
#   ben   benefit of clinical trial, if success
#   d     true difference in means
#   t     a vector, time(s) at interim(s)
#   c     a vector, critical value(s) for B-value(s)
#   z     critical value for final analysis
#   N0    sample size of group 0
#   r     sample size ratio, n0 : n1 = N0 : N1 = 1 : r
#         DEFAULT = 1
#   s0    standard deviation of group 0
#         DEFAULT = 1
#   s1    standard deviation of group 1
#         DEFAULT = 1
#   show  whether to show results
# Output
#   u     expected utility
#   pS    a vector, probability of each scenario
#   nS0   a vector, comsumed sample size (group 0) of each scenario
#   nS1   a vector, comsumed sample size (group 1) of each scenario
#   uS    a vector, utility of each scenario
###############################################################################
UtilityWrap <- function(ben, d, t, c, z, N0, r = 1, s0 = 1, s1 = 1, show = TRUE) {
  nInterim <- NInterim(t, N0, r, s0, s1)
  N1 <- nInterim$N1
  n0 <- nInterim$n0
  n1 <- nInterim$n1
  uAnalysis <- Utility(ben, d, c, z, n0, n1, N0, N1, s0, s1)
  if (show) {
    cat("Benefit =", ben, "\n")
    cat("Diff.   =", d, "\n")
    cat("N0 : N1 = 1 :", r, "\n")
    cat("Sigma 0 =", s0, "\n")
    cat("Sigma 1 =", s1, "\n")
    cat("Drift   =", Drift(d, N0, N1, s0, s1), "\n")
    cat("           ", format(c(seq_along(t), "Fail", "Success"), width = 12, justify = "right"), "\n")
    cat("Information", format(c(t, 1, 1), width = 12), "\n")
    cat("N0         ", format(uAnalysis$nS0, width = 12), "\n")
    cat("N1         ", format(uAnalysis$nS1, width = 12), "\n")
    if (length(t) > 0) {
      cat("B-value    ", format(c(c, z, z), digits = 3, nsmall = 2, width = 12), "\n")
      cat("Z-value    ", format(c(c / sqrt(t), z, z), digits = 3, nsmall = 2, width = 12), "\n")
    }
    cat("Probability", paste(format(100 * uAnalysis$pS, digits = 3, nsmall = 2, width = 11), "%", sep = ""), "\n")
    cat("Utility    ", format(uAnalysis$uS, digits = 3, width = 12), "\n")
    cat("Expected Utility =", uAnalysis$u, "\n")
  }
  return(uAnalysis)
}

# UtilityWrap(100, 0.2, NULL, NULL, 1.6, 450)

###############################################################################
# Design
# Input
#   dA    difference in means under alternative hypothesis
#   t     a vector, time(s) at interim(s)
#   c     a vector, critical value(s) for B-value(s)
#   z     critical value for final analysis
#   N0    sample size of group 0
#   r     sample size ratio, n0 : n1 = N0 : N1 = 1 : r
#         DEFAULT = 1
#   s0    standard deviation of group 0
#         DEFAULT = 1
#   s1    standard deviation of group 1
#         DEFAULT = 1
###############################################################################
Design <- function(dA, t, c, z, N0, r = 1, s0 = 1, s1 = 1) {
  nInterim <- NInterim(t, N0, r, s0, s1)
  N1 <- nInterim$N1
  n0 <- nInterim$n0
  n1 <- nInterim$n1
  muA <- Drift(dA, N0, N1, s0, s1)
  cp <- CondPower(muA, c, t, z)
  muA <- Drift(dA, N0, N1, s0, s1)
  k <- length(t)
  prob <- matrix(NA, 2, k + 2)
  for (i in seq_len(k)) {
    prob[1, i] <- ProbPass(c[seq_len(i)], t[seq_len(i)], 0, c(rep(TRUE, i - 1), FALSE))
    prob[2, i] <- ProbPass(c[seq_len(i)], t[seq_len(i)], muA, c(rep(TRUE, i - 1), FALSE))
  }
  prob[1, k + 1] <- ProbPass(c(c, z), c(t, 1), 0, c(rep(TRUE, k), FALSE))
  prob[2, k + 1] <- ProbPass(c(c, z), c(t, 1), muA, c(rep(TRUE, k), FALSE))
  prob[1, k + 2] <- ProbPass(c(c, z), c(t, 1), 0, rep(TRUE, k + 1))
  prob[2, k + 2] <- ProbPass(c(c, z), c(t, 1), muA, rep(TRUE, k + 1))
  cat("Drift (std. diff.) =", muA, "\n")
  cat("           ", format(c(seq_along(t), "Fail", "Success"), width = 12, justify = "right"), "\n")
  cat("Information", format(c(t, 1, 1), width = 12), "\n")
  cat("N0         ", format(c(n0, N0, N0), width = 12), "\n")
  cat("N1         ", format(c(n1, N1, N1), width = 12), "\n")
  if (k > 0) {
    cat("B-value    ", format(c(c, z, z), digits = 3, nsmall = 2, width = 12), "\n")
    cat("Z-value    ", format(c(c / sqrt(t), z, z), digits = 3, nsmall = 2, width = 12), "\n")
    cat("C.P.(alt.) ", format(cp$alt, digits = 3, nsmall = 2, width = 12), "\n")
    cat("C.P.(trend)", format(cp$trend, digits = 3, nsmall = 2, width = 12), "\n")
    cat("C.P.(pred.)", format(cp$pred, digits = 3, nsmall = 2, width = 12), "\n")
  }
  cat("Prob.(H0)  ", paste(format(100 * prob[1, ], digits = 3, nsmall = 2, width = 11), "%", sep = ""), "\n")
  cat("Prob.(HA)  ", paste(format(100 * prob[2, ], digits = 3, nsmall = 2, width = 11), "%", sep = ""), "\n")
  return()
}

# Design(0.2, NULL, NULL, qnorm(0.05, lower.tail = F), 500)

###############################################################################
# Utility Analysis
# Input
#   ben   benefit of clinical trial, if success
#   dA    difference under alternative hypothesis
#   d     difference(s) of interest
#         DEFAULT = seq(0, dA, length.out = 11)
#   t     a vector, time(s) at interim(s)
#   c     a vector, critical value(s) for B-value(s)
#   z     critical value for final analysis
#         if z = NULL (default) then z is set by alpha
#   N0    sample size of group 0
#         if N0 = NULL (default) then N0 is set by power
#   r     sample size ratio, n0 : n1 = N0 : N1 = 1 : r
#         DEFAULT = 1
#   s0    standard deviation of group 0
#         DEFAULT = 1
#   s1    standard deviation of group 1
#         DEFAULT = 1
#   alpha designed type I error
#         DEFAULT = 0.05
#   power overall power
#         DEFAULT = 0.8
#   minpower
#         minimum overall power
#         DEFAULT = power
#   maxpower
#         maximum overall power
#         DEFAULT = 0.95
#   heat  draw heatmap instead of curve(s), only works when length(d) = 1
#         DEFAULT = FALSE
###############################################################################
UtilityAnalysis <- function(ben, dA, d = seq(0, dA, length.out = 7),
                            t, c, z = NULL, N0 = NULL,
                            r = 1, s0 = 1, s1 = 1,
                            alpha = 0.05, power = 0.8,
                            minpower = power, maxpower = 0.95,
                            heat = FALSE) {
  if (is.null(z)) z <- qnorm(alpha, lower.tail = FALSE)
  if (is.null(N0)) {
    Fun <- function(N0) {
      N1 <- ceiling(N0 * r)
      muA <- Drift(dA, N0, N1, s0, s1)
      p <- ProbPass(c(c, z), c(t, 1), muA, rep(TRUE, length(c) + 1))
      return(p - power)
    }
    N0 <- ceiling(uniroot(Fun, c(1, 1e4), extendInt = "upX")$root)
  }
  Design(dA, t, c, z, N0, r, s0, s1)
  u <- rep(NA, length(d))
  for (j in seq_along(d)) {
    u[j] <- UtilityWrap(ben, d[j], t, c, z, N0, r, s0, s1, show = FALSE)$u
  }
  cat("Diff.  ", format(d, digits = 2, width = 6), "\n")
  cat("Utility", format(u, digits = 2, width = 6), "\n")
  
  N0Seq <- ceiling(seq(max(min(0.5 * N0, N0 - 200), 0), max(2.0 * N0, N0 + 200), 10))
  if (length(t) == 0) {
    par(mfrow = c(2, 2), mar = c(4, 4, 0.5, 3))
    pSuccessMat <- uMat <- matrix(NA, length(N0Seq), length(d))
    for (i in seq_along(N0Seq)) {
      for (j in seq_along(d)) {
        uTmp <- UtilityWrap(ben, d[j], t, c, z, N0Seq[i], r, s0, s1, FALSE)
        pSuccessMat[i, j] <- uTmp$pS[2]
        uMat[i, j] <- uTmp$u
      }
    }
    matplot(N0Seq, pSuccessMat, type = "l", lwd = 2, lty = 1, ylim = c(0, 1),
            col = topo.colors(length(d)), xlab = "N0", ylab = "P(success)")
    abline(v = N0, lwd = 2, lty = 2)
    abline(h = c(alpha, power), lwd = 2, lty = 3, col = "red")
    axis(4, pSuccessMat[length(N0Seq), ], format(d, digits = 1), las = 1)
    matplot(N0Seq, uMat, type = "l", lwd = 2, lty = 1,
            col = topo.colors(length(d)), xlab = "N0", ylab = "Utility")
    abline(v = N0, lwd = 2, lty = 2)
    axis(4, uMat[length(N0Seq), ], format(d, digits = 1), las = 1)
  } else if (length(t) == 1) {
    tSeq <- seq(0.1, 0.9, 0.02)
    cSeq <- seq(-z, z, 0.05)
    if ((length(d) == 1) && heat) {
      par(mfrow = c(1, 2), mar = c(4, 4, 0.5, 0.5))
      # pI/p/u ~ (c, t)
      powerMatN0 <- uMatN0 <- matrix(NA, length(cSeq), length(tSeq))
      for (i in seq_along(cSeq)) {
        for (j in seq_along(tSeq)) {
          uTmp <- UtilityWrap(ben, d, tSeq[j], cSeq[i], z, N0, r, s0, s1, show = FALSE)
          uMatN0[i, j] <- uTmp$u
          muA <- Drift(dA, N0, ceiling(N0 * r), s0, s1)
          powerMatN0[i, j] <- ProbPass(c(cSeq[i], z), c(tSeq[j], 1), muA, c(T, T))
        }
      }
      # pI/p/u ~ (N0, t)
      powerMatC <- uMatC <- matrix(NA, length(N0Seq), length(tSeq))
      for (i in seq_along(N0Seq)) {
        for (j in seq_along(tSeq)) {
          uTmp <- UtilityWrap(ben, d, tSeq[j], c, z, N0Seq[i], r, s0, s1, show = FALSE)
          uMatC[i, j] <- uTmp$u
          muA <- Drift(dA, N0Seq[i], ceiling(N0Seq[i] * r), s0, s1)
          powerMatC[i, j] <- ProbPass(c(c, z), c(tSeq[j], 1), muA, c(T, T))
        }
      }
      # pI/p/u ~ (N0, c)
      powerMatT <- uMatT <- matrix(NA, length(N0Seq), length(cSeq))
      for (i in seq_along(N0Seq)) {
        for (j in seq_along(cSeq)) {
          uTmp <- UtilityWrap(ben, d, t, cSeq[j], z, N0Seq[i], r, s0, s1, show = FALSE)
          uMatT[i, j] <- uTmp$u
          muA <- Drift(dA, N0Seq[i], ceiling(N0Seq[i] * r), s0, s1)
          powerMatT[i, j] <- ProbPass(c(cSeq[j], z), c(t, 1), muA, c(T, T))
        }
      }
      zRange <- range(uMatN0, uMatC, uMatT)
      #
      image(cSeq, tSeq, uMatN0, col = heat.colors(32),
            zlim = zRange,
            xlab = "B-value", ylab = "Information")
      contour(cSeq, tSeq, powerMatN0,
              levels = sort(unique(c(power, minpower, maxpower, seq(0, 1, 0.2)))),
              lty = 3, lwd = 2, col = "blue", add = TRUE)
      indMat <- (powerMatN0 >= minpower) & (powerMatN0 <= maxpower)
      ind <- which(uMatN0 == max(uMatN0[indMat]),
                   arr.ind = TRUE)
      points(cSeq[ind[1]], tSeq[ind[2]], pch = 4, cex = 2)
      abline(h = t, v = c, lty = 2, lwd = 2)
      #
      image(N0Seq, tSeq, uMatC, col = heat.colors(32),
            zlim = zRange,
            xlab = "N0", ylab = "Information")
      contour(N0Seq, tSeq, powerMatC,
              levels = sort(unique(c(power, minpower, maxpower, seq(0, 1, 0.2)))),
              lty = 3, lwd = 2, col = "blue", add = TRUE)
      indMat <- (powerMatC >= minpower) & (powerMatC <= maxpower)
      ind <- which(uMatC == max(uMatC[indMat]),
                   arr.ind = TRUE)
      points(N0Seq[ind[1]], tSeq[ind[2]], pch = 4, cex = 2)
      abline(h = t, v = N0, lty = 2, lwd = 2)
      #
      image(cSeq, N0Seq, t(uMatT), col = heat.colors(32),
            zlim = zRange,
            ylab = "N0", xlab = "B-value")
      contour(cSeq, N0Seq, t(powerMatT),
              levels = sort(unique(c(power, minpower, maxpower, seq(0, 1, 0.2)))),
              lty = 3, lwd = 2, col = "blue", add = TRUE)
      indMat <- (powerMatT >= minpower) & (powerMatT <= maxpower)
      ind <- which(uMatT == max(uMatT[indMat]),
                   arr.ind = TRUE)
      points(cSeq[ind[2]], N0Seq[ind[1]], pch = 4, cex = 2)
      abline(v = c, h = N0, lty = 2, lwd = 2)
      color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
        scale = (length(lut)-1)/(max-min)
        plot(c(0,50), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
        axis(2, format(ticks, digits = 3, width = 5), las=1)
        for (i in 1:(length(lut)-1)) {
          y = (i-1)/scale + min
          rect(0,y,10,y+1/scale, col=lut[i], border=NA)
        }
      }
      color.bar(heat.colors(32), min(zRange), max(zRange))
#       plot(N0, N0, type = "p", pch = 4, cex = 2, xlim = range(N0Seq), ylim = range(N0Seq))
#       abline(h = N0, v = N0, lty = 2, lwd = 2)
#       Fun <- function(N0) {
#         N1 <- ceiling(N0 * r)
#         p <- ProbPass(c(c, z), c(t, 1), Drift(dA, N0, N1, s0, s1), c(T, T))
#         return(p)
#       }
#       powerSeq <- sapply(N0Seq, Fun)
#       points(N0Seq, N0Seq, type = "p", pch = 20, cex = 2, col = "blue")
#       text(N0Seq, N0Seq, format(powerSeq, digits = 2), pos = 4)
    } else {
      par(mfcol = c(3, 3), mar = c(4, 4, 0.5, 4))
      # pI/p/u ~ N0
      p1Mat <- pSuccessMat <- uMat <- matrix(NA, length(N0Seq), length(d))
      for (i in seq_along(N0Seq)) {
        for (j in seq_along(d)) {
          uTmp <- UtilityWrap(ben, d[j], t, c, z, N0Seq[i], r, s0, s1, show = FALSE)
          p1Mat[i, j] <- sum(uTmp$pS[-1])
          pSuccessMat[i, j] <- uTmp$pS[3]
          uMat[i, j] <- uTmp$u
        }
      }
      matplot(N0Seq, p1Mat, type = "l", lwd = 2, lty = 1, ylim = c(0, 1),
              col = topo.colors(length(d)), xlab = "N0", ylab = "P(pass)")
      abline(v = N0, lwd = 2, lty = 2)
      axis(4, p1Mat[length(N0Seq), ], format(d, digits = 1), las = 1)
      matplot(N0Seq, pSuccessMat, type = "l", lwd = 2, lty = 1, ylim = c(0, 1),
              col = topo.colors(length(d)), xlab = "N0", ylab = "P(success)")
      abline(v = N0, lwd = 2, lty = 2)
      abline(h = c(alpha, power), lwd = 2, lty = 3, col = "red")
      axis(4, pSuccessMat[length(N0Seq), ], format(d, digits = 1), las = 1)
      matplot(N0Seq, uMat, type = "l", lwd = 2, lty = 1,
              col = topo.colors(length(d)), xlab = "N0", ylab = "Utility")
      abline(v = N0, lwd = 2, lty = 2)
      axis(4, uMat[length(N0Seq), ], format(d, digits = 1), las = 1)
      # pI/p/u ~ c
      p1Mat <- pSuccessMat <- uMat <- matrix(NA, length(cSeq), length(d))
      for (i in seq_along(cSeq)) {
        for (j in seq_along(d)) {
          uTmp <- UtilityWrap(ben, d[j], t, cSeq[i], z, N0, r, s0, s1, show = FALSE)
          p1Mat[i, j] <- sum(uTmp$pS[-1])
          pSuccessMat[i, j] <- uTmp$pS[3]
          uMat[i, j] <- uTmp$u
        }
      }
      matplot(cSeq, p1Mat, type = "l", lwd = 2, lty = 1, ylim = c(0, 1),
              col = topo.colors(length(d)), xlab = "B-value", ylab = "P(pass)")
      abline(v = c, lwd = 2, lty = 2)
      axis(4, p1Mat[length(cSeq), ], format(d, digits = 1), las = 1)
      matplot(cSeq, pSuccessMat, type = "l", lwd = 2, lty = 1, ylim = c(0, 1),
              col = topo.colors(length(d)), xlab = "B-value", ylab = "P(success)")
      abline(v = c, lwd = 2, lty = 2)
      abline(h = c(alpha, power), lwd = 2, lty = 3, col = "red")
      axis(4, pSuccessMat[length(cSeq), ], format(d, digits = 1), las = 1)
      matplot(cSeq, uMat, type = "l", lwd = 2, lty = 1,
              col = topo.colors(length(d)), xlab = "B-value", ylab = "Utility")
      abline(v = c, lwd = 2, lty = 2)
      axis(4, uMat[length(cSeq), ], format(d, digits = 1), las = 1)
      # pI/p/u ~ t
      p1Mat <- pSuccessMat <- uMat <- matrix(NA, length(tSeq), length(d))
      for (i in seq_along(tSeq)) {
        for (j in seq_along(d)) {
          uTmp <- UtilityWrap(ben, d[j], tSeq[i], c, z, N0, r, s0, s1, show = FALSE)
          p1Mat[i, j] <- sum(uTmp$pS[-1])
          pSuccessMat[i, j] <- uTmp$pS[3]
          uMat[i, j] <- uTmp$u
        }
      }
      matplot(tSeq, p1Mat, type = "l", lwd = 2, lty = 1, xlim = c(0, 1), ylim = c(0, 1),
              col = topo.colors(length(d)), xlab = "Information", ylab = "P(pass)")
      abline(v = t, lwd = 2, lty = 2)
      axis(4, p1Mat[length(tSeq), ], format(d, digits = 1), las = 1)
      matplot(tSeq, pSuccessMat, type = "l", lwd = 2, lty = 1, xlim = c(0, 1), ylim = c(0, 1),
              col = topo.colors(length(d)), xlab = "Information", ylab = "P(success)")
      abline(v = t, lwd = 2, lty = 2)
      abline(h = c(alpha, power), lwd = 2, lty = 3, col = "red")
      axis(4, pSuccessMat[length(tSeq), ], format(d, digits = 1), las = 1)
      matplot(tSeq, uMat, type = "l", lwd = 2, lty = 1, xlim = c(0, 1),
              col = topo.colors(length(d)), xlab = "Information", ylab = "Utility")
      abline(v = t, lwd = 2, lty = 2)
      axis(4, uMat[length(tSeq), ], format(d, digits = 1), las = 1)
    }
  }
  return()
}

# UtilityAnalysis(3000, 0.2, d = 0.1, t = 0.6, c = 0.3, heat = T)


###############################################################################
# Utility Optimization
# Input
#   ben   benefit of clinical trial, if success
#   dA    difference under alternative hypothesis
#   d     difference of interest (scenario to be optimized)
#   k     number of interim(s)
#         DEFAULT = 1
#   z     critical value for final analysis
#         if z = NULL (default) then z is set by alpha
#   r     sample size ratio, n0 : n1 = N0 : N1 = 1 : r
#         DEFAULT = 1
#   s0    standard deviation of group 0
#         DEFAULT = 1
#   s1    standard deviation of group 1
#         DEFAULT = 1
#   alpha designed type I error
#         DEFAULT = 0.05
#   minpower
#         minimum overall power
#         DEFAULT = 0.8
#   maxpower
#         maximum overall power
#         DEFAULT = 0.95
###############################################################################
UtilityOptim <- function(ben, dA, d, k = 1, z = NULL,
                         r = 1, s0 = 1, s1 = 1,
                         alpha = 0.05, minpower = 0.8, maxpower = 0.95) {
  if (is.null(z)) z <- qnorm(alpha, lower.tail = FALSE)
  MinN0 <- function(c, t) {
    Fun <- function(N0) {
      N1 <- round(N0 * r)
      muA <- Drift(dA, N0, N1, s0, s1)
      p <- ProbPass(c(c, z), c(t, 1), muA, rep(TRUE, length(c) + 1))
      return(p - minpower)
    }
    N0 <- ceiling(uniroot(Fun, c(1, 1e4), extendInt = "upX")$root)
    return(N0)
  }
  MaxN0 <- function(c, t) {
    Fun <- function(N0) {
      N1 <- round(N0 * r)
      muA <- Drift(dA, N0, N1, s0, s1)
      p <- ProbPass(c(c, z), c(t, 1), muA, rep(TRUE, length(c) + 1))
      return(p - maxpower)
    }
    N0 <- floor(uniroot(Fun, c(1, 1e4), extendInt = "upX")$root)
    return(N0)
  }
  xIni <- rep(0.5, k)
  cIni <- rep(0, k)
  parIni <- c(cIni, xIni)
  NUByCX <- function(par) {
    c <- par[seq_len(k)]
    x <- par[seq_len(k) + k]
    t <- rev(cumprod(rev(x)))
    N0Min <- MinN0(c, t)
    N0Max <- MaxN0(c, t)
    UByCTN <- function(N0) {
      N0 <- round(N0)
      Fun <- function(d) {
        u <- UtilityWrap(ben, d, t, c, z, N0, r, s0, s1, show = FALSE)$u
        return(u)
      }
      u <- sum(sapply(d, Fun))
      return(u)
    }
    res <- optimize(UByCTN, interval = c(N0Min, N0Max), maximum = TRUE)
    return(list(N0 = round(res$maximum), u = res$objective))
  }
  UByCX <- function(par) {
    u <- NUByCX(par)$u
    return(u)
  }
  res <- optimx(parIni, UByCX,
                lower = c(rep(- z, k), rep(0.1, k)),
                upper = c(rep(z, k), rep(0.9, k)),
                method = "bobyqa",
                control = list(maximize = TRUE))
  c <- as.numeric(res[, seq_len(k)])
  x <- as.numeric(res[, seq_len(k) + k])
  t <- rev(cumprod(rev(x)))
  N0 <- NUByCX(c(c, x))$N0
  u <- UtilityWrap(ben, d, t, c, z, N0, r, s0, s1, show = FALSE)$u
  print(res$convcode)
  cat("Information:", t, "\n")
  cat("B-value:    ", c, "\n")
  cat("Z-value:    ", c / sqrt(t), "\n")
  cat("N0:         ", N0, "\n")
  cat("Utility:    ", u, "\n")
  return(list(c = c, t = t, N0 = N0, u = u))
}



###############################################################################
# Test Area
###############################################################################
ben <- 3000
dA <- 0.2
alpha <- 0.05
minpower <- power <- 0.8
maxpower <- 0.95
d <- 0.1
dSeq <- seq(0, dA, length.out = 21)

# res <- UtilityOptim(ben = ben, dA = dA, d = d, k = 1,
#                     minpower = minpower, maxpower = maxpower)
# Design(dA = dA, t = res$t, c = res$c, z = qnorm(alpha, lower.tail = F), N0 = res$N0)
# UtilityAnalysis(ben = ben, dA = dA, d = dSeq,
#                 t = res$t, c = res$c, N0 = res$N0,
#                 alpha = alpha, power = power, minpower = minpower, maxpower = maxpower,
#                 heat = T)
# Design(dA = dA, t = NULL, c = NULL, z = qnorm(alpha, lower.tail = F), N0 = res$N0)

tSeq <- cSeq <- N0Seq <- uSeq <- uSeq1 <- uSeq2 <- uSeq3 <- rep(NA, length(dSeq))
res <- UtilityOptim(ben = ben, dA = dA, d = d, k = 1, alpha = alpha, minpower = minpower, maxpower = maxpower)
t2 <- res$t
c2 <- res$c
N02 <- res$N0
res <- UtilityOptim(ben = ben, dA = dA, d = 0.5*d, k = 1, alpha = alpha, minpower = minpower, maxpower = maxpower)
t1 <- res$t
c1 <- res$c
N01 <- res$N0
res <- UtilityOptim(ben = ben, dA = dA, d = 1.5*d, k = 1, alpha = alpha, minpower = minpower, maxpower = maxpower)
t3 <- res$t
c3 <- res$c
N03 <- res$N0
for (i in seq_along(dSeq)) {
  res <- UtilityOptim(ben = ben, dA = dA, d = dSeq[i], k = 1, alpha = alpha, minpower = minpower, maxpower = maxpower)
  tSeq[i] <- res$t
  cSeq[i] <- res$c
  N0Seq[i] <- res$N0
  uSeq[i] <- res$u
  uSeq2[i] <- UtilityWrap(ben = ben, d = dSeq[i], t = t2, c = c2, z = qnorm(alpha, lower.tail = FALSE), N0 = N02, show = FALSE)$u
  uSeq1[i] <- UtilityWrap(ben = ben, d = dSeq[i], t = t1, c = c1, z = qnorm(alpha, lower.tail = FALSE), N0 = N01, show = FALSE)$u
  uSeq3[i] <- UtilityWrap(ben = ben, d = dSeq[i], t = t3, c = c3, z = qnorm(alpha, lower.tail = FALSE), N0 = N03, show = FALSE)$u
}

##############################################################################
# TO DO
# Slides: intro -> theory -> 3 examples: no futility, normal
# (t, c, N0)* ~ d

# pick a bc ratio s.t. optimal is not at boundary.
# 1: d = 0, 0.1, 0.2, t-c heatmap ONLY
# 2: rainbow (u ~ c,t; last row)
# 3: fixed (t, c, N0), draw the u~d (ONLY), other 3 in backup
# conservative (c = 0), aggresive (c = 1), recommended
##############################################################################
# par(mfrow = c(2, 2), mar = c(4, 4, 1, 1))
plot(dSeq, tSeq)
plot(dSeq, cSeq)
plot(dSeq, N0Seq)
plot(dSeq, uSeq, type = "l", lwd = 2, xlab = "d", ylab = "U",
     ylim = range(uSeq, uSeq1, uSeq2, uSeq3))
lines(dSeq, uSeq1, col = "red", lwd = 2, lty = 2)
lines(dSeq, uSeq2, col = "green", lwd = 2, lty = 2)
lines(dSeq, uSeq3, col = "blue", lwd = 2, lty = 2)
points(dSeq[c(6, 11, 16)], uSeq[c(6, 11, 16)], pch = 20, cex = 2, col = c("red", "green", "blue"))
# 
# i <- 14
# UtilityAnalysis(ben = ben, dA = dA, d = dSeq[i],
#                 t = tSeq[i], c = cSeq[i], N0 = N0Seq[i],
#                 alpha = alpha, power = power, minpower = minpower, maxpower = maxpower,
#                 heat = TRUE)

# alpha <- 0.05
# power <- 0.8
# muA <- qnorm(power, lower.tail = TRUE) + qnorm(alpha, lower.tail = FALSE)
# mu <- 1
# z <- qnorm(alpha, lower.tail = FALSE)
# t <- 0.6
# tSeq <- seq(0, t, length.out = 500)
# bSeq <- cumsum(c(0, rnorm(length(tSeq) - 1)) * c(0, sqrt(diff(tSeq)))) + mu * tSeq
# b <- bSeq[length(bSeq)]
# cp <- CondPower(muA, b, t, z)
# muTrend <- b / t
# zTrend <- b + (1 - t) * muTrend
# zA <- b + (1 - t) * muA
# plot(tSeq, bSeq, type = "l",
#      xlim = c(0, 1), ylim = range(bSeq, z, zTrend, zA),
#      xlab = "t", ylab = "B(t)", bty = "n", xaxs = "i", yaxs = "i")
# points(t, b, pch = 20, cex = 2, col = "blue")
# abline(h = 0, v = c(0, 1), lty = 1)
# abline(h = z, lty = 2)
# lines(c(0, t, 1), c(0, b, zTrend), lty = 2, lwd = 2, col = "blue")
# lines(c(t, 1), c(b, zA), lty = 2, lwd = 2, col = "red")
