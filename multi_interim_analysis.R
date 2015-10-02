###############################################################################
# Info
#
# Input
#   n0    number(s) of enrolled subjects in Group 0 at interim(s)
#   n1    number(s) of enrolled subjects in Group 1 at interim(s)
#   N0    total sample size of Group 0
#   N1    total sample size of Group 1
#   s0    standard deviation of Group 0, default = 1
#   s1    standard deviation of Group 1, default = 1
#
# Output
#   t     information fraction(s) at interim(s)
###############################################################################
Info <- function(n0, n1, N0, N1, s0 = 1, s1 = 1) {
  v0 <- s0 ^ 2
  v1 <- s1 ^ 2
  t <- (v0 / N0 + v1 / N1) / (v0 / n0 + v1 / n1)
  return(t)
}

###############################################################################
# NInterim
#
# Input
#   t     information fraction(s) at interim(s)
#   N0    total sample size of Group 0
#   r     ratio of n0 : n1 = N0 : N1 = 1 : r, default = 1
#   s0    standard deviation of Group 0, default = 1
#   s1    standard deviation of Group 1, default = 1
#
# Output
#   n0    number(s) of enrolled subjects in Group 0 at interim(s)
#   n1    number(s) of enrolled subjects in Group 1 at interim(s)
#   N0    total sample size of Group 0
#   N1    total sample size of Group 1
###############################################################################
NInterim <- function(t, N0, r = 1, s0 = 1, s1 = 1) {
  v0 <- s0 ^ 2
  v1 <- s1 ^ 2
  N1 <- N0 * r
  n0 <- t * (v0 + v1 / r) / (v0 / N0 + v1 / N1)
  n1 <- n0 * r
  return(list(n0 = round(n0), n1 = round(n1), N0 = round(N0), N1 = round(N1)))
}

###############################################################################
# Drift
#
# Input
#   d     difference in means between 2 groups (Group 1 - Group 0)
#   N0    total sample size of Group 0
#   N1    total sample size of Group 1
#   s0    standard deviation of Group 0, default = 1
#   s1    standard deviation of Group 1, default = 1
#
# Output
#   mu    drift of Brownian motion
###############################################################################
Drift <- function(d, N0, N1, s0 = 1, s1 = 1) {
  mu <- d / sqrt(s0 ^ 2 / N0 + s1 ^ 2 / N1)
  return(mu)
}

###############################################################################
# ProbState
#
# Input
#   mu    drift of Brownian motion
#   t     information fraction(s) at interim(s)
#   lower lower bound(s) of B-value(s) at interim(s), if B < lower then stop
#         for futility
#   upper upper bound(s) of B-value(s) at interim(s), if B > upper then stop
#         for efficacy
#   state state(s) at interim(s), 1 indicates B > upper (efficacy),
#         -1 indicates B < lower (futility), 0 indicates lower <= B <= upper
#
# Output
#   p     probability of such status occurs
#
# Note
#   recursive algorithm
###############################################################################
ProbState <- function(mu, t, lower, upper, state) {
  if (length(t) == 1) {
    if (state == 1) {
      p <- pnorm(q = upper, mean = mu * t, sd = sqrt(t), lower.tail = FALSE)
    } else if (state == -1) {
      p <- pnorm(q = lower, mean = mu * t, sd = sqrt(t), lower.tail = TRUE)
    } else {
      p <- pnorm(q = upper, mean = mu * t, sd = sqrt(t), lower.tail = TRUE) -
        pnorm(q = lower, mean = mu * t, sd = sqrt(t), lower.tail = TRUE)
    }
  } else {
    FunSingle <- function(bSingle) {
      f <- ProbState(mu, t[-1] - t[1],
                     lower[-1] - bSingle, upper[-1] - bSingle, state[-1]) *
        dnorm(x = bSingle, mean = mu * t[1], sd = sqrt(t[1]))
      return(f)
    }
    FunVector <- function(b) {
      f <- sapply(X = b, FUN = FunSingle)
      return(f)
    }
    if (state[1] == 1) {
      p <- integrate(f = FunVector, lower = upper[1], upper = Inf)$value
    } else if (state[1] == -1) {
      p <- integrate(f = FunVector, lower = -Inf, upper = lower[1])$value
    } else {
      p <- integrate(f = FunVector, lower = lower[1], upper = upper[1])$value
    }
  }
  return(p)
}

###############################################################################
# CondPower
#
# Input
#   muA   drift of Brownian motion, under alternative hypothesis
#   t     information fraction(s) at interim(s)
#   lower lower bound(s) of B-value(s) at interim(s), if B < lower then stop
#         for futility
#   upper upper bound(s) of B-value(s) at interim(s), if B > upper then stop
#         for efficacy
#   z     critical B/Z-value in final analysis (this is where lower = upper)
#   b     B-value(s) at interim(s) to calculate conditional power(s),
#         default = lower
#
# Output
#   null  conditional power(s) at interim(s), under null hypothesis
#   alt   conditional power(s) at interim(s), under alternative hypothesis
#   trend conditional power(s) at interim(s), estimating trend by B-value
#   pred  conditional power(s) at interim(s), estimating trend by posterior
###############################################################################
CondPower <- function(muA, t, lower, upper, z, b = lower) {
  n <- length(t)
  null <- alt <- trend <- pred <- numeric(n)
  tExt <- c(t, 1)
  lowerExt <- c(lower, z)
  upperExt <- c(upper, z)
  for (i in 1 : n) {
    for (j in (i + 1) : (n + 1)) {
      null[i] <- null[i] + ProbState(mu = 0,
                                     t = tExt[(i + 1) : j] - t[i],
                                     lower = lowerExt[(i + 1) : j] - b[i],
                                     upper = upperExt[(i + 1) : j] - b[i],
                                     state = c(rep(0, j - i - 1), 1))
      alt[i] <- alt[i] + ProbState(mu = muA,
                                   t = tExt[(i + 1) : j] - t[i],
                                   lower = lowerExt[(i + 1) : j] - b[i],
                                   upper = upperExt[(i + 1) : j] - b[i],
                                   state = c(rep(0, j - i - 1), 1))
      trend[i] <- trend[i] + ProbState(mu = b[i] / t[i],
                                       t = tExt[(i + 1) : j] - t[i],
                                       lower = lowerExt[(i + 1) : j] - b[i],
                                       upper = upperExt[(i + 1) : j] - b[i],
                                       state = c(rep(0, j - i - 1), 1))
      FunSingle <- function(muSingle) {
        f <- ProbState(mu = muSingle,
                       t = tExt[(i + 1) : j] - t[i],
                       lower = lowerExt[(i + 1) : j] - b[i],
                       upper = upperExt[(i + 1) : j] - b[i],
                       state = c(rep(0, j - i - 1), 1)) *
          dnorm(muSingle, b[i] / t[i], 1 / sqrt(t[i]))
        return(f)
      }
      FunVector <- function(mu) {
        f <- sapply(X = mu, FUN = FunSingle)
        return(f)
      }
      pred[i] <- pred[i] + integrate(f = FunVector,
                                     lower = -Inf, upper = Inf)$value
    }
  }
  return(list(null = null, alt = alt, trend = trend, pred = pred))
}