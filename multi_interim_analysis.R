###############################################################################
# Info
#
# Input
#   n0      number(s) of enrolled subjects in Group 0 at interim(s)
#   n1      number(s) of enrolled subjects in Group 1 at interim(s)
#   n0.tot  total sample size of Group 0
#   n1.tot  total sample size of Group 1
#   s0      standard deviation of Group 0, default = 1
#   s1      standard deviation of Group 1, default = 1
#
# Output
#   t     information fraction(s) at interim(s)
###############################################################################
Info <- function(n0, n1, n0.tot, n1.tot, s0 = 1, s1 = 1) {
  v0 <- s0 ^ 2
  v1 <- s1 ^ 2
  t <- (v0 / n0.tot + v1 / n1.tot) / (v0 / n0 + v1 / n1)
  return(t)
}

###############################################################################
# NInterim
#
# Input
#   t       information fraction(s) at interim(s)
#   n.tot   total sample size
#   r       ratio of n0 : n1 = n0.tot : n1.tot = 1 : r, default = 1
#   s0      standard deviation of Group 0, default = 1
#   s1      standard deviation of Group 1, default = 1
#
# Output
#   n0      number(s) of enrolled subjects in Group 0 at interim(s)
#   n1      number(s) of enrolled subjects in Group 1 at interim(s)
#   n0.tot  total sample size of Group 0
#   n1.tot  total sample size of Group 1
###############################################################################
NInterim <- function(t, n.tot, r = 1, s0 = 1, s1 = 1) {
  v0 <- s0 ^ 2
  v1 <- s1 ^ 2
  n0.tot <- n.tot / (1 + r)
  n1.tot <- n0.tot * r
  n0 <- t * (v0 + v1 / r) / (v0 / n0.tot + v1 / n1.tot)
  n1 <- n0 * r
  return(list(n0 = round(n0), n1 = round(n1), n0.tot = round(n0.tot), n1.tot = round(n1.tot)))
}

###############################################################################
# Drift
#
# Input
#   d       difference in means between 2 groups (Group 1 - Group 0)
#   n0.tot  total sample size of Group 0
#   n1.tot  total sample size of Group 1
#   s0      standard deviation of Group 0, default = 1
#   s1      standard deviation of Group 1, default = 1
#
# Output
#   mu      drift of Brownian motion
###############################################################################
Drift <- function(d, n0.tot, n1.tot, s0 = 1, s1 = 1) {
  mu <- d / sqrt(s0 ^ 2 / n0.tot + s1 ^ 2 / n1.tot)
  return(mu)
}

###############################################################################
# ProbState
#
# Input
#   mu      drift of Brownian motion
#   t       information fraction(s) at interim(s)
#   lower   lower bound(s) of B-value(s) at interim(s), if B < lower then stop
#           for futility
#   upper   upper bound(s) of B-value(s) at interim(s), if B > upper then stop
#           for efficacy
#   state   state(s) at interim(s), 1 indicates B > upper (efficacy),
#           -1 indicates B < lower (futility), 0 indicates lower <= B <= upper
#
# Output
#   p       probability of such status occurs
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
    FunSingle <- function(b.single) {
      f <- ProbState(mu, t[-1] - t[1],
                     lower[-1] - b.single, upper[-1] - b.single, state[-1]) *
        dnorm(x = b.single, mean = mu * t[1], sd = sqrt(t[1]))
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
#   mu.alt  drift of Brownian motion, under alternative hypothesis
#   t       information fraction(s) at interim(s)
#   lower   lower bound(s) of B-value(s) at interim(s), if B < lower then stop
#           for futility
#   upper   upper bound(s) of B-value(s) at interim(s), if B > upper then stop
#           for efficacy
#   z       critical B/Z-value in final analysis (this is where lower = upper)
#   b       B-value(s) at interim(s) to calculate conditional power(s),
#           default = lower
#
# Output
#   conditional.power
#           a data.frame consists of conditional powers
###############################################################################
CondPower <- function(mu.alt, t, lower, upper, z, b = lower) {
  n <- length(t)
  null <- alt <- trend <- pred <- numeric(n)
  t.ext <- c(t, 1)
  lower.ext <- c(lower, z)
  upper.ext <- c(upper, z)
  for (i in 1 : n) {
    for (j in (i + 1) : (n + 1)) {
      null[i] <- null[i] + ProbState(mu = 0,
                                     t = t.ext[(i + 1) : j] - t[i],
                                     lower = lower.ext[(i + 1) : j] - b[i],
                                     upper = upper.ext[(i + 1) : j] - b[i],
                                     state = c(rep(0, j - i - 1), 1))
      alt[i] <- alt[i] + ProbState(mu = mu.alt,
                                   t = t.ext[(i + 1) : j] - t[i],
                                   lower = lower.ext[(i + 1) : j] - b[i],
                                   upper = upper.ext[(i + 1) : j] - b[i],
                                   state = c(rep(0, j - i - 1), 1))
      trend[i] <- trend[i] + ProbState(mu = b[i] / t[i],
                                       t = t.ext[(i + 1) : j] - t[i],
                                       lower = lower.ext[(i + 1) : j] - b[i],
                                       upper = upper.ext[(i + 1) : j] - b[i],
                                       state = c(rep(0, j - i - 1), 1))
      FunSingle <- function(mu.single) {
        f <- ProbState(mu = mu.single,
                       t = t.ext[(i + 1) : j] - t[i],
                       lower = lower.ext[(i + 1) : j] - b[i],
                       upper = upper.ext[(i + 1) : j] - b[i],
                       state = c(rep(0, j - i - 1), 1)) *
          dnorm(mu.single, b[i] / t[i], 1 / sqrt(t[i]))
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
  conditional.power <- data.frame(info.fraction = t,
                                  futility.bound = lower,
                                  efficacy.bound = upper,
                                  b.value = b,
                                  cp.null = null,
                                  cp.alt = alt,
                                  cp.trend = trend,
                                  cp.pred = pred)
  return(conditional.power)
}

###############################################################################
# B2C
#
# Input
#   b       B-value(s) at interim(s)
#   t       information fraction(s) at interim(s)
#
# Output
#   z       Z-value(s) at interim(s)
###############################################################################
B2C <- function(b, t) {
  z <- b / sqrt(t)
  z[t == 0] <- 0
  return(z)
}

###############################################################################
# Design
#
# Input
#   d.alt   difference in means between 2 groups, under alternative hypoethesis
#   d.eval  difference in means as scenarios for performance evaluation,
#           default = c(0, d.alt)
#   t       information fraction(s) at interim(s)
#   lower   lower bound(s) of B-value(s) at interim(s), if B < lower then stop
#           for futility
#   upper   upper bound(s) of B-value(s) at interim(s), if B > upper then stop
#           for efficacy
#   z       critical B/Z-value in final analysis (this is where lower = upper)
#   n.tot   total sample size
#   r       ratio of n0 : n1 = n0.tot : n1.tot = 1 : r, default = 1
#   s0      standard deviation of Group 0, default = 1
#   s1      standard deviation of Group 1, default = 1
#   Cost    cost function, default = NULL
#             input: t, n.tot
#             output: cost(s) at interim(s)
#   Benefit benefit function, default = NULL
#             input: t, n.tot
#             output: benefit(s) at interim(s)
#
# Output
###############################################################################
Design <- function(d.alt, d.eval = c(0, d.alt), t, lower, upper, z, n.tot,
                   r = 1, s0 = 1, s1 = 1, Cost = NULL, Benefit = NULL) {
  
}