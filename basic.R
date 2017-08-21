###############################################################################
# Information
#
# Argument
#   n0, n1: scalar/vector, sample size(s) of group 0/1
#   s0, s1: scalar, sdandard deviation of group 0/1, default = 1
#
# Details
#   Information is defined as the inverse of variance (of estimated difference):
#   I(n0, n1) = 1 / s(n0, n1)^2, where
#   s(n0, n1)^2 = s0^2/n0 + s1^2/n1
#
# Value
#   scalar/vector, information

Information <- function(n0, n1, s0 = 1, s1 = 1)
{
  info <- 1 / (s0^2 / n0 + s1^2 / n1)
  return(info)
}

###############################################################################
# InformationFraction
#
# Argument
#   n0, n1: scalar/vector, sample size(s) of group 0/1 at interim(s)
#   N0, N1: scalar, sample size of gorup 0/1 at final
#   s0, s1: scalar, sdandard deviation of group 0/1, default = 1
#
# Details
#   Information fraction is defined as:
#   t(n0, n1) = I(n0, n1) / I(N0, N1)
#
# Value
#   scalar/vector, information fraction at interim(s)

InformationFraction <- function(n0, n1, N0, N1, s0 = 1, s1 = 1)
{
  t <- Information(n0, n1, s0, s1) / Information(N0, N1, s0, s1)
  return(t)
}

###############################################################################
# Drift
#
# Argument
#   d: scalar/vector, difference(s) in means (group 1 - group 0)
#   N0, N1: scalar, sample size of gorup 0/1 at final
#   s0, s1: scalar, sdandard deviation of group 0/1, default = 1
#
# Details
#   Drift (mu) is computed by the following formula:
#   mu = d / sqrt(s0^2/N0 + s1^2/N1)
#
# Value
#   scalar/vector, drift(s) of Brownian motion

Drift <- function(d, N0, N1, s0 = 1, s1 = 1)
{
  mu <- d / sqrt(s0^2 / N0 + s1^2 / N1)
  return(mu)
}

###############################################################################
# ProbState
#
# Argument
#   mu: scalar, drift of Brownian motion
#   t: information fraction(s) at interim(s)/final
#   lower, upper: scalar/vector, lower/upper B-bound(s) at interim(s)/final,
#     exceeding bounds means stopping for futility/efficacy
#   state: scalar/vector, state(s) at interim(s), takes value -1/0/+1, meaning
#     futility/continue/efficacy respectively
#
# Details
#   Recursive algorithm.
#
# Value
#   scalar, probability of state configuration

ProbState <- function(mu, t, lower, upper, state)
{
  if (length(t) == 1)
  {
    if (state == 1)
    {
      p <- 1 - pnorm(upper, mean = mu * t, sd = sqrt(t))
    }
    else if (state == 0)
    {
      p <- pnorm(upper, mean = mu * t, sd = sqrt(t)) - pnorm(lower, mean = mu * t, sd = sqrt(t))
    }
    else if (state == -1)
    {
      p <- pnorm(lower, mean = mu * t, sd = sqrt(t))
    }
    else
    {
      warning("state not recognized")
      p <- NA
    }
    return(p)
  }
  else
  {
    DensityFirstStep <- function(b) # b must be scalar
    {
      f <- ProbState(mu, t = t[-1] - t[1], lower = lower[-1] - b, upper = upper[-1] - b, state = state[-1]) *
        dnorm(b, mean = mu * t[1], sd = sqrt(t[1]))
      return(f)
    }
    VecDensityFirstStep <- function(b) # vectorized version
    {
      f <- sapply(b, FUN = DensityFirstStep)
      return(f)
    }
    if (state[1] == 1)
    {
      p <- integrate(VecDensityFirstStep, lower = upper[1], upper = Inf)$value
    }
    else if (state[1] == 0)
    {
      p <- integrate(VecDensityFirstStep, lower = lower[1], upper = upper[1])$value
    }
    else if (state[1] == -1)
    {
      p <- integrate(VecDensityFirstStep, lower = -Inf, upper = lower[1])$value
    }
    else
    {
      warning("state not recognized")
      p <- NA
    }
  }
  return(p)
}

###############################################################################
# Utility
#
# Argument
#   d: scalar, true difference(s) in means (group 1 - group 0)
#   n0, n1: scalar/vector, sample size(s) of group 0/1 at interim(s)
#   N0, N1: scalar, sample size of gorup 0/1 at final
#   lower, upper: scalar/vector, lower/upper bound(s) (B-value) at interim(s),
#     exceeding bounds means stopping for futility/efficacy
#   z: critical B/Z-value at final (where lower = upper)
#   UFun: function, returns utility given sample sizes and state 
#   s0, s1: scalar, sdandard deviation of group 0/1
#
# Details
#   Compute expected utility.
#
# Value (list)
#   u: scalar, utility
#   u.state: matrix, utility at each state
#   p.state: matrix, probability at each state

Utility <- function(d, n0, n1, N0, N1, lower, upper, z, UFun, s0 = 1, s1 = 1)
{
  # design parameter
  mu <- Drift(d, N0, N1, s0, s1)
  t <- InformationFraction(n0, n1, N0, N1, s0, s1)
  k <- length(t)
  
  # u.state
  u.state.success <- UFun(c(n0, N0), c(n1, N1), 1)
  u.state.fail <- UFun(c(n0, N0), c(n1, N1), -1)
  u.state <- rbind(u.state.success, u.state.fail)
  
  # p.state
  p.state <- matrix(NA, 2, k + 1)
  for (i in 1 : k)
  {
    p.state[1, i] <- ProbState(mu, t[1 : i], lower[1 : i], upper[1 : i], c(numeric(i - 1), 1))
    p.state[2, i] <- ProbState(mu, t[1 : i], lower[1 : i], upper[1 : i], c(numeric(i - 1), -1))
  }
  p.state[1, k + 1] <- ProbState(mu, c(t, 1), c(lower, z), c(upper, z), c(numeric(k), 1))
  p.state[2, k + 1] <- ProbState(mu, c(t, 1), c(lower, z), c(upper, z), c(numeric(k), -1))
  
  # u
  u <- sum(u.state * p.state)
  
  return(list(u = u, u.state = u.state, p.state = p.state))
}


# ###############################################################################
# # Test Area
# 
# n0 <- c(200, 322)
# n1 <- c(200, 322)
# 
# N0 <- 450
# N1 <- 450
# 
# d <- 0.10
# 
# lower <- c(0.30, 0.60)
# upper <- c(Inf, Inf)
# 
# z <- 1.96
# 
# Cost <- function(n0, n1)
# {
#   return(0.1 * (n0 + n1))
# }
# 
# Benefit <- function(n0, n1)
# {
#   return(300)
# }
# 
# UFun <- function(n0, n1, state)
# {
#   cost <- Cost(n0, n1)
#   benefit <- Benefit(n0, n1)
#   u <- benefit * (state == 1) - cost
#   return(u)
# }
# 
# z <- 1.96
# 
# out <- Utility(d, n0, n1, N0, N1, lower, upper, z, UFun)