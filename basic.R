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
#   mu: scalar, difference(s) in means (group 1 - group 0)
#   t: information fraction(s) at interim(s)
#   lower, upper: scalar/vector, lower/upper bound(s) (B-value) at interim(s),
#     exceeding bounds means stopping for futility/efficacy
#   state: scalar/vector, state(s) at interim(s), takes value -1/0/+1, meaning
#     futility/continue/efficacy respectively
#   method: "exact" (default) or "simulation", method to compute probability,
#     exact method implements a recursive algorithm, therefore simulation is
#     recommended for large number of interims
#
# Details
#   Exact method implements a recursive algorithm, therefore simulation is
#   recommended for large number of interims.
#
# Value
#   scalar, probability of state configuration

ProbState <- function(mu, t, lower, upper, state, method = "exact")
{
  
}


###############################################################################
# Test Area

n0 <- c(200, 300)
n1 <- c(200, 300)

N0 <- 500
N1 <- 500

t <- InformationFraction(n0, n1, N0, N1)

d <- 0.3

mu <- Drift(d, N0, N1)