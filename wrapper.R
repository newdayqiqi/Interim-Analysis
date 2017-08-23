source("basic.R")

Cost <- function(n0, n1)
{
  return(0.1 * (n0 + n1))
}

Benefit <- function(n0, n1)
{
  return(300)
}

UFun <- function(n0, n1, state)
{
  cost <- Cost(n0, n1)
  benefit <- Benefit(n0, n1)
  u <- benefit * (state == 1) - cost
  return(u)
}

UtilityWrapper <- function(d, n, boundary, N = 450, z = qnorm(0.025, lower.tail = FALSE))
{
  if (length(d) == 1)
  {
    out <- Utility(d, n0 = n, n1 = n, N0 = N, N1 = N, lower = boundary, upper = Inf, z = z, UFun = UFun)
    u <- out$u
    power <- out$power
  }
  else
  {
    out <- sapply(d, Utility, n0 = n, n1 = n, N0 = N, N1 = N, lower = boundary, upper = Inf, z = z, UFun = UFun)
    u <- mean(unlist(out[1, ]))
    power <- mean(unlist(out[2, ]))
  }
  return(list(u = u, power = power))
}

UtilityGridSearch <- function(d, n.grid, boundary.grid, N = 450, z = qnorm(0.025, lower.tail = FALSE))
{
  u <- power <- matrix(NA, length(n.grid), length(boundary.grid))
  for (i in seq_along(n.grid))
  {
    cat(i, "\t")
    n <- n.grid[i]
    for (j in seq_along(boundary.grid))
    {
      boundary <- boundary.grid[j]
      out <- UtilityWrapper(d, n, boundary, N, z)
      u[i, j] <- out$u
      power[i, j] <- out$power
      cat(".")
    }
    cat("\n")
  }
  return(list(u = u, power = power))
}