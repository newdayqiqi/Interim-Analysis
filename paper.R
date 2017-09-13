source("wrapper.R")
require(ggplot2)
require(scales)

###############################################################################
# set parameters
N <- 450
n.grid <- seq(45, 360, 3) # t from 0.1 to 0.8
boundary.grid <- seq(-0.5, 1.0, 0.01)
d.seq <- seq(0.00, 0.30, 0.05)
t.grid <- n.grid / N

# ###############################################################################
# # run grid search
# for (d in d.seq)
# {
#   out <- UtilityGridSearch(d, n.grid, boundary.grid)
#   file <- paste("results/grid_d_", d, ".RData", sep = "")
#   save(out, file = file)
#   cat("Finished: d =", d, "\n")
# }

###############################################################################
# Table 1

table1 <- NULL
for (d in d.seq)
{
  file <- paste("results/grid_d_", d, ".RData", sep = "")
  load(file)
  u.opt <- max(out$u)
  ind <- which(out$u == u.opt, arr.ind = TRUE)
  t.opt <- t.grid[ind[1]]
  boundary.opt <- boundary.grid[ind[2]]
  z.opt <- boundary.opt / sqrt(t.opt)
  table1 <- rbind(table1,
                  data.frame(d = d, t = t.opt, b = boundary.opt, z = z.opt, u = u.opt))
}
print(table1)

###############################################################################
# Figure 1

d.prior <- seq(-0.3, 0.7, 0.005)
d.weight <- dnorm(d.prior, 0.2, 1 / sqrt(50))
plot(d.prior, d.weight)

# out <- UtilityGridSearch(rbind(d.prior, d.weight), n.grid, boundary.grid)
# file <- "results/grid_d_prior.RData"
# save(out, file = file)
# cat("Finished: d = prior", "\n")

file <- "results/grid_d_prior.RData"
load(file)
u.opt <- max(out$u)
ind <- which(out$u == u.opt, arr.ind = TRUE)
t.opt <- t.grid[ind[1]]
boundary.opt <- boundary.grid[ind[2]]

dat.figure1 <- data.frame(boundary = rep(boundary.grid, each = length(t.grid)),
                          t = rep(t.grid, times = length(boundary.grid)),
                          u = as.vector(out$u))

file <- "results/grid_d_0.2.RData"
load(file)

dat.figure1$power <- as.vector(out$power)

dat.figure1 <- subset(dat.figure1, subset = (boundary >= -0.25) & (boundary <= 0.75) & (t >= 0.2))
dat.figure1.text <- data.frame(x = c(0.45, 0.43, 0.39, 0.29),
                               y = c(0.30, 0.35, 0.42, 0.71),
                               label = c("0.70", "0.75", "0.80", "0.85"))

Trans1 <- function(x)
{
  y <- 1.1 - log(117.5 - x)
  return(y)
}
InvTrans1 <- function(y)
{
  x <- 117.5 - exp(1.1 - y)
  return(x)
}
trans1_trans <- function() trans_new(name = "trans1", transform = "Trans1", inverse = "InvTrans1")

hist(Trans1(dat.figure1$u), breaks = 50)

figure1 <- ggplot(data = dat.figure1, mapping = aes(x = boundary, y = t)) +
  geom_raster(mapping = aes(fill = u), interpolate = TRUE) +
  scale_fill_gradientn(colours = c("#400000", "#800000", "#FF0000", "#FF4000", "#FF8000", "#FFFF00", "#FFFFFF"),
                       trans = "trans1",
                       breaks = c(80, 90, 100, 110, 111, 112, 113, 114, 115, 116, 117),
                       guide = guide_colorbar(barheight = 25, title = "Utility\n ($M)")) +
  geom_hline(yintercept = t.opt, linetype = 2, size = 1) +
  geom_vline(xintercept = boundary.opt, linetype = 2, size = 1) +
  geom_contour(mapping = aes(z = power), breaks = seq(0.70, 0.85, 0.05), colour = "blue", linetype = 2, size = 1) +
  geom_text(data = dat.figure1.text, mapping = aes(x = x, y = y, label = label), colour = "blue") +
  labs(x = "B-value", y = "Information Fraction") +
  theme_minimal()
  # geom_contour(mapping = aes(z = z), breaks = seq(-1, 1, 0.1), colour = "green", linetype = 4, size = 1)

###############################################################################
# Figure 2

file <- "results/grid_d_prior.RData"
load(file)
ind <- apply(out$u, 1, which.max)
ind <- (ind - 1) * length(t.grid) + seq_along(t.grid)

file <- "results/grid_d_0.RData"
load(file)

dat.figure2.alpha <- data.frame(t = t.grid, p = 0.5 + 20 * out$power[ind], type = "alpha") # axis1 = 0.5 + 20 * axis2; axis2 = -0.025 + 0.05 * axis1

file <- "results/grid_d_0.2.RData"
load(file)

dat.figure2.power <- data.frame(t = t.grid, p = out$power[ind], type = "power")

dat.figure2 <- rbind(dat.figure2.alpha, dat.figure2.power)

figure2 <- ggplot(data = dat.figure2, mapping = aes(x = t, y = p, group = type, color = type)) +
  geom_vline(xintercept = t.opt, linetype = 2, size = 1) +
  geom_line(size = 1) +
  scale_colour_manual(name = NULL,values =c("blue", "red"),
                      labels = c("Overall power at optimal futility boundary", "Type I error at optimal futility boundary")) +
  geom_hline(yintercept = 0.85, linetype = 2, size = 1, colour = "red") +
  geom_hline(yintercept = 1, linetype = 2, size = 1, colour = "blue") +
  coord_cartesian(ylim = c(0.50, 1.00)) +
  scale_y_continuous(sec.axis = sec_axis(~ . * 0.05 - 0.025, name = "Type I Error Rate")) +
  labs(x = "Information Fraction", y = "Overall Power") +
  theme_minimal() %+replace%
  theme(legend.position = c(0.8, 0.18))

###############################################################################
# Figure 5a

d.prior2 <- seq(0, 0.3, 0.002)

# out <- UtilityGridSearch(d.prior2, n.grid, boundary.grid)
# file <- "results/grid_d_prior2.RData"
# save(out, file = file)
# cat("Finished: d = prior2", "\n")

file <- "results/grid_d_prior2.RData"
load(file)
u.opt <- max(out$u)
ind <- which(out$u == u.opt, arr.ind = TRUE)
t.opt <- t.grid[ind[1]]
boundary.opt <- boundary.grid[ind[2]]

dat.figure5a <- data.frame(boundary = rep(boundary.grid, each = length(t.grid)),
                           t = rep(t.grid, times = length(boundary.grid)),
                           u = as.vector(out$u))

file <- "results/grid_d_0.2.RData"
load(file)

dat.figure5a$power <- as.vector(out$power)

dat.figure5a <- subset(dat.figure5a, subset = (boundary >= -0.25) & (boundary <= 0.75) & (t >= 0.2))
dat.figure5a.text <- data.frame(x = c(0.90, 0.85, 0.76, 0.42),
                                y = c(0.47, 0.51, 0.57, 0.75),
                                label = c("0.70", "0.75", "0.80", "0.85"))

Trans5a <- function(x)
{
  y <- 1.6 - log(86 - x)
  return(y)
}
InvTrans5a <- function(y)
{
  x <- 86 - exp(1.6 - y)
  return(x)
}
trans5a_trans <- function() trans_new(name = "trans5a", transform = "Trans5a", inverse = "InvTrans5a")

hist(Trans5a(dat.figure5a$u), breaks = 50)

figure5a <- ggplot(data = dat.figure5a, mapping = aes(x = boundary, y = t)) +
  geom_raster(mapping = aes(fill = u), interpolate = TRUE) +
  scale_fill_gradientn(colours = c("#400000", "#800000", "#FF0000", "#FF4000", "#FF8000", "#FFFF00", "#FFFFFF"),
                       trans = "trans5a",
                       breaks = c(40, 50, 60, 70, 80, 81, 82, 83, 84, 85),
                       guide = guide_colorbar(barheight = 25, title = "Utility\n ($M)")) +
  geom_hline(yintercept = t.opt, linetype = 2, size = 1) +
  geom_vline(xintercept = boundary.opt, linetype = 2, size = 1) +
  geom_contour(mapping = aes(z = power), breaks = seq(0.70, 0.85, 0.05), linetype = 2, size = 1) +
  geom_text(data = dat.figure1.text, mapping = aes(x = x, y = y, label = label), colour = "blue") +
  labs(x = "B-value", y = "Information Fraction") +
  theme_minimal()

###############################################################################
# Figure 5b

file <- "results/grid_d_prior2.RData"
load(file)
ind <- apply(out$u, 1, which.max)
ind <- (ind - 1) * length(t.grid) + seq_along(t.grid)

file <- "results/grid_d_0.RData"
load(file)

dat.figure5b.alpha <- data.frame(t = t.grid, p = 0.5 + 20 * out$power[ind], type = "alpha") # axis1 = 0.5 + 20 * axis2; axis2 = -0.025 + 0.05 * axis1

file <- "results/grid_d_0.2.RData"
load(file)

dat.figure5b.power <- data.frame(t = t.grid, p = out$power[ind], type = "power")

dat.figure5b <- rbind(dat.figure5b.alpha, dat.figure5b.power)

figure5b <- ggplot(data = dat.figure5b, mapping = aes(x = t, y = p, group = type, color = type)) +
  geom_vline(xintercept = t.opt, linetype = 2, size = 1) +
  geom_line(size = 1) +
  scale_colour_manual(name = NULL,values =c("blue", "red"),
                      labels = c("Overall power at optimal futility boundary", "Type I error at optimal futility boundary")) +
  geom_hline(yintercept = 0.85, linetype = 2, size = 1, colour = "red") +
  geom_hline(yintercept = 1, linetype = 2, size = 1, colour = "blue") +
  coord_cartesian(ylim = c(0.50, 1.00)) +
  scale_y_continuous(sec.axis = sec_axis(~ . * 0.05 - 0.025, name = "Type I Error Rate")) +
  labs(x = "Information Fraction", y = "Overall Power") +
  theme_minimal() %+replace%
  theme(legend.position = c(0.8, 0.18))

###############################################################################
# Figure 3

dat.figure3 <- NULL
for (scenario in d.seq)
{
  # Bayes design
  file <- "results/grid_d_prior.RData"
  load(file)
  u.opt <- max(out$u)
  ind <- which(out$u == u.opt, arr.ind = TRUE)
  n.opt <- n.grid[ind[1]]
  boundary.opt <- boundary.grid[ind[2]]
  out <- UtilityWrapper(scenario, n.opt, boundary.opt)
  dat.figure3 <- rbind(dat.figure3,
                       data.frame(design = paste("Optimal design under d ~ N(0.2,", round(sqrt(1/50), 3) ,")", sep = ""), true = scenario, u = out$u))
  
  # Fixed design
  out <- Utility(scenario, numeric(0), numeric(0), N, N, numeric(0), numeric(0), qnorm(0.025, lower.tail = FALSE), UFun)
  dat.figure3 <- rbind(dat.figure3,
                       data.frame(design = "Fixed design", true = scenario, u = out$u))
  
  # Different design
  d.design <- c(0, 0.05, 0.10, 0.15, 0.30)
  for (d in d.design)
  {
    file <- paste("results/grid_d_", d, ".RData", sep = "")
    load(file)
    u.opt <- max(out$u)
    ind <- which(out$u == u.opt, arr.ind = TRUE)
    n.opt <- n.grid[ind[1]]
    boundary.opt <- boundary.grid[ind[2]]
    out <- UtilityWrapper(scenario, n.opt, boundary.opt)
    dat.figure3 <- rbind(dat.figure3,
                         data.frame(design = paste("Optimal design under d =", d), true = scenario, u = out$u))
  }
  
  # Local optimal
  u <- max(subset(dat.figure3, subset = (true == scenario))$u)
  dat.figure3 <- rbind(dat.figure3,
                       data.frame(design = "Local optimal", true = scenario, u = u))
}
dat.figure3.text1 <- data.frame(x = 0.1 * c(-0.05, -0.05, -0.05, 3.03, 3.03),
                               y = c(-42, -60, -70, 0, 38),
                               text = c(3, 4, 7, 1, 2))
dat.figure3.text2 <- data.frame(x = 0.1 * 0.03,
                                y = seq(144, 95, length.out = 5),
                                text = c(1, 2, 3, 4, 7))

figure3 <- ggplot(data = dat.figure3, mapping = aes(x = true, y = u)) +
  geom_line(mapping = aes(group = design, color = design, linetype = design, size = design)) +
  scale_colour_manual(name = NULL,
                      values = c("#00FF00", "#0000FF", rep("#BBBBBB", length(d.design)), "#FF0000")) +
  scale_linetype_manual(name = NULL,
                        values = c(rep(1, length(d.design) + 2), 2)) +
  scale_size_manual(name = NULL,
                    values = c(1.2, 1.2, rep(1, length(d.design)), 1.2)) +
  labs(x = "True treatment difference", y = "Expected net benefit ($M)") +
  geom_text(data = dat.figure3.text1, mapping = aes(x = x, y = y, label = text)) +
  geom_text(data = dat.figure3.text2, mapping = aes(x = x, y = y, label = text)) +
  theme_minimal() %+replace%
  theme(legend.position = c(0.22, 0.7))

###############################################################################
# Figure 4 (benefit)

d.prior <- seq(-0.3, 0.7, 0.005)
d.weight <- dnorm(d.prior, 0.2, 1 / sqrt(50))
plot(d.prior, d.weight)

# Benefit <- function(n0, n1)
# {
#   return(150)
# }
# 
# out <- UtilityGridSearch(rbind(d.prior, d.weight), n.grid, boundary.grid)
# file <- "results/grid_d_prior_benefit_150.RData"
# save(out, file = file)
# cat("Finished: d = prior, benefit = 150", "\n")

file <- "results/grid_d_prior_benefit_150.RData"
load(file)
u.opt <- max(out$u)
ind <- which(out$u == u.opt, arr.ind = TRUE)
t.opt <- t.grid[ind[1]]
boundary.opt <- boundary.grid[ind[2]]

dat.figure4a <- data.frame(boundary = rep(boundary.grid, each = length(t.grid)),
                          t = rep(t.grid, times = length(boundary.grid)),
                          u = as.vector(out$u))

file <- "results/grid_d_0.2.RData"
load(file)

dat.figure4a$power <- as.vector(out$power)

dat.figure4a <- subset(dat.figure4a, subset = (boundary >= -0.25) & (boundary <= 0.75))
dat.figure4a.text <- data.frame(x = c(0.45, 0.43, 0.39, 0.29),
                               y = c(0.30, 0.35, 0.42, 0.71),
                               label = c("0.70", "0.75", "0.80", "0.85"))

Trans4a <- function(x)
{
  y <- 0.85 - log(20 - x)
  return(y)
}
InvTrans4a <- function(y)
{
  x <- 20 - exp(0.85 - y)
  return(x)
}
trans4a_trans <- function() trans_new(name = "trans4a", transform = "Trans4a", inverse = "InvTrans4a")

hist(Trans4a(dat.figure4a$u), breaks = 50)

figure4a <- ggplot(data = dat.figure4a, mapping = aes(x = boundary, y = t)) +
  geom_raster(mapping = aes(fill = u), interpolate = TRUE) +
  scale_fill_gradientn(colours = c("#400000", "#800000", "#FF0000", "#FF4000", "#FF8000", "#FFFF00", "#FFFFFF"),
                       trans = "trans4a",
                       breaks = c(5, 10, 15, 16, 17, 18, 19),
                       guide = guide_colorbar(barheight = 25, title = "Utility\n ($M)")) +
  geom_hline(yintercept = t.opt, linetype = 2, size = 1) +
  geom_vline(xintercept = boundary.opt, linetype = 2, size = 1) +
  geom_contour(mapping = aes(z = power), breaks = seq(0.70, 0.85, 0.05), linetype = 2, size = 1) +
  geom_text(data = dat.figure1.text, mapping = aes(x = x, y = y, label = label), colour = "blue") +
  labs(x = "B-value", y = "Information Fraction") +
  theme_minimal()

# Benefit <- function(n0, n1)
# {
#   return(1200)
# }
# 
# out <- UtilityGridSearch(rbind(d.prior, d.weight), n.grid, boundary.grid)
# file <- "results/grid_d_prior_benefit_1200.RData"
# save(out, file = file)
# cat("Finished: d = prior, benefit = 1200", "\n")

file <- "results/grid_d_prior_benefit_1200.RData"
load(file)
u.opt <- max(out$u)
ind <- which(out$u == u.opt, arr.ind = TRUE)
t.opt <- t.grid[ind[1]]
boundary.opt <- boundary.grid[ind[2]]

dat.figure4b <- data.frame(boundary = rep(boundary.grid, each = length(t.grid)),
                           t = rep(t.grid, times = length(boundary.grid)),
                           u = as.vector(out$u))

file <- "results/grid_d_0.2.RData"
load(file)

dat.figure4b$power <- as.vector(out$power)

dat.figure4b <- subset(dat.figure4b, subset = (boundary >= -0.25) & (boundary <= 0.75))
dat.figure4b.text <- data.frame(x = c(0.45, 0.43, 0.39, 0.29),
                                y = c(0.30, 0.35, 0.42, 0.71),
                                label = c("0.70", "0.75", "0.80", "0.85"))

Trans4b <- function(x)
{
  y <- 2 - log(719.5 - x)
  return(y)
}
InvTrans4b <- function(y)
{
  x <- 719.5 - exp(2 - y)
  return(x)
}
trans4b_trans <- function() trans_new(name = "trans4b", transform = "Trans4b", inverse = "InvTrans4b")

hist(Trans4b(dat.figure4b$u), breaks = 50)

figure4b <- ggplot(data = dat.figure4b, mapping = aes(x = boundary, y = t)) +
  geom_raster(mapping = aes(fill = u), interpolate = TRUE) +
  scale_fill_gradientn(colours = c("#400000", "#800000", "#FF0000", "#FF4000", "#FF8000", "#FFFF00", "#FFFFFF"),
                       trans = "trans4b",
                       breaks = c(300, 400, 500, 600, 700, 710, 715, 716, 717, 718, 719),
                       guide = guide_colorbar(barheight = 25, title = "Utility\n ($M)")) +
  geom_hline(yintercept = t.opt, linetype = 2, size = 1) +
  geom_vline(xintercept = boundary.opt, linetype = 2, size = 1) +
  geom_contour(mapping = aes(z = power), breaks = seq(0.70, 0.85, 0.05), linetype = 2, size = 1) +
  geom_text(data = dat.figure1.text, mapping = aes(x = x, y = y, label = label), colour = "blue") +
  labs(x = "B-value", y = "Information Fraction") +
  theme_minimal()

###############################################################################
# Figure 6 (cost)

Benefit <- function(n0, n1)
{
  return(300)
}

Cost <- function(n0, n1)
{
  return(0.1 * (n0 + n1))
}

n.seq <- 0 : 450
cost.seq <- Cost(n.seq, n.seq)
dat.figure6a.linear <- data.frame(size = n.seq * 2, cost = cost.seq, type = "Linear")

Cost <- function(n0, n1)
{
  # black magic
  x <- 90 / (log(450 * 2 + 100) - log(100))
  y <- - x * log(100)
  cost <- log(n0 + n1 + 100) * x + y
  return(cost)
}

# out <- UtilityGridSearch(rbind(d.prior, d.weight), n.grid, boundary.grid)
# file <- "results/grid_d_prior_cost2.RData"
# save(out, file = file)
# cat("Finished: d = prior, concave cost", "\n")

cost.seq2 <- Cost(n.seq, n.seq)
dat.figure6a.concave <- data.frame(size = n.seq * 2, cost = cost.seq2, type = "Concave")

dat.figure6a <- rbind(dat.figure6a.linear, dat.figure6a.concave)

figure6a <- ggplot(data = dat.figure6a, mapping = aes(x = size, y = cost, group = type, color = type)) +
  geom_line(size = 1) +
  scale_colour_manual(name = NULL,values =c("blue", "red"),
                      labels = c("Linear cost-spending function", "Concave cost-spending function")) +
  labs(x = "Sample size available at interim analysis", y = "Cost spent ($M)") +
  theme_minimal() %+replace%
  theme(legend.position = c(0.8, 0.18))

file <- "results/grid_d_prior_cost2.RData"
load(file)
u.opt <- max(out$u)
ind <- which(out$u == u.opt, arr.ind = TRUE)
t.opt <- t.grid[ind[1]]
boundary.opt <- boundary.grid[ind[2]]

dat.figure6b <- data.frame(boundary = rep(boundary.grid, each = length(t.grid)),
                           t = rep(t.grid, times = length(boundary.grid)),
                           u = as.vector(out$u))

file <- "results/grid_d_0.2.RData"
load(file)

dat.figure6b$power <- as.vector(out$power)

dat.figure6b <- subset(dat.figure6b, subset = (boundary >= -0.25) & (boundary <= 0.75))
dat.figure6b.text <- data.frame(x = c(0.45, 0.43, 0.39, 0.29),
                                y = c(0.30, 0.35, 0.42, 0.71),
                                label = c("0.70", "0.75", "0.80", "0.85"))

Trans6b <- function(x)
{
  y <- 1.1 - log(114 - x)
  return(y)
}
InvTrans6b <- function(y)
{
  x <- 114 - exp(1.1 - y)
  return(x)
}
trans6b_trans <- function() trans_new(name = "trans6b", transform = "Trans6b", inverse = "InvTrans6b")

hist(Trans6b(dat.figure6b$u), breaks = 50)

figure6b <- ggplot(data = dat.figure6b, mapping = aes(x = boundary, y = t)) +
  geom_raster(mapping = aes(fill = u), interpolate = TRUE) +
  scale_fill_gradientn(colours = c("#400000", "#800000", "#FF0000", "#FF4000", "#FF8000", "#FFFF00", "#FFFFFF"),
                       trans = "trans6b",
                       breaks = c(10, 50, 60, 70, 80, 90, 100, 110, 111, 112, 113),
                       guide = guide_colorbar(barheight = 25, title = "Utility\n ($M)")) +
  geom_hline(yintercept = t.opt, linetype = 2, size = 1) +
  geom_vline(xintercept = boundary.opt, linetype = 2, size = 1) +
  geom_contour(mapping = aes(z = power), breaks = seq(0.70, 0.85, 0.05), linetype = 2, size = 1) +
  geom_text(data = dat.figure1.text, mapping = aes(x = x, y = y, label = label), colour = "blue") +
  labs(x = "B-value", y = "Information Fraction") +
  theme_minimal()