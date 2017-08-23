source("wrapper.R")
require(ggplot2)
require(scales)

###############################################################################
# set parameters
N <- 450
n.grid <- seq(45, 360, 3) # t from 0.1 to 0.8
boundary.grid <- seq(-0.5, 1.5, 0.02)
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

# set.seed(12138)
# n.sim <- 1000
# d.prior <- rnorm(n.sim, 0.2, 1 / sqrt(50))
# out <- UtilityGridSearch(d.prior, n.grid, boundary.grid)
# file <- paste("results/grid_d_prior", ".RData", sep = "")
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

dat.figure1 <- subset(dat.figure1, subset = (boundary <= 1) & (t >= 0.2))
dat.figure1.text <- data.frame(x = c(0.90, 0.85, 0.76, 0.42),
                               y = c(0.47, 0.51, 0.57, 0.75),
                               label = c("0.70", "0.75", "0.80", "0.85"))

Trans1 <- function(x)
{
  y <- 1.2 - log(114 - x)
  return(y)
}
InvTrans1 <- function(y)
{
  x <- 114 - exp(1.2 - y)
  return(x)
}
trans1_trans <- function() trans_new(name = "trans1", transform = "Trans1", inverse = "InvTrans1")

hist(Trans1(dat.figure1$u), breaks = 50)

figure1 <- ggplot(data = dat.figure1, mapping = aes(x = boundary, y = t)) +
  geom_raster(mapping = aes(fill = u)) +
  scale_fill_gradientn(colours = c("#400000", "#800000", "#FF0000", "#FF4000", "#FF8000", "#FFFF00", "#FFFFFF"),
                       trans = "trans1",
                       breaks = c(50, 60, 70, 80, 90, 100, 110, 111, 112, 113),
                       guide = guide_colorbar(barheight = 30, title = "Utility\n ($M)")) +
  geom_hline(yintercept = t.opt) +
  geom_vline(xintercept = boundary.opt) +
  geom_contour(mapping = aes(z = power), breaks = seq(0.70, 0.85, 0.05), linetype = 2, size = 1.2) +
  geom_text(data = dat.figure1.text, mapping = aes(x = x, y = y, label = label), colour = "blue") +
  labs(x = "B-value", y = "Information Fraction") +
  theme_minimal()

###############################################################################
# Figure 2

file <- "results/grid_d_prior.RData"
load(file)
ind <- apply(out$u, 1, which.max)
ind <- (ind - 1) * length(t.grid) + seq_along(t.grid)

# axis2 = -0.025 + 0.05 * axis1
# axis1 = 0.5 + 20 * axis2

file <- "results/grid_d_0.RData"
load(file)

dat.figure2.alpha <- data.frame(t = t.grid, p = 0.5 + 20 * out$power[ind], type = "alpha")

file <- "results/grid_d_0.2.RData"
load(file)

dat.figure2.power <- data.frame(t = t.grid, p = out$power[ind], type = "power")

dat.figure2 <- rbind(dat.figure2.alpha, dat.figure2.power)

figure2 <- ggplot(data = dat.figure2, mapping = aes(x = t, y = p, group = type, color = type)) +
  geom_vline(xintercept = t.opt, linetype = 2, size = 1) +
  geom_line(size = 1.2) +
  scale_colour_manual(name = NULL,values =c("blue", "red"),
                      labels = c("Overall power at optimal futility boundary", "Type I error at optimal futility boundary")) +
  geom_hline(yintercept = 0.85, linetype = 2, size = 1, colour = "red") +
  geom_hline(yintercept = 1, linetype = 2, size = 1, colour = "blue") +
  coord_cartesian(ylim = c(0.50, 1.00)) +
  scale_y_continuous(sec.axis = sec_axis(~ . * 0.05 - 0.025, name = "Type I Error Rate")) +
  labs(x = "Information Fraction", y = "Overall Power") +
  theme_minimal() %+replace%
  theme(legend.position = c(0.8, 0.2))

###############################################################################
# Figure 3