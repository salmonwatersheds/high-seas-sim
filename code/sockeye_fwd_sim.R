###############################################################################
# Sockeye salmon forward simulation: hierarchical Ricker (AR1 + covariates)
# plus simple harvest
#
# Organization and naming loosely follow fwd_sim.R from the Yukon Chinook
# stock assessment (Pacific-salmon-assess): setup -> inputs -> simulation
# settings -> main loop -> outputs (S_fwd_*.csv, H_fwd_*.csv style files).
# Unlike that script, everything happens sequentially here (no process()
# helper, no parallelism, no HCR comparison, no forecast error / outcome
# uncertainty yet).
#
# Parameterization matches hbm_stat_inter.stan (all on the log scale):
#   y[t]   = ln(R[t] / S[t])
#   x1[t]  = (S[t] - S.mean) / S.sd                       (standardized spawners)
#   y[t]   = alpha + beta*x1 + gamma*x2 + kappa*x3 + chi*x2*x3 + eps[t]
#   eps[t] = phi * eps[t-1] + e[t],      e[t] ~ Normal(0, sigma)
#            (sigma is the Stan "sigma": innovation SD, already AR1-corrected)
#
# Population dynamics (per stock)
#   Recruits by brood year:  R[t] = S[t] * exp(y[t])
#   Run in return year y:    N[y] = sum_a R[y - a] * p[a]
#   Harvest:                 H[y] = ER * N[y]
#   Spawners:                S[y] = N[y] - H[y]
###############################################################################

library(here)
library(tidyverse)

set.seed(2)

# 1. Scenario and file paths ---------------------------------------------

k         <- "base"                                    # scenario label used in output names
fit.file  <- here("output", "stan_fit.rds")            # stanfit object from hbm_stat_inter.stan
out.dir   <- here("output")
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

# 2. Simulation settings ---------------------------------------------------

num.sims     <- 1000      # number of simulation replicates (one posterior draw each)
ny           <- 34        # number of projection years
last.obs.yr  <- 2025      # last observed (brood) year; projections start the year after
ER           <- 0.20      # fixed exploitation rate

# Stocks (series), in the same order as the series index in the Stan model - PLACEHOLDERS
stock.names <- c("Stock A", "Stock B", "Stock C")
n.stocks    <- length(stock.names)

# Standardization used for x1 (spawners), one value per stock - PLACEHOLDERS
# Needed to convert projected spawners to x1 = (S - S.mean) / S.sd
S.mean <- c(100000, 60000, 30000)
S.sd   <- c(50000, 30000, 15000)

# Row index of the last observed year of each series in the Stan data (y_end)
# Used to pull the last fitted residual as the starting AR1 state - PLACEHOLDERS
y.end <- c(30, 60, 90)

# Age-at-return proportions by total age, one row per stock - PLACEHOLDERS
ages    <- c(4, 5, 6)
max.age <- max(ages)
p <- matrix(rep(c(0.3, 0.6, 0.1), each = n.stocks), nrow = n.stocks,
            dimnames = list(stock.names, ages))
stopifnot(abs(rowSums(p) - 1) < 1e-8)

# Future covariates: flat at 0 (= historical mean, since they were standardized).
# Replace with a vector of length ny to project a different path.
x2.proj <- rep(0, ny)
x3.proj <- rep(0, ny)

# Initial spawners as a fraction of each draw's equilibrium (where ln(R/S) = 0)
init.frac <- 0.5

# 3. Read posterior draws -----------------------------------------------------
# Expected list `post` as returned by rstan::extract(), with draws in rows:
#   alpha, beta, gamma, kappa, chi, sigma : [draws x n_series]
#   phi                                   : [draws]
#   epsilon                               : [draws x N]  (residuals, indexed by data row)

if (file.exists(fit.file)) {
  fit  <- readRDS(fit.file)
  post <- rstan::extract(fit, pars = c("alpha", "beta", "gamma", "kappa", "chi",
                                       "sigma", "phi", "epsilon"))
} else {
  message("'", fit.file, "' not found: using made-up placeholder draws.")
  n.mock <- 4000
  mock <- function(mu, sd) matrix(rnorm(n.mock * n.stocks, mu, sd), n.mock, n.stocks)
  post <- list(
    alpha   = mock(1.5, 0.15),
    beta    = -exp(mock(log(0.5), 0.15)),            # negative: density dependence
    gamma   = mock(0, 0.1),
    kappa   = mock(0, 0.1),
    chi     = mock(0, 0.05),
    sigma   = exp(mock(log(0.5), 0.1)),
    phi     = pmin(pmax(rnorm(n.mock, 0.4, 0.1), -0.95), 0.95),
    epsilon = matrix(rnorm(n.mock * max(y.end), 0, 0.5), n.mock, max(y.end))
  )
}
stopifnot(ncol(post$alpha) == n.stocks)

# Residual in the last observed year of each series = starting AR1 state
eps.last <- post$epsilon[, y.end, drop = FALSE]      # [draws x n_stocks]

# Take one joint draw for each simulation. The same draw is used for all stocks,
# which keeps the hierarchical structure and the shared phi intact.
n.draws  <- nrow(post$alpha)
draw.idx <- sample(n.draws, num.sims, replace = num.sims > n.draws)

if (any(post$beta[draw.idx, ] >= 0))
  warning("Some draws have beta >= 0 (no density dependence); equilibrium is set to 0 for those.")

# 4. Set up arrays and initial conditions ---------------------------------------------
# Brood years run from (1 - max.age) to ny so the first returns are seeded.

brood.yrs <- (1 - max.age):ny
R.fwd <- array(NA, c(num.sims, n.stocks, length(brood.yrs)),
               dimnames = list(NULL, stock.names, brood.yrs))   # recruits by brood year
N.fwd <- array(NA, c(num.sims, n.stocks, ny), dimnames = list(NULL, stock.names, 1:ny))  # run
S.fwd <- N.fwd    # spawners
H.fwd <- N.fwd    # harvest (catch)

# 5. Main simulation loop --------------------------------------------------------

for (i in 1:num.sims) { # Iteration loop
  
  d   <- draw.idx[i]
  phi <- post$phi[d]
  
  for (s in 1:n.stocks) { # Stock loop
    
    alpha <- post$alpha[d, s]
    beta  <- post$beta[d, s]
    gamma <- post$gamma[d, s]
    kappa <- post$kappa[d, s]
    chi   <- post$chi[d, s]
    sigma <- post$sigma[d, s]
    eps   <- eps.last[d, s]       # AR1 state: last fitted residual
    
    # Seed brood years before the projection with equilibrium-fraction spawners
    # (x2 = x3 = 0 and eps = 0 here, so equilibrium is where alpha + beta*x1 = 0)
    S.eq   <- max(S.mean[s] - S.sd[s] * alpha / beta, 0)
    S.init <- init.frac * S.eq
    x1     <- (S.init - S.mean[s]) / S.sd[s]
    R.fwd[i, s, as.character((1 - max.age):0)] <- S.init * exp(alpha + beta * x1)
    
    for (y in 1:ny) { # Year loop
      
      # Run returning this year, from recruits of brood years y - age
      N.fwd[i, s, y] <- sum(R.fwd[i, s, as.character(y - ages)] * p[s, ])
      
      # Harvest and escapement
      H.fwd[i, s, y] <- ER * N.fwd[i, s, y]
      S.fwd[i, s, y] <- N.fwd[i, s, y] - H.fwd[i, s, y]
      
      # Standardize spawners to match the x1 used in the model
      x1 <- (S.fwd[i, s, y] - S.mean[s]) / S.sd[s]
      
      # AR1 process error: no bias correction, matches the mean-zero normal in the Stan model
      eps <- phi * eps + rnorm(1, 0, sigma)
      
      # Recruits of brood year y produced by this year's spawners
      R.fwd[i, s, as.character(y)] <- S.fwd[i, s, y] *
        exp(alpha + beta * x1 + gamma * x2.proj[y] + kappa * x3.proj[y] +
              chi * x2.proj[y] * x3.proj[y] + eps)
    }
  }
}

# 6. Outputs -----------------------------------------------------------------------
# One row per simulation x stock x year - recruits are by brood year

yrs <- last.obs.yr + 1:ny

sim.out <- tibble(
  sim      = rep(1:num.sims, times = n.stocks * ny),
  stock    = rep(rep(stock.names, each = num.sims), times = ny),
  year     = rep(yrs, each = num.sims * n.stocks),
  recruits = as.vector(R.fwd[, , as.character(1:ny), drop = FALSE]),
  run      = as.vector(N.fwd),
  spawners = as.vector(S.fwd),
  catch    = as.vector(H.fwd)
)

# Summaries by stock and year: median and quartiles
summarize_fwd <- function(df, var) {
  df %>%
    group_by(stock, year) %>%
    summarise(
      med = median({{ var }}),
      lwr = quantile({{ var }}, 0.25),
      upr = quantile({{ var }}, 0.75),
      .groups = "drop"
    )
}

S.fwd.summary <- summarize_fwd(sim.out, spawners)
H.fwd.summary <- summarize_fwd(sim.out, catch)
R.fwd.summary <- summarize_fwd(sim.out, recruits)

write_csv(sim.out,       file.path(out.dir, paste0("sim_out_", k, ".csv")))
write_csv(S.fwd.summary, file.path(out.dir, paste0("S_fwd_", k, ".csv")))
write_csv(H.fwd.summary, file.path(out.dir, paste0("H_fwd_", k, ".csv")))
write_csv(R.fwd.summary, file.path(out.dir, paste0("R_fwd_", k, ".csv")))

# Simulation parameter quantiles (parameters of the draws actually used)
pars <- tibble(
  stock = rep(stock.names, each = num.sims),
  alpha = as.vector(post$alpha[draw.idx, , drop = FALSE]),
  beta  = as.vector(post$beta[draw.idx, , drop = FALSE]),
  gamma = as.vector(post$gamma[draw.idx, , drop = FALSE]),
  kappa = as.vector(post$kappa[draw.idx, , drop = FALSE]),
  chi   = as.vector(post$chi[draw.idx, , drop = FALSE]),
  sigma = as.vector(post$sigma[draw.idx, , drop = FALSE]),
  phi   = rep(post$phi[draw.idx], times = n.stocks)
)

pars.out <- pars %>%
  group_by(stock) %>%
  summarise(across(everything(), list(
    q10 = ~ quantile(.x, 0.10), med = ~ median(.x), q90 = ~ quantile(.x, 0.90)
  )), .groups = "drop")
write_csv(pars.out, file.path(out.dir, paste0("sim_pars_", k, ".csv")))

# Quick look plot
bind_rows(
  mutate(S.fwd.summary, metric = "Spawners"),
  mutate(R.fwd.summary, metric = "Recruits"),
  mutate(H.fwd.summary, metric = "Catch")
) %>%
  ggplot(aes(year, med)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
  geom_line() +
  facet_grid(metric ~ stock, scales = "free_y") +
  labs(x = "Year", y = "Fish (median, 50% interval)") +
  theme_bw()