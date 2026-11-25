###############################################################################
# Sockeye salmon forward simulation: Ricker recruitment + simple harvest
#
# Organization and naming loosely follow fwd_sim.R from the Yukon Chinook
# stock assessment (Pacific-salmon-assess): setup -> inputs -> simulation
# settings -> main loop -> outputs (S_fwd_*.csv, H_fwd_*.csv style files).
# Unlike that script, everything happens sequentially here (no process()
# helper, no parallelism, no HCR comparison, no AR1/TVA productivity, no
# forecast error / outcome uncertainty yet).
#
# Dynamics
#   Recruits by brood year:  R[t] = S[t] * exp(alpha - beta * S[t] + eps[t])
#                            eps[t] ~ Normal(-sigma^2/2, sigma)
#   Run in return year y:    N[y] = sum_a R[y - a] * p[a]
#   Harvest:                 H[y] = ER * N[y]
#   Spawners:                S[y] = N[y] - H[y]
###############################################################################

library(here)
library(tidyverse)

set.seed(2)

# 1. Scenario and file paths ---------------------------------------------

k          <- "base"                                 # scenario label used in output names
draws.file <- here("output", "posterior_draws.csv")  # alpha, beta, sigma
out.dir    <- here("output")

# 2. Simulation settings ---------------------------------------------------

num.sims     <- 1000      # number of simulation replicates (one posterior draw each)
ny           <- 34        # number of projection years
last.obs.yr  <- 2025      # last observed year; projections start the year after
ER           <- 0.20      # fixed exploitation rate
bias.correct <- FALSE      # mean-one lognormal process error on the recruit scale

# Age-at-return proportions - PLACEHOLDERS
p        <- c("4" = 0.3, "5" = 0.6, "6" = 0.1)
ages     <- as.integer(names(p))
max.age  <- max(ages)
stopifnot(abs(sum(p) - 1) < 1e-8)

# Initial spawners as a fraction of each draw's equilibrium (alpha / beta)
init.frac <- 0.5

# 3. Read posterior draws -----------------------------------------------------
# Needs columns alpha, beta, sigma - may need to rename

# Use fake parameter draws for now:
n.mock <- 4000
draws <- data.frame(
  alpha = rnorm(n.mock, 1.5, 0.15),
  beta  = rlnorm(n.mock, log(1 / 200000), 0.15),   # Smax ~ 200,000
  sigma = rlnorm(n.mock, log(0.6), 0.1))


# Take one joint draw for each simulation from alpha, beta, and sigma
pars <- draws[sample(nrow(draws), num.sims, replace = F), ]

# 4. Set up arrays and initial conditions ---------------------------------------------
# Brood years run from (1 - max.age) to ny so the first returns are seeded.

brood.yrs <- (1 - max.age):ny
R.fwd <- matrix(NA, num.sims, length(brood.yrs), dimnames = list(NULL, brood.yrs))
N.fwd <- matrix(NA, num.sims, ny)    # run (total return)
S.fwd <- matrix(NA, num.sims, ny)    # spawners
H.fwd <- matrix(NA, num.sims, ny)    # harvest (catch)

# 5. Main simulation loop --------------------------------------------------------

for (i in 1:num.sims) { # Iteration loop

  alpha <- pars$alpha[i]
  beta  <- pars$beta[i]
  sigma <- pars$sigma[i]

  # Seed brood years before the projection with equilibrium-fraction spawners
  S.init <- init.frac * alpha / beta
  R.fwd[i, as.character((1 - max.age):0)] <- S.init * exp(alpha - beta * S.init)

  for (y in 1:ny) { # Year loop

    # Run returning this year, from recruits of brood years y - age
    N.fwd[i, y] <- sum(R.fwd[i, as.character(y - ages)] * p)

    # Harvest and escapement
    H.fwd[i, y] <- ER * N.fwd[i, y]
    S.fwd[i, y] <- N.fwd[i, y] - H.fwd[i, y]

    # Process error, epsilon with or without bias correction (should match what's done in s-r models)
    eps <- rnorm(1, mean = if (bias.correct) (-sigma^2)/2 else 0, sd = sigma)
    
    # Recruits of brood year y produced by this year's spawners
    R.fwd[i, as.character(y)] <- S.fwd[i, y] * exp(alpha - beta * S.fwd[i, y] + eps) # Need to add gamma and kappa
  }
}

# 6. Outputs -----------------------------------------------------------------------
# One row per simulation x year - by brood year 

yrs <- last.obs.yr + 1:ny

sim.out <- tibble(
  sim      = rep(1:num.sims, times = ny),
  year     = rep(yrs, each = num.sims),
  recruits = as.vector(R.fwd[, as.character(1:ny)]),
  run      = as.vector(N.fwd),
  spawners = as.vector(S.fwd),
  catch    = as.vector(H.fwd)
)

# Summaries by year: median and quartiles
summarize_fwd <- function(df, var) {
  df %>%
    group_by(year) %>%
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

# Simulation parameter quantiles
pars.out <- pars %>%
  summarise(across(everything(), list(
    q10 = ~ quantile(.x, 0.10), med = ~ median(.x), q90 = ~ quantile(.x, 0.90)
  )))
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
  facet_wrap(~ metric, scales = "free_y") +
  labs(x = "Year", y = "Fish (median, 50% interval)") +
  theme_bw()
