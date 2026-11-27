## Load data and model fits ## ----------------------------------------------

# Point to salmon warming crowded oceans repository
wc.dir <- "C:/Users/hhunter/Documents/Salmon-warming-crowded-oceans"

# Load brood table and associated covariate values
sock.brood.covar <- read.csv(here(wc.dir, "data/sockeye/master_sockeye_brood_table_covar.csv"))

# Raw brood tables - For (some of the) age data
sock.raw <- read.csv(here(wc.dir, "data-downloaded/salmon_productivity_compilation2026-08-10.csv")) %>%
  filter(species == "Sockeye", stock.id %in% sock.brood.covar$Stock.ID)

# Load Bayesian model fits

# Stationary - with interaction term - pink+chum+sockeye competitor index #NOTE: wanting to replace with fits to custom pink+chum index
load(here(wc.dir, "output/models/stat/sockeye/stat_inter.RData")) # Name: stat_inter, use verbose=T to check name of loaded object

dim(rstan::extract(stat_inter, pars='alpha')$alpha)

# Eras / regime shift - without interaction term - pink+chum+sockeye competitor index #NOTE: wanting to replace with fits to custom pink+chum index
load(here(wc.dir, "output/models/dyn/sockeye/hbm_era_2c.RData")) # Name: era.2c, use verbose=T to check name of loaded object
dim(rstan::extract(era.2c, pars='alpha')$alpha)


# Load Age data

# We need not only age at ocean entry, but also total age proportions

# Get total age info
rx.cols <- grep("^r[[:digit:]]$", names(sock.raw), value = TRUE)
rx.x.cols <- grep("^r[[:digit:]]\\.[[:digit:]]", names(sock.raw), value = TRUE)

# Case 1: add up r.x columns
sock.raw %>% mutate(r.sum = rowSums(.[rx.cols], na.rm=T)) %>% 
  mutate(r.diff = recruits - r.sum) %>% View()

# case 2: add up rx.x columns

# age at ocean entry will go into competitor index -> which will ACCOUNT FOR time lag in the same way the original project does, 
# using a simple weighted average matched to brood year 
# Age at ocean entry in "ocean_x" column in sock.brood.covar