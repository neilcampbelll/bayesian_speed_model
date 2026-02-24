rm(list = ls())

Sys.setenv(TMPDIR = "C:/Users/Neil.Campbell/stan_tmp")
dir.create("C:/Users/Neil.Campbell/stan_tmp", showWarnings = FALSE)


install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))

library(brms)
library(dplyr)
library(tidyr)
library(tidybayes)
library(ggplot2)
library(purrr)
library(patchwork)
library(cmdstanr)

set_cmdstan_path("C:/Users/Neil.Campbell/.cmdstan/cmdstan-2.38.0") ## this will need to be set to your own values

# =============================================================================
# 1. Load and prepare data
# =============================================================================

load("Data/tacsatEflalo2024.RData")

vms_df <- tacsatEflalo %>%
  select(SI_SP, LE_GEAR, LE_L5MET) %>%
  mutate(
    LE_GEAR = as.factor(LE_GEAR),
    LE_L5MET = as.factor(LE_L5MET)
  )

# Check nesting: each L5MET should belong to only one GEAR
nesting_check <- vms_df %>%
  distinct(LE_GEAR, LE_L5MET) %>%
  group_by(LE_L5MET) %>%
  filter(n() > 1)

if (nrow(nesting_check) > 0) {
  warning("Some LE_L5MET values appear in multiple LE_GEAR groups. Fixing via L5 prefix.")
  print(nesting_check)
  vms_df$LE_GEAR <- as.factor(
    sapply(strsplit(as.character(vms_df$LE_L5MET), split = "_"), `[`, 1)
  )
}

## clean mismatching gear types at l4 and l5
vms_df$LE_L5MET[vms_df$LE_L5MET == "OTM_SPF" & vms_df$LE_GEAR == "PS"] <- "PS_SPF"
vms_df$LE_L5MET[vms_df$LE_L5MET == "OTB_DEF" & vms_df$LE_GEAR == "PTM"] <- "PTM_DEF"
vms_df$LE_L5MET[vms_df$LE_L5MET == "SSC_DEF" & vms_df$LE_GEAR == "GTR"] <- "GTR_DEF"
vms_df$LE_L5MET[vms_df$LE_L5MET == "PTM_DEF" & vms_df$LE_GEAR == "OTB"] <- "OTB_DEF"




# Add small constant to avoid zero speeds (gamma requires > 0)
vms_df$SI_SP <- vms_df$SI_SP + 0.001

# ---- Optional: subsample for development/testing ----
# Uncomment the following to use a random subset while tuning the model.
# Once satisfied, comment it out and run on full data.
# set.seed(42)
# vms_df <- vms_df %>% slice_sample(prop = 0.2)

cat("Data dimensions:", nrow(vms_df), "rows,",
    nlevels(vms_df$LE_GEAR), "gear groups,",
    nlevels(vms_df$LE_L5MET), "metiers\n")

# =============================================================================
# 2. Model specification
# =============================================================================
# Key simplifications for speed:
#   - Random effects on mu only (not shape) — shape is hard to estimate per group
#   - Single grouping level (LE_L5MET) rather than nested (LE_GEAR/LE_L5MET)
#     since LE_L5MET already identifies the metier uniquely
#   - K = 4 components retained; reduce to 3 if appropriate for your data

K <- 4

mix_gamma_k <- brms::mixture("gamma", nmix = K, order = TRUE)

bform <- bf(
  SI_SP ~ 1,
  family = mix_gamma_k,
  # Random effects on component means only
  mu1 ~ 1 + (1 | LE_L5MET),
  mu2 ~ 1 + (1 | LE_L5MET),
  mu3 ~ 1 + (1 | LE_L5MET),
  mu4 ~ 1 + (1 | LE_L5MET),
  # Fixed (intercept-only) shape parameters — no random effects
  shape1 ~ 1,
  shape2 ~ 1,
  shape3 ~ 1,
  shape4 ~ 1
)

# =============================================================================
# 3. Priors
# =============================================================================

priors <- c(
  # Global intercepts for component means (on log scale for gamma)
  prior(normal(log(1), 1.5),  class = Intercept, dpar = mu1),
  prior(normal(log(3), 1.5),  class = Intercept, dpar = mu2),
  prior(normal(log(6), 1.5),  class = Intercept, dpar = mu3),
  prior(normal(log(10), 1.5), class = Intercept, dpar = mu4),
  
  # Shape parameter intercepts
  prior(normal(log(2), 1), class = Intercept, dpar = shape1),
  prior(normal(log(2), 1), class = Intercept, dpar = shape2),
  prior(normal(log(2), 1), class = Intercept, dpar = shape3),
  prior(normal(log(2), 1), class = Intercept, dpar = shape4),
  
  # SD of random effects for mu parameters at LE_L5MET level
  prior(exponential(1), class = sd, group = LE_L5MET, dpar = mu1),
  prior(exponential(1), class = sd, group = LE_L5MET, dpar = mu2),
  prior(exponential(1), class = sd, group = LE_L5MET, dpar = mu3),
  prior(exponential(1), class = sd, group = LE_L5MET, dpar = mu4)
)

# =============================================================================
# 4. Fit the model
# =============================================================================


fit <- brm(
  formula = bform,
  data    = vms_df,
  prior   = priors,
  iter    = 2000,
  warmup  = 1000,
  chains  = 4,
  cores   = 4,
  backend = "cmdstanr",
  threads = threading(2),  # within-chain parallelism (needs cmdstanr)
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  seed    = 123,
  refresh = 100,
  file    = "brms_vms_mix_k4_gamma_simplified"
)

# =============================================================================
# 5. Diagnostics
# =============================================================================

summary(fit)
plot(fit)
pp_check(fit, ndraws = 20)

# =============================================================================
# 6. Extract component parameters per metier
# =============================================================================

# Get unique metiers for prediction
new_data <- distinct(vms_df, LE_L5MET, LE_GEAR)

# Extract posterior draws for each parameter manually
# This is more reliable than posterior_epred for mixture models

draws <- as_draws_df(fit)

# Helper: get the random effect offset for a given metier and parameter
get_metier_params <- function(draws_df, new_data_df) {
  
  # Global intercepts (on link/log scale)
  globals <- tibble(
    k = 1:K,
    b_mu    = map_dbl(1:K, ~ median(draws_df[[paste0("b_mu", .x, "_Intercept")]])),
    b_shape = map_dbl(1:K, ~ median(draws_df[[paste0("b_shape", .x, "_Intercept")]]))
  )
  
  # For each metier, get the random effect offset for each mu_k
  metier_list <- new_data_df$LE_L5MET
  
  results <- map_dfr(seq_along(metier_list), function(i) {
    met <- as.character(metier_list[i])
    gear <- as.character(new_data_df$LE_GEAR[i])
    
    map_dfr(1:K, function(k) {
      # Random effect column name pattern: r_LE_L5MET__mu{k}[{metier},Intercept]
      re_col <- paste0("r_LE_L5MET__mu", k, "[", met, ",Intercept]")
      
      if (re_col %in% colnames(draws_df)) {
        re_draws <- draws_df[[re_col]]
        mu_draws <- exp(draws_df[[paste0("b_mu", k, "_Intercept")]] + re_draws)
      } else {
        mu_draws <- exp(draws_df[[paste0("b_mu", k, "_Intercept")]])
      }
      
      shape_draws <- exp(draws_df[[paste0("b_shape", k, "_Intercept")]])
      
      # Mixing weights (theta) — brms uses theta1..theta(K-1) on log-odds scale
      # For ordered mixtures, the thetas are on the log scale
      # Extract from draws directly
      tibble(
        LE_L5MET = met,
        LE_GEAR  = gear,
        k        = k,
        mu       = median(mu_draws),
        mu_lower = quantile(mu_draws, 0.1),
        mu_upper = quantile(mu_draws, 0.9),
        shape    = median(shape_draws)
      )
    })
  })
  
  results
}

metier_params <- get_metier_params(draws, new_data)
print(metier_params)

# =============================================================================
# 7. Find intersection thresholds between adjacent components
# =============================================================================

find_intersection <- function(mu1, shape1, lambda1, mu2, shape2, lambda2,
                              speed_range = c(0.01, 25), tol = 0.001) {
  # Find speed where weighted density of component 1 = weighted density of component 2
  f <- function(x) {
    d1 <- lambda1 * dgamma(x, shape = shape1, rate = shape1 / mu1)
    d2 <- lambda2 * dgamma(x, shape = shape2, rate = shape2 / mu2)
    d1 - d2
  }
  
  # Check if there's a sign change
  x_seq <- seq(speed_range[1], speed_range[2], length.out = 500)
  f_vals <- sapply(x_seq, f)
  sign_changes <- which(diff(sign(f_vals)) != 0)
  
  if (length(sign_changes) == 0) return(NA_real_)
  
  # Use the last sign change (between the two component peaks)
  idx <- sign_changes[length(sign_changes)]
  result <- tryCatch(
    uniroot(f, interval = c(x_seq[idx], x_seq[idx + 1]), tol = tol)$root,
    error = function(e) NA_real_
  )
  result
}

# For mixing weights, we need to extract theta from the model
# brms stores theta1, theta2, theta3 for K=4 (K-1 parameters)
# These are on the log-ratio scale for ordered mixtures

get_mixing_weights <- function(draws_df, K) {
  # For ordered mixtures, brms uses a stick-breaking parameterisation
  # Extract theta columns
  theta_cols <- paste0("theta", 1:(K - 1))
  available <- theta_cols[theta_cols %in% colnames(draws_df)]
  
  if (length(available) == 0) {
    # Equal weights as fallback
    warning("Could not find theta parameters. Using equal weights.")
    return(rep(1 / K, K))
  }
  
  # Median of each theta draw, then convert to simplex
  # For ordered mixture, thetas are already mixture probabilities
  thetas <- map_dbl(available, ~ median(draws_df[[.x]]))
  
  # The Kth weight is 1 - sum of others
  lambdas <- c(thetas, 1 - sum(thetas))
  lambdas <- pmax(lambdas, 0)  # ensure non-negative
  lambdas / sum(lambdas)        # normalise
}

lambdas <- get_mixing_weights(draws, K)
cat("Estimated mixing weights:", round(lambdas, 3), "\n")

# Calculate thresholds per metier
threshold_results <- map_dfr(seq_len(nrow(new_data)), function(i) {
  met <- as.character(new_data$LE_L5MET[i])
  gear <- as.character(new_data$LE_GEAR[i])
  
  params <- metier_params %>% filter(LE_L5MET == met) %>% arrange(k)
  
  if (nrow(params) < 2) return(tibble())
  
  map_dfr(1:(nrow(params) - 1), function(j) {
    threshold <- find_intersection(
      mu1 = params$mu[j],     shape1 = params$shape[j],     lambda1 = lambdas[j],
      mu2 = params$mu[j + 1], shape2 = params$shape[j + 1], lambda2 = lambdas[j + 1]
    )
    
    tibble(
      LE_GEAR   = gear,
      LE_L5MET  = met,
      threshold = paste0("T", j, "_", j + 1),
      speed     = threshold
    )
  })
})

print(threshold_results)

# Key thresholds (adjust these labels to match your interpretation)
# T1_2 = boundary between "not moving" and "slow/fishing"
# T3_4 = boundary between "fishing" and "steaming"
key_thresholds <- threshold_results %>%
  filter(threshold %in% c("T1_2", "T3_4"))

print(key_thresholds)

# =============================================================================
# 8. Plotting
# =============================================================================

speed_max <- quantile(vms_df$SI_SP, 0.99)
speed_grid <- seq(0.01, speed_max, length.out = 300)

# Generate density lines per metier
density_data <- map_dfr(seq_len(nrow(new_data)), function(i) {
  met <- as.character(new_data$LE_L5MET[i])
  gear <- as.character(new_data$LE_GEAR[i])
  params <- metier_params %>% filter(LE_L5MET == met) %>% arrange(k)
  
  map_dfr(1:nrow(params), function(j) {
    mu_j    <- params$mu[j]
    shape_j <- params$shape[j]
    lam_j   <- lambdas[j]
    
    tibble(
      LE_L5MET = met,
      LE_GEAR  = gear,
      k        = j,
      speed    = speed_grid,
      density  = lam_j * dgamma(speed_grid, shape = shape_j, rate = shape_j / mu_j)
    )
  })
})

# Plot per gear group
plot_list <- map(levels(vms_df$LE_GEAR), function(gear_group) {
  metiers_in_gear <- new_data %>%
    filter(LE_GEAR == gear_group) %>%
    pull(LE_L5MET) %>%
    unique()
  
  sub_plots <- map(metiers_in_gear, function(met) {
    met_data      <- filter(vms_df, LE_L5MET == met)
    met_density   <- filter(density_data, LE_L5MET == met)
    met_threshold <- filter(key_thresholds, LE_L5MET == met)
    
    ggplot(met_data, aes(x = SI_SP)) +
      geom_histogram(aes(y = after_stat(density)),
                     bins = 50, alpha = 0.4, fill = "grey") +
      geom_line(data = met_density,
                aes(x = speed, y = density, colour = factor(k)),
                linewidth = 1) +
      geom_vline(data = met_threshold,
                 aes(xintercept = speed, linetype = threshold),
                 colour = "black", linewidth = 1) +
      scale_colour_viridis_d(name = "Component") +
      scale_linetype_manual(
        name = "Threshold",
        values = c("T1_2" = "dashed", "T3_4" = "dotted")
      ) +
      ggtitle(paste0("Metier: ", met)) +
      coord_cartesian(xlim = c(0, speed_max)) +
      theme_bw() +
      theme(legend.position = "bottom")
  })
  
  wrap_plots(sub_plots) +
    plot_annotation(title = paste("Gear Group:", gear_group))
})

# Display or save
# plot_list[[1]]
walk2(plot_list, levels(vms_df$LE_GEAR), ~ {
  ggsave(paste0("Threshold_Plot_", .y, ".png"), .x, width = 14, height = 10)
})
