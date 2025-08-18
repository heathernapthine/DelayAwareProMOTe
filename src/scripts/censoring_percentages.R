# Purpose: Summarise censoring patterns (left/right) per cluster in the synthetic onset-time dataset.

library(dplyr)
library(tidyr)
library(readr)

data <- readRDS("data/generated_onset_promote_style.rds")

z        <- data$z      
d        <- data$d       
t_mat    <- data$t       
rho      <- data$rho      
tau      <- data$tau      
iota     <- data$iota     
N        <- data$N
K        <- data$K


# Censoring per cluster

# Initialise results
cluster_summary <- data.frame(
  Cluster = 1:K,
  N_Patients = integer(K),
  Diagnosis_Events = integer(K),
  Left_Censored = integer(K),
  Right_Censored = integer(K),
  Left_Censored_Perc = numeric(K),
  Right_Censored_Perc = numeric(K)
)

for (k in 1:K) {
  idx <- which(z == k)
  n_k <- length(idx)

  d_k <- d[idx, , drop = FALSE]
  t_k <- t_mat[idx, , drop = FALSE]
  rho_k <- rho[idx]
  tau_k <- tau[idx]
  iota_k <- iota[idx]

  # Repeat rho/tau/iota to match (N_k x M)
  rho_rep <- matrix(rho_k, nrow = n_k, ncol = ncol(d_k))
  tau_rep <- matrix(tau_k, nrow = n_k, ncol = ncol(d_k))
  iota_rep <- matrix(iota_k, nrow = n_k, ncol = ncol(d_k))

  # Only where disease is present
  present_mask <- d_k == 1

  # Left-censored: t < rho
  left_mask <- (t_k < rho_rep) & present_mask

  # Right-censored: t > tau AND iota == 0 (not dead)
  right_mask <- (t_k > tau_rep) & (iota_rep == 0) & present_mask

  total_events <- sum(present_mask, na.rm = TRUE)
  left_censored <- sum(left_mask, na.rm = TRUE)
  right_censored <- sum(right_mask, na.rm = TRUE)

  cluster_summary[k, ] <- list(
    Cluster = k,
    N_Patients = n_k,
    Diagnosis_Events = total_events,
    Left_Censored = left_censored,
    Right_Censored = right_censored,
    Left_Censored_Perc = round(100 * left_censored / total_events, 2),
    Right_Censored_Perc = round(100 * right_censored / total_events, 2)
  )
}

print(cluster_summary)
