# Does the same as generate_mixed_delays.R but with an additional 
# flag for early vs late on the bimodal conditions.
set.seed(499)

library(extraDistr)
library(MASS)
library(truncnorm)

# Parameters.
N <- 200000  # number of individuals
K <- 10      # number of clusters
cond_list <- c(
  "Allergic and chronic rhinitis", "Asbestosis", "Asthma", "Atrial fibrillation",
  "Bronchiectasis", "Cardiomyopathy", "Chronic liver disease", "Chronic renal disease",
  "Coeliac disease", "Conduction disorders and other arrhythmias", "Coronary heart disease",
  "Cystic Fibrosis", "Diabetes NOS", "Diverticular disease of intestine", "Erectile dysfunction",
  "Fatty Liver", "Gastro-oesophageal reflux, gastritis and similar", "Gout", "Heart failure",
  "Heart valve disorders", "Hyperplasia of prostate", "Hypertension", "Hypo or hyperthyroidism",
  "Inflammatory arthritis and other inflammatory conditions", "Inflammatory bowel disease",
  "Irritable bowel syndrome", "Non-acute cystitis", "Osteoporosis and vertebral crush fractures",
  "Peptic ulcer disease", "Peripheral arterial disease", "Primary pulmonary hypertension",
  "Sleep apnoea", "Spinal stenosis", "Stroke", "Transient ischaemic attack", "Type 1 diabetes",
  "Type 2 diabetes", "Urinary Incontinence", "Addisons disease", "Alcohol Problems",
  "Anorexia and bulimia nervosa", "Anxiety disorders", "Autism and Asperger’s syndrome",
  "Benign neoplasm of brain and other parts of CNS", "Bipolar affective disorder and mania",
  "Cerebral Palsy", "Dementia", "Down’s syndrome", "Epilepsy", "Hearing loss",
  "HIV", "Immunodeficiencies", "Intellectual disability",
  "Iron and vitamin deficiency anaemia conditions", "Macular degeneration", "Migraine",
  "Meniere disease", "Motor neuron disease", "Multiple sclerosis", "Myasthenia gravis",
  "Non-melanoma skin malignancies", "Obsessive-compulsive disorder",
  "Other psychoactive substance misuse", "Parkinson’s disease", "Peripheral or autonomic neuropathy",
  "Post-traumatic stress disorder", "Postviral fatigue syndrome, neurasthenia and fibromyalgia",
  "Psoriasis", "Sarcoidosis", "Schizophrenia, schizotypal and delusional disorders",
  "Sickel-cell anaemia", "Solid organ malignancies", "Thalassaemia",
  "Tuberculosis", "Visual impairment and blindness"
)

M <- length(cond_list)

# Map condition names to indices.
idx_of <- function(names_vec) {
  idx <- match(names_vec, cond_list)
  if (anyNA(idx)) {
    warning("These condition names were not found in cond_list: ",
            paste(names_vec[is.na(idx)], collapse = ", "))
  }
  idx[!is.na(idx)]
}

# Pick non-overlapping groups from the pool.
pool <- cond_list
pick_group <- function(pool, n) {
  chosen <- sample(pool, size = n, replace = FALSE)
  pool <- setdiff(pool, chosen)
  list(group = chosen, pool = pool)
}

# Six groups of five conditions.
res <- pick_group(pool, 5);  grp_uniform_asymptomatic <- res$group; pool <- res$pool
res <- pick_group(pool, 5);  grp_early_narrow         <- res$group; pool <- res$pool
res <- pick_group(pool, 5);  grp_late_wide            <- res$group; pool <- res$pool
res <- pick_group(pool, 5);  grp_psoriasis_medium     <- res$group; pool <- res$pool
res <- pick_group(pool, 5);  grp_mix_two_modes        <- res$group; pool <- res$pool
res <- pick_group(pool, 5);  grp_near_zero            <- res$group; pool <- res$pool

# Delay family & params per condition.
delay_dist_cond <- rep("gaussian", M)   # "gaussian", "uniform", "mixture2"
delay_group     <- rep("gaussian_other", M)
delay_mu_cond   <- rep(5, M)
delay_sd_cond   <- rep(1.5, M)
unif_a <- rep(NA_real_, M); unif_b <- rep(NA_real_, M)
mix_w1 <- rep(NA_real_, M)
mix_mu1 <- rep(NA_real_, M); mix_sd1 <- rep(NA_real_, M)
mix_mu2 <- rep(NA_real_, M); mix_sd2 <- rep(NA_real_, M)

# Near zero gaussian.
ix <- idx_of(grp_near_zero)
delay_dist_cond[ix] <- "gaussian"; delay_mu_cond[ix] <- 0.05; delay_sd_cond[ix] <- 0.05
delay_group[ix]     <- "gaussian_near_zero"

# Uniform.
ix <- idx_of(grp_uniform_asymptomatic)
delay_dist_cond[ix] <- "uniform"; unif_a[ix] <- 0; unif_b[ix] <- 15
delay_mu_cond[ix] <- NA_real_; delay_sd_cond[ix] <- NA_real_
delay_group[ix]   <- "uniform"

# Early tight gaussian.
ix <- idx_of(grp_early_narrow)
delay_dist_cond[ix] <- "gaussian"; delay_mu_cond[ix] <- 3; delay_sd_cond[ix] <- 0.75
delay_group[ix]     <- "gaussian_early"

# Late wide gaussian.
ix <- idx_of(grp_late_wide)
delay_dist_cond[ix] <- "gaussian"; delay_mu_cond[ix] <- 8; delay_sd_cond[ix] <- 3
delay_group[ix]     <- "gaussian_late"

# Medium gaussian.
ix <- idx_of(grp_psoriasis_medium)
delay_dist_cond[ix] <- "gaussian"; delay_mu_cond[ix] <- 5; delay_sd_cond[ix] <- 1.5
delay_group[ix]     <- "gaussian_medium"

# Two-mode mixture.
ix <- idx_of(grp_mix_two_modes)
delay_dist_cond[ix] <- "mixture2"
mix_w1[ix]   <- 0.5
mix_mu1[ix]  <- 1.5
mix_sd1[ix]  <- 0.5
mix_mu2[ix]  <- 7.0
mix_sd2[ix]  <- 2.0
delay_mu_cond[ix] <- NA_real_; delay_sd_cond[ix] <- NA_real_
delay_group[ix]   <- "mixture2"

# Bimodal indicator vector (needed by analysis).
is_bimodal <- as.integer(delay_dist_cond == "mixture2")

# Clusters.
cluster_sizes <- runif(K, 0.03, 0.15); cluster_sizes <- cluster_sizes / sum(cluster_sizes)
z <- sample(1:K, N, replace = TRUE, prob = cluster_sizes)

# Presence independent of clusters.
pi_m  <- rbeta(M, 0.07, 0.49)
pi_mk <- matrix(pi_m, nrow = M, ncol = K)

d <- matrix(0L, N, M)
for (n in 1:N) d[n, ] <- rbinom(M, 1, prob = pi_mk[, z[n]])

# Onset parameters (NIG).
mu_mk <- matrix(NA_real_, M, K)
sigma2_mk <- matrix(NA_real_, M, K)
for (m in 1:M) {
  for (k in 1:K) {
    sigma2 <- rinvgamma(1, 5, 300)
    mu <- rnorm(1, mean = 50, sd = sqrt(sigma2 / 0.3))
    mu_mk[m, k]   <- mu
    sigma2_mk[m,k] <- sigma2
  }
}

# Generate onset, delay, diagnosis.
onset <- matrix(0, N, M)
delay <- matrix(0, N, M)
t     <- matrix(0, N, M)

# NEW: per-observation early/late labels (only for mixture2 conditions).
delay_component <- matrix(NA_character_, N, M)

for (n in 1:N) {
  zn <- z[n]
  for (m in 1:M) {
    if (d[n, m] == 1L) {
      mu    <- mu_mk[m, zn]
      sigma <- sqrt(sigma2_mk[m, zn])
      onset[n, m] <- max(rnorm(1, mean = mu, sd = sigma), 0)

      distm <- delay_dist_cond[m]
      if (distm == "gaussian") {
        delay[n, m] <- truncnorm::rtruncnorm(1, a = 0, mean = delay_mu_cond[m], sd = delay_sd_cond[m])
        # leave delay_component[n, m] as NA for unimodal
      } else if (distm == "uniform") {
        delay[n, m] <- runif(1, unif_a[m], unif_b[m])
        # leave label NA (unimodal-like)
      } else if (distm == "mixture2") {
        comp <- rbinom(1, 1, prob = mix_w1[m])  # 1 = early, 0 = late
        if (comp == 1) {
          delay[n, m] <- truncnorm::rtruncnorm(1, a = 0, mean = mix_mu1[m], sd = mix_sd1[m])
          delay_component[n, m] <- "early"
        } else {
          delay[n, m] <- truncnorm::rtruncnorm(1, a = 0, mean = mix_mu2[m], sd = mix_sd2[m])
          delay_component[n, m] <- "late"
        }
      } else {
        stop("Unknown delay_dist_cond: ", distm)
      }

      t[n, m] <- onset[n, m] + delay[n, m]
    }
  }
}

# Individual censoring.
rho  <- runif(N, 20, 60)
tau  <- rho + 30
iota <- rbinom(N, 1, 0.8)

# Delay prior frame for inference.
delay_prior_df <- data.frame(
  condition  = cond_list,
  delay_dist = delay_dist_cond,
  delay_group = factor(delay_group, levels = c(
    "gaussian_near_zero","gaussian_early","gaussian_late","gaussian_medium",
    "gaussian_other","uniform","mixture2"
  )),
  delay_mu = delay_mu_cond,
  delay_sd = delay_sd_cond,
  unif_a   = unif_a,
  unif_b   = unif_b,
  mix_w1   = mix_w1,
  mix_mu1  = mix_mu1,
  mix_sd1  = mix_sd1,
  mix_mu2  = mix_mu2,
  mix_sd2  = mix_sd2,
  stringsAsFactors = FALSE
)

prior_groups <- list(
  uniform_asymptomatic = grp_uniform_asymptomatic,
  early_narrow         = grp_early_narrow,
  late_wide            = grp_late_wide,
  psoriasis_medium     = grp_psoriasis_medium,
  mix_two_modes        = grp_mix_two_modes,
  near_zero            = grp_near_zero
)

# Package dataset.
data_list <- list(
  z = z, d = d, t = t, onset = onset, delay = delay,
  rho = rho, tau = tau, iota = iota,
  N = N, K = K, M = M, cond_list = cond_list,
  pi_mk = pi_mk, mu_mk = mu_mk, sigma2_mk = sigma2_mk,
  delay_dist_cond = delay_dist_cond,
  delay_mu_cond   = delay_mu_cond,
  delay_sd_cond   = delay_sd_cond,
  delay_prior_df  = delay_prior_df,
  prior_groups    = prior_groups,
  # NEW:
  is_bimodal = is_bimodal,                 # 1 where mixture2
  delay_component = delay_component        # "early"/"late" only for mixture2; NA otherwise
)

# Quick shapes.
cat("Shapes of key items in data_list:\n")
for (name in c("z","d","t","onset","delay","rho","tau","iota",
               "pi_mk","mu_mk","sigma2_mk","delay_prior_df","delay_component","is_bimodal")) {
  item <- data_list[[name]]
  item_shape <- if (is.matrix(item)) {
    paste0("(", paste(dim(item), collapse = ", "), ")")
  } else if (is.data.frame(item)) {
    paste0("data.frame: (", nrow(item), ", ", ncol(item), ")")
  } else {
    paste0("Length: ", length(item))
  }
  cat(sprintf("%-18s : %s\n", name, item_shape))
}

# Save RDS.
saveRDS(data_list, "data/generated_promote_style_mixed_delays.rds")
