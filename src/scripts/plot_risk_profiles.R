# Purpose: Visualise cluster-specific trajectories (mean onset/diagnosis errors and onset distributions by condition) 
# from the delay-aware posterior. 

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(stringr)
library(readr)
library(forcats)
library(grid)



target_cluster <- 1

data_rds_path      <- "data/generated_promote_style_mixed_delays.rds"
posterior_delay_rds<- "src/resultsmixeddata/posterior_val_delay_train.rds"
directory            <- "src/clustertrajectoryplots/mixeddata"

# Load data and posterior


data_all        <- readRDS(data_rds_path)
posterior_train <- readRDS(posterior_delay_rds)

pp <- posterior_train$posterior.parameters

expected_t_full <- pp$expected_t      
gap_mu_full     <- pp$gap_mu_star     
C_star_full     <- pp$C_star         

# Recreate the same train/test split (seed = 42, 80/20)
set.seed(42)
n_total   <- nrow(data_all$d)
train_prop <- 0.8
train_idx <- sample(1:n_total, size = floor(train_prop * n_total))
test_idx  <- setdiff(1:n_total, train_idx)

N_train <- length(train_idx)
M       <- ncol(data_all$d)


# Cluster assignment (use which.min as you specified)

cluster_assign_train <- apply(C_star_full, 1, which.min)
idx_in_cluster_train <- which(cluster_assign_train == target_cluster)
idx_global           <- train_idx[idx_in_cluster_train]

# Pull TRUE data for those patients 
onset_true  <- data_all$onset[idx_global, , drop = FALSE]  # True onset (if available)
diag_true   <- data_all$t[idx_global,     , drop = FALSE]  # True diagnosis ages
d_true      <- data_all$d[idx_global,     , drop = FALSE]  # Presence indicator

# Predicted (from delay-aware posterior)
pred_onset   <- expected_t_full[idx_in_cluster_train, , drop = FALSE]
pred_diag    <- pred_onset + gap_mu_full[idx_in_cluster_train, , drop = FALSE]

# Build df with Columns: 
# Condition (index), Age, Event_Type (Onset/Diagnosis), Type (True/Predicted)
to_long <- function(mat, event_type, type_label) {
  if (is.null(mat)) return(tibble())
  as_tibble(mat) |>
    mutate(Patient = row_number()) |>
    pivot_longer(cols = -Patient, names_to = "ConditionIdx", values_to = "Age") |>
    mutate(
      Condition = as.integer(gsub("^V", "", ConditionIdx)),  # from default col names
      Event_Type = event_type,
      Type = type_label
    ) |>
    select(Condition, Age, Event_Type, Type) |>
    filter(!is.na(Age) & is.finite(Age))
}

# True onset: keep only finite values (if present)
df_true_onset <- to_long(onset_true, "Onset", "True")

# True diagnosis: only where event actually occurred (d_true == 1)
diag_mask <- d_true == 1 & is.finite(diag_true)
diag_true_masked <- diag_true
diag_true_masked[!diag_mask] <- NA_real_
df_true_diag <- to_long(diag_true_masked, "Diagnosis", "True")

# Pred onset/diag
df_pred_onset <- to_long(pred_onset, "Onset", "Predicted")
df_pred_diag  <- to_long(pred_diag,  "Diagnosis", "Predicted")

df_all <- bind_rows(df_true_onset, df_true_diag, df_pred_onset, df_pred_diag)

# Top 5 conditions by frequency (true diagnosis count) within the cluster
condition_counts <- colSums(d_true == 1, na.rm = TRUE)

cond_list <- data_all$cond_list


top_conds <- order(condition_counts, decreasing = TRUE)
top_conds <- top_conds[!is.na(top_conds) & condition_counts[top_conds] > 0]
top_conds <- top_conds[top_conds <= M]
top_conds <- head(top_conds, 5)



# Filter df_all to the selected conditions
df_plot <- df_all |>
  filter(Condition %in% top_conds) |>
  mutate(
    ConditionName = factor(cond_list[Condition], levels = cond_list[top_conds])
  )

# Shorten labels to <= 20 characters for display
short_labels <- function(x) {
  ifelse(nchar(x) > 20, paste0(substr(x, 1, 17), "..."), x)
}
levels(df_plot$ConditionName) <- short_labels(levels(df_plot$ConditionName))

AGE_MIN <- 0
AGE_MAX <- 110
ABS_ERR_CAP <- 30
TRIM_FRAC <- 0.10

tmean <- function(x, trim = TRIM_FRAC) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x, trim = trim)
}

M <- ncol(d_true)

diag_obs_mask  <- (d_true == 1) & is.finite(diag_true) & is.finite(pred_diag)
onset_obs_mask <- is.finite(onset_true) & is.finite(pred_onset)
valid_age <- function(A) (A >= AGE_MIN & A <= AGE_MAX)

diag_obs_mask  <- diag_obs_mask & valid_age(diag_true)  & valid_age(pred_diag)
onset_obs_mask <- onset_obs_mask & valid_age(onset_true) & valid_age(pred_onset)

diag_err  <- pred_diag - diag_true
onset_err <- pred_onset - onset_true
diag_err[!diag_obs_mask] <- NA_real_
onset_err[!onset_obs_mask] <- NA_real_

diag_err[is.finite(diag_err) & abs(diag_err) > ABS_ERR_CAP] <- NA_real_
onset_err[is.finite(onset_err) & abs(onset_err) > ABS_ERR_CAP] <- NA_real_

mean_diag_err  <- apply(diag_err,  2, tmean)
mean_onset_err <- apply(onset_err, 2, tmean)

diag_true_clean <- diag_true
diag_true_clean[!diag_obs_mask] <- NA_real_
diag_true_clean[is.finite(diag_true_clean) & !valid_age(diag_true_clean)] <- NA_real_
mean_diag_true <- apply(diag_true_clean, 2, tmean, trim = TRIM_FRAC)

condition_counts <- colSums(d_true == 1, na.rm = TRUE)
top_conds <- head(order(condition_counts, decreasing = TRUE)[condition_counts > 0], 10)

if (length(top_conds) > 0) {
  df_plot <- tibble(
    ConditionIdx    = top_conds,
    ConditionName   = cond_list[top_conds],
    Diagnosis_Error = mean_diag_err[top_conds],
    Onset_Error     = mean_onset_err[top_conds],
    Start_Age       = mean_diag_true[top_conds]
  ) %>%
    filter(if_any(c(Start_Age, Diagnosis_Error, Onset_Error), is.finite))

  df_plot$ConditionName <- factor(df_plot$ConditionName, levels = rev(df_plot$ConditionName))
  levels(df_plot$ConditionName) <- ifelse(
    nchar(levels(df_plot$ConditionName)) > 25,
    paste0(substr(levels(df_plot$ConditionName), 1, 22), "..."),
    levels(df_plot$ConditionName)
  )

  light_pink <- "#eb6fc2ff"
  light_blue <- "#7cc1e6ff"

  arrow_plot <- ggplot(df_plot, aes(y = ConditionName)) +
    geom_segment(aes(x = Start_Age, xend = Start_Age + Diagnosis_Error, yend = ConditionName),
                 arrow = arrow(length = unit(0.2, "cm")), color = light_pink, linewidth = 1.1, na.rm = TRUE) +
    geom_segment(aes(x = Start_Age, xend = Start_Age + Onset_Error,
                     y = as.numeric(ConditionName) - 0.25,
                     yend = as.numeric(ConditionName) - 0.25),
                 arrow = arrow(length = unit(0.2, "cm")), color = light_blue, linewidth = 1.1, na.rm = TRUE) +
    labs(title = paste("Mean Patient-Level Errors by Condition – Cluster", target_cluster),
         x = "Age (years)", y = "Condition") +
    theme_classic(base_size = 14) +
    theme(panel.grid.major.x = element_line(color = "grey90", linetype = "dashed")) +
    coord_cartesian(xlim = c(AGE_MIN, AGE_MAX))

  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(directory, sprintf("cluster_%d_arrow_error_plot.png", target_cluster)),
         plot = arrow_plot, width = 10, height = 6, dpi = 300)
  print(arrow_plot)
} else {
  warning("No conditions with positive diagnosis count in this cluster. Skipping arrow plot.")
}

# Onset: per-condition coloured histogram + density lines (overlay)
if (length(top_conds) > 0) {
  df_onset_true_sel <- df_true_onset %>%
    filter(Condition %in% top_conds) %>%
    mutate(ConditionName = factor(cond_list[Condition], levels = cond_list[top_conds]))

  df_onset_pred_sel <- df_pred_onset %>%
    filter(Condition %in% top_conds) %>%
    mutate(ConditionName = factor(cond_list[Condition], levels = cond_list[top_conds]))

  # shorten labels 
  levels(df_onset_true_sel$ConditionName) <- short_labels(levels(df_onset_true_sel$ConditionName))
  levels(df_onset_pred_sel$ConditionName) <- short_labels(levels(df_onset_pred_sel$ConditionName))

  # palette shared by fill & color so legend merges
  n_conds <- length(top_conds)
  pal <- scales::hue_pal()(n_conds)
  names(pal) <- levels(df_onset_true_sel$ConditionName)

  BINWIDTH <- 5

  p_onset_by_cond <- ggplot() +
    # Overlayed histograms, one per condition
    geom_histogram(
      data = df_onset_true_sel,
      aes(x = Age,
          y = after_stat(density),  # density so it matches the lines
          fill = ConditionName,
          group = ConditionName),          
      binwidth = BINWIDTH, boundary = 0, closed = "left",
      position = "identity", alpha = 0.35, color = NA
    ) +
    # Density lines per condition (predicted onsets)
    geom_density(
      data = df_onset_pred_sel,
      aes(x = Age, color = ConditionName),
      linewidth = 1.1, na.rm = TRUE
    ) +
    scale_fill_manual(values = pal, name = "Condition") +
    scale_color_manual(values = pal, name = "Condition") +
    labs(
      title = paste0("Cluster ", target_cluster),
      x = "Age of Onset", y = "Density"
    ) +
    coord_cartesian(xlim = c(0, 100)) +
    theme_classic(base_size = 14) +
    theme(legend.position = "right")

  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(directory, sprintf("cluster_%d_onset_hist_lines_by_condition.png", target_cluster)),
         p_onset_by_cond, width = 10, height = 6, dpi = 300)
} 