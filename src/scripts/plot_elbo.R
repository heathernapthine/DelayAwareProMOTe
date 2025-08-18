library(ggplot2)

# Thesis colors 
light_pink <- "#fbaee1"
light_blue <- "#A6CEE3"
purple     <- "#bd5bd3"

path_delay   <- "src/paramresults/posterior_val_delay_train.rds"
path_nodelay <- "src/paramresults/posterior_val_no_delay_train.rds"
obj_delay   <- readRDS(path_delay)
obj_nodelay <- readRDS(path_nodelay)

elbo_delay <- (as.numeric(obj_delay$elbo))[3:69]   
elbo_nodelay <- as.numeric(obj_nodelay$elbo)[3:65]
param_diff_delay <- (as.numeric(obj_delay$param_diffs))[3:69]   
param_diff_nodelay <- as.numeric(obj_nodelay$param_diffs)[3:65]

print(param_diff_nodelay)

df_list <- list()
if (length(elbo_delay) > 0) {
  df_list[[length(df_list) + 1]] <- data.frame(
    iter  = seq_along(elbo_delay),
    elbo  = elbo_delay,
    model = "Delay-aware",
    stringsAsFactors = FALSE
  )
}
if (length(elbo_nodelay) > 0) {
  df_list[[length(df_list) + 1]] <- data.frame(
    iter  = seq_along(elbo_nodelay),
    elbo  = elbo_nodelay,
    model = "No-delay",
    stringsAsFactors = FALSE
  )
}

df <- if (length(df_list) == 1) df_list[[1]] else do.call(rbind, df_list)

p <- ggplot(df, aes(x = iter, y = elbo, colour = model)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = c("Delay-aware" = purple, "No-delay" = light_blue)) +
  labs(
    title = "ELBO over VB iterations",
    x = "Iteration",
    y = "Evidence Lower Bound (ELBO)",
    colour = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )

ggsave("src/paramresults/fig_elbo_comparison.png", p, width = 8, height = 4.5, units = "in", dpi = 300)

library(ggplot2)

# put into one data frame
df <- data.frame(
  iter       = seq_along(elbo_delay),
  elbo       = elbo_delay,
  param_diff = param_diff_delay[seq_along(elbo_delay)]
)

# Map param_diff linearly into the ELBO range (match min/max, not just scale)
rng_elbo <- range(df$elbo, na.rm = TRUE)
rng_par  <- range(df$param_diff, na.rm = TRUE)

scaleFactor <- diff(rng_elbo) / diff(rng_par)
offset      <- rng_elbo[1] - rng_par[1] * scaleFactor

p <- ggplot(df, aes(iter)) +
  geom_line(aes(y = elbo, colour = "ELBO"), linewidth = 1.1) +
  geom_line(aes(y = param_diff * scaleFactor + offset, colour = "Param diff"), linewidth = 1.1) +
  scale_y_continuous(
    name = "ELBO",
    sec.axis = sec_axis(~ (. - offset) / scaleFactor, name = "Parameter difference")
  ) +
  scale_color_manual(values = c("ELBO" = purple, "Param diff" = light_blue)) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top", plot.title = element_text(face = "bold")) +
  labs(title = "ELBO and Parameter Differences over VB iterations", x = "Iteration", colour = NULL)


ggsave("src/plots/fig_elbo_paramdiff.png", p, width = 8, height = 4.5, units = "in", dpi = 300)
