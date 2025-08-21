library(ggplot2)

# Thesis colors 
light_pink <- "#fbaee1"
light_blue <- "#A6CEE3"
purple     <- "#bd5bd3"

path_delay   <- "src/paramresults/posterior_val_delay_train.rds"
path_nodelay <- "src/paramresults/posterior_val_no_delay_train.rds"
obj_delay   <- readRDS(path_delay)
obj_nodelay <- readRDS(path_nodelay)

elbo_delay <- (as.numeric(obj_delay$elbo))[4:40]   
elbo_nodelay <- as.numeric(obj_nodelay$elbo)[4:40]
param_diff_delay <- (as.numeric(obj_delay$param_diffs))[4:40]   
param_diff_nodelay <- as.numeric(obj_nodelay$param_diffs)[4:40]


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
    x = "Iteration",
    y = "ELBO",
    colour = NULL
  ) +
  theme_classic(base_size = 20) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    legend.text = element_text(size = 29)
  )

ggsave("src/plots/fig_elbo_comparison.png", p, width = 8, height = 4.5, units = "in", dpi = 300)

library(ggplot2)

# Put into one data frame.
df <- data.frame(
  iter       = seq_along(elbo_delay),
  elbo       = elbo_delay,
  param_diff = param_diff_delay[seq_along(elbo_delay)]
)
# ELBO + Param-diff dual-axis plot 
# Align lengths explicitly
L <- min(length(elbo_delay), length(param_diff_delay))
df <- data.frame(
  iter       = seq_len(L),
  elbo       = elbo_delay[seq_len(L)],
  param_diff = param_diff_delay[seq_len(L)]
)

# Ranges
rng_elbo <- range(df$elbo, na.rm = TRUE)
rng_par  <- range(df$param_diff, na.rm = TRUE)

eps <- .Machine$double.eps * 100
if (diff(rng_par) < eps) {
  a <- 1.0
  b <- mean(rng_elbo) - mean(rng_par) * a
} else {
  a <- diff(rng_elbo) / diff(rng_par)                      
  b <- rng_elbo[1] - a * rng_par[1]                        
}

df$param_diff_mapped <- df$param_diff * a + b

p <- ggplot(df, aes(iter)) +
  geom_line(aes(y = elbo, colour = "ELBO"), linewidth = 1.1) +
  geom_line(aes(y = param_diff_mapped, colour = "Param diff"), linewidth = 1.1) +
  scale_y_continuous(
    name = "ELBO",
    sec.axis = sec_axis(~ if (abs(a) < eps) . else ((. - b) / a),
                        name = "Parameter difference")
  ) +
  scale_color_manual(values = c("ELBO" = "#bd5bd3", "Param diff" = "#A6CEE3")) +
  labs(x = "Iteration", colour = NULL) +
  theme_classic(base_size = 20) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    legend.text = element_text(size = 29)
  )

ggsave("src/plots/fig_elbo_paramdiff.png", p, width = 8, height = 4, units = "in", dpi = 300)
