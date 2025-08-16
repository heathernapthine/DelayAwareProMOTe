library(ggplot2)

# Thesis colors 
light_pink <- "#fbaee1"
light_blue <- "#A6CEE3"
purple     <- "#bd5bd3"

path_delay   <- "src/resultsmixeddata/posterior_val_delay_train.rds"
path_nodelay <- "src/resultsmixeddata/posterior_val_no_delay_train.rds"
obj_delay   <- readRDS(path_delay)
obj_nodelay <- readRDS(path_nodelay)

elbo_delay <- as.numeric(obj_delay$elbo)   
elbo_nodelay <- as.numeric(obj_nodelay$elbo) 

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

ggsave("src/plots/fig_elbo_comparison.png", p, width = 8, height = 4.5, units = "in", dpi = 300)

