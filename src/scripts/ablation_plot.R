library(ggplot2)
library(dplyr)
library(readr)
library(scales)
library(patchwork)
library(cowplot)

# Thesis color scheme 
light_pink <- "#fbaee1"
light_blue <- "#A6CEE3"
purple     <- "#bd5bd3"


directory <- "src/ablationresultsnew"
tops_path    <- file.path(directory, "ablation_prior_strength_top.csv")
shrink_path  <- file.path(directory, "ablation_prior_shrinkage_condition.csv")

tops   <- read.csv(tops_path, check.names = FALSE)
shrink <- read.csv(shrink_path, check.names = FALSE)
tops$var_scale   <- as.numeric(tops$var_scale)
shrink$var_scale <- as.numeric(shrink$var_scale)

# Output directory for figures
directory_plots <- file.path(directory, "plots")

# for clean x-axis labels
scale_x_s <- scale_x_continuous(breaks = c(0.25, 0.5, 1, 2, 4), labels = c("0.25", "0.5", "1", "2", "4"))

# figure A: Global sensitivity (ARI & NMI vs prior scale)
pA <- ggplot(tops, aes(x = var_scale)) +
  geom_line(aes(y = cluster_ari, color = "ARI"), linewidth = 1) +
  geom_point(aes(y = cluster_ari, color = "ARI"), size = 2) +
  geom_line(aes(y = cluster_nmi, color = "NMI"), linewidth = 1) +
  geom_point(aes(y = cluster_nmi, color = "NMI"), size = 2) +
  scale_color_manual(values = c("ARI" = purple, "NMI" = light_blue)) +
  scale_y_continuous(limits = c(0.6, 1), breaks = seq(0.6, 1, by = 0.05)) +
  scale_x_s +
  labs(x = "Prior variance scale (s)", y = "Cluster recovery (ARI / NMI)",
       color = NULL,
       title = "Effect of prior strength on cluster recovery (Dataset 3)") +
  theme_bw() +
  theme(legend.position = "top")

ggsave(file.path(directory_plots, "figA_prior_scale_vs_cluster_recovery.png"), pA, width = 7, height = 4.5, dpi = 300)


# figure B: SR_w by scale
shrink <- shrink %>%
  mutate(family = factor(family, levels = c("gaussian", "mixture2", "uniform")))

pal <- c(
  gaussian = light_pink,
  mixture2 = light_blue,
  uniform  = purple
)

# Common y-limits so all panels share the same log scale
ylims <- range(shrink$SR_w, na.rm = TRUE)

# Base plot (no caption will put in spare cell)
base <- ggplot(shrink, aes(x = factor(var_scale), y = SR_w, fill = family)) +
  geom_boxplot(outlier.alpha = 0.15) +
  geom_hline(yintercept = 1, linetype = 2, color = "gray40") +
  scale_y_log10() +
  coord_cartesian(ylim = ylims) +
  annotation_logticks(sides = "l") +
  scale_fill_manual(values = pal, drop = FALSE) +   # keep all keys in legend
  labs(x = "Prior variance scale (s)",
       y = "SR_w (log scale)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

# Three single-family panels
p_g  <- base %+% filter(shrink, family == "gaussian") + ggtitle("Gaussian")
p_m2 <- base %+% filter(shrink, family == "mixture2") + ggtitle("Mixture (2)")
p_u  <- base %+% filter(shrink, family == "uniform")  + ggtitle("Uniform")

# Legend extracted from a plot that contains all families
p_for_leg <- base + theme(legend.position = "right") +
  labs(fill = "Delay family")
leg <- cowplot::get_legend(p_for_leg + geom_boxplot())
cap_text <- "(SR_w): < 1; posterior variance smaller than prior variance,\n              = 1; posterior variance equals prior variance,\n              > 1; posterior variance larger than prior variance"
cap <- ggdraw() + draw_text(cap_text, x = 0, y = 1, hjust = 0, vjust = 1, size = 10)

# Stack legend over caption in the spare cell
leg_cap <- plot_grid(ggdraw(leg), cap, ncol = 1, rel_heights = c(2, 1))

# 2×2 grid: top row = Gaussian | Mixture(2); bottom row = Uniform | Legend+Caption
p_grid <- plot_grid(
  p_g, p_m2,
  p_u, leg_cap,
  ncol = 2, align = "hv", axis = "tblr"
)

p_final <- plot_grid(p_grid, ncol = 1, rel_heights = c(0.08, 1))

ggsave(file.path(directory_plots, "figB_shrinkage_SRw_by_scale_2x2.pdf"),
       p_final, width = 10, height = 7, device = cairo_pdf)
ggsave(file.path(directory_plots, "figB_shrinkage_SRw_by_scale_2x2.png"),
       p_final, width = 10, height = 7, dpi = 300)
