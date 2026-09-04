# ---- 1. Libraries & Config ----
library(tidyverse)
library(ggpubr)

cancers <- c("HNSC", "LUSC", "LUAD", "Kidney", "LIHC", "COAD", "BRCA", 
             "GBM",  "PAAD", "CESC", "OV", "SKCM")

elife_cols <- c("xCell" = "#88CCEE", "CellTypist" = "#44AA99", "JadeCell" = "#117733")

# ---- 2. Load & Combine Data ----
methods <- c("xCell", "CellTypist", "JadeCell")

combined_df <- map_dfr(cancers, function(cancer) {
  map_dfr(methods, function(m) {
    read.csv(paste0(m, "_", cancer, "_scores.csv")) |>
      setNames(c("CellType", "Score")) |>
      mutate(Method = m, Cancer = cancer)
  })
}) |>
  mutate(
    Method = factor(Method, levels = c("xCell", "CellTypist", "JadeCell")),
    Cancer = factor(Cancer, levels = cancers)
  )

# ---- 3. Combined Multi-Cancer Mean Bar Plot ----
combined_summary <- combined_df |>
  group_by(Cancer, Method) |>
  summarize(max_top = mean(Score, na.rm = TRUE) + (sd(Score, na.rm = TRUE) / sqrt(n())), .groups = "drop")

start_y_combined <- max(combined_summary$max_top, na.rm = TRUE) * 0.85

p_bar_combined <- ggplot(combined_df, aes(x = Cancer, y = Score, fill = Method)) +
  stat_summary(
    fun = "mean", geom = "bar",
    position = position_dodge(width = 0.65),
    width = 0.55, alpha = 0.85, color = "black", linewidth = 0.4
  ) +
  stat_summary(
    fun.data = "mean_se", geom = "errorbar",
    position = position_dodge(width = 0.65),
    width = 0.12, linewidth = 0.6
  ) +
  geom_pwc(
    aes(group = Method),
    method = "t.test",
    label = "p.signif",
    ref.group = 3,
    y.position = start_y_combined,
    step.increase = 0.04,
    tip.length = 0.005,
    size = 0.4,
    vjust = 0.5,
    label.size = 3,
    hide.ns = TRUE
  ) +
  scale_fill_manual(values = elife_cols) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  annotate("segment", x = 0.5, xend = 7.4, y = -0.17, color = "#4477AA", linewidth = 4, alpha = 0.8) +
  annotate("text", x = 4.0, y = -0.23, label = "In Atlas", fontface = "bold", size = 3.5, color = "#4477AA") +
  annotate("segment", x = 7.5, xend = 12.4, y = -0.17, color = "#EE6677", linewidth = 4, alpha = 0.8) +
  annotate("text", x = 10.0, y = -0.23, label = "Not In Atlas", fontface = "bold", size = 3.5, color = "#EE6677") +
  labs(x = "Cancer Type", y = "Mean AUPRC") +
  theme_classic(base_size = 11) +
  theme(
    axis.title.x = element_text(margin = margin(t = 30)),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.title = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 15, l = 10)
  )

p_bar_combined

ggsave("bar_plot_combined_cancers_AUPRC.pdf", plot = p_bar_combined,
       units = "in", dpi = 300, height = 5, width = 6)
