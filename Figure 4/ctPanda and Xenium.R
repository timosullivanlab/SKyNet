# ---- 1. Libraries ----
library(tidyverse)
library(readxl)
library(ggvenn)
library(ggplot2)

# ---- 2. Load Data ----
xenium_p_cutoff  <- 0.10
ctpanda_p_cutoff <- 0.10

skynet_df  <- read.csv("PAAD_prioritized_permuted.csv", row.names = 1)
ctpanda_df <- read_excel("paad_ctpanda.xlsx")
xenium_df  <- read.csv("all_spatialdm_lr_pairs_summary.csv")

endothelial <- c("Capillary.endothelial", "Vein.endothelial", "Artery.endothelial")

# ---- 3. Join Keys & Subsets ----
skynet_df <- skynet_df |>
  mutate(
    interaction_id = paste(sender, receiver, ligand, receptor, sep = "_"),
    sl_key = paste(sender, ligand, sep = "_"),
    lr_key = paste(ligand, receptor, sep = "_")
  ) |>
  filter(sender %in% endothelial | receiver %in% endothelial)

ctpanda_valid_sl <- ctpanda_df |>
  filter(pval < ctpanda_p_cutoff, `Direction (SKyNet)` == `Direction (ctPanda)`) |>
  mutate(sl_key = paste(Sender, Ligand, sep = "_")) |>
  pull(sl_key) |>
  unique()

patient_pval_cols <- colnames(xenium_df)[grep("pval", colnames(xenium_df))]
xenium_valid_lr <- xenium_df |>
  filter(if_any(all_of(patient_pval_cols), ~ .x < xenium_p_cutoff)) |>
  mutate(lr_key = paste(ligand, receptor, sep = "_")) |>
  pull(lr_key) |>
  unique()

ctpanda_validated_ids <- skynet_df |> filter(sl_key %in% ctpanda_valid_sl) |> pull(interaction_id)
xenium_validated_ids  <- skynet_df |> filter(lr_key %in% xenium_valid_lr)  |> pull(interaction_id)

total_interactions <- nrow(skynet_df)
validated_union    <- length(union(ctpanda_validated_ids, xenium_validated_ids))
unvalidated_count  <- total_interactions - validated_union

# ---- 4. Formatted Labels & Palette ----
venn_sets <- list(
  "ctPanda" = ctpanda_validated_ids,
  "Xenium"  = xenium_validated_ids
)

elife_cols <- c(
  "ctPanda" = "#77AADD",
  "Xenium"  = "#AA4499"
)

# ---- 5. Generate Framed Plot ----
p <- ggvenn(
  venn_sets,
  fill_color      = elife_cols,
  fill_alpha      = 0.45,
  stroke_color    = "black",
  stroke_size     = 1.0,
  set_name_size   = 5.5,
  text_size       = 6.0,
  show_percentage = FALSE
) +

  annotate(
    "rect", xmin = -2.2, xmax = 2.2, ymin = -1.6, ymax = 1.6,
    fill = "#F9F9F9", color = "#CCCCCC", linewidth = 0.8, linetype = "solid"
  ) +

  ggvenn(
    venn_sets,
    fill_color      = elife_cols,
    fill_alpha      = 0.45,
    stroke_color    = "white",
    stroke_size     = 1.0,
    set_name_size   = 5.5,
    text_size       = 6.0,
    show_percentage = FALSE
  )$layers +
  
  annotate(
    "rect", xmin = -2.1, xmax = -1.05, ymin = -1.55, ymax = -1.13,
    fill = "#FFFFFF", color = "#E0E0E0", linewidth = 0.5
  ) +
  annotate(
    "text", x = -2.0, y = -1.26, 
    label = paste0("Total Ints: ", total_interactions), 
    hjust = 0, fontface = "bold", size = 3.1, color = "#222222"
  ) +
  annotate(
    "text", x = -2.0, y = -1.42, 
    label = paste0("Validated: ", validated_union), 
    hjust = 0, fontface = "plain", size = 3.1, color = "#117733"
  ) +
  
  theme_void(base_size = 14) +
  theme(
    plot.margin = margin(20)
  )

print(p)

ggsave("paad_ctpanda_xenium.pdf", plot = p, width = 6, height = 4, dpi = 300)


# ---- 6. Donut Plot ----
plot_validation_donut <- function(df, key_col, valid_keys, method_name, fill_color, unit_label) {
  all_keys    <- unique(df[[key_col]])
  total_cnt   <- length(all_keys)
  valid_cnt   <- sum(all_keys %in% valid_keys)
  unvalid_cnt <- total_cnt - valid_cnt
  pct_str     <- sprintf("%.1f%%", (valid_cnt / total_cnt) * 100)
  
  unvalid_color <- "#ECECEC"
  
  donut_df <- data.frame(
    category = factor(c("Validated", "Unvalidated"), levels = c("Validated", "Unvalidated")),
    count    = c(valid_cnt, unvalid_cnt)
  )
  
  label_df <- data.frame(
    x_pos = c(2, 2),
    y_pos = c(total_cnt * 0.75, total_cnt * 0.25),
    count = c(valid_cnt, unvalid_cnt)
  )
  
  ggplot(donut_df, aes(x = 2, y = count, fill = category)) +
    geom_col(width = 0.7, color = "white", linewidth = 1.0) +
    geom_text(
      data        = label_df,
      aes(x = x_pos, y = y_pos, label = count),
      inherit.aes = FALSE,
      color       = "#222222",
      size        = 4.5,
      fontface    = "bold"
    ) +
    coord_polar(theta = "y") +
    xlim(c(0.5, 2.5)) +

    scale_fill_manual(
      values = c("Validated" = fill_color, "Unvalidated" = unvalid_color),
      labels = c("Validated", "Unvalidated")
    ) +

    annotate("text", x = 0.5, y = 0, label = pct_str, size = 7.0, fontface = "bold", color = "#222222") +
  
    labs(
      title    = paste(method_name, "Validation of Endothelial Interactions"),
      subtitle = paste0("Total Unique Endothelial ", unit_label, ": ", total_cnt)
    ) +
    theme_void(base_size = 12, base_family = "sans") +
    theme(
      legend.position = "bottom",
      legend.title    = element_blank(),
      legend.text     = element_text(size = 10, color = "#333333"),
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 13, color = "#222222", margin = margin(t = 8, b = 2)),
      plot.subtitle   = element_text(hjust = 0.5, face = "plain", size = 10, color = "#555555", margin = margin(b = 6)),
      plot.margin     = margin(12, 12, 12, 12)
    )
}

# ---- 7. Generate Plots ----
p_ctpanda <- plot_validation_donut(skynet_df, "sl_key", ctpanda_valid_sl, "ctPanda", elife_cols[["ctPanda"]], "S-L Pairs")
p_xenium  <- plot_validation_donut(skynet_df, "lr_key", xenium_valid_lr,  "Xenium",  elife_cols[["Xenium"]],  "L-R Pairs")

print(p_ctpanda)
print(p_xenium)

ggsave("ctpanda_unique_sl_validation.pdf", plot = p_ctpanda, width = 4.5, height = 4.8, dpi = 300)
ggsave("xenium_unique_lr_validation.pdf",  plot = p_xenium,  width = 4.5, height = 4.8, dpi = 300)

