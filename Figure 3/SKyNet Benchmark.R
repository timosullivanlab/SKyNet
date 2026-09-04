# ==== Libraries & Config ====
library(tidyverse)
library(ggvenn)
library(patchwork)
library(CellChat)

elife_cols <- c(
  "CellChat" = "#EE99AA",
  "NicheNet" = "#E69F00",
  "SKyNet"   = "#999933"
)

elife_colors_im_strom <- c(
  "Immune->Immune"   = "#377EB8",
  "Immune->Stromal"  = "#4DAF4A",
  "Stromal->Immune"  = "#FF7F00",
  "Stromal->Stromal" = "#E41A1C"
)

cancers_list <- c("CESC", "COAD", "GBM", "HNSC", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "OV", "PAAD")

# ---- 1. Interaction Family Dot Plot & Clustering ----
# SKyNet
data_df <- read.csv("interactions_SKyNet_annotated.csv", row.names = 1) |>
  filter(interaction_family != "Other") |>
  mutate(direction = as.character(sign(interaction_score)))

data_df$X <- NULL

cancers_order <- c("LUAD", "LUSC", "PAAD", "KIRP", "KIRC", "CESC", "OV", "COAD", "HNSC", "LIHC", "GBM")

interaction_families <- c(
  "Chemokines", "Interleukins & cytokines", "TGF-beta superfamily",
  "TNF superfamily", "Growth factors",
  "Immune receptors & checkpoints", "MHC & MHC-like",
  "Galectins & lectins", "Selectins", "Cell adhesion",
  "Semaphorins", "Ephrins", "Notch ligands", "Apolipoproteins",
  "Collagens", "Laminins", "ECM glycoproteins",
  "Proteases", "Protease inhibitors", "Complement"
)

data_df <- data_df |>
  mutate(cancer = factor(cancer, levels = cancers_order),
         interaction_family = factor(interaction_family, levels = interaction_families))

p_main <- ggplot(data_df, aes(x = interaction_family, y = cancer)) +
  geom_jitter(
    aes(color = direction, size = abs(interaction_score)),
    width = 0.25, height = 0.25, alpha = 0.8, shape = 16
  ) +
  scale_y_discrete(limits = rev) +
  scale_x_discrete(limits = interaction_families) +
  scale_color_manual(
    name = "Direction",
    values = c("-1" = "#0072B2", "1" = "#D55E00"),
    labels = c("-1" = "Decreased survival", "1" = "Increased survival")
  ) +
  scale_size_continuous(name = "Interaction score", range = c(2.0, 5.0)) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(size = 9, color = "#222222", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9, color = "#222222"),
    axis.title.y = element_text(color = "#222222", margin = margin(r = 8)),
    panel.grid.major = element_line(color = "#F0F0F0", linewidth = 0.4),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    plot.margin = margin(10)
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 4)),
    size = guide_legend(override.aes = list(color = "grey40"))
  )

p_main

ggsave("LR_Individual_Interactions_DotPlot_SKyNet.pdf", plot = p_main, width = 10, height = 6, dpi = 300)

# CellChat
data_df <- read.csv("interactions_CellChat_annotated.csv", row.names = 1) |>
  filter(interaction_family != "Other") |>
  mutate(direction = ifelse(regulation == "Higher in Tumor", "-1", "1"))

data_df <- data_df |>
  mutate(cancer = factor(cancer, levels = cancers_order),
         interaction_family = factor(interaction_family, levels = interaction_families))

p_main <- ggplot(data_df, aes(x = interaction_family, y = cancer)) +
  geom_jitter(
    aes(color = direction),
    width = 0.25, height = 0.25, alpha = 0.8, shape = 16, size = 3.5
  ) +
  scale_y_discrete(limits = rev) +
  scale_x_discrete(limits = interaction_families) +
  scale_color_manual(
    name = "Direction",
    values = c("-1" = "#0072B2", "1" = "#D55E00"),
    labels = c("-1" = "Higher in Tumor", "1" = "Lower in Tumor")
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(size = 9, color = "#222222", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9, color = "#222222"),
    axis.title.y = element_text(color = "#222222", margin = margin(r = 8)),
    panel.grid.major = element_line(color = "#F0F0F0", linewidth = 0.4),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    plot.margin = margin(10)
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 4)),
    size = guide_legend(override.aes = list(color = "grey40"))
  )

p_main

ggsave("LR_Individual_Interactions_DotPlot_CellChat.pdf", plot = p_main, width = 10, height = 6, dpi = 300)

# NicheNet
data_df <- read.csv("interactions_NicheNet_annotated.csv", row.names = 1) |>
  filter(interaction_family != "Other") |>
  mutate(direction = ifelse(driven_state == "Tumor", "-1", "1"))

data_df <- data_df |>
  mutate(cancer = factor(cancer, levels = cancers_order),
         interaction_family = factor(interaction_family, levels = interaction_families))

p_main <- ggplot(data_df, aes(x = interaction_family, y = cancer)) +
  geom_jitter(
    aes(color = direction),
    width = 0.25, height = 0.25, alpha = 0.8, shape = 16, size = 3.5
  ) +
  scale_y_discrete(limits = rev) +
  scale_x_discrete(limits = interaction_families) +
  scale_color_manual(
    name = "Direction",
    values = c("-1" = "#0072B2", "1" = "#D55E00"),
    labels = c("-1" = "Higher in Tumor", "1" = "Lower in Tumor")
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(size = 9, color = "#222222", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9, color = "#222222"),
    axis.title.y = element_text(color = "#222222", margin = margin(r = 8)),
    panel.grid.major = element_line(color = "#F0F0F0", linewidth = 0.4),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    plot.margin = margin(10)
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 4)),
    size = guide_legend(override.aes = list(color = "grey40"))
  )

p_main

ggsave("LR_Individual_Interactions_DotPlot_NicheNet.pdf", plot = p_main, width = 10, height = 6, dpi = 300)


# ---- 2. Number of Interactions Across Networks ----
nints <- lapply(cancers_list, function(cancer) {
  skynet   <- read.csv(paste0(cancer, "_prioritized_permuted.csv"), row.names = 1) |> nrow()
  nichenet <- read.csv(paste0(cancer, "_NicheNet_Top_Interactions.csv")) |>
    filter(aupr_zscore > 2.58) |>
    nrow()
  cellchat <- read.csv(paste0(cancer, "_CellChat_Top_Interactions.csv")) |> nrow()
  
  c(SKyNet = skynet, NicheNet = nichenet, CellChat = cellchat)
}) |> setNames(cancers_list)

nints_long <- do.call(rbind, nints) |>
  as.data.frame() |>
  rownames_to_column(var = "Cancer") |>
  pivot_longer(
    cols = -Cancer,
    names_to = "Network",
    values_to = "Interactions"
  )

ggplot(nints_long, aes(x = Cancer, y = Interactions, fill = Network)) +
  geom_col(position = position_dodge(width = 0.65), width = 0.55, alpha = 0.85, color = "black", linewidth = 0.4) +
  labs(x = "Cancer", y = "Number of Interactions") +
  scale_fill_manual(values = elife_cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  theme_classic(base_size = 10) +
  theme(
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

ggsave("num_ints_by_network_by_cancer.pdf", units = "in", dpi = 300, height = 3, width = 5)

# ---- 3. Immune vs Stromal Interactions ----
all_ints <- lapply(cancers_list, function(cancer) {
  skynet <- read.csv(paste0(cancer, "_prioritized_permuted.csv"), row.names = 1) |>
    select(sender, receiver) |>
    mutate(network = "SKyNet")
  
  nichenet <- read.csv(paste0(cancer, "_NicheNet_Top_Interactions.csv")) |>
    filter(aupr_zscore > 2.58) |>
    select(sender, receiver) |>
    mutate(network = "NicheNet")
  
  cellchat <- read.csv(paste0(cancer, "_CellChat_Top_Interactions.csv")) |>
    select(sender = source, receiver = target) |>
    mutate(network = "CellChat")
  
  bind_rows(skynet, nichenet, cellchat) |> mutate(cancer = cancer)
})

all_ints_df <- bind_rows(all_ints)

immune_cells <- c(
  "IFN.primed.macrophage", "PVM", "Activated.CD8..T.cell", "Mature.NK", "Mast.cell", 
  "trNK", "cDC2b", "Treg", "Neutrophil", "Interstitial.macrophage", "LAM", "cDC1", 
  "BCell", "nCD8", "nCD4", "Exhausted.T.cell", "Plasma", "Inflammatory.macrophage", 
  "cMo", "Th17", "pDC", "Migratory.DC1", "SPP1..LAM", "ncMo", "iMo", "Microglia", 
  "CD8..Trm", "CXCL9.10.11..LAM", "IL22RA2..Langerhans", "Migratory.DC2", "Kupffer.cell", 
  "Immature.NK", "Alveolar.macrophage", "CD4..effector.T.cell", "IFN.primed.T.cell", 
  "ILC.progenitor", "Th1", "Langerin..cDC2b", "ILC3", "ILC2", "MAIT", "Activated.cDC1", 
  "Langerhans.cell"
)

all_ints_df <- all_ints_df |>
  mutate(
    sender_type = if_else(sender %in% immune_cells, "Immune", "Stromal"),
    receiver_type = if_else(receiver %in% immune_cells, "Immune", "Stromal"),
    int_type = paste0(sender_type, "->", receiver_type)
  )

# ---- 4. Immune vs Stromal Total ----
stacked_total <- all_ints_df |>
  count(network, int_type) |>
  group_by(network) |>
  mutate(
    prop = n / sum(n),
    percent_label = if_else(prop >= 0.05, sprintf("%.0f%%", prop * 100), "")
  ) |>
  ungroup() |>
  mutate(
    network = factor(network, levels = c("CellChat", "NicheNet", "SKyNet")),
    int_type = factor(int_type, levels = names(elife_colors_im_strom))
  )

ggplot(stacked_total, aes(x = network, y = prop, fill = int_type)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.4) +
  geom_text(
    aes(label = percent_label),
    position = position_fill(vjust = 0.5),
    color = "white",
    size = 2.6,
    fontface = "bold"
  ) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = elife_colors_im_strom, name = "Interaction Type") +
  labs(x = NULL, y = "Proportion of Interactions") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8.5, color = "#333333"),
    axis.text.y = element_text(size = 8.5, color = "#333333"),
    axis.title.y = element_text(size = 10, margin = margin(r = 8)),
    axis.line = element_line(linewidth = 0.4, color = "#444444"),
    axis.ticks = element_line(linewidth = 0.4, color = "#444444"),
    legend.position = "bottom",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8.5),
    legend.key.size = unit(0.4, "cm")
  )

ggsave("immune_stromal_all.pdf", units = "in", dpi = 300, height = 4, width = 5)

# ---- 5. Immune vs Stromal by Cancer ----
stacked_cancer <- all_ints_df |>
  count(cancer, network, int_type) |>
  group_by(cancer, network) |>
  mutate(
    prop = n / sum(n),
    percent_label = if_else(prop >= 0.05, sprintf("%.0f%%", prop * 100), "")
  ) |>
  ungroup() |>
  mutate(
    network = factor(network, levels = c("CellChat", "NicheNet", "SKyNet")),
    int_type = factor(int_type, levels = names(elife_colors_im_strom))
  )

ggplot(stacked_cancer, aes(x = network, y = prop, fill = int_type)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.4) +
  geom_text(
    aes(label = percent_label),
    position = position_fill(vjust = 0.5),
    color = "white",
    size = 2.6,
    fontface = "bold"
  ) +
  facet_wrap(~ cancer, nrow = 2, scales = "free_x") +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = elife_colors_im_strom, name = "Interaction Type") +
  labs(x = NULL, y = "Proportion of Interactions") +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "#F0F0F0", color = NA),
    strip.text = element_text(face = "bold", size = 10, color = "#222222"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8.5, color = "#333333"),
    axis.text.y = element_text(size = 8.5, color = "#333333"),
    axis.title.y = element_text(size = 10, margin = margin(r = 8)),
    axis.line = element_line(linewidth = 0.4, color = "#444444"),
    axis.ticks = element_line(linewidth = 0.4, color = "#444444"),
    legend.position = "bottom",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8.5),
    legend.key.size = unit(0.4, "cm"),
    panel.spacing = unit(0.8, "lines")
  )

ggsave("immune_stromal_by_cancer.pdf", units = "in", width = 11, height = 6.5, dpi = 300)

# ---- 6. Overlap Venn Diagrams ----
all_ints_lr <- lapply(cancers_list, function(cancer) {
  skynet <- read.csv(paste0(cancer, "_prioritized_permuted.csv"), row.names = 1) |>
    distinct(ligand, receptor) |>
    mutate(network = "SKyNet")
  
  nichenet <- read.csv(paste0(cancer, "_NicheNet_Top_Interactions.csv")) |>
    filter(aupr_zscore > 2.58) |>
    distinct(ligand, receptor) |>
    mutate(network = "NicheNet")
  
  cellchat <- read.csv(paste0(cancer, "_CellChat_Top_Interactions.csv")) |>
    distinct(ligand, receptor) |>
    mutate(network = "CellChat")
  
  bind_rows(skynet, nichenet, cellchat) |>
    mutate(cancer = cancer, lr_pair = paste0(ligand, "_", receptor))
})

all_ints_lr_df <- bind_rows(all_ints_lr)

# ---- 7. Venn Diagrams By Cancer ----
venn_plots <- lapply(cancers_list, function(cancer_id) {
  cancer_sub <- filter(all_ints_lr_df, cancer == cancer_id)
  venn_sets <- split(cancer_sub$lr_pair, cancer_sub$network)
  
  ggvenn(
    venn_sets,
    fill_color = elife_cols, 
    stroke_size = 0.6,
    set_name_size = 3.5,
    text_size = 5,
    show_percentage = FALSE
  ) +
    labs(title = cancer_id) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5, margin = margin(b = 6)),
      plot.margin = margin(10)
    )
})

combined_plot <- (wrap_plots(venn_plots, ncol = 4) & 
                    theme(
                      panel.spacing.x = unit(1.5, "cm"), 
                      panel.spacing.y = unit(1.2, "cm")
                    )) +
  plot_annotation(
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, color = "#111111", hjust = 0.5, margin = margin(b = 15))
    )
  )

combined_plot

ggsave("venn_diagram_by_cancer.pdf", plot = combined_plot, width = 10, height = 8, dpi = 300)

# ---- 8. Venn Diagrams Total ----
all_ints_lr_df <- mutate(all_ints_lr_df, cancer_lr = paste0(cancer, "_", lr_pair))

venn_sets_global <- split(all_ints_lr_df$cancer_lr, all_ints_lr_df$network)

global_venn <- ggvenn(
  venn_sets_global,
  fill_color = elife_cols,
  stroke_size = 1.0,
  set_name_size = 5.5,
  text_size = 6,
  show_percentage = FALSE
) +
  theme_void(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, color = "#111111", margin = margin(b = 20)),
    plot.margin = margin(20)
  )

global_venn

ggsave("venn_diagram_pred_ints.pdf", plot = global_venn, width = 6, height = 6, dpi = 300)

# ---- 9. Database Overlap Venn Diagram ----
venn_sets_db <- list(
  SKyNet = readRDS("lr_scdiffcom_v2.rds") |>
    unite("interaction_name", 2:6, sep = "_", na.rm = TRUE) |>
    pull(interaction_name) |>
    unique(),
  
  CellChat = CellChatDB.human$interaction |>
    pull(interaction_name) |>
    unique(),
  
  NicheNet = readRDS(url("https://zenodo.org/records/10229222/files/lr_network_human_allInfo_30112033.rds")) |>
    unite("interaction_name", 1:2, sep = "_", na.rm = TRUE) |>
    pull(interaction_name) |>
    unique()
)

db_venn <- ggvenn(
  venn_sets_db,
  fill_color = elife_cols,
  stroke_size = 1.0,      
  set_name_size = 5.5,    
  text_size = 6,        
  show_percentage = FALSE   
) +
  theme_void(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, color = "#111111", margin = margin(b = 20)),
    plot.margin = margin(20)
  )

ggsave("venn_diagram_lr_database.pdf", plot = db_venn, width = 6, height = 6, dpi = 300)