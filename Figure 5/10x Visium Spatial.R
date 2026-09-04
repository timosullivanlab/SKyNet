# 0. Libraries and Functions ===================================================
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggsignif)
library(ggrepel)

SeuratAnalysis1 <- function(SeuratObject){
  SeuratObject <- NormalizeData(SeuratObject)
  SeuratObject <- FindVariableFeatures(SeuratObject)
  SeuratObject <- ScaleData(SeuratObject)
  SeuratObject <- RunPCA(SeuratObject, features = VariableFeatures(SeuratObject))
  return(SeuratObject)
}

SeuratAnalysis2 <- function(SeuratObject, res = 0.5, n){
  SeuratObject <- FindNeighbors(SeuratObject, dims = 1:n)
  SeuratObject <- FindClusters(SeuratObject, resolution = res)
  SeuratObject <- RunUMAP(SeuratObject, dims = 1:n)
  return(SeuratObject)
}

compute_jaccard_gene <- function(obj, mod1, mod2) {
  obj <- AddModuleScore(obj, mod1, name = "Module1")
  obj <- AddModuleScore(obj, mod2, name = "Module2")
  
  counts     <- obj@meta.data[, c("Module11", "Module21")]
  counts$co  <- (counts$Module11 > 0) + (counts$Module21 > 0)
  
  sum(counts$co == 2) / sum(counts$co > 0)
}

# 1. Load Data =================================================================
files <- list.files(recursive = TRUE)
lineage <- read.csv("LUAD_Lineage_Markers.csv", header = T)
lineage <- sapply(lineage, function(x) str_replace_all(x, "\\s", ""))
lineage <- as.data.frame(lineage)
data <- read.csv('LUAD_prioritized_permuted.csv', row.names = 1)

luad <- list(c("P10_T1", "P10_B1"), c("P15_T1", "P15_B1"),
             c("P16_T1", "P10_B1"), c("P24_T1", "P24_B1"),
             c("P25_T1", "P25_B1"))

# 2. Load Seurat Objects =======================================================
seurat_list <- map(luad, function(pair) {
  map(pair, function(pat) {
    h5_file      <- files[str_detect(files, paste0("^", pat, ".*\\.h5$"))]
    short_h5_file <- str_split_i(h5_file, "/", 2)
    
    obj <- Load10X_Spatial(data.dir = pat,
                           filename  = short_h5_file,
                           assay     = "Spatial",
                           slice     = pat)
    
    obj@project.name <- pat
    obj$orig.ident   <- str_split_1(pat, "_")[[1]]
    
    # QC
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    
    obj <- subset(obj, subset = percent.mt < 20 & nCount_RNA >= 100 & 
                    abs(nCount_RNA - mean(nCount_RNA)) <= 2 * sd(nCount_RNA) & 
                    abs(nFeature_RNA - mean(nFeature_RNA)) <= 2 * sd(nFeature_RNA))
    
    obj <- SeuratAnalysis1(obj)
    obj <- SeuratAnalysis2(obj, n = 5)
    
    for (cell in colnames(lineage)) {
      genes <- lineage[,cell]
      genes <- genes[genes %in% row.names(obj)]
      obj <- AddModuleScore(obj, list(genes), name = cell)
      
      # Rename columns
      colnames(obj@meta.data)[
        colnames(obj@meta.data) == paste0(cell, "1")
      ] <- cell
    }
    
    obj
  })
})

# 3. Jaccard Similarity Index - All Genes Slopegraphs ==========================
interactions <- data |>
  distinct(ligand, receptor, .keep_all = T) |>
  mutate(ligand = str_split(ligand, pattern = "_"),
         receptor = str_split(receptor, pattern = "_"))

# --- Step 1: Loop over interactions only ---
for (i in 1:nrow(interactions)) {
  mod1 <- interactions$ligand[i]
  mod2 <- interactions$receptor[i]
  
  d <- map(seurat_list, function(pair) {
    c(
      compute_jaccard_gene(pair[[1]], mod1, mod2),  # Tumor
      compute_jaccard_gene(pair[[2]], mod1, mod2)   # Adjacent
    )
  })
  
  name <- str_c(
    str_c(unlist(mod1), collapse = "_"),
    str_c(unlist(mod2), collapse = "_"),
    sep = "+"
  )
  
  d_long <- as.data.frame(d) |>
    setNames(c("P10", "P15", "P16", "P24", "P25")) |>
    mutate(Type = c("Tumor", "Adjacent")) |>
    pivot_longer(-Type, names_to = "Patient", values_to = "Value") |>
    mutate(
      Patient = factor(Patient, levels = c("P10", "P15", "P16", "P24", "P25")),
      Type    = factor(Type,    levels = c("Adjacent", "Tumor"))
    )
  
  a <- do.call(rbind, d) |> as.data.frame()
  colnames(a) <- c('Tumor', 'Adjacent')
  
  pval <- t.test(
    a$Tumor,
    a$Adjacent,
    paired = TRUE
  )$p.value
  
  ggplot(d_long, aes(x = Type, y = Value, group = Patient)) +
    geom_line(color = 'black', linewidth = 0.5) +
    geom_point(aes(fill = Type), stroke = 0.5, color = 'black', shape = 21, size = 2.8) +
    scale_fill_manual(values = c("Adjacent" = "black", "Tumor" = "#FF321A")) +
    scale_y_continuous(limits = c(0, max(d_long$Value)*1.15), expand = c(0, 0)) +
    scale_x_discrete(expand = expansion(add = c(0.55, 0.55))) +
    geom_signif(
      comparisons = list(c("Tumor", "Adjacent")),
      y_position = max(d_long$Value) * 1.05,
      annotations = paste0("p = ", signif(pval, 3)),
      textsize = 3
    ) +
    labs(title = paste0(name, ' Axis'), x = NULL, y = 'Jaccard Similarity', legend = NULL) +
    theme_classic(base_size = 10, base_family = 'Arial') +
    theme(legend.title = element_blank(),
          axis.ticks.length = unit(2, 'mm'),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'))
  
  ggsave(
    filename = paste0("../LUAD Figures/", name, "_Jaccard.pdf"),
    device = cairo_pdf,
    dpi = 300, 
    height = 5, 
    width = 3, 
    units = "in"
  )
}

# 4. Example Individual Colocalization - Gene ==================================
patient <- 'P25'

patient.1 <- paste0(patient, "_T1")
patient.2 <- paste0(patient, "_B1")

# Pat 1
h5_file <- files[str_detect(files, pattern = paste0("^", patient.1, ".*\\.h5$"))]
short_h5_file <- str_split_i(h5_file, pattern = "/", 2)
pat.1 <- Load10X_Spatial(data.dir = patient.1,
                         filename = short_h5_file,
                         assay = 'Spatial',
                         slice = patient.1)

pat.1@project.name <- patient.1
pat.1$orig.ident <- str_split_1(patient.1, pattern = "_")[[1]]

# Pat 2
h5_file <- files[str_detect(files, pattern = paste0("^", patient.2, ".*\\.h5$"))]
short_h5_file <- str_split_i(h5_file, pattern = "/", 2)
pat.2 <- Load10X_Spatial(data.dir = patient.2,
                         filename = short_h5_file,
                         assay = 'Spatial',
                         slice = patient.2)

pat.2@project.name <- patient.2
pat.2$orig.ident <- str_split_1(patient.2, pattern = "_")[[1]]

merged.pats <- merge(pat.1, pat.2, add.cell.ids = c(pat.1@project.name, pat.2@project.name))
merged.pats <- JoinLayers(merged.pats)

merged.pats <- SeuratAnalysis1(merged.pats)
merged.pats <- SeuratAnalysis2(merged.pats, n = 5)

# Module Scores
mod1 <- list('COL5A2')
mod2 <- list('ITGA1', 'ITGB1')

name1 <- paste0(mod1, collapse = "+")
name2 <- paste0(mod2, collapse = "+")

merged.pats <- AddModuleScore(merged.pats, mod1, name = "Module1")
merged.pats <- AddModuleScore(merged.pats, mod2, name = "Module2")

counts <- merged.pats@meta.data[, c("Module11", "Module21")]

counts$coexpress <- "No"
counts$coexpress[counts[,"Module11"] > 0] <- name1
counts$coexpress[counts[,"Module21"] > 0] <- name2
counts$coexpress[counts[,"Module11"] > 0 & counts[,"Module21"] > 0] <- "Both"

merged.pats$coexpress <- counts$coexpress

SpatialPlot(merged.pats, group.by = "coexpress", ncol = 1, pt.size.factor = 2.5) &
  scale_fill_manual(values = setNames(
    c("black", "red", "yellow", "green"),
    c("No", name1, name2, "Both")
  ),
  breaks = c("No", name1, name2, "Both")) &
  theme(
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

ggsave(paste0('../', name1, '_', name2, '_', patient, '.pdf'), dpi = 300, height = 6, width = 4, units = 'in')

# 5. Jaccard Similarity Scatter Plot - Genes ===================================
interactions <- data |>
  distinct(ligand, receptor, .keep_all = T) |>
  mutate(ligand = str_split(ligand, pattern = "_"),
         receptor = str_split(receptor, pattern = "_"))

# --- Step 1: Loop over interactions ---
ints <- list()

for (i in 1:nrow(interactions)) {
  mod1 <- interactions$ligand[i]
  mod2 <- interactions$receptor[i]
  
  d <- map(seurat_list, function(pair) {
    c(
      compute_jaccard_gene(pair[[1]], mod1, mod2),  # Tumor
      compute_jaccard_gene(pair[[2]], mod1, mod2)   # Adjacent
    )
  })
  
  name <- str_c(
    str_c(unlist(mod1), collapse = "_"),
    str_c(unlist(mod2), collapse = "_"),
    sep = "+"
  )
  
  diffs <- unlist(lapply(d, function(x) x[1] - x[2]))
  
  p <- t.test(vapply(d, function(x) x[1], numeric(1)), vapply(d, function(x) x[2], numeric(1)), paired = T)[3]
  
  # Conditions for three categories
  if (mean(unlist(d)) < 0.1) {
    cat <- 3 # Not expressed
  } else if (p < 0.05) {
    cat <- 1 # Tumor high
  } else {
    cat <- 2 # Both high
  }
  
  ints[[i]] <- c(name, mean(vapply(d, function(x) x[1], numeric(1))), mean(vapply(d, function(x) x[2], numeric(1))), cat, p, sign(interactions$interaction_score[i]), interactions$interaction_score[i])
}

ints <- do.call(rbind, ints)
ints <- as.data.frame(ints)
colnames(ints) <- c('Interaction', 'Tumor', 'Adjacent', 'Category', 'pval', 'Direction', 'Score')
ints$Category <- factor(ints$Category, levels = c(1, 2, 3))
ints$Interaction <- vapply(ints$Interaction, function(x) x[[1]], character(1))
ints$Interaction <- ifelse
class(ints$Tumor) <- 'numeric'
class(ints$Adjacent) <- 'numeric'
class(ints$pval) <- 'numeric'
ints$Direction <- ifelse(ints$Direction == 1, "+", "-")
class(ints$Score) <- 'numeric'

ggplot(ints, aes(x = Tumor, y = Adjacent)) +

  geom_point(aes(color = Category, shape = Direction, size = abs(Score)), alpha = 0.85) +
  stat_ellipse(aes(group = Category, color = Category), level = 0.90, alpha = 0.4, linewidth = 0.4) +
  
  geom_text_repel(
    aes(label = Interaction), 
    size = 4,
    max.overlaps = 5
  ) +

  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey50", linewidth = 0.35) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "grey50", linewidth = 0.35) +
  
  theme_classic(base_size = 8, base_family = "sans") +
  
  theme(
  
    axis.line = element_line(linewidth = 0.4, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    axis.title = element_text(size = 12, face = "plain", color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.35, "cm"),
    legend.spacing.y = unit(0.1, "cm"),
    legend.background = element_blank(),
    
    plot.margin = margin(t = 8, r = 8, b = 8, l = 8, unit = "pt")
  ) +
  
  labs(
    x = "Mean Tumor Jaccard",
    y = "Mean Adjacent Jaccard"
  ) +

  scale_color_manual(
    name = "Category",
    values = c(
      "1" = "#D55E00", # Color for Tumor high
      "2" = "#0072B2", # Color for Both high
      "3" = "#999999"  # Color for Not expressed
    ),
    labels = c(
      "1" = "Tumor high",
      "2" = "Both high",
      "3" = "Not expressed"
    )
  ) +
  
  scale_size_continuous(
    name = "Interaction Score", 
    breaks = c(0, 1000, 2000, 3000, 4000), 
    limits = c(0, NA),
    range = c(1, 4)
  ) +
  
  scale_shape_manual(
    name = "Interaction Direction",
    labels = c(
      "+" = "Positive",
      "-" = "Negative"
    ),
    values = c("+" = 15, "-" = 16)
  ) +
  
  annotate("text", x = Inf, y = -Inf, label = "Higher in Tumor", 
           hjust = 1.1, vjust = -0.5, size = 3, fontface = "italic", color = "grey30") +
  annotate("text", x = -Inf, y = Inf, label = "Higher in Adjacent", 
           hjust = -0.1, vjust = -0.5, angle = 270, size = 3, fontface = "italic", color = "grey30")

ggsave('LUAD_SKyNet_Jaccard_Scatter.pdf', units = 'in', dpi = 300, height = 8, width = 12)
