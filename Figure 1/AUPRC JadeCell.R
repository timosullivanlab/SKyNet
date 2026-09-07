# ---- Libraries ----
library(Seurat)
library(ROCR)
library(tidyverse)
library(harmony)
library(patchwork)

# ---- Seurat Functions ----
SeuratAnalysis1 <- function(SeuratObject){
  SeuratObject <- NormalizeData(SeuratObject)
  SeuratObject <- FindVariableFeatures(SeuratObject)
  SeuratObject <- ScaleData(SeuratObject)
  SeuratObject <- RunPCA(SeuratObject, features = VariableFeatures(SeuratObject))
  return(SeuratObject)
}

SeuratHarmony <- function(SeuratObject, harmony.by, res = 0.5, n){
  SeuratObject <- RunHarmony(SeuratObject, group.by.vars = harmony.by)
  SeuratObject <- FindNeighbors(SeuratObject, dims = 1:n, reduction = "harmony")
  SeuratObject <- FindClusters(SeuratObject, resolution = res)
  SeuratObject <- RunUMAP(SeuratObject, dims = 1:n, reduction = "harmony")
  return(SeuratObject)
}

# ---- 1. Load in Seurat RDS ----
cancer <- 'PAAD'

seurat.obj <- readRDS(paste0(cancer, '_final.rds'))
seurat.obj <- SeuratAnalysis1(seurat.obj)
ElbowPlot(seurat.obj)
seurat.obj <- SeuratHarmony(seurat.obj, harmony.by = 'Patient', n = 10)

# ---- 2. Load and Pre-process JadeCell Gene Sets ----
jadecell_modules <- read.csv('jadecell.csv', row.names = 1)

intersect(colnames(jadecell_modules), unique(seurat.obj$celltype))

# ---- 3. Add Module Scores & Calculate Scalar AUPRCs ----

# SKCM
# shared_cells <- c("pDC", "BCell", "Plasma", "nCD4", "Activated.CD8..T.cell",
#                   "CD4..effector.T.cell", "Treg", "nCD8", "LAM",
#                   "IFN.primed.macrophage", "Interstitial.macrophage",
#                   "Mast.cell", "Vein.endothelial", "Fibro", "cDC2b", "iMo",
#                   "CXCL9.10.11..LAM", "SPP1..LAM", "cDC1", "Basal",
#                   "Epithelial", "Lymphatic.endothelial", "trNK",
#                   "Exhausted.T.cell", "Tcm")

# # HNSC
# shared_cells <- c("pDC", "Plasma", "BCell", "Activated.CD8..T.cell", "CD8..Trm",
#                   "Exhausted.T.cell", "nCD4", "Treg", "PVM", "cDC2b",
#                   "Vein.endothelial", "Artery.endothelial", "Lymphatic.endothelial",
#                   "Mast.cell", "Epithelial", "Fibro", "Mature.NK", "trNK")

# # BRCA
# shared_cells <- c("pDC", "PVM", "BCell", "Plasma", "nCD4", "Activated.CD8..T.cell",
#                   "gdT", "Treg", "Immature.NK", "nCD8", "LAM", "IFN.primed.macrophage",
#                   "Interstitial.macrophage", "Mast.cell", "Vein.endothelial", "Fibro",
#                   "cDC2b", "IL22RA2..Langerhans", "iMo", "CXCL9.10.11..LAM", "cDC1",
#                   "Migratory.DC1", "Basal", "Lymphatic.endothelial", "Artery.endothelial",
#                   "Capillary.endothelial", "Luminal.epithelial.hr",
#                   "Luminal.epithelial.secretory", "Mature.NK", "trNK", "Inflammatory.macrophage")

# # CESC
# shared_cells <- c("PVM", "BCell", "Plasma", "nCD4", "Treg", "nCD8", "LAM", "cMo",
#                   "Mast.cell", "Vein.endothelial", "Fibro", "cDC2b", "Neutrophil",
#                   "cDC1", "Basal", "Epithelial", "Lymphatic.endothelial", "Ciliated.cell",
#                   "Artery.endothelial", "Capillary.endothelial", "Mature.NK",
#                   "trNK", "Exhausted.T.cell", "Inflammatory.macrophage")

# # COAD
# shared_cells <- c("pDC", "PVM", "BCell", "Plasma", "nCD4", "Activated.CD8..T.cell",
#                   "gdT", "CD4..effector.T.cell", "Treg", "LAM", "IFN.primed.macrophage",
#                   "cMo", "ncMo", "Interstitial.macrophage", "Mast.cell",
#                   "Vein.endothelial", "Fibro", "cDC2b", "cDC1", "Migratory.DC1",
#                   "Basal", "Epithelial", "Myofibroblast", "Lymphatic.endothelial",
#                   "Artery.endothelial", "Capillary.endothelial", "Mature.NK", "trNK",
#                   "Exhausted.T.cell", "Inflammatory.macrophage", "Tcm", "Th17",
#                   "Intestine.goblet.cell", "Distal.enterocyte", "Intestinal.stem.cell")

# # GBM
# shared_cells <- c("BCell", "nCD4", "Activated.CD8..T.cell", "Treg", "nCD8",
#                   "LAM", "cMo", "ncMo", "Vein.endothelial", "Fibro", "cDC2b",
#                   "iMo", "SPP1..LAM", "cDC1", "Artery.endothelial",
#                   "Capillary.endothelial", "Mature.NK", "trNK")

# # Kidney
# shared_cells <- c("pDC", "BCell", "Plasma", "nCD4", "Treg", "Fibro", "cDC2b", "iMo", "cDC1", "trNK")

# # LIHC
# shared_cells <- c("pDC", "PVM", "BCell", "Plasma", "nCD4", "Activated.CD8..T.cell",
#                   "Treg", "nCD8", "LAM", "cMo", "ncMo", "Mast.cell", "Fibro", "cDC2b",
#                   "Kupffer.cell", "cDC1", "Migratory.DC2", "Lymphatic.endothelial",
#                   "Macrovascular.endothelial", "Sinusoidal.endothelial",
#                   "Cholangiocyte", "Mature.NK", "trNK")

# # LUAD
# shared_cells <- c("pDC", "PVM", "BCell", "Plasma", "nCD4", "Activated.CD8..T.cell",
#                   "CD8..Trm", "IFN.primed.T.cell", "gdT", "Th1", "CD4..effector.T.cell",
#                   "Treg", "Immature.NK", "nCD8", "LAM", "IFN.primed.macrophage",
#                   "cMo", "ncMo", "Interstitial.macrophage", "Mast.cell", "Club.cell",
#                   "Alveolar.cell.type.2", "Vein.endothelial", "Fibro", "Alveolar.macrophage",
#                   "cDC2b", "IL22RA2..Langerhans", "Langerin..cDC2b", "iMo", "Neutrophil",
#                   "CXCL9.10.11..LAM", "SPP1..LAM", "cDC1", "Migratory.DC1", "Migratory.DC2",
#                   "Squamous.epithelial", "Lung.goblet.cell", "Basal", "Epithelial",
#                   "TFF1.S100P..tumor.cell", "Inflammatory.epithelial", "Myofibroblast",
#                   "Lymphatic.endothelial", "Alveolar.cell.type.1", "Ciliated.cell",
#                   "Mature.NK", "Exhausted.T.cell")

# # LUSC
# shared_cells <- c("pDC", "PVM", "BCell", "Plasma", "nCD4", "Activated.CD8..T.cell",
#                   "gdT", "Th1", "Treg", "Immature.NK", "nCD8", "LAM", "cMo", "ncMo",
#                   "Interstitial.macrophage", "Mast.cell", "Club.cell", "Alveolar.cell.type.2",
#                   "Vein.endothelial", "Fibro", "Alveolar.macrophage", "cDC2b", "iMo",
#                   "CXCL9.10.11..LAM", "SPP1..LAM", "cDC1", "Lung.goblet.cell",
#                   "Myofibroblast", "Lymphatic.endothelial", "Alveolar.cell.type.1",
#                   "Ciliated.cell", "Mature.NK", "trNK", "Exhausted.T.cell", "Inflammatory.macrophage")

# # OV
# shared_cells <- c("pDC", "PVM", "BCell", "Activated.CD8..T.cell", "nCD8", "LAM",
#                   "IFN.primed.macrophage", "Mast.cell", "Vein.endothelial", "Fibro",
#                   "cDC2b", "iMo", "SPP1..LAM", "Epithelial", "Lymphatic.endothelial",
#                   "Artery.endothelial", "Capillary.endothelial", "Mature.NK", "trNK")

# # PAAD
# shared_cells <- c("PVM", "BCell", "Plasma", "nCD4", "Activated.CD8..T.cell", "Treg",
#                   "nCD8", "LAM", "IFN.primed.macrophage", "cMo", "Interstitial.macrophage",
#                   "Mast.cell", "Vein.endothelial", "Fibro", "cDC2b", "Neutrophil",
#                   "Artery.endothelial", "Capillary.endothelial", "Erythrocyte",
#                   "Mature.NK", "trNK")

# 4. ---- Run AUPRC ----
auprc_scores <- numeric(0)
plot_list <- list()

for (cell in shared_cells) {
  message(paste0('Analyzing: ', cell))
  
  # Load gene module - take top 20 genes present in scRNA-seq
  gene_list <- jadecell_modules[, cell]
  gene_list <- gene_list[gene_list %in% rownames(seurat.obj)]
  gene_list <- head(gene_list, 20)
  
  # Module score
  seurat.obj <- AddModuleScore(seurat.obj, features = list(gene_list), name = cell)
  seurat.obj[[cell]] <- seurat.obj[[paste0(cell, '1')]]
  seurat.obj[[paste0(cell, '1')]] <- NULL
  
  # Ground Truth vs Prediction
  truth <- ifelse(seurat.obj$celltype == cell, 1, 0)
  truth <- factor(truth, levels = c(0, 1))
  cell_module <- seurat.obj@meta.data[, cell, drop = F]
  
  # Run Area Under Precision-Recall Curve
  pred <- prediction(predictions = cell_module, labels = truth)
  perf_pr <- performance(pred, measure = "prec", x.measure = "rec")
  perf_auc <- performance(pred, measure = "aucpr")
  
  # Curve Values
  pr_df <- data.frame(
    Recall = perf_pr@x.values[[1]],
    Precision = perf_pr@y.values[[1]]
  )
  pr_df$Precision[is.nan(pr_df$Precision)] <- 1
  
  # Extract area under curve
  auc_val <- perf_auc@y.values[[1]]
  auprc_scores[[cell]] <- auc_val
  cat("  -> AUPRC:", round(auc_val, 3), "\n")
  
  # Save Plot
  plot_list[[cell]] <- ggplot(pr_df, aes(x = Recall, y = Precision)) +
    geom_line(color = "#1f77b4", linewidth = 1) +
    scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0)) +
    scale_y_continuous(limits = c(0, 1.05), expand = c(0, 0.01)) +
    theme_classic(base_size = 10) +
    labs(
      title = cell,
      subtitle = paste0("AUPRC = ", round(auc_val, 3)),
      x = "Recall",
      y = "Precision"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 8, color = "gray30"),
      axis.title = element_text(size = 9)
    )
}

message("JadeCell Average AUPRC: ", round(mean(auprc_scores, na.rm = TRUE), 3))

pr_grid <- wrap_plots(plot_list, ncol = 4) +
  plot_annotation(
    title = paste0("JadeCell Lineage Markers Precision-Recall Curves: ", cancer),
    subtitle = paste0("Average AUPRC: ", round(mean(auprc_scores, na.rm = TRUE), 3)),
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 13, color = "gray20")
    )
  )

pr_grid

# ---- 5. Save Grid to PDF ----
ggsave(paste0("JadeCell_All_PR_Curves_Grid_", cancer, ".pdf"), plot = pr_grid, width = 14, height = 10)

# ---- 6. Save Raw Numbers ----
write.csv(file = paste0('JadeCell_', cancer, '_scores.csv'), auprc_scores)
