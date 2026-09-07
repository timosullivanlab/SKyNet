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
cancer <- 'SKCM'
seurat.obj <- readRDS(paste0(cancer, '_final.rds'))
seurat.obj <- SeuratAnalysis1(seurat.obj)
seurat.obj <- SeuratHarmony(seurat.obj, harmony.by = 'Patient', n = 9)

# ---- 2. Load and Pre-process xCell Gene Sets ----
xcell_raw <- read.csv('xcell.csv', stringsAsFactors = FALSE)

# Strip dataset source suffix to get canonical base cell type
xcell_raw$Base_CellType <- sub("_(HPCA|IRIS|FANTOM|ENCODE|NOVERS|BACON|BLUEPRINT|LM22|VALI).*$", "", xcell_raw$Celltype_Source_ID)
gene_cols <- setdiff(names(xcell_raw), c("Celltype_Source_ID", "X..of.genes", "# of genes", "Base_CellType"))

xcell_gene_list <- list()
for (ct in unique(xcell_raw$Base_CellType)) {
  genes_matrix <- xcell_raw[xcell_raw$Base_CellType == ct, gene_cols]
  genes_vec <- unique(na.omit(as.character(unlist(genes_matrix))))
  xcell_gene_list[[ct]] <- genes_vec[genes_vec != ""]
}


# ---- 3. Map Seurat Ground Truth -> xCell Base Cell Types ----

# # SKCM
# xcell_map <- c(
#   "pDC"                     = "pDC",
#   "BCell"                   = "B-cells",
#   "Plasma"                  = "Plasma cells",
#   "nCD4"                    = "CD4+ naive T-cells",
#   "Activated.CD8..T.cell"   = "CD8+ Tem",
#   "CD4..effector.T.cell"    = "CD4+ Tem",
#   "Treg"                    = "Tregs",
#   "nCD8"                    = "CD8+ naive T-cells",
#   "LAM"                     = "Macrophages",
#   "IFN.primed.macrophage"   = "Macrophages M1",
#   "Interstitial.macrophage" = "Macrophages",
#   "Mast.cell"               = "Mast cells",
#   "Vein.endothelial"        = "Endothelial cells",
#   "Fibro"                   = "Fibroblasts",
#   "cDC2b"                   = "cDC",
#   "iMo"                     = "Monocytes",
#   "CXCL9.10.11..LAM"        = "Macrophages M1",
#   "SPP1..LAM"               = "Macrophages M2",
#   "cDC1"                    = "cDC",
#   "Basal"                   = "Epithelial cells",
#   "Epithelial"              = "Epithelial cells",
#   "Lymphatic.endothelial"   = "ly Endothelial cells",
#   "trNK"                    = "NK cells",
#   "Exhausted.T.cell"        = "CD8+ T-cells",
#   "Tcm"                     = "CD8+ Tcm"
# )

# HNSC
# xcell_map <- c(
#   "pDC"                   = "pDC",
#   "Plasma"                = "Plasma cells",
#   "BCell"                 = "B-cells",
#   "Activated.CD8..T.cell" = "CD8+ Tem",
#   "CD8..Trm"              = "CD8+ Tem",
#   "Exhausted.T.cell"      = "CD8+ T-cells",
#   "nCD4"                  = "CD4+ naive T-cells",
#   "Treg"                  = "Tregs",
#   "PVM"                   = "Macrophages",
#   "cDC2b"                 = "cDC",
#   "Vein.endothelial"      = "Endothelial cells",
#   "Artery.endothelial"    = "Endothelial cells",
#   "Lymphatic.endothelial" = "ly Endothelial cells",
#   "Mast.cell"             = "Mast cells",
#   "Epithelial"            = "Epithelial cells",
#   "Fibro"                 = "Fibroblasts",
#   "Mature.NK"             = "NK cells",
#   "trNK"                  = "NK cells"
# )

# # BRCA
# xcell_map <- c(
#   "pDC" = "pDC",
#   "PVM" = "Macrophages",
#   "BCell" = "B-cells",
#   "Plasma" = "Plasma cells",
#   "nCD4" = "CD4+ naive T-cells",
#   "Activated.CD8..T.cell" = "CD8+ Tem",
#   "gdT" = "Tgd cells",
#   "Treg" = "Tregs",
#   "Immature.NK" = "NK cells",
#   "nCD8" = "CD8+ naive T-cells",
#   "LAM" = "Macrophages",
#   "IFN.primed.macrophage" = "Macrophages M1",
#   "Interstitial.macrophage" = "Macrophages",
#   "Mast.cell" = "Mast cells",
#   "Vein.endothelial" = "Endothelial cells",
#   "Fibro" = "Fibroblasts",
#   "cDC2b" = "cDC",
#   "IL22RA2..Langerhans" = "DC",
#   "iMo" = "Monocytes",
#   "CXCL9.10.11..LAM" = "Macrophages",
#   "cDC1" = "cDC",
#   "Migratory.DC1" = "cDC",
#   "Basal" = "Epithelial cells",
#   "Lymphatic.endothelial" = "ly Endothelial cells",
#   "Artery.endothelial" = "Endothelial cells",
#   "Capillary.endothelial" = "mv Endothelial cells",
#   "Luminal.epithelial.hr" = "Epithelial cells",
#   "Luminal.epithelial.secretory" = "Epithelial cells",
#   "Mature.NK" = "NK cells",
#   "trNK" = "NK cells",
#   "Inflammatory.macrophage" = "Macrophages M1"
# )

# # CESC
# xcell_map <- c(
#   "PVM" = "Macrophages",
#   "BCell" = "B-cells",
#   "Plasma" = "Plasma cells",
#   "nCD4" = "CD4+ naive T-cells",
#   "Treg" = "Tregs",
#   "nCD8" = "CD8+ naive T-cells",
#   "LAM" = "Macrophages",
#   "cMo" = "Monocytes",
#   "Mast.cell" = "Mast cells",
#   "Vein.endothelial" = "Endothelial cells",
#   "Fibro" = "Fibroblasts",
#   "cDC2b" = "cDC",
#   "Neutrophil" = "Neutrophils",
#   "cDC1" = "cDC",
#   "Basal" = "Epithelial cells",
#   "Epithelial" = "Epithelial cells",
#   "Lymphatic.endothelial" = "ly Endothelial cells",
#   "Ciliated.cell" = "Epithelial cells",
#   "Artery.endothelial" = "Endothelial cells",
#   "Capillary.endothelial" = "mv Endothelial cells",
#   "Mature.NK" = "NK cells",
#   "trNK" = "NK cells",
#   "Exhausted.T.cell" = "CD8+ Tem",
#   "Inflammatory.macrophage" = "Macrophages M1"
# )

# # COAD
# xcell_map <- c(
#   "pDC" = "pDC",
#   "PVM" = "Macrophages",
#   "BCell" = "B-cells",
#   "Plasma" = "Plasma cells",
#   "nCD4" = "CD4+ naive T-cells",
#   "Activated.CD8..T.cell" = "CD8+ T-cells",
#   "gdT" = "Tgd cells",
#   "CD4..effector.T.cell" = "CD4+ Tem",
#   "Treg" = "Tregs",
#   "LAM" = "Macrophages",
#   "IFN.primed.macrophage" = "Macrophages M1",
#   "cMo" = "Monocytes",
#   "ncMo" = "Monocytes",
#   "Interstitial.macrophage" = "Macrophages",
#   "Mast.cell" = "Mast cells",
#   "Vein.endothelial" = "Endothelial cells",
#   "Fibro" = "Fibroblasts",
#   "cDC2b" = "cDC",
#   "cDC1" = "cDC",
#   "Migratory.DC1" = "cDC",
#   "Basal" = "Epithelial cells",
#   "Epithelial" = "Epithelial cells",
#   "Myofibroblast" = "Fibroblasts",
#   "Lymphatic.endothelial" = "ly Endothelial cells",
#   "Artery.endothelial" = "Endothelial cells",
#   "Capillary.endothelial" = "mv Endothelial cells",
#   "Mature.NK" = "NK cells",
#   "trNK" = "NK cells",
#   "Exhausted.T.cell" = "CD8+ T-cells",
#   "Inflammatory.macrophage" = "Macrophages M1",
#   "Tcm" = "CD4+ Tcm",
#   "Th17" = "CD4+ T-cells",
#   "Intestine.goblet.cell" = "Epithelial cells",
#   "Distal.enterocyte" = "Epithelial cells",
#   "Intestinal.stem.cell" = "Epithelial cells"
# )

# # GBM
# xcell_map <- c(
#   "BCell" = "B-cells",
#   "nCD4" = "CD4+ naive T-cells",
#   "Activated.CD8..T.cell" = "CD8+ T-cells",
#   "Treg" = "Tregs",
#   "nCD8" = "CD8+ naive T-cells",
#   "LAM" = "Macrophages",
#   "cMo" = "Monocytes",
#   "ncMo" = "Monocytes",
#   "Vein.endothelial" = "Endothelial cells",
#   "Fibro" = "Fibroblasts",
#   "cDC2b" = "cDC",
#   "iMo" = "Monocytes",
#   "SPP1..LAM" = "Macrophages",
#   "cDC1" = "cDC",
#   "Artery.endothelial" = "Endothelial cells",
#   "Capillary.endothelial" = "mv Endothelial cells",
#   "Mature.NK" = "NK cells",
#   "trNK" = "NK cells"
# )

# # Kidney
# xcell_map <- c(
#   "pDC" = "pDC",
#   "BCell" = "B-cells",
#   "Plasma" = "Plasma cells",
#   "nCD4" = "CD4+ naive T-cells",
#   "Treg" = "Tregs",
#   "Fibro" = "Fibroblasts",
#   "cDC2b" = "cDC",
#   "iMo" = "Monocytes",
#   "cDC1" = "cDC",
#   "trNK" = "NK cells"
# )

# # LIHC
# xcell_map <- c(
#   "pDC" = "pDC",
#   "PVM" = "Macrophages",
#   "BCell" = "B-cells",
#   "Plasma" = "Plasma cells",
#   "nCD4" = "CD4+ naive T-cells",
#   "Activated.CD8..T.cell" = "CD8+ T-cells",
#   "Treg" = "Tregs",
#   "nCD8" = "CD8+ naive T-cells",
#   "LAM" = "Macrophages",
#   "cMo" = "Monocytes",
#   "ncMo" = "Monocytes",
#   "Mast.cell" = "Mast cells",
#   "Fibro" = "Fibroblasts",
#   "cDC2b" = "cDC",
#   "Kupffer.cell" = "Macrophages",
#   "cDC1" = "cDC",
#   "Migratory.DC2" = "cDC",
#   "Lymphatic.endothelial" = "ly Endothelial cells",
#   "Macrovascular.endothelial" = "mv Endothelial cells",
#   "Sinusoidal.endothelial" = "Endothelial cells",
#   "Cholangiocyte" = "Epithelial cells",
#   "Mature.NK" = "NK cells",
#   "trNK" = "NK cells"
# )

# # LUAD
# xcell_map <- c(
#   "pDC" = "pDC",
#   "PVM" = "Macrophages",
#   "BCell" = "B-cells",
#   "Plasma" = "Plasma cells",
#   "nCD4" = "CD4+ naive T-cells",
#   "Activated.CD8..T.cell" = "CD8+ Tem",
#   "CD8..Trm" = "CD8+ Tem",
#   "IFN.primed.T.cell" = "CD4+ T-cells",
#   "gdT" = "Tgd cells",
#   "Th1" = "Th1 cells",
#   "CD4..effector.T.cell" = "CD4+ Tem",
#   "Treg" = "Tregs",
#   "Immature.NK" = "NK cells",
#   "nCD8" = "CD8+ naive T-cells",
#   "LAM" = "Macrophages",
#   "IFN.primed.macrophage" = "Macrophages M1",
#   "cMo" = "Monocytes",
#   "ncMo" = "Monocytes",
#   "Interstitial.macrophage" = "Macrophages",
#   "Mast.cell" = "Mast cells",
#   "Club.cell" = "Epithelial cells",
#   "Alveolar.cell.type.2" = "Epithelial cells",
#   "Vein.endothelial" = "Endothelial cells",
#   "Fibro" = "Fibroblasts",
#   "Alveolar.macrophage" = "Macrophages",
#   "cDC2b" = "cDC",
#   "IL22RA2..Langerhans" = "DC",
#   "Langerin..cDC2b" = "cDC",
#   "iMo" = "Monocytes",
#   "Neutrophil" = "Neutrophils",
#   "CXCL9.10.11..LAM" = "Macrophages",
#   "SPP1..LAM" = "Macrophages",
#   "cDC1" = "cDC",
#   "Migratory.DC1" = "cDC",
#   "Migratory.DC2" = "cDC",
#   "Squamous.epithelial" = "Epithelial cells",
#   "Lung.goblet.cell" = "Epithelial cells",
#   "Basal" = "Epithelial cells",
#   "Epithelial" = "Epithelial cells",
#   "TFF1.S100P..tumor.cell" = "Epithelial cells",
#   "Inflammatory.epithelial" = "Epithelial cells",
#   "Myofibroblast" = "Fibroblasts",
#   "Lymphatic.endothelial" = "ly Endothelial cells",
#   "Alveolar.cell.type.1" = "Epithelial cells",
#   "Ciliated.cell" = "Epithelial cells",
#   "Mature.NK" = "NK cells",
#   "Exhausted.T.cell" = "CD8+ Tem"
# )

# # LUSC
# xcell_map <- c(
#   "pDC" = "pDC",
#   "PVM" = "Macrophages",
#   "BCell" = "B-cells",
#   "Plasma" = "Plasma cells",
#   "nCD4" = "CD4+ naive T-cells",
#   "Activated.CD8..T.cell" = "CD8+ Tem",
#   "gdT" = "Tgd cells",
#   "Th1" = "Th1 cells",
#   "Treg" = "Tregs",
#   "Immature.NK" = "NK cells",
#   "nCD8" = "CD8+ naive T-cells",
#   "LAM" = "Macrophages",
#   "cMo" = "Monocytes",
#   "ncMo" = "Monocytes",
#   "Interstitial.macrophage" = "Macrophages",
#   "Mast.cell" = "Mast cells",
#   "Club.cell" = "Epithelial cells",
#   "Alveolar.cell.type.2" = "Epithelial cells",
#   "Vein.endothelial" = "Endothelial cells",
#   "Fibro" = "Fibroblasts",
#   "Alveolar.macrophage" = "Macrophages",
#   "cDC2b" = "cDC",
#   "iMo" = "Monocytes",
#   "CXCL9.10.11..LAM" = "Macrophages",
#   "SPP1..LAM" = "Macrophages",
#   "cDC1" = "cDC",
#   "Lung.goblet.cell" = "Epithelial cells",
#   "Myofibroblast" = "Fibroblasts",
#   "Lymphatic.endothelial" = "ly Endothelial cells",
#   "Alveolar.cell.type.1" = "Epithelial cells",
#   "Ciliated.cell" = "Epithelial cells",
#   "Mature.NK" = "NK cells",
#   "trNK" = "NK cells",
#   "Exhausted.T.cell" = "CD8+ T-cells",
#   "Inflammatory.macrophage" = "Macrophages M1"
# )

# # OV
# xcell_map <- c(
#   "pDC" = "pDC",
#   "PVM" = "Macrophages",
#   "BCell" = "B-cells",
#   "Activated.CD8..T.cell" = "CD8+ Tem",
#   "nCD8" = "CD8+ naive T-cells",
#   "LAM" = "Macrophages",
#   "IFN.primed.macrophage" = "Macrophages M1",
#   "Mast.cell" = "Mast cells",
#   "Vein.endothelial" = "Endothelial cells",
#   "Fibro" = "Fibroblasts",
#   "cDC2b" = "cDC",
#   "iMo" = "Monocytes",
#   "SPP1..LAM" = "Macrophages",
#   "Epithelial" = "Epithelial cells",
#   "Lymphatic.endothelial" = "ly Endothelial cells",
#   "Artery.endothelial" = "Endothelial cells",
#   "Capillary.endothelial" = "mv Endothelial cells",
#   "Mature.NK" = "NK cells",
#   "trNK" = "NK cells"
# )

# # PAAD
# xcell_map <- c(
#   "PVM" = "Macrophages",
#   "BCell" = "B-cells",
#   "Plasma" = "Plasma cells",
#   "nCD4" = "CD4+ naive T-cells",
#   "Activated.CD8..T.cell" = "CD8+ Tem",
#   "Treg" = "Tregs",
#   "nCD8" = "CD8+ naive T-cells",
#   "LAM" = "Macrophages",
#   "IFN.primed.macrophage" = "Macrophages M1",
#   "cMo" = "Monocytes",
#   "Interstitial.macrophage" = "Macrophages",
#   "Mast.cell" = "Mast cells",
#   "Vein.endothelial" = "Endothelial cells",
#   "Fibro" = "Fibroblasts",
#   "cDC2b" = "cDC",
#   "Neutrophil" = "Neutrophils",
#   "Artery.endothelial" = "Endothelial cells",
#   "Capillary.endothelial" = "Endothelial cells",
#   "Erythrocyte" = "Erythrocytes",
#   "Mature.NK" = "NK cells",
#   "trNK" = "NK cells"
# )

# ---- 4. Add Module Scores & Calculate Scalar AUPRCs ----
auprc_scores <- numeric(0)
plot_list <- list()

for (cell in names(xcell_map)) {
  message(paste0('Analyzing: ', cell))
  
  # Load gene module - take top 20 genes present in scRNA-seq
  xc_name <- xcell_map[[cell]]
  gene_list <- xcell_gene_list[[xc_name]]
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
      subtitle = paste0("AUPRC = ", round(auc_val, 3), "\n(xCell: ", xc_name, ")"),
      x = "Recall",
      y = "Precision"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 8, color = "gray30"),
      axis.title = element_text(size = 9)
    )
}

message("XCell Average AUPRC: ", round(mean(auprc_scores, na.rm = TRUE), 3))

pr_grid <- wrap_plots(plot_list, ncol = 4) +
  plot_annotation(
    title = paste0("XCell Lineage Markers Precision-Recall Curves: ", cancer),
    subtitle = paste0("Average AUPRC: ", round(mean(auprc_scores, na.rm = TRUE), 3)),
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 13, color = "gray20")
    )
  )

pr_grid

# 5. Save Grid to PDF ----
ggsave(paste0("xCell_All_PR_Curves_Grid_", cancer, ".pdf"), plot = pr_grid, width = 20, height = 24)

# 6. Save Raw Numbers ----
write.csv(file = paste0('xCell_', cancer, '_scores.csv'), auprc_scores)
