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
cancer <- 'HNSC'
seurat.obj <- readRDS(paste0(cancer, '_final.rds'))
seurat.obj <- SeuratAnalysis1(seurat.obj)
seurat.obj <- SeuratHarmony(seurat.obj, harmony.by = 'Patient', n = 8)

# ---- 2. Load and Pre-process CellTypist Gene Sets ----
celltypist_modules <- read.csv('celltypist.csv', row.names = 1)

# ---- 3. Map Seurat Ground Truth -> CellTypist Base Cell Types ----

# SKCM
# celltypist_map <- c(
#   "pDC"                     = "pDC",
#   "BCell"                   = "B.cells",
#   "Plasma"                  = "Plasma.cells",
#   "nCD4"                    = "Tcm.Naive.helper.T.cells",
#   "Activated.CD8..T.cell"   = "Tem.Temra.cytotoxic.T.cells",
#   "CD4..effector.T.cell"    = "Tem.Effector.helper.T.cells",
#   "Treg"                    = "Regulatory.T.cells",
#   "nCD8"                    = "Tcm.Naive.cytotoxic.T.cells",
#   "LAM"                     = "Macrophages",
#   "IFN.primed.macrophage"   = "Macrophages",
#   "Interstitial.macrophage" = "Macrophages",
#   "Mast.cell"               = "Mast.cells",
#   "Vein.endothelial"        = "Endothelial.cells",
#   "Fibro"                   = "Fibroblasts",
#   "cDC2b"                   = "DC2",
#   "iMo"                     = "Monocytes",
#   "CXCL9.10.11..LAM"        = "Macrophages",
#   "SPP1..LAM"               = "Macrophages",
#   "cDC1"                    = "DC1",
#   "Basal"                   = "Epithelial.cells",
#   "Epithelial"              = "Epithelial.cells",
#   "Lymphatic.endothelial"   = "Endothelial.cells",
#   "trNK"                    = "NK.cells",
#   "Exhausted.T.cell"        = "Tem.Trm.cytotoxic.T.cells",
#   "Tcm"                     = "Tcm.Naive.helper.T.cells"
# )

# # HNSC
# celltypist_map <- c(
#   "pDC"                   = "pDC",
#   "Plasma"                = "Plasma.cells",
#   "BCell"                 = "B.cells",
#   "Activated.CD8..T.cell" = "Tem.Temra.cytotoxic.T.cells",
#   "CD8..Trm"              = "Trm.cytotoxic.T.cells",
#   "Exhausted.T.cell"      = "Tem.Trm.cytotoxic.T.cells",
#   "nCD4"                  = "Tcm.Naive.helper.T.cells",
#   "Treg"                  = "Regulatory.T.cells",
#   "PVM"                   = "Macrophages",
#   "cDC2b"                 = "DC2",
#   "Vein.endothelial"      = "Endothelial.cells",
#   "Artery.endothelial"    = "Endothelial.cells",
#   "Lymphatic.endothelial" = "Endothelial.cells",
#   "Mast.cell"             = "Mast.cells",
#   "Epithelial"            = "Epithelial.cells",
#   "Fibro"                 = "Fibroblasts",
#   "Mature.NK"             = "CD16..NK.cells",
#   "trNK"                  = "NK.cells"
# )

# # BRCA
# celltypist_map <- c(
#   "pDC" = "pDC",
#   "PVM" = "Macrophages",
#   "BCell" = "B cells",
#   "Plasma" = "Plasma cells",
#   "nCD4" = "Tcm/Naive helper T cells",
#   "Activated.CD8..T.cell" = "Tem/Temra cytotoxic T cells",
#   "gdT" = "gamma-delta T cells",
#   "Treg" = "Regulatory T cells",
#   "Immature.NK" = "NK cells",
#   "nCD8" = "Tcm/Naive cytotoxic T cells",
#   "LAM" = "Macrophages",
#   "IFN.primed.macrophage" = "Macrophages",
#   "Interstitial.macrophage" = "Macrophages",
#   "Mast.cell" = "Mast cells",
#   "Vein.endothelial" = "Endothelial cells",
#   "Fibro" = "Fibroblasts",
#   "cDC2b" = "DC2",
#   "IL22RA2..Langerhans" = "DC",
#   "iMo" = "Monocytes",
#   "CXCL9.10.11..LAM" = "Macrophages",
#   "cDC1" = "DC1",
#   "Migratory.DC1" = "Migratory DCs",
#   "Basal" = "Epithelial cells",
#   "Lymphatic.endothelial" = "Endothelial cells",
#   "Artery.endothelial" = "Endothelial cells",
#   "Capillary.endothelial" = "Endothelial cells",
#   "Luminal.epithelial.hr" = "Epithelial cells",
#   "Luminal.epithelial.secretory" = "Epithelial cells",
#   "Mature.NK" = "NK cells",
#   "trNK" = "NK cells",
#   "Inflammatory.macrophage" = "Macrophages"
# )

# # CESC
# celltypist_map <- c(
#   "PVM" = "Macrophages",
#   "BCell" = "B cells",
#   "Plasma" = "Plasma cells",
#   "nCD4" = "Tcm/Naive helper T cells",
#   "Treg" = "Regulatory T cells",
#   "nCD8" = "Tcm/Naive cytotoxic T cells",
#   "LAM" = "Macrophages",
#   "cMo" = "Classical monocytes",
#   "Mast.cell" = "Mast cells",
#   "Vein.endothelial" = "Endothelial cells",
#   "Fibro" = "Fibroblasts",
#   "cDC2b" = "DC2",
#   "Neutrophil" = "Neutrophils",
#   "cDC1" = "DC1",
#   "Basal" = "Epithelial cells",
#   "Epithelial" = "Epithelial cells",
#   "Lymphatic.endothelial" = "Endothelial cells",
#   "Ciliated.cell" = "Epithelial cells",
#   "Artery.endothelial" = "Endothelial cells",
#   "Capillary.endothelial" = "Endothelial cells",
#   "Mature.NK" = "NK cells",
#   "trNK" = "Transitional NK",
#   "Exhausted.T.cell" = "Tem/Effector helper T cells PD1+",
#   "Inflammatory.macrophage" = "Macrophages"
# )

# # COAD
# celltypist_map <- c(
#   "pDC" = "pDC",
#   "PVM" = "Macrophages",
#   "BCell" = "B cells",
#   "Plasma" = "Plasma cells",
#   "nCD4" = "Tcm/Naive helper T cells",
#   "Activated.CD8..T.cell" = "Tem/Temra cytotoxic T cells",
#   "gdT" = "gamma-delta T cells",
#   "CD4..effector.T.cell" = "Tem/Effector helper T cells",
#   "Treg" = "Regulatory T cells",
#   "LAM" = "Alveolar macrophages",
#   "IFN.primed.macrophage" = "Macrophages",
#   "cMo" = "Classical monocytes",
#   "ncMo" = "Non-classical monocytes",
#   "Interstitial.macrophage" = "Macrophages",
#   "Mast.cell" = "Mast cells",
#   "Vein.endothelial" = "Endothelial cells",
#   "Fibro" = "Fibroblasts",
#   "cDC2b" = "DC2",
#   "cDC1" = "DC1",
#   "Migratory.DC1" = "Migratory DCs",
#   "Basal" = "Epithelial cells",
#   "Epithelial" = "Epithelial cells",
#   "Myofibroblast" = "Fibroblasts",
#   "Lymphatic.endothelial" = "Endothelial cells",
#   "Artery.endothelial" = "Endothelial cells",
#   "Capillary.endothelial" = "Endothelial cells",
#   "Mature.NK" = "NK cells",
#   "trNK" = "NK cells",
#   "Exhausted.T.cell" = "Tem/Temra cytotoxic T cells",
#   "Inflammatory.macrophage" = "Macrophages",
#   "Tcm" = "Tcm/Naive helper T cells",
#   "Th17" = "Type 17 helper T cells",
#   "Intestine.goblet.cell" = "Epithelial cells",
#   "Distal.enterocyte" = "Epithelial cells",
#   "Intestinal.stem.cell" = "Epithelial cells"
# )

# # GBM
# celltypist_map <- c(
#   "BCell" = "B cells",
#   "nCD4" = "Tcm/Naive helper T cells",
#   "Activated.CD8..T.cell" = "Tem/Temra cytotoxic T cells",
#   "Treg" = "Regulatory T cells",
#   "nCD8" = "Tcm/Naive cytotoxic T cells",
#   "LAM" = "Macrophages",
#   "cMo" = "Classical monocytes",
#   "ncMo" = "Non-classical monocytes",
#   "Vein.endothelial" = "Endothelial cells",
#   "Fibro" = "Fibroblasts",
#   "cDC2b" = "DC2",
#   "iMo" = "Monocytes",
#   "SPP1..LAM" = "Macrophages",
#   "cDC1" = "DC1",
#   "Artery.endothelial" = "Endothelial cells",
#   "Capillary.endothelial" = "Endothelial cells",
#   "Mature.NK" = "NK cells",
#   "trNK" = "Transitional NK"
# )

# # Kidney
# celltypist_map <- c(
#   "pDC" = "pDC",
#   "BCell" = "B cells",
#   "Plasma" = "Plasma cells",
#   "nCD4" = "Tcm/Naive helper T cells",
#   "Treg" = "Regulatory T cells",
#   "Fibro" = "Fibroblasts",
#   "cDC2b" = "DC2",
#   "iMo" = "Monocytes",
#   "cDC1" = "DC1",
#   "trNK" = "NK cells"
# )

# # LIHC
# celltypist_map <- c(
#   "pDC" = "pDC",
#   "PVM" = "Macrophages",
#   "BCell" = "B cells",
#   "Plasma" = "Plasma cells",
#   "nCD4" = "Tcm/Naive helper T cells",
#   "Activated.CD8..T.cell" = "Tem/Temra cytotoxic T cells",
#   "Treg" = "Regulatory T cells",
#   "nCD8" = "Tcm/Naive cytotoxic T cells",
#   "LAM" = "Alveolar macrophages",
#   "cMo" = "Classical monocytes",
#   "ncMo" = "Non-classical monocytes",
#   "Mast.cell" = "Mast cells",
#   "Fibro" = "Fibroblasts",
#   "cDC2b" = "DC2",
#   "Kupffer.cell" = "Kupffer cells",
#   "cDC1" = "DC1",
#   "Migratory.DC2" = "Migratory DCs",
#   "Lymphatic.endothelial" = "Endothelial cells",
#   "Macrovascular.endothelial" = "Endothelial cells",
#   "Sinusoidal.endothelial" = "Endothelial cells",
#   "Cholangiocyte" = "Epithelial cells",
#   "Mature.NK" = "NK cells",
#   "trNK" = "CD16- NK cells"
# )

# # LUAD
# celltypist_map <- c(
#   "pDC" = "pDC",
#   "PVM" = "Macrophages",
#   "BCell" = "B cells",
#   "Plasma" = "Plasma cells",
#   "nCD4" = "Tcm/Naive helper T cells",
#   "Activated.CD8..T.cell" = "Tem/Temra cytotoxic T cells",
#   "CD8..Trm" = "Trm cytotoxic T cells",
#   "IFN.primed.T.cell" = "Tem/Effector helper T cells",
#   "gdT" = "gamma-delta T cells",
#   "Th1" = "Type 1 helper T cells",
#   "CD4..effector.T.cell" = "Tem/Effector helper T cells",
#   "Treg" = "Regulatory T cells",
#   "Immature.NK" = "Transitional NK",
#   "nCD8" = "Tcm/Naive cytotoxic T cells",
#   "LAM" = "Macrophages",
#   "IFN.primed.macrophage" = "Macrophages",
#   "cMo" = "Classical monocytes",
#   "ncMo" = "Non-classical monocytes",
#   "Interstitial.macrophage" = "Macrophages",
#   "Mast.cell" = "Mast cells",
#   "Club.cell" = "Epithelial cells",
#   "Alveolar.cell.type.2" = "Epithelial cells",
#   "Vein.endothelial" = "Endothelial cells",
#   "Fibro" = "Fibroblasts",
#   "Alveolar.macrophage" = "Alveolar macrophages",
#   "cDC2b" = "DC2",
#   "IL22RA2..Langerhans" = "Migratory DCs",
#   "Langerin..cDC2b" = "DC2",
#   "iMo" = "Intermediate macrophages",
#   "Neutrophil" = "Neutrophils",
#   "CXCL9.10.11..LAM" = "Macrophages",
#   "SPP1..LAM" = "Macrophages",
#   "cDC1" = "DC1",
#   "Migratory.DC1" = "Migratory DCs",
#   "Migratory.DC2" = "Migratory DCs",
#   "Squamous.epithelial" = "Epithelial cells",
#   "Lung.goblet.cell" = "Epithelial cells",
#   "Basal" = "Epithelial cells",
#   "Epithelial" = "Epithelial cells",
#   "TFF1.S100P..tumor.cell" = "Epithelial cells",
#   "Inflammatory.epithelial" = "Epithelial cells",
#   "Myofibroblast" = "Fibroblasts",
#   "Lymphatic.endothelial" = "Endothelial cells",
#   "Alveolar.cell.type.1" = "Epithelial cells",
#   "Ciliated.cell" = "Epithelial cells",
#   "Mature.NK" = "NK cells",
#   "Exhausted.T.cell" = "Tem/Temra cytotoxic T cells"
# )

# # LUSC
# celltypist_map <- c(
#   "pDC" = "pDC",
#   "PVM" = "Macrophages",
#   "BCell" = "B cells",
#   "Plasma" = "Plasma cells",
#   "nCD4" = "Tcm/Naive helper T cells",
#   "Activated.CD8..T.cell" = "Tem/Temra cytotoxic T cells",
#   "gdT" = "gamma-delta T cells",
#   "Th1" = "Type 1 helper T cells",
#   "Treg" = "Regulatory T cells",
#   "Immature.NK" = "Transitional NK",
#   "nCD8" = "Tcm/Naive cytotoxic T cells",
#   "LAM" = "Macrophages",
#   "cMo" = "Classical monocytes",
#   "ncMo" = "Non-classical monocytes",
#   "Interstitial.macrophage" = "Macrophages",
#   "Mast.cell" = "Mast cells",
#   "Club.cell" = "Epithelial cells",
#   "Alveolar.cell.type.2" = "Epithelial cells",
#   "Vein.endothelial" = "Endothelial cells",
#   "Fibro" = "Fibroblasts",
#   "Alveolar.macrophage" = "Alveolar macrophages",
#   "cDC2b" = "DC2",
#   "iMo" = "Monocytes",
#   "CXCL9.10.11..LAM" = "Macrophages",
#   "SPP1..LAM" = "Macrophages",
#   "cDC1" = "DC1",
#   "Lung.goblet.cell" = "Epithelial cells",
#   "Myofibroblast" = "Fibroblasts",
#   "Lymphatic.endothelial" = "Endothelial cells",
#   "Alveolar.cell.type.1" = "Epithelial cells",
#   "Ciliated.cell" = "Epithelial cells",
#   "Mature.NK" = "NK cells",
#   "trNK" = "NK cells",
#   "Exhausted.T.cell" = "Tem/Effector helper T cells PD1+",
#   "Inflammatory.macrophage" = "Macrophages"
# )

# # OV
# celltypist_map <- c(
#   "pDC" = "pDC",
#   "PVM" = "Macrophages",
#   "BCell" = "B cells",
#   "Activated.CD8..T.cell" = "Tem/Temra cytotoxic T cells",
#   "nCD8" = "Tcm/Naive cytotoxic T cells",
#   "LAM" = "Macrophages",
#   "IFN.primed.macrophage" = "Macrophages",
#   "Mast.cell" = "Mast cells",
#   "Vein.endothelial" = "Endothelial cells",
#   "Fibro" = "Fibroblasts",
#   "cDC2b" = "DC2",
#   "iMo" = "Monocytes",
#   "SPP1..LAM" = "Macrophages",
#   "Epithelial" = "Epithelial cells",
#   "Lymphatic.endothelial" = "Endothelial cells",
#   "Artery.endothelial" = "Endothelial cells",
#   "Capillary.endothelial" = "Endothelial cells",
#   "Mature.NK" = "NK cells",
#   "trNK" = "NK cells"
# )

# # PAAD
# celltypist_map <- c(
#   "PVM" = "Macrophages",
#   "BCell" = "B cells",
#   "Plasma" = "Plasma cells",
#   "nCD4" = "Tcm/Naive helper T cells",
#   "Activated.CD8..T.cell" = "Tem/Temra cytotoxic T cells",
#   "Treg" = "Regulatory T cells",
#   "nCD8" = "Tcm/Naive cytotoxic T cells",
#   "LAM" = "Macrophages",
#   "IFN.primed.macrophage" = "Macrophages",
#   "cMo" = "Classical monocytes",
#   "Interstitial.macrophage" = "Macrophages",
#   "Mast.cell" = "Mast cells",
#   "Vein.endothelial" = "Endothelial cells",
#   "Fibro" = "Fibroblasts",
#   "cDC2b" = "DC2",
#   "Neutrophil" = "Neutrophils",
#   "Artery.endothelial" = "Endothelial cells",
#   "Capillary.endothelial" = "Endothelial cells",
#   "Erythrocyte" = "Erythrocytes",
#   "Mature.NK" = "NK cells",
#   "trNK" = "NK cells"
# )

# ---- 4. Add Module Scores & Calculate Scalar AUPRCs ----
celltypist_map <- gsub("\\W", ".", celltypist_map)

auprc_scores <- numeric(0)
plot_list <- list()

for (cell in names(celltypist_map)) {
  message(paste0('Analyzing: ', cell))
  
  # Load gene module - take top 20 genes present in scRNA-seq
  celltypist_name <- celltypist_map[[cell]]
  gene_list <- celltypist_modules[, celltypist_name]
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
      subtitle = paste0("AUPRC = ", round(auc_val, 3), "\n(CellTypist: ", celltypist_name, ")"),
      x = "Recall",
      y = "Precision"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 8, color = "gray30"),
      axis.title = element_text(size = 9)
    )
}

message("CellTypist Average AUPRC: ", round(mean(auprc_scores, na.rm = TRUE), 3))

pr_grid <- wrap_plots(plot_list, ncol = 4) +
  plot_annotation(
    title = paste0("CellTypist Lineage Markers Precision-Recall Curves: ", cancer),
    subtitle = paste0("Average AUPRC: ", round(mean(auprc_scores, na.rm = TRUE), 3)),
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 13, color = "gray20")
    )
  )

pr_grid

# ---- 5. Save Grid to PDF ----
ggsave(paste0("CellTypist_All_PR_Curves_Grid_", cancer, ".pdf"), plot = pr_grid, width = 20, height = 24)

# ---- 6. Save Raw Numbers ----
write.csv(file = paste0('CellTypist_', cancer, '_scores.csv'), auprc_scores)

