# ---- 1. Libraries and Functions ----
library(tidyverse)
library(harmony)
library(Seurat)
library(nichenetr)
library(CellChat)

preprocess_seurat <- function(rds_path) {
  readRDS(rds_path) |>
    NormalizeData() |>
    FindVariableFeatures() |>
    ScaleData() |>
    RunPCA(verbose = FALSE) |>
    RunHarmony(group.by.vars = "Patient", verbose = FALSE) |>
    FindNeighbors(dims = 1:10, reduction = "harmony", verbose = FALSE) |>
    FindClusters(resolution = 0.5, verbose = FALSE) |>
    RunUMAP(dims = 1:10, reduction = "harmony", verbose = FALSE)
}

# ---- 2. NicheNet Pipeline ----
score_receiver_ligands <- function(receiver_cell, rds, lr_network, global_expressed_genes_sender, 
                                   ligand_target_matrix, weighted_networks, global_sender_mapping) {
  Idents(rds) <- "celltype"
  expressed_genes_receiver <- get_expressed_genes(receiver_cell, rds, assay_oi = "RNA", pct = 0.05)
  
  expressed_ligands <- intersect(unique(lr_network$from), global_expressed_genes_sender)
  expressed_receptors <- intersect(unique(lr_network$to), expressed_genes_receiver)
  background_expressed_genes <- intersect(expressed_genes_receiver, rownames(ligand_target_matrix))
  
  receiver.rds <- subset(rds, idents = receiver_cell)
  Idents(receiver.rds) <- "condition"
  
  DE_table <- FindMarkers(receiver.rds, ident.1 = "Tumor", ident.2 = "Adjacent", min.pct = 0, verbose = FALSE)
  
  tumor_genes <- DE_table |> filter(p_val_adj < 0.05 & avg_log2FC > 0.25) |> rownames() |> intersect(rownames(ligand_target_matrix))
  adjacent_genes <- DE_table |> filter(p_val_adj < 0.05 & avg_log2FC < -0.25) |> rownames() |> intersect(rownames(ligand_target_matrix))
  
  potential_ligands <- lr_network |> filter(from %in% expressed_ligands & to %in% expressed_receptors) |> pull(from) |> unique()
  if (length(potential_ligands) == 0) return(NULL)
  
  run_nichenet_for_state <- function(geneset_oi, state_name) {
    if (length(geneset_oi) == 0) return(NULL)
    
    predict_ligand_activities(geneset_oi, background_expressed_genes, ligand_target_matrix, potential_ligands) |>
      select(test_ligand, aupr_corrected) |>
      inner_join(
        weighted_networks |> 
          filter(from %in% potential_ligands & to %in% expressed_receptors) |> 
          semi_join(lr_network, by = c("from", "to")),
        by = c("test_ligand" = "from")
      ) |>
      inner_join(global_sender_mapping, by = c("test_ligand" = "ligand")) |>
      mutate(driven_state = state_name, receiver = receiver_cell) |>
      select(driven_state, sender, receiver, ligand = test_ligand, receptor = to, ligand_aupr = aupr_corrected, prior_weight = weight)
  }
  
  bind_rows(
    run_nichenet_for_state(tumor_genes, "Tumor"),
    run_nichenet_for_state(adjacent_genes, "Adjacent")
  )
}

run_NicheNet <- function(rds_obj, cancer_name, lr_network, ligand_target_matrix, weighted_networks, zscore_threshold = 2.5) {
  Idents(rds_obj) <- "celltype"
  all_celltypes <- unique(rds_obj$celltype)
  global_expressed_genes_sender <- get_expressed_genes(levels(rds_obj), rds_obj, assay_oi = "RNA", pct = 0.05) |> unlist() |> unique()
  
  potential_global_ligands <- intersect(global_expressed_genes_sender, unique(lr_network$from))
  
  expr_matrix <- GetAssayData(rds_obj, assay = "RNA", layer = "data")[potential_global_ligands, , drop = FALSE]
  pct_expr_matrix <- sapply(all_celltypes, function(ct) {
    cells <- WhichCells(rds_obj, idents = ct)
    rowMeans(expr_matrix[, cells, drop = FALSE] > 0)
  })
  
  global_sender_mapping <- data.frame(
    ligand = potential_global_ligands,
    sender = colnames(pct_expr_matrix)[apply(pct_expr_matrix, 1, which.max)]
  )
  
  z_scored_interactions <- lapply(all_celltypes, function(ct) {
    score_receiver_ligands(ct, rds_obj, lr_network, global_expressed_genes_sender, 
                           ligand_target_matrix, weighted_networks, global_sender_mapping)
  }) |> 
    bind_rows() |>
    group_by(receiver, driven_state) |> 
    mutate(aupr_zscore = as.numeric(scale(ligand_aupr))) |>
    ungroup() |>
    filter(aupr_zscore >= zscore_threshold) |> 
    arrange(driven_state, desc(aupr_zscore))
  
  write.csv(z_scored_interactions, file = paste0(cancer_name, "_NicheNet_Top_Interactions.csv"), row.names = FALSE)
  return(z_scored_interactions)
}

# ---- 3. CellChat Pipeline ----
run_CellChat <- function(rds_obj, cancer_name, min_log2FC = log2(100), min_base_prob = 0.10) {
  seurat_list <- SplitObject(rds_obj, split.by = "condition")
  
  process_cellchat <- function(seurat_sub) {
    meta <- data.frame(celltype = droplevels(as.factor(seurat_sub$celltype)), row.names = colnames(seurat_sub))
    
    createCellChat(object = GetAssayData(seurat_sub, assay = "RNA", layer = "data"), meta = meta, group.by = "celltype") |>
      (`@<-`)("DB", CellChatDB.human) |>
      subsetData() |>
      identifyOverExpressedGenes() |>
      identifyOverExpressedInteractions() |>
      computeCommunProb(type = "triMean", raw.use = TRUE) |>
      filterCommunication(min.cells = 10) |>
      computeCommunProbPathway() |>
      aggregateNet()
  }
  
  cellchat.tumor <- process_cellchat(seurat_list[["Tumor"]])
  cellchat.adjacent <- process_cellchat(seurat_list[["Adjacent"]])
  
  df_diff <- bind_rows(
    Adjacent = subsetCommunication(cellchat.adjacent, thresh = 1),
    Tumor    = subsetCommunication(cellchat.tumor, thresh = 1),
    .id      = "dataset"
  ) |>
    select(dataset, source, target, interaction_name_2, pathway_name, prob, pval) |>
    pivot_wider(
      names_from  = dataset,
      values_from = c(prob, pval),
      values_fill = list(prob = 0, pval = 1)
    ) |>
    filter(pval_Adjacent < 0.05 | pval_Tumor < 0.05) |>
    mutate(
      log2FC     = log2((prob_Tumor + 1e-5) / (prob_Adjacent + 1e-5)),
      abs_log2FC = abs(log2FC),
      regulation = ifelse(log2FC > 0, "Higher in Tumor", "Lower in Tumor"),
      ligand     = str_trim(sapply(str_split(interaction_name_2, " - "), `[`, 1)),
      receptor   = str_trim(gsub("[()]", "", sapply(str_split(interaction_name_2, " - "), `[`, 2)))
    )
  
  top_diff_interactions <- df_diff |>
    filter(source != target, ligand != receptor) |>
    filter(prob_Tumor >= min_base_prob | prob_Adjacent >= min_base_prob) |>
    filter(abs_log2FC >= min_log2FC) |>
    group_by(target, ligand, receptor) |>
    slice_max(order_by = abs_log2FC, n = 1, with_ties = FALSE) |>
    ungroup() |>
    arrange(desc(abs_log2FC)) |>
    select(source, target, interaction_name_2, ligand, receptor, pathway_name, regulation, log2FC, prob_Tumor, prob_Adjacent, pval_Tumor, pval_Adjacent)
  
  write.csv(top_diff_interactions, file = paste0(cancer_name, "_CellChat_Top_Interactions.csv"), row.names = FALSE)
  return(top_diff_interactions)
}

# ---- 4. Global Data & Execution ----
lr_network <- readRDS('lr_network.rds') |>
  distinct(from, to, .keep_all = TRUE) |>
  filter(!database %in% c("ppi_prediction_go", "ppi_prediction"))

ligand_target_matrix <- readRDS('ligand_target_matrix.rds')
weighted_networks <- readRDS('weighted_networks.rds')$lr_sig |> 
  inner_join(lr_network, by = c("from", "to"))

cancers_list <- c('HNSC', 'CESC', 'COAD', 'GBM', 'KIRC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', 'OV', 'PAAD')

for (cancer in cancers_list) {
  message("\n===> Processing: ", cancer)
  
  seurat_obj <- preprocess_seurat(paste0(cancer, '_final.RDS'))
  
  run_NicheNet(seurat_obj, cancer, lr_network, ligand_target_matrix, weighted_networks)
  run_CellChat(seurat_obj, cancer)
  
  rm(seurat_obj)
  gc()
}

message("\nAll cohorts processed successfully!")