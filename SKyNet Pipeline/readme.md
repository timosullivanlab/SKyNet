Our full survival analysis has been split into six annotated .Rmd files.

1. Import TCGA.Rmd

This file extracts the TCGA data (clinical, counts, and gene metadata) for each cancer. The gene metadata is only used to convert ENSEMBL gene IDs to gene names; as such, BLCA gene metadata will be used as the example later in the pipeline. 

2. scRNA Expression.Rmd

This file extracts lineage markers, percent expression by cell type, and average expression by cell type from a scRNA-seq object. For cell types found across cancer types, the code references the merged atlas to pick up more generalizable cell type lineage markers; otherwise, it will identify the top 20 genes ranked by positive log2FC.

The file also extracts expression values for all ligands, receptors, and receptor complexes found in scDiffCom. For receptor complexes, all components must be present for the cell to express the receptor complex (in which case expression is the average of individual gene expression).

3. Deconvolution.Rmd

This file conducts single and double lineage deconvolution of cell types. Based on the lineage markers derived from scRNA Expression.Rmd and TCGA data from Import TCGA.Rmd, the code evaluates the degree to which each cell is associated with patient survival outcomes.

4. SKyNet.Rmd

This file runs the communication analysis network. Built upon the framework from Deconvolution.Rmd, this code evaluates the degree to which each ligand and receptor within a given cell type is associated with patient survival outcomes. It filters through the list of ligand-receptor interactions from scDiffCom and infers the most likely sender and receiver cell types based on survival implications, resulting in a list of the ligand-receptor interactions associated with significant survival outcomes.

5. Visualization.Rmd

This file generates summary plots for the interactions for each cancer. It can also zoom in on specific cell types and the interactionsthat cell type are involved in.