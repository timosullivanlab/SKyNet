library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(tidyverse)
library(rhdf5)
library(harmony)

#function -----
LoadUnfiltered_10x <- function(data.dir, project){
  rawdata <- Read10X(data.dir)
  rawdata_object = CreateSeuratObject(rawdata, project = project, min.cells = 3)
  rawdata_object[["percent.mt"]] <- PercentageFeatureSet(rawdata_object, pattern = "^MT-")
  rawdata_object[['percent.rp']] = PercentageFeatureSet(rawdata_object, pattern = "^RP[SL]")
  rawdata_object$Type = 'SeuratObject'
  return(rawdata_object)

}


SeuratAnalysis1 <- function(SeuratObject){
  SeuratObject <- JoinLayers(SeuratObject) #can't run this after converting sce back to seurat 
  SeuratObject <- NormalizeData(SeuratObject, normalization.method = "LogNormalize", scale.factor = 10000)
  SeuratObject <- FindVariableFeatures(object = SeuratObject, selection.method = 'vst',nfeatures = 5000) #choose 5k features?
  SeuratObject <- ScaleData(SeuratObject)
  SeuratObject <- RunPCA(SeuratObject, features = VariableFeatures({SeuratObject}))
}

SeuratAnalysis2 <- function(SeuratObject, harmony, harmony_by, n){
  if (harmony){
    SeuratObject <- RunHarmony(SeuratObject, harmony_by)
    SeuratObject <- FindNeighbors(SeuratObject, dims = 1:n, reduction = "harmony")
    SeuratObject <- FindClusters(SeuratObject, resolution = 0.8)
    SeuratObject <- RunUMAP(SeuratObject, reduction = "harmony", label = TRUE, dims = 1:n)
  }
  else {
    SeuratObject <- FindNeighbors(SeuratObject, dims = 1:n)
    SeuratObject <- FindClusters(SeuratObject, resolution = 0.8)
    SeuratObject <- RunUMAP(SeuratObject, reduction = "pca", label = TRUE, dims = 1:n)
  }
}


#breast ----- 
#load raw counts -- a list of raw H5 files 
untar("GSE195665_RAW.tar", exdir = '~/jliang/Breast_H5_ready')
h5.fileList <- list.files("~/jliang/Human tissue harmonization/Breast_H5_rawdata", pattern = "\\.h5$")
load_h5_file <- function(h5_file) {
  raw_matrix = Read10X_h5(h5_file) 
  seu = CreateSeuratObject(counts = raw_matrix, min.cells = 3,project='Breast')
  return(seu)
}
all_h5_data <- lapply(c(h5.fileList[105:112], h5.fileList[115]), load_h5_file) #RETURN A LIST
names(all_h5_data) <-c('c01','c02','c03','c04','c05','c08','c11') #basename(h5.fileList[105:109]) #NAME EACH COMPONENT OF THE LIST 

p1= all_h5_data[[1]] #737280 cells
p1[["percent.mt"]] <- PercentageFeatureSet(p1, pattern = "^MT-")
p1$percent.mt[is.na(p1$percent.mt)] <- 0 #most have NaN or 0 values, so force to 0 
p1$Type = 'SeuratObject'
plot = VlnPlot(p1,  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, raster = F, group.by = 'Type')
p3 = all_h5_data[[3]] #737280 cells
p3[["percent.mt"]] <- PercentageFeatureSet(p3, pattern = "^MT-")
p3$percent.mt[is.na(p3$percent.mt)] <- 0 #most have NaN or 0 values, so force to 0 
 
p3$Type = 'SeuratObject'

p2$nFeature_RNA %>% fivenum %>% unname
p2$nCount_RNA %>% fivenum %>% unname
p2$percent.mt %>% fivenum %>% unname


QC= function(rawdata_object){
  rawdata_object[["percent.mt"]] <- PercentageFeatureSet(rawdata_object, pattern = "^MT-")
  rawdata_object$percent.mt[is.na(rawdata_object$percent.mt)] <- 0 
  rawdata_object <- subset(rawdata_object, subset = nFeature_RNA <= 5000 & nCount_RNA > 100 &
                             nCount_RNA <= 25000 &
                             percent.mt <= 15)
  return(rawdata_object)
}

QC_h5_data <- lapply(all_h5_data, QC)

for (i in 1:5) {
  QC_h5_data[[i]]$Patient <- names(QC_h5_data)[i]
}
 
Breast_all = Reduce(function(x, y) merge(x, y), QC_h5_data) 

#merging ---- 
kidney_normal = readRDS("normalkidney_QC.rds")
kidney_cancer = readRDS("ccRCC_QC.rds")
hnscc_normal = readRDS("hnscc_normal_QC.rds")
hnscc_tumor = readRDS("hnscc_tumor_QC.rds")
lusc_normal = readRDS("LUSC_normal_QC.rds")
lusc_tumor = readRDS("LUSC_tumor_QC.rds")
luad= readRDS("LUAD_Combined_Updated_1020.rds")
HCC = readRDS("FinalRDS.rds")

#uniform metadata
Breast_all$Patient = paste0("Breast", "_", Breast_all$Patient) 
Breast_all$condition = 'Tumor'
unique(kidney_cancer$donor_id)
#ss 2022, ss 2017 is bening adjacent ; ss_2017
kidney_cancer$Patient = paste0("Kidney", "_", kidney_cancer$donor_id) 
kidney_normal$Patient = paste0("Kidney", "_", kidney_normal$donor_id) 
kidney_cancer$orig.ident = 'Kidney'
kidney_normal$orig.ident = 'Kidney'
kidney_cancer$condition = 'Tumor'
kidney_normal$condition = 'Adjacent'

hnscc_normal$Patient = paste0("hnscc", "_", hnscc_normal$patient.id) 
hnscc_tumor$Patient = paste0("hnscc", "_", hnscc_tumor$patient.id)  
hnscc_tumor$condition = 'Tumor'
hnscc_normal$condition = 'Adjacent'

lusc_normal$condition = 'Adjacent'
lusc_normal$Patient = paste0("LUSC", "_", lusc_normal$Patient) 
lusc_normal$orig.ident = 'LUSC'
lusc_tumor$condition = 'Tumor'
lusc_tumor$Patient = paste0("LUSC", "_", lusc_tumor$Patient) 
lusc_tumor$orig.ident = 'LUSC'

HCC$condition = HCC$site 
HCC$condition =gsub('Normal','Adjacent', HCC$condition)
HCC$orig.ident ='HCC'
HCC$Patient = paste0("HCC", "_", HCC$patient)

luad = JoinLayers(luad)
luad$condition = luad$orig.ident
luad$condition =gsub('normal','Adjacent', luad$condition)
luad$condition =gsub('LUAD','Tumor', luad$condition)
luad$orig.ident = 'LUAD'
luad$Patient = paste0("LUAD", "_", luad$Patient)

breast_metadata = Breast_all@meta.data[, c("orig.ident",'Patient','condition','nCount_RNA','nFeature_RNA', 'percent.mt')]
kidney_normal_metadata = kidney_normal@meta.data[, c("orig.ident",'Patient', 'condition','nCount_RNA','nFeature_RNA', 'percent.mt')]
kidney_cancer_metadata = kidney_cancer@meta.data[, c("orig.ident",'Patient', 'condition','nCount_RNA','nFeature_RNA', 'percent.mt')]
hnscc_normal_metadata = hnscc_normal@meta.data[,  c("orig.ident",'Patient', 'condition','nCount_RNA','nFeature_RNA', 'percent.mt')]
hnscc_tumor_metadata = hnscc_tumor@meta.data[,  c("orig.ident",'Patient', 'condition','nCount_RNA','nFeature_RNA', 'percent.mt')]
LUSC_normal_metadata = lusc_normal@meta.data[,  c("orig.ident",'Patient', 'condition','nCount_RNA','nFeature_RNA', 'percent.mt')]
LUSC_tumor_metadata = lusc_tumor@meta.data[,  c("orig.ident",'Patient', 'condition','nCount_RNA','nFeature_RNA', 'percent.mt')]
HCC_metadata = HCC@meta.data[,  c("orig.ident",'Patient', 'condition','nCount_RNA','nFeature_RNA', 'percent.mt')]
luad_metadata = luad@meta.data[,  c("orig.ident",'Patient', 'condition','nCount_RNA','nFeature_RNA', 'percent.mt')]


total_metadata = rbind(breast_metadata, kidney_cancer_metadata, kidney_normal_metadata,
                       hnscc_normal_metadata, hnscc_tumor_metadata,
                       LUSC_normal_metadata, LUSC_tumor_metadata, luad_metadata,
                       HCC_metadata)
 
breast_export = CreateSeuratObject(counts = LayerData(Breast_all, assay = "RNA", layer = "counts"), meta.data= breast_metadata)
lusc_normal_export = CreateSeuratObject(counts = LayerData(lusc_normal, assay = "RNA", layer = "counts"), meta.data= LUSC_normal_metadata)
lusc_tumor_export = CreateSeuratObject(counts = LayerData(lusc_tumor, assay = "RNA", layer = "counts"), meta.data= LUSC_tumor_metadata)
HCC_export = CreateSeuratObject(counts = LayerData(HCC, assay = "RNA", layer = "counts"), meta.data= HCC_metadata)
luad_export = CreateSeuratObject(counts = LayerData(luad, assay = "RNA", layer = "counts"), meta.data= luad_metadata)

total = merge(x = breast_export, y = c(kidney_cancer, kidney_normal, hnscc_normal, hnscc_tumor,
                                       lusc_normal_export, lusc_tumor_export, luad_export, HCC_export
                                          ))
Layers(total)
total@meta.data = total_metadata

total <- NormalizeData(total)
total <- FindVariableFeatures(total)  
total <- ScaleData(total)
total <- RunPCA(total)

total <- IntegrateLayers(
  object = total, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

ElbowPlot(total, ndims = 50)

total <- FindNeighbors(total, reduction = "harmony", dims = 1:20)

total <- FindClusters(total, resolution = 0.5, cluster.name = "harmony_clusters") #lower resolution (default) to make it faster 

total <- RunUMAP(total, reduction = "harmony", dims = 1:20, reduction.name = "umap.harmony")


p1= DimPlot(
  total,
  reduction = "umap.harmony",
  group.by = c("orig.ident", 'harmony_clusters'), raster = FALSE, label = TRUE, combine = FALSE
)

ggsave('integrate_harmony_clusters.pdf', plot =p1, width = 20, height = 20)

DimPlot(
  total,
  reduction = "umap.harmony",
  group.by = 'seurat_clusters', raster = FALSE, label = TRUE, combine = FALSE
)

 
total = JoinLayers(total)

#find general clusters: lymphoid, myeloid, structural ---------------
T_NK_cell_genes = c('CD3D','CD3E','CD40LG','CD8A','CD8B','GZMB','NKG7','KLRF1','PRF1','CTSW')
B_cell_genes =  c("CD19",'CD79A','BANK1','MS4A1','JCHAIN','IGHA1','IGHG1')
Th1_genes = c('CD3D','CD3E','CD40LG','IFNG', 'TBX21','TNF','RUNX3','ANXA3','FASLG')
Tem_genes = c('KLRG1', 'KLF2', 'GZMA', 'GZMB', 'GZMH', 'PRDM1','TCF7',
              'TNF', 'RORA', 'ID2', 'ZEB2')
ILC1_genes = c('AREG', 'SOX4', 'KIT', 'TNFRSF18', 'TBX21', 
               'ZNF683', 'CXCR6', 'AQP3', 'GEM', 'CEBPB', 'S100A10', 'RBPJ'
)
Macro_mono_genes = c('CD14', 'CD163', 'CD68', 'MARCO', 'MRC1',
                     'C1QA', 'C1QB','IL1B', 'FCN1', 'CXCL8')
DC_genes = c("CLEC10A", 'CD1C','ITGAM','XCR1','CLEC9A')
Epi_genes = 'EPCAM'
Endo_genes = c( "PECAM1", "NOTCH4", "TEK",'VWF','PLVAP', 'RAMP3' )
fibro_genes = c('LUM', 'DCN', 'PI16','PDGFRA','COL1A1')
SMC_genes = c("TPM1",'TPM2','TAGLN','PDGFRB','RGS5')

Idents(total)= 'harmony_clusters'

p3 = FeaturePlot(total, 
                 features = B_cell_genes,
                 raster=F, label = T,  cols = c('grey','red'))
ggsave('FeaturePlot_BCell_Plasma.pdf', plot =p3, width = 20, height = 20)
#6,part 7,22; 21=plasma 

p5 = FeaturePlot(total, 
                 features =T_NK_cell_genes, 
                 raster=F, label = T,  cols = c('grey','red'))

#1,2,3,20,10,13,14; part 9  
ggsave('FeaturePlot_T_NK_general.pdf', plot =p5, width = 20, height = 20)


p6 = FeaturePlot(total, 
                 features =Macro_mono_genes, 
                 raster=F, label = T,  cols = c('grey','red'))

# 
ggsave('FeaturePlot_macrophage_monocyte_general.pdf', plot =p6, width = 20, height = 20)

p7 = FeaturePlot(total, 
                 features =DC_genes, 
                 raster=F, label = T,  cols = c('grey','red'))

#tip of 12 = CDC1; some 4 and 12 and 7 = cdc2 
ggsave('FeaturePlot_DC_general.pdf', plot =p7, width = 20, height = 20)

p8 = FeaturePlot(total, 
                 features = 'EPCAM', 
                 raster=F, label = T,  cols = c('grey','red'))

#part 18 other than CD14 expression + 23,11,16,34,33 + a little mixed in 31,37,38
ggsave('FeaturePlot_epi_general.pdf', plot =p8, width = 20, height = 20)

p9 = FeaturePlot(total, 
                 features = Endo_genes, 
                 raster=F, label = T,  cols = c('grey','red'))

#
ggsave('FeaturePlot_endo_general.pdf', plot =p9, width = 20, height = 20)

p10 = FeaturePlot(total, 
                 features = fibro_genes, 
                 raster=F, label = T,  cols = c('grey','red'))

#
ggsave('FeaturePlot_fibro_general.pdf', plot =p10, width = 20, height = 20)

p11 = FeaturePlot(total, 
                  features = SMC_genes, 
                  raster=F, label = T,  cols = c('grey','red'))

#
ggsave('FeaturePlot_SMC_general.pdf', plot =p11, width = 20, height = 20)

p12= DotPlot(total, group.by = 'harmony_clusters', features = c(T_NK_cell_genes, B_cell_genes,
                                                                Macro_mono_genes, DC_genes, Endo_genes, SMC_genes, fibro_genes, 'EPCAM')) + coord_flip()+
  theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8))

ggsave('DotPlot_all_markers_harmony_clusters.pdf', plot =p12, width = 20, height = 20)

#subcluster level 1 ----- 
Lymphoid = subset(total, idents = c(6,7,22,21,1,2,3,20,10,13,14,9 ))
Lymphoid <- NormalizeData(Lymphoid)
Lymphoid <- FindVariableFeatures(Lymphoid)  
Lymphoid <- ScaleData(Lymphoid)
Lymphoid <- RunPCA(Lymphoid)
ElbowPlot(Lymphoid, ndims = 50)

Lymphoid = RunHarmony(Lymphoid, 'orig.ident')
Lymphoid <- FindNeighbors(Lymphoid, reduction = "harmony", dims = 1:20)
Lymphoid <- FindClusters(Lymphoid, resolution = 1)
Lymphoid <- RunUMAP(Lymphoid, reduction = "harmony", dims = 1:20)

DimPlot(Lymphoid,
  reduction = "umap",label = T, 
  raster = FALSE,
  group.by = 'seurat_clusters' 
)

saveRDS(Lymphoid, file = 'Lymphoid_raw.rds')

p13=FeaturePlot(Lymphoid, 
            features =B_cell_genes, 
            raster=F, label = T, repel = T, cols = c('grey','red'))

ggsave('Lymphoid_BCell_general.pdf', plot =p13, width = 20, height = 20)
#BCell = 0,20; plasma = 19,34,18,9
#20 has some CD3

p14= FeaturePlot(Lymphoid, 
            features =T_NK_cell_genes, 
            raster=F, label = T, repel = T, cols = c('grey','red'))
ggsave('Lymphoid_T_NK_general.pdf', plot =p14, width = 20, height = 20)

#based on feature plots but may not be precise: 
#T cells =6,3,31,15,27,40,20,22,2,41,7,23,24,1,11,12,5,4,8,39,25
#cd4+: 27,40,31,3; dotplot 2,3,15
#cd8+: 2,22,23,24,1,8; dotplot 1,8,11,17,39
#NK= 32,38,5,29 (need to subcluster to distinguish from T)

DotPlot(Lymphoid, features = c('CD3D','CD3E','CD40LG','IFNG', 'TBX21','TNF',
                               'CD8A','CD8B','SELL','CCR7',
                               'KLRG1', 'KLF2', 'GZMA', 'GZMB', 'GZMH','TCF7',
                               'ZNF683', 'TRDC', 'TRDV1', 'TRGC1','KIR2DL4'), cols = 'RdYlBu')+
  coord_flip() + theme(axis.text.x = element_text(size = 8),
                       axis.text.y = element_text(size = 8))

#3: uniquely high TCF7; 
#5,11,12,32,33 = Th1, TBX21+ (can be NK too)
#gdt= 4,5,12,32,38: seems 4,12 the most confident 
#not T cell: 0,5,9,10,14,16,18,19,20,21,25,26,28,29,20,32,33,34,35,37,38,
p15= FeaturePlot(Lymphoid, 
                 features = c('CD40LG','IFNG', 'TBX21','TNF',
                              'CD8A','CD8B','SELL','CCR7',
                              'KLRG1', 'KLF2', 'GZMA', 'GZMB', 'GZMH','TCF7',
                              'ZNF683', 'TRDC', 'TRDV1', 'TRGC1','KIR2DL4'), 
                 raster=F, label = T, repel = T, cols = c('grey','red'))
ggsave('Lymphoid_specific_1.pdf', plot =p15, width = 20, height = 20)

#B and NK specific genes: 
#NK: 4 has high xcl1 and cd247 
#5, 12,32,33 --> 4,12 is ncam1+fcgr3a- cd56brightcd16dim immature 
# 5,32,33 maybe cd56dim and cd16bright, but need to verify 
#check 4, verify immature or mature 
#so subcluster 5,12,32,33 

DotPlot(Lymphoid, features = c( #"CD19",'CD79A','BANK1','MS4A1','JCHAIN',
                                'KLRF1','PRF1','GZMB','CXCR6','CD160','XCL1','CD247', #MATURE 
                                'NCAM1','FCGR3A','ICAM1','KLRK1','NCR1', #IMMATURE 
                                'KLRC2','CD2' #ADAPTIVE
                                ), cols = 'RdYlBu')+
  coord_flip() + theme(axis.text.x = element_text(size = 8),
                       axis.text.y = element_text(size = 8))

#how I subcluster inside Lymphoid: ------
#NK: 4,5,12,32,33 
#CD4+: 2,3,15 
#CD8+: 1,8,11,17,39
#foxp3: 6
#bcell: 0,30,19,34,18,9


#subcluster NK:------
NK = subset(Lymphoid, idents = c(4,5,12,32,33 ))
NK_final = SeuratAnalysis1(NK_final)
ElbowPlot(NK_final, ndims = 50)
NK_final = SeuratAnalysis2(NK_final, harmony = T, harmony_by = 'orig.ident', n = 20)

DimPlot(NK_final,
        reduction = "umap",label = T, 
        raster = FALSE, group.by = 'celltype'
)


DotPlot(NK_final, features = c(,
                               'IFNG', 'TBX21','TNF',
                               'KLRG1', 'KLF2','TCF7','RORA','RORC','CCR7',#Tem
                               'SELL', 'CD44', 'CARS1', 'FGFBP2', #Tcm 
                               'ITGA1','ITGAE','CD69','GZMK','EOMES', #Trm
                               'ZNF683', 'KLRC2', 'TRDC', 'TRDV1',  'TRGC1', 'TRGV9', #GDT 
                               'KLRF1','PRF1',
                               'GZMB','CD160','XCL1','CD247', 
                               'NCAM1','FCGR3A', #IMMATURE 
                               'ADGRG1','FCER1G','LAIR2','SPON2' #SPECIFIC MATURE NK GENES 
                               
), group.by ='celltype', cols = 'RdYlBu')+
  coord_flip() + theme(axis.text.x = element_text(size = 8),
                       axis.text.y = element_text(size = 8))

VlnPlot(NK_final, features = c('CD40LG','IFNG', 'TBX21','TNF',#TH1 
                              'GZMH','PRDM1','KLRG1', 'KLF2',#Tem
                              'SELL', 'TCF7','CCR7', #Tcm 
                              'ITGA1','ITGAE','TOB2'), group.by = 'celltype') #TRM

DotPlot(NK_final, features = c('CD40LG', 'CD8A','CD8B','CCR7', 'SELL', 'TCF7',  #tcm 
                          'CXCR6','ITGAE','ITGA1', #trm 
                          'RORA','RORC','CX3CR1','TNF','PRDM1','GZMK')) #tem 

DotPlot(NK, features = c('CD3D','CD3E',
  'KLRF1','PRF1','GZMB','CXCR6','CD160','XCL1','CD247', #MATURE 
  'NCAM1','FCGR3A','ICAM1','KLRK1','NCR1', #IMMATURE 
  'KLRC2','CD2' #ADAPTIVE
), cols = 'RdYlBu')+
  coord_flip() + theme(axis.text.x = element_text(size = 8),
                       axis.text.y = element_text(size = 8))
#NK:0,11, 5, all mature ;  
#9 is trNK; 

#look into mature genes
DotPlot(NK, features = c('NCAM1','FCGR3A','CD160','CD247',
                          'ADGRG1','FCER1G','LAIR2','SPON2',
                         'CLIC3','CHST2','IGFBP7'), cols = 'RdYlBu')
cluster2 = FindMarkers(NK, ident.1= 2, min.pct = 0.25)
#6 is a lot of MT genes, subset out 

#NK = subset(NK, idents = c(0,1,2,3,4,5,6,7,8,10,12,13,14,15))
#NK = subset(NK, idents = 0:14)
#NK = subset(NK, idents = c(0,1,2,3,5,6,7,8,9,10,11,12,13,14))
NK = subset(NK, idents = c(0,1,2,3,5,6,8,9,10,11,12,13))

#check T cells: 1,2,3,4,6,7,10,12,13 
  
DotPlot(NK, features = c('CD40LG','IFNG', 'TBX21','TNF',
                         'CD8A','CD8B', 
                         'KLRG1', 'KLF2', 'GZMA', 'GZMB', 'GZMH','TCF7',
                         'ZNF683', 'TRDC', 'TRDV1', 'TRGC1','KIR2DL4','TRGV4','KLRC1'), cols = 'RdYlBu') 
#cd8+: 2,3,10 
#cd4+: 13,12; 13 = th1 
#gdt = DEF 1 , 7
#8 = TCM 
 
cluster.ids = c("Mature NK",'gdT','CytoCD8','CytoCD8','Unknown',
                'Mature NK','Trm','gdT','Tcm','trNK','CytoCD8',
                'Mature NK','nCD4','Th1')
names(cluster.ids) <- levels(NK)
NK <- RenameIdents(NK, cluster.ids)
DimPlot(NK, reduction = "umap", label = TRUE, repel = TRUE)  

NK$celltype = Idents(NK)
NK@meta.data$celltype %>% table



VlnPlot(NK, Th1_genes)
VlnPlot(NK, features = c('GZMH', 'KLRG1', 'KLF2', 'GZMA', 
                         'GZMB', 'CX3CR1', 'PRDM1', 'TNF',
                         'RORC', 'RORA', 'ID2', 'ZEB2', 'TCF7'))  #Tem genes
VlnPlot(NK, features = c('KLRB1','IFNGR1','SLC4A10','TRAV1-2','TRBV6-2',
                         'IL17A','IL22','IL21','CD40LG')) #MAIT & Th17 gene
VlnPlot(NK, features=c('CD3D', 'CD3E','CCR7', 'SELL', 'CD44', 'TCF7',
                       'CARS1', 'FGFBP2')) #Tcm 
VlnPlot(NK, features = c( 'KLRF1','PRF1','GZMB','CXCR6','CD160','XCL1','XCL2','CD247',
                          'IRF8','TOX2')) #TRNK 
VlnPlot(NK, features = c('CD3D','CD3E','CD8A','CD8B',
                         'ITGA1','ITGAE','CD69','GZMK','EOMES'))
#6= TRM
VlnPlot(NK,  features = c('NCAM1','FCGR3A', 
                          'ICAM1', 'KLRK1' ,'NCR1',
                          'KIT', 'KLRC1', 'IL2RB', 'TCF7', 
                          'GPR183', 'SELL', 'LTB', 'FOS', 'IFITM3', 'COTL1', 'PIK3R1', 'ZFP36L2'
                          ))


VlnPlot(NK, features = c(#'LAG3','PDCD1','HAVCR2',
                         'CD27','CD28','LEF1','CD8A','CD8B','GZMB','GZMA')) #NAIVE, BUT STILL HAVE  HIGH GZM, SO NOT NAIVE 


VlnPlot(NK, features = c('CD3D','CD3E','CD8A','CD8B',
                         'ZNF683', 'KLRC2', 'TRDC', 'TRDV1',
                         'TRGC1', 'TRGV9','TRG-AS1','KIR2DL4', 'KLRD1', 'KLRG1')) 

#subcluster 4 
cluster4 = subset(NK, idents = c(4,7))
cluster4 =SeuratAnalysis1(cluster4)
ElbowPlot(cluster4, ndims = 50)
cluster4= SeuratAnalysis2(cluster4, harmony = T, harmony_by = 'orig.ident', n = 20)

DimPlot(cluster4,
        reduction = "umap",label = T, 
        raster = FALSE
)

DotPlot(cluster4, features = c('CD3D','CD3E',
                         'KLRF1','PRF1','GZMB','CXCR6','CD160','XCL1','CD247', #MATURE 
                         'NCAM1','FCGR3A','ICAM1','KLRK1','NCR1', #IMMATURE 
                         'KLRC2','CD2' #ADAPTIVE
), cols = 'RdYlBu')

DotPlot(cluster4, features = c('CD40LG','IFNG', 'TBX21','TNF',
                         'CD8A','CD8B', 
                         'KLRG1', 'KLF2', 'GZMA', 'GZMB', 'GZMH','TCF7',
                         'ZNF683', 'TRDC', 'TRDV1', 'TRGC1','KIR2DL4','TRGV4','KLRC1'), cols = 'RdYlBu')

#2,3, 5 =cytocd8; 
#0,1 WITH SELL, FGFBP2, KLF2, KLRG1 
#4,6 NK 


VlnPlot(cluster4, features = c('CD8A','CD8B', 'GZMH', 'KLRG1', 'KLF2', 'GZMA', 'GZMB', 'CX3CR1',
                               'PRDM1', 'TNF', 'RORC', 'RORA', 'ID2', 'ZEB2', 'TCF7','CCR7')) #Tem genes 
VlnPlot(cluster4, features = c('CD3D','CD3E','CD40LG','CD8A','CD8B',
                               'CCR7', 'SELL', 'CD44', 'TCF7', 'CARS1', 'FGFBP2'))
VlnPlot(cluster4, features = c('ZNF683', 'KLRC2', 'TRDC', 'TRDV1',
                              'TRGC1', 'TRGV9','TRG-AS1','KIR2DL4', 'KLRD1', 'KLRG1'))

cluster0 = FindMarkers(cluster4, ident.1 = 0, min.pct = 0.25)

#cluster4 = subset(cluster4, idents = c(1,2,3,5))

cluster.ids = c("Tem",'Tem','Unknown','Unknown','Unknown','CytoCD8','trNK')
names(cluster.ids) <- levels(cluster4)
cluster4 <- RenameIdents(cluster4, cluster.ids)
DimPlot(cluster4, reduction = "umap", label = TRUE, repel = TRUE)  

cluster4$celltype = Idents(cluster4)
cluster4@meta.data$celltype %>% table

cluster4= subset(cluster4, idents = c('CytoCD8','trNK', 'Tem'))

#subcluster cd4+ T cell------
CD4 = subset(Lymphoid, idents = c(2,3,15 ))
CD4_final = SeuratAnalysis1(CD4_final)
ElbowPlot(CD4_final, ndims = 50)
CD4_final = SeuratAnalysis2(CD4_final, harmony = T, harmony_by = 'orig.ident', n = 20)
DimPlot(CD4_final,
        reduction = "umap",label = T, group.by = 'condition',
        raster = FALSE
)

#13,14 is not trash 
cluster9 = FindMarkers(CD4, ident.1 = 9, min.pct = 0.25)
VlnPlot(CD4, features = c(#"IGKC",'IGHA1','JCHAIN',
                          'CD3D','CD3E','CD40LG',
                          'CD8A','CD8B',
                          'ISG15','IFIT1','IFIT2', 'IFITM2'))

DotPlot(CD4_final, features = c('CD3D','CD3E','CD40LG', 'SELL', 'CCR7', 'LEF1', #ncd4 
                          'CD8A','CD8B',
                          'CD69','GZMA','PRDM1','CXCR6', #EFFECTOR 
                          'IFNG', 'TBX21','TNF',#Th1  
                          'GATA3', 'STAT6', 'STAT5A', 'IRF4', 'IL4R','IL5','IL13','IL4', #Th2 
                          'ISG15','IFIT1','IFIT2', 'IFITM2', #ifn-primed 
                          "LAG3",'PDCD1','CTLA4','CXCL13', 'HAVCR2','TOX'#exhausted
                          
                          ), group.by='celltype') + coord_flip()+theme(axis.text.x = element_text(size = 8),
                                                                        axis.text.y = element_text(size = 8))
#cd8+:7
#cd4: 0,1,3,4,5,6,8,9,10,13,14
#Th1: 0,1,6,8,9,14  ; yes 0,1,8,9 def effector 
#IFN PRIMED: 11,12 
#2 may have nCD8 
#EXHAUSTED: 4,9; 9 definitely exhausted 


cluster.ids = c('Th1','Th1','Unknown','nCD4', 'Unknown', 'nCD4','Th1',
                'CytoCD8', 'Th1', 'Exhausted  T cell',
                'Unknown','Unknown','IFN-primed',
                'Unknown','Th1')
names(cluster.ids) <- levels(CD4)
CD4 <- RenameIdents(CD4, cluster.ids)
DimPlot(CD4, reduction = "umap", label = TRUE, repel = TRUE)  

CD4$celltype = Idents(CD4)
CD4@meta.data$celltype %>% table
CD4=subset(CD4, idents = c('Th1','nCD4', 'CytoCD8'))



#LOOK INTO GENES: TOX, TOX2, CD27, CD28 , CCR7, TCF7, 

VlnPlot(CD4, features = c('IFNG', 'TBX21','TNF','RUNX3','ANXA3', 'FASLG' ))#th1 
VlnPlot(CD4, features = c("LAG3",'PDCD1','CTLA4','CXCL13',
                          'HAVCR2','TOX','TOX2'))
#exhausted: 4,9,10,13 --> 9 IS REAL EXHAUSTED, OTHERS NEED TO CHECK 


DotPlot(CD4, features = c( 'CD40LG', 'SELL', 'CCR7', 'LEF1',# naive 
                           'CD69','GZMA','PRDM1','CXCR6', #EFFECTOR 
                           'TOX', 'ICOS', 'PDCD1', 'BCL6', 'CXCR5', #TFH
                           'CXCL13', 'SLAMF6','CD38','IL4'),cols = 'RdYlBu')#TFH
                          # 'CD27', 'CD28', 'SELL', 'TCF7', 'LEF1',
                           #'CCR7', 'SELL', 'CD44') )  
                          
#3,4,5,6,13 is the subset with tcf7+cd27+cd28+sell+ccr7+; Tcm?
#SOMEHOW EXHAUSTED ONE IS TFH? ;4 SEEMS LIKELY TO BE  TFH 
#naive: def 4,5,13 = naive, 4,9,12 may have tfh 

VlnPlot(CD4, features = c( 'FOXP3', 'TIGIT', 'ICOS', 'CTLA4', 'TNFRSF9',
                           'TNFRSF4', 'IL2RA', 'IKZF2', 'RTKN2' )) #nCD8 
VlnPlot(CD4, features = c('IL17A', 'IL22', 'IL21', 'MAST4','RORC','STAT3',
                          'GATA3', 'STAT6', 'STAT5A', 'IRF4', 'IL4R','LIMA1','IL5','IL13', 'SPINT2','TRIB3','CST7')) #TH2 



#FIRST GET RID OF 6,9,13 
#CD4= subset(CD4, idents = c(0,1,2,3,4,5,7,8,10,11,12,14,15))
#CD4= subset(CD4, idents = c(1,2,3,4,5,6,7,9,10))
#CD4= subset(CD4, idents = c(0,1,2,3,4,5,6,7,8,9,10,11,12,14,16))
#CD4= subset(CD4, idents = c(0,1,2,3,4,5,6,8,9,10,11,12,13,14,15))


#FURTHER SUBCLUSTER 2,10,11,4,13, 
sub_CD4 = subset(CD4, idents = c(2,10,11,4,13 )) #10 and 11 have useful cells other than trash so subcluster
sub_CD4 = SeuratAnalysis1(sub_CD4)
ElbowPlot(sub_CD4, ndims = 50)
sub_CD4 = SeuratAnalysis2(sub_CD4, harmony = T, harmony_by = 'orig.ident', n = 20)
DimPlot(sub_CD4,
        reduction = "umap",label = T, 
        raster = FALSE
)


cluster5 = FindMarkers(sub_CD4, ident.1 = 5, min.pct = 0.25)
VlnPlot(sub_CD4, features = c("IGKC",'IGHA1','JCHAIN',
  'CD3D','CD3E','CD40LG',
  'CD8A','CD8B',
  'ISG15','IFIT1','IFIT2', 'IFITM2'))

DotPlot(sub_CD4, features = c('CD3D','CD3E','CD40LG','CCR6','CXCR6',#ncd4 
                          'CD8A','CD8B',
                          'IFNG', 'TBX21','TNF',
                          'FOXP3','TNFRSF9' ,
                          'ISG15','IFIT1','IFIT2', 'IFITM2',
                          "LAG3",'PDCD1','CTLA4','CXCL13', 'HAVCR2','TOX','TOX2'), cols = 'RdYlBu')
VlnPlot(sub_CD4, features = c( 'CD69','GZMA','PRDM1','CXCR6',
                               'GATA3', 'STAT6', 'STAT5A', 'IRF4', 'IL4R','LIMA1','IL5','IL13', 'SPINT2','TRIB3','CST7'))
VlnPlot(sub_CD4, features = c("LAG3",'PDCD1','CTLA4','CXCL13',
                          'HAVCR2','TOX','TOX2'))
VlnPlot(sub_CD4, features = c('IFNG', 'TBX21','TNF','RUNX3','ANXA3', 'FASLG' ))#th1 
DotPlot(sub_CD4, features = c( 'CD40LG', 'SELL', 'CCR7', 'LEF1',# naive 
                           'CD69','GZMA','PRDM1','CXCR6', #EFFECTOR 
                           'IL17A','IFNG','TBX21','TNF', #SPECIFIC EFFECTOR
                           'TOX', 'ICOS', 'PDCD1', 'BCL6', 'CXCR5', #TFH
                           'CXCL13', 'SLAMF6','CD38','IL4'),cols = 'RdYlBu')#TFH
#2 = cytocd8 
#1,3  = ncd4
#5 exhausted
#0  Th1 ; 4 IS GENERAL CD4 T EFFECTOR 

cluster.ids = c('Th1','nCD4','CytoCD8','nCD4','CD4+ effector T cell',
                'Unknown','Unknown')
names(cluster.ids) <- levels(sub_CD4)
sub_CD4 <- RenameIdents(sub_CD4, cluster.ids)
DimPlot(sub_CD4, reduction = "umap", label = TRUE, repel = TRUE)  

sub_CD4$celltype = Idents(sub_CD4)
sub_CD4@meta.data$celltype %>% table

sub_CD4 = subset(sub_CD4, idents = c('Th1','nCD4','CytoCD8','nCD4','CD4+ effector T cell'))

sub_CD4_1= sub_CD4

sub_CD4 = subset(CD4, idents = c(9,12)) 
sub_CD4 = SeuratAnalysis1(sub_CD4)
ElbowPlot(sub_CD4, ndims = 50)
sub_CD4 = SeuratAnalysis2(sub_CD4, harmony = T, harmony_by = 'orig.ident', n = 20)
DimPlot(sub_CD4,
        reduction = "umap",label = T, 
        raster = FALSE
)

DotPlot(sub_CD4, features = c('CD3D','CD3E','CD40LG','CCR6','CXCR6',#ncd4 
                          'CD8A','CD8B',
                          'IFNG', 'TBX21','TNF', 
                          'GATA3', 'STAT6', 'STAT5A', 'IRF4', 'IL4R','IL5','IL13','IL4',
                          'ISG15','IFIT1','IFIT2', 'IFITM2',
                          "LAG3",'PDCD1','CTLA4','CXCL13', 'HAVCR2','TOX'), cols = 'RdYlBu')

cluster0 = FindMarkers(sub_CD4, ident.1 =0, min.pct = 0.25)
#5=cytocd8, 1= ifn primed ; 2 = th1 ; 3= TH2 
#get rid of 4 

VlnPlot(sub_CD4, features = c( 'FOXP3', 'TIGIT', 'ICOS', 'CTLA4', 'TNFRSF9',
                           'TNFRSF4', 'IL2RA', 'IKZF2', 'RTKN2' )) #nCD8 nothing!
VlnPlot(sub_CD4, features = c('CD3D','CD3E','CD40LG', 'IFNG', 'TBX21','TNF', 'ISG15','IFIT1','IFIT2', 'IFITM2',
                               "LAG3",'PDCD1','CTLA4','CXCL13', 'HAVCR2','TOX'))

cluster.ids = c('Exhausted T cell','IFN-primed T cell','Th1',
                'Th2', 'Unknown','CytoCD8')
names(cluster.ids) <- levels(sub_CD4)
sub_CD4 <- RenameIdents(sub_CD4, cluster.ids)
DimPlot(sub_CD4, reduction = "umap", label = TRUE, repel = TRUE)  

sub_CD4$celltype = Idents(sub_CD4)
sub_CD4@meta.data$celltype %>% table

sub_CD4 = subset(sub_CD4, idents=c('Exhausted T cell','IFN-primed T cell','Th1',
                                   'Th2','CytoCD8'))
sub_CD4_2= sub_CD4


#sub_CD4 = subset(sub_CD4, idents = c(0,1,2,4,5,6,7))

#subcluster cd8+ T cell ------
CD8 = subset(Lymphoid, idents = c(1,8,11,17,39))
CD8_Final = SeuratAnalysis1(CD8_Final)
ElbowPlot(CD8_Final, ndims = 50)
CD8_Final = SeuratAnalysis2(CD8_Final, harmony = T, harmony_by = 'orig.ident', n = 20)

DimPlot(CD8_Final,
        reduction = "umap",label = T, 
        raster = FALSE, group.by = 'celltype'
)

table(CD8_Final$seurat_clusters)
saveRDS(CD8_Final, file = 'CD8_final_0126.rds')

#11 is hepatocyte but has some weird genes in it --> subcluster
cluster9 = FindMarkers(CD8_Final, ident.1 = 9, min.pct = 0.25)
NKT_genes = FindMarkers(CD8_Final, ident.1 = 'NCAM1+ NKT cell', min.pct=0.25)
IFN_genes =FindMarkers(CD8_Final, ident.1 = 'IFN-primed T cell', min.pct = 0.25)

DotPlot(CD8_Final, features = c('CD3D','CD3E','CD40LG',
                                'CD8A','CD8B','IFNG',    'ITGB1','GZMB','GZMK',
                                'CCR7', 'SELL', 'TCF7',  #tcm 
                                'CXCR6','ITGAE','ITGA1', #trm 
                                "LAG3",'PDCD1','HAVCR2',
                                'ISG15','IFIT1','IFIT2', 'IFITM2',
                                'IFNGR1','SLC4A10','TRAV1-2','TRBV6-2',
                                'NCAM1','FCGR3A','NCR1','ICAM1','KLRK1','GNLY','PRF1','KLRF1'
                                ), cols = 'RdYlBu')  + coord_flip()
DotPlot(CD8_Final, features = c('ISG15','IFIT1','IFIT2', 'IFITM2'), cols = 'RdYlBu')  + coord_flip()

VlnPlot(CD8_Final, features = c(#'CD3D','CD3E',
                                'CD40LG',
                                'CD8A','CD8B','IFNG', 'ITGB1',#'GZMB','GZMK',
                                'CCR7', 'SELL', 'TCF7',  #tcm 
                                'CXCR6','ITGAE','ITGA1', #trm 
                                "LAG3",'PDCD1','HAVCR2',
                                'ISG15','IFIT1','IFIT2', 'IFITM2'
                                #'IFNGR1','SLC4A10','TRAV1-2','TRBV6-2',
                               # 'NCAM1','FCGR3A','NCR1','ICAM1','KLRK1','GNLY','PRF1','KLRF1'
))

#look into trtex 
DotPlot(CD8_Final, features = TRTex_genes, cols = 'RdYlBu')+ coord_flip() + theme(axis.text.x = element_text(size = 8),
                                                                                  axis.text.y = element_text(size = 8))

FeaturePlot(CD8, features = c(#"IGKC",'IGHA1','JCHAIN','MKI67','TOP2A',
                          'CD3D','CD3E','CD40LG',
                          'CD8A','CD8B','IFNG','IFNGR1',
                          'ITGB1','GZMB','GZMK'), label = T, reduction ='umap')

#cd4:5,8,13 --> subcluster, have more others ; 8 NAIVE, 13 with il18ra, il7r, ifngr1, cd28,cd69
#cd8: all others 
#6 and 12 are ifn-primed, subcluster 6 
#1,5,8,10 CX3CR1 

#new:
#5 with xcl1, xcl2, ccl3, ccl3l1, ccl4
#2 with cx3cr1 

#first round: get rid of 3,9,12, 13,14  
#able to find foxp3+ nCD8 

#CD8 = subset(CD8, idents = c(0,1,2,4,5,6,7,8,10,11))
#CD8 = subset(CD8, idents = c(0,1,2,3,4,5,6,7,9,10,12))
#CD8 = subset(CD8, idents = c(0,1,2,3,4,5,6,7,8,9,11))
#CD8 = subset(CD8, idents = c(0,1,2,4,5,6,8,9,10))
#CD8=subset(CD8, idents = c(0,1,2,3,4,7,9,10,11))
#CD8= subset(CD8, idents = c(0,1,2,3,4,5,6,7,8,9))
#CD8= subset(CD8,idents = c(1,2,3,4,5,6,7,8))

DotPlot(CD8_Final, features = c('CCR7', 'SELL', 'CD44', 'TCF7',  #tcm 
                           'CXCR6','ITGAE','ITGA1', #trm 
                          'RORA','RORC','CX3CR1','TNF','PRDM1','GZMK')) #tem 

VlnPlot(CD8, features = c( 'FOXP3', 'TIGIT', 'ICOS', 'CTLA4', 'TNFRSF9',
                           'TNFRSF4', 'IL2RA', 'IKZF2', 'RTKN2' ))
VlnPlot(CD8, features =c("IFNG",'TBX21','ITGB1', 'GZMB','CD27','CD28' ))
VlnPlot(CD8, features = c(  'CXCR6','ITGAE','ITGA1','TNF','CD69','BCL6',
                            'GZMB','ITGB1','SELL'))
DotPlot(CD8, features = c( "LAG3",'PDCD1','CTLA4','CXCL13', 'HAVCR2','TOX', 'SLAMF6','CX3CR1'))

#some cd4 in 5 ??
#2= TCM (CCR7+TCF7+),  0= TEM , 1 = trm ; 3 = EXHAUSTED 
#Clusters with itgb1, ifng, gzm are not naive 
#others all CYTOCD8 ; no naive

cluster.ids = c('CytoCD8','CD8+ Trm','CD8+ Tcm','Exhausted T cell','CytoCD8',
                'CytoCD8', 'CytoCD8', 'CytoCD8', 'CytoCD8', 'CytoCD8')

names(cluster.ids) <- levels(CD8)
CD8 <- RenameIdents(CD8, cluster.ids)
DimPlot(CD8, reduction = "umap", label = TRUE, repel = TRUE)  

CD8$celltype = Idents(CD8)
CD8@meta.data$celltype %>% table



DotPlot(CD8, features = c('ISG15','IFIT1','IFIT2', 'IFITM2'))
DotPlot(CD8_Final,features = c('KLRC1', 'TRDC', 'TRDV1', 'TRGC1', 'TRGV9','TRG-AS1','KIR2DL4'))

#hepatocyte genes: 11 but be careful to validate it again
VlnPlot(CD8, features = c('HPR', 'F12', 'SLC13A5', 'HGFAC', 'PROZ', 
                          'AFM', 'AGXT', 'SAA4')) 




sub_CD8_4 = subset(CD8, idents = c(5,8,13))
sub_CD8_4 = SeuratAnalysis1(sub_CD8_4)
ElbowPlot(sub_CD8_4, ndims = 50)
sub_CD8_4 = SeuratAnalysis2(sub_CD8_4, harmony = T, harmony_by = 'orig.ident', n = 20)
DimPlot(sub_CD8_4,
        reduction = "umap",label = T, 
        raster = FALSE
)

#1 DEF C-NAIVE, 
#5= TNF+ GZMK+; 4 
#1,3 cd8+, 3 is MAIT!!! 
#2 HAS NCAM1, nk

VlnPlot(sub_CD8_4, features = c('CD3D','CD3E','CD40LG', 'SELL', 'CCR7', 'LEF1',
                                'TRAV1-2','IFNGR1','SLC4A10','TRBV6-2'))
                                #"IFNG",'TBX21','TNF','IL17A'))
DotPlot(sub_CD8_4, features = c("IL7R",'CD69','GZMA','PRDM1','CXCR6',
                                'ITGA1','ITGAE', 'XCL1','CD40LG','TNF','KLF4'))#TRM or T effector
#trm (cxcr6, itgae, itga1) 

VlnPlot(sub_CD8_4, features = c('TOX', 'ICOS', 'PDCD1', 'BCL6', 'CXCR5','CXCL13'))
#6 is not t cell, 

VlnPlot(sub_CD8_4, features = c('NCAM1','FCGR3A','NCR1','ICAM1','KLRK1',
                                'GNLY','GZMB','PRF1','KLRF1'))

cluster0= FindMarkers(sub_CD8_4, ident.1 =0, min.pct = 0.25)

cluster.ids = c('nCD4','CD4+CD8+ T cell','NCAM1+ NKT cell','MAIT','nCD4','nCD4')
names(cluster.ids) <- levels(sub_CD8_4)
sub_CD8_4 <- RenameIdents(sub_CD8_4, cluster.ids)
DimPlot(sub_CD8_4, reduction = "umap", label = TRUE, repel = TRUE)  

sub_CD8_4$celltype = Idents(sub_CD8_4)
sub_CD8_4@meta.data$celltype %>% table


#DOUBLE POSITIVE CONFIRMED 
#sub_CD8_4= subset(sub_CD8_4, idents =c(0,1,2,3,4,5,6))
cluster0_1 = subset(CD8, idents = 0)
cluster0_1 = SeuratAnalysis1(cluster0_1)
ElbowPlot(cluster0_1, ndims = 50)
cluster0_1 = SeuratAnalysis2(cluster0_1, harmony =F, n = 20)
DimPlot(cluster0_1,
        reduction = "umap",label = T, 
        raster = FALSE
)

cluster2 = FindMarkers(cluster0_1, ident.1 = 2, min.pct = 0.25)
VlnPlot(cluster0_1, features = c('CD3D','CD3E','CD40LG', 'SELL',  'LEF1',
                                 'CD8A','CD8B'))

DotPlot(cluster0_1, features = c('CCR7', 'SELL', 'CD44', 'TCF7',  #tcm 
                          'CXCR6','ITGAE','ITGA1', #trm 
                          'RORA','RORC','CX3CR1','TNF','PRDM1','GZMK')) #tem 

VlnPlot(cluster0_1, features = c( 'FOXP3', 'TIGIT', 'ICOS', 'CTLA4', 'TNFRSF9',
                           'TNFRSF4', 'IL2RA', 'IKZF2', 'RTKN2' ))
VlnPlot(cluster0_1, features =c('CD8A','CD8B',"IFNG",'TBX21',
                                    'ITGB1', 'GZMB','GZMK','CD27','CD28',
                                'SELL','LEF1'))

cluster0_1$celltype = 'CytoCD8'
#cluster0_1= subset(cluster0_1, idents = c(0,1,2,5,4))


IFN_CD8 = subset(CD8, idents = c(6,12))
IFN_CD8 = SeuratAnalysis1(IFN_CD8)
ElbowPlot(IFN_CD8, ndims = 50)
IFN_CD8 = SeuratAnalysis2(IFN_CD8, harmony = T, harmony_by = 'orig.ident', n = 20)
DimPlot(IFN_CD8,
        reduction = "umap",label = T, 
        raster = FALSE
)
cluster2 = FindMarkers(IFN_CD8, ident.1 = 2, min.pct = 0.25)

#get rid of 2,7,0,5
#ALL CD8+; 1= IFN ; 8 = MAIT, 4 HAS HIGH TCF7 ; 3 = TRM ; 4 = TCM 
#6 IS FGFBP2+FCGR3A+KLF, ITGB1

#0 = ifn
#5=mait, 1= tcm 
#3,4= trm 
#6 effector 

cluster.ids = c('IFN-primed T cell', 'CD8+ Tcm','CytoCD8','CD8+ Trm')
names(cluster.ids) <- levels(IFN_CD8)
IFN_CD8 <- RenameIdents(IFN_CD8, cluster.ids)
DimPlot(IFN_CD8, reduction = "umap", label = TRUE, repel = TRUE)  

IFN_CD8$celltype = Idents(IFN_CD8)
IFN_CD8@meta.data$celltype %>% table


VlnPlot(IFN_CD8,features = c( 'CD3D','CD3E','CD40LG','CD8A','CD8B',
                              'ISG15','IFIT1','IFIT2', 'IFITM2'))
DotPlot(IFN_CD8, features = c("IL7R",'CD69','GZMA','PRDM1','CXCR6',
                                'ITGA1','ITGAE', 
                              'TRAV1-2','IFNGR1','SLC4A10','TRBV6-2',
                              'TCF7','CCR7'
                              ))
VlnPlot(IFN_CD8, features = c('CCR7', 'SELL', 'CD44', 'TCF7', 'TRAV1-2','IFNGR1','SLC4A10','TRBV6-2',
                        'CXCR6','ITGAE','ITGA1',
                        'RORA','RORC','CX3CR1','TNF','PRDM1')) #tEM 


DotPlot(IFN_CD8, features = c( 'ZNF683', 'KLRC2', 'TRDC', 'TRDV1',  'TRGC1', 'TRGV9'))
#IFN_CD8=subset(IFN_CD8, idents = c(1,3,4,6,8))
IFN_CD8=subset(IFN_CD8, idents = c(0,1,3,4))
#3 CAN BE TRM ; 7,8 EFFECTOR 

#does not find anything in NKT subclustering
#the ifn-primed cluster seems to be mostly junk and doublet with IG genes--> check again 

TRTex = subset(CD8_Final, idents =c("Exhausted T cell",'IFN-primed T cell'))
TRTex=subset(TRTex, idents = c(3,5,7))
TRTex = SeuratAnalysis1(TRTex)
ElbowPlot(TRTex, ndims = 50)
TRTex = SeuratAnalysis2(TRTex, harmony = T,harmony_by = 'orig.ident', n = 20)
DimPlot(TRTex,
        reduction = "umap",label = T, 
        raster = FALSE
) 

FeaturePlot(TRTex, features=c('CD3D','CD3E',
                              'LAG3','PDCD1','CTLA4','CXCL13',
                              'CD8A','CD8B','CD40LG'), reduction = 'umap', label = T,
            cols = c('gray',"red"))
#0,1 NOT EXHAUSTED BUT CYTOTOXIC T CELLS 
#2,3 TRTEX 


FeaturePlot(TRTex, features=c('CD8A','CD8B','ITGB1','GZMB','IFNG', 'IFNGR1'), reduction = 'umap', label = T, cols = c('gray',"red"))

DotPlot(TRTex, features = TRTex_genes, cols = 'RdYlBu')+ coord_flip() + theme(axis.text.x = element_text(size = 8),
                                                                              axis.text.y = element_text(size = 8))

DotPlot(TRTex, features=c('CD3D','CD3E', 'LAG3','PDCD1','CTLA4','CXCL13'))

#TRTex=subset(TRTex, idents= c(0,1,2,4,5,6))
#7 SEEMS TO BE TRTEX 

cluster0 = FindMarkers(TRTex, ident.1 = 0, min.pct = 0.25)

cluster.ids = c("CytoCD8",'CytoCD8','TRTex','TRTex')
names(cluster.ids) <- levels(TRTex)
TRTex <- RenameIdents(TRTex, cluster.ids)
DimPlot(TRTex, reduction = "umap", label = TRUE, repel = TRUE)  

TRTex$celltype = Idents(TRTex)
TRTex@meta.data$celltype %>% table


#rename 
CD8_Final$barcodes = rownames(CD8_Final@meta.data)
Idents(CD8_Final)='barcodes'
CD8_rest = subset(CD8_Final, idents =CD8_Final$barcodes[!CD8_Final$barcodes %in% rownames(TRTex@meta.data)] )

CD8_Final= merge(CD8_rest, TRTex)
Idents(CD8_Final)='celltype'
CD8_Final=RenameIdents(CD8_Final, 'NCAM1+ NKT cell'='gdT')
CD8_Final$celltype=Idents(CD8_Final)
table(CD8_Final$celltype)


#concluding T and NK-----
#T_NK_final_0127 = merge(x=CD8_Final, y = c(NK_final, CD4_Final))
T_NK_final_0128= merge(T_NK_final_0127, Treg_Final)
T_NK_final_0128=  SeuratAnalysis1(T_NK_final_0128)
ElbowPlot(T_NK_final_0128, ndims = 50)
T_NK_final_0128 = SeuratAnalysis2(T_NK_final_0128, harmony = T, 
                             harmony_by = 'orig.ident', n =25)

DimPlot(T_NK_final_0128,
        reduction = "umap",label = T, repel=T,
        raster = FALSE , group.by = 'celltype',label.size = 3
)  

table(T_NK_final_0128$celltype)


saveRDS(T_NK_final_0128, file = 'T_NK_combined_0203.rds')

#check trNK, Cd4+ Tcm, CD4+ Tem (compare with cd4+ effector T)
Idents(T_NK_final_0128)='celltype'
unique(T_NK_final_0128$celltype)
T_NK_final_0128 = subset(T_NK_final_0128, idents = unique(T_NK_final_0128$celltype)[1:18] )


trNK_genes = FindMarkers(T_NK_final, ident.1 = 'trNK', min.pct = 0.25)
#these are trNK genes: indeed highly expressed: 
DotPlot(T_NK_final, features = c('KLRF1','PRF1','CD160',
                                 'XCL1','XCL2','CD247','GZMB','CXCR6'), group.by = 'orig.ident',cols = 'RdYlBu')

#tcm does have high tcf7 and sell ; tem has high rora, klrg1, tbx21,'GZMH','GZMB'
DotPlot(T_NK_final, features = c('CD40LG', 'IFNG', 'TBX21','TNF',#TH1 
                                 'KLRG1', 'KLF2','RORA','RORC','CCR7','CX3CR1','PRDM1','GZMH','GZMB','GZMA','ID2', #Tem
                                 'SELL', 'CD44', 'TCF7', #tcm 
                                 'CXCR6','ITGAE','ITGA1') %>% unique(), group.by = 'celltype', cols = 'RdYlBu')
 
#find foxp3 and naive cd8 
VlnPlot(T_NK_final, features = c(#'FOXP3','TIGIT','ICOS', 
                                 'CD8A', 'CD8B','CD40LG','ITGB1','IFNGR1','CD27','CD28','SELL'), group.by = 'seurat_clusters')
DotPlot(Lymphoid, features = c('CD3D','CD3E','FOXP3','TIGIT','ICOS', 
  'CD8A', 'CD8B','CD40LG','ITGB1','IFNGR1','CD27','CD28','SELL'), cols='RdYlBu')
DotPlot(Lymphoid, features = c('CD3D','CD3E','FOXP3','TIGIT','ICOS'), cols='RdYlBu', group.by='orig.ident')


#maybe 6,21,26 have nCD8? ; 0 has the highest foxp3 though not significant  

#still no foxp3 
cluster0_0 = subset(T_NK_final, idents =0) 
cluster0_0 = SeuratAnalysis1(cluster0_0)
ElbowPlot(cluster0_0, ndims = 50)
cluster0_0 = SeuratAnalysis2(cluster0_0, harmony = F,n = 20)
DimPlot(cluster0_0,
        reduction = "umap",label = T,  
        raster = FALSE, group.by='celltype'
) 
VlnPlot(cluster0_0, features = c('FOXP3','TIGIT','ICOS')) 


nCD8 = subset(T_NK_final, idents = c(6,26,21)) 
nCD8 = SeuratAnalysis1(nCD8)
ElbowPlot(nCD8, ndims = 50)
nCD8 = SeuratAnalysis2(nCD8, harmony = T,harmony_by = 'orig.ident', n = 20)
DimPlot(nCD8,
        reduction = "umap",label = T, 
        raster = FALSE, group.by = 'celltype'
) 
DotPlot(nCD8, features = c('CD8A', 'CD8B','CD40LG','ITGB1','IFNGR1','CD27','CD28','SELL'), cols = 'RdYlBu') 
#5,6 maybe 
VlnPlot(nCD8, features = c('CD8A', 'CD8B','CD40LG','ITGB1','IFNGR1',
                               'CD27','CD28','SELL')) 
cluster0 = FindMarkers(nCD8, ident.1 = 0, min.pct = 0.25 )
#clean up 7
#still not nCD8: the candidates are either junk (7), or solid IFN-primed (5), or solid GZM expressing, TCF7 expressing, or exhausted


T_NK_final$barcodes = rownames(T_NK_final@meta.data)
Idents(T_NK_final)='barcodes'
T_NK_final_0127 = subset(T_NK_final, idents = T_NK_final$barcodes[!T_NK_final$barcodes %in% rownames(nCD8@meta.data[nCD8$seurat_clusters == '7',])])

#subcluster Treg ----- 
Treg= subset(Lymphoid, idents =6)
Treg_Final= SeuratAnalysis1(Treg_Final)
ElbowPlot(Treg_Final, ndims = 50)
Treg_Final = SeuratAnalysis2(Treg_Final, harmony = T,harmony_by = 'orig.ident', n = 20)
DimPlot(Treg_Final,
        reduction = "umap",label = T, repel=T,
        raster = FALSE , group.by = 'celltype'
)

table(Treg$seurat_clusters)

FeaturePlot(Treg_Final, features = c('CD3D','CD3E',
                                     'FOXP3','TIGIT','ICOS','CTLA4',
                                     'TNFRSF9', 'TNFRSF4', 'IL2RA' 
                                     ), reduction='umap', label = F, cols = c("gray",'red'))
DotPlot(Treg_Final, features = c('CD3D','CD3E','CD40LG','CD8A','CD8B',
                           'TRDC','TRDV1','KLRD1','KLRC1',
                           'FOXP3','TIGIT','ICOS',
                           #"HAVCR2",'PDCD1','LAG3','CXCL13','CTLA4',
                           'ISG15','IFIT1','IFIT2',
                           'JCHAIN','MS4A1','BANK1','CD79A',
                           'CCR7','LEF1','TCF7'
                           ), cols='RdYlBu', group.by = 'celltype') 

VlnPlot(Treg, features = c('CD3D','CD3E','CD40LG','CD8A','CD8B',
                          #'SELL','CCR7','LEF1','CD27','CD28','ITGB1','GZMB',
                           'FOXP3','TIGIT','ICOS'
                           ))
VlnPlot(Treg, features = c("HAVCR2",'PDCD1','LAG3','CXCL13','CTLA4'))

#8, 10, 4, 7, 15,16 not T cell 
#treg: 14, 1, 3, 
#14 IS IFN-PRIMED 
#16=NCD4; 15=BCELL ; 10=HEATSHOCK (GET RID OF IT ) 
#2,12= CYTOCD8; 5 IS LIKE A NAIVE T WITH NO CLEAR CD4 OR CD8 
#GET RID OF 8,10,0,4 ,7 , 13 (CYCLING), 11 DOUBLET WITH PLASMA

#NEW:
#6= BCELL, 9= NCD4, 2,3 =CYTOCD8, 7 =IFN-PRIMED T CELL,10= GDT 
cluster.ids = c("Treg","Treg","CytoCD8","CytoCD8","Treg","Treg",
                'BCell','IFN-primed T cell','Treg','nCD4','gdT')
names(cluster.ids) <- levels(Treg)
Treg <- RenameIdents(Treg, cluster.ids)
DimPlot(Treg, reduction = "umap", label = TRUE, repel = TRUE)  

Treg$celltype = Idents(Treg)
Treg@meta.data$celltype %>% table


#14 IS THE KERATINOCYTE/EPI THING ; 4 = EXHAUSTED? 
#GET RID OF 0,11

#these are subclustering: 
#8 is a doublet with myeloid(maybe mono or neutrophil); 12 is a doublet with a structural cell
#2 has hcc tissue genes


cluster10 =FindMarkers(Treg, ident.1 = 10, min.pct = 0.25)
#Treg= subset(Treg, idents = c(1,2,3,5,6,7,9,12,14,15,16))
#Treg = subset(Treg, idents = c(0,3,4,5,7,9,10,11,13,14))
#Treg= subset(Treg, idents = c(1,2,3,5,6,7,8,9,10,12,13,14))
#Treg= subset(Treg, idents = c(0,1,2,3,4,5,7,8,9,11))


suB_myeloidreg = subset(Treg, idents = c(1,2,8,12,6))
suB_myeloidreg = SeuratAnalysis1(suB_myeloidreg)
ElbowPlot(suB_myeloidreg, ndims = 50)
suB_myeloidreg = SeuratAnalysis2(suB_myeloidreg, harmony = T,harmony_by = 'orig.ident', n = 20)
DimPlot(suB_myeloidreg,
        reduction = "umap",label = T,
        raster = FALSE , group.by = 'orig.ident'
) 

FeaturePlot(suB_myeloidreg, features = c('CD3D','CD3E','CD40LG','CD8A','CD8B',
                               'SELL','CCR7','LEF1','FOXP3','TIGIT','ICOS'
), reduction='umap', label = T)

cluster7 =FindMarkers(suB_myeloidreg, ident.1 = 7, min.pct = 0.25)
VlnPlot(suB_myeloidreg, features = c("PI3",'LY6D','KRT6B', 'KRT16',
                               'SPRR2E',
                               'KRT1','KRT10','KRT5','KRT14',
                               'ISG15','IFIT1','IFIT2'), group.by = 'seurat_clusters')

VlnPlot(suB_myeloidreg, features = c('ALDOB','HPR','APOA1','APOA2', 'F12',  'HGFAC', 'AGXT', 'SAA4','ARG1'))
table(suB_myeloidreg$seurat_clusters)

#8 = terminally differentiated keratinocytes 
#9 HAS KRT GENES UPREGULATED BUT ALSO CXCR6, IL21R, IL1R2  AND S100A8/9/14
#9= inflammatory keratinocyte 

#new: 6,7 = inflammatory keratinocyte; 3 is just naive t 

cluster.ids = c("Treg","Treg","Treg","Treg","Treg","Treg",
                'Terminally differentiated inflammatory keratinocyte','Inflammatory keratinocyte')
names(cluster.ids) <- levels(suB_myeloidreg)
suB_myeloidreg <- RenameIdents(suB_myeloidreg, cluster.ids)
DimPlot(suB_myeloidreg, reduction = "umap", label = TRUE, repel = TRUE)  

suB_myeloidreg$celltype = Idents(suB_myeloidreg)
suB_myeloidreg@meta.data$celltype %>% table


#suB_myeloidreg= subset(suB_myeloidreg, idents = c(0,2,3,4,5,7))

#subcluster B and plasma -------
BCell = subset(Lymphoid, idents = c(0,30,19,34,18,9))  
BCell=  SeuratAnalysis1(BCell)
ElbowPlot(BCell, ndims = 50)
BCell = SeuratAnalysis2(BCell, harmony = T, 
                                  harmony_by = 'orig.ident', n =30)

DimPlot(BCell,
        reduction = "umap",label = T,
        raster = FALSE ,label.size = 7
) 

table(BCell$seurat_clusters)

DotPlot(BCell, features = c( 'JCHAIN','IGHA1','IGHG1','IGHM','IGKC',
                             'MS4A1','BANK1','CD79A','CD79B',"CD19"), cols = 'RdYlBu')
#NON B CELL: 2,3,5,18, 20, 22
#PLASMA: 7,8,13,24,19 ; 14 and 15 only jchain pos not the ig gene 
#b cell: 0,1,4,6,16,17,12 (maybe 9,10,21) 

#LOOK AT THE LUNG/BREAST PAPER FOR SPECIFIC B CELL MARKERS 
DotPlot(BCell, features = c('CD3D','CD3E',
                            'CD40LG','CD8A','CD8B','KLRF1','PRF1'), cols = 'RdYlBu')
#2,3,20,21,22,5,9,10 = T CELLS , 10 = CD8 
#MAYBE 2,3,11 also t cell 
#18 is not junk, but it contains IFN genes, CLEC4A... 


FeaturePlot(BCell, features = c( 'JCHAIN','IGHA1','IGHG1','IGHM',
                                 'MS4A1','BANK1','CD79A','CD79B',"CD19"
                                 #'CD3D','CD3E',
                                 #'CD40LG','CD8A','CD8B'
                                 ), reduction='umap', label = T, cols = c('gray','red'))
cluster18 = FindMarkers(BCell, ident.1 = 18, min.pct = 0.25)


sub_B= subset(BCell, idents = c(0,4,16,9,6,1,23,17,12))
sub_B=  SeuratAnalysis1(sub_B)
ElbowPlot(sub_B, ndims = 50)
sub_B = SeuratAnalysis2(sub_B, harmony = T, 
                        harmony_by = 'orig.ident', n =30)
DimPlot(sub_B,
        reduction = "umap",label = T,
        raster = FALSE ,label.size = 5
) 

FeaturePlot(sub_B, features = c( 'JCHAIN','IGHA1','IGHG1','IGHM',
                                 'MS4A1','BANK1','CD79A','CD79B',"CD19",
                                 'CD3D','CD3E',
                                 'CD40LG','CD8A','CD8B'
), reduction='umap', label = T, cols = c('gray','red'))

DotPlot(sub_B, features = c( 'JCHAIN','IGHA1','IGHG1','SDC1',
                             'IGHE',
                             #'MZB1','TOP2A','CDK1','MKI67', #PLASMABLAST
                             'IGHM','TCL1A','FCER2', # NAIVE 
                             'IGHD','TNFRSF13B', #MEMORY 
                             'MS4A1','BANK1','CD79A','CD79B',"CD19",
                             'CD3D','CD3E',
                             'CD40LG'), cols = 'RdYlBu')+coord_flip()

DotPlot(sub_B, features = c( 'CD3D','CD3E',
                             'CD40LG','SELL','CCR6','CXCR6',
                             'CD8A','CD8B'))
VlnPlot(sub_B, features = c("TOP2A",'MKI67','IGHE'))

cluster0 = FindMarkers(sub_B, ident.1 = 0, min.pct = 0.25)


cluster.ids = c("BCell", "Plasma", "BCell",'naive BCell', "BCell","nCD4","BCell","BCell","BCell","BCell",
                "BCell","BCell",'Plasma',"nCD4","BCell","IgE+ BCell")

names(cluster.ids) <- levels(sub_B)
sub_B <- RenameIdents(sub_B, cluster.ids)
DimPlot(sub_B, reduction = "umap", label = TRUE, repel = TRUE)  

sub_B$celltype = Idents(sub_B)
sub_B@meta.data$celltype %>% table


#new: 
#3= naive b
#t cell = 5,13 (naive cd4)
#plasma = 1,12
#20= keratinocyte features ; 16= IGHE 

#GET RID OF 18,16,5
#sub_B=subset(sub_B, idents= c(0:4, 6:15, 17,19:22))
#sub_B= subset(sub_B, idents = c(0:4, 6:15))
#sub_B= subset(sub_B, idents = c(0:15, 17,18,19,20))
#sub_B = subset(sub_B, idents = c(0:13,15,16))
saveRDS(sub_B, file = 'BCell_subcluster.rds')

plasma = subset(BCell, idents = c(7,8,13,24,19,14,15))
plasma=  SeuratAnalysis1(plasma)
ElbowPlot(plasma, ndims = 50)
plasma = SeuratAnalysis2(plasma, harmony = T, 
                        harmony_by = 'orig.ident', n =15)
DimPlot(plasma,
        reduction = "umap",label = T,
        raster = FALSE ,label.size = 5, group.by = 'orig.ident'
) 

cluster2 = FindMarkers(plasma, ident.1 = 2, min.pct =0.25)


DotPlot(plasma, features = c( 'JCHAIN','IGHA1','IGHG1','SDC1',
                             'IGHE',
                             'MZB1','TOP2A','CDK1','MKI67', #PLASMABLAST
                             'IGHM','TCL1A','FCER2', # NAIVE 
                             'IGHD','TNFRSF13B', #MEMORY 
                             'MS4A1','BANK1','CD79A','CD79B',"CD19",
                             'CD3D','CD3E',
                             'CD40LG'), cols = 'RdYlBu')+coord_flip()
DotPlot(plasma, features = c( 'CD3D','CD3E',
                             'CD40LG','SELL','CCR6','CXCR6',
                             'CD8A','CD8B'))
VlnPlot(plasma, features = c("TOP2A",'MKI67','IGHE'))
VlnPlot(plasma, features=c('CLEC4C', 'CD1C','PTGDS',  'IL3RA', 'CLIC3') ) 
VlnPlot(plasma, features = c( 'CD3D','CD3E','TRDC','TRDV1','TRGV9','KLRC2','KLRD1'))
FeaturePlot(plasma, features=c('JCHAIN','IGHA1','IGHG1', 'MS4A1','BANK1','CD79A','CD79B',"CD19") , reduction='umap', label = T) 

cluster.ids = c('pDC','Plasma','Plasma','pDC','Plasma','Plasma','pDC',
                'Plasma','Plasma','Plasmablast','CytoCD8','Plasma','Plasma','gdT')
names(cluster.ids) <- levels(plasma)
plasma <- RenameIdents(plasma, cluster.ids)
DimPlot(plasma, reduction = "umap", label = TRUE, repel = TRUE)  

plasma$celltype = Idents(plasma)
plasma@meta.data$celltype %>% table


#0,3,6,13 = pdc 
#10,13 = t cell; 10 = cytocd8 
#9 =plasmablast  

  
#plasma= subset(plasma, idents = c(0:6, 8:18))
#plasma = subset(plasma, idents = c(0:15, 17))
#plasma = subset(plasma, idents = c(0:11,13,14,16,17))
#plasma = subset(plasma, idents = 1:14)
saveRDS(plasma, 'Plasma.rds')


other_B = subset(BCell, idents = c(2,3,11,10,20,21,22,5,18))
other_B=  SeuratAnalysis1(other_B)
ElbowPlot(other_B, ndims = 50)
other_B = SeuratAnalysis2(other_B, harmony = T, 
                         harmony_by = 'orig.ident', n = 20)
DimPlot(other_B,
        reduction = "umap",label = T,
        raster = FALSE ,label.size = 5)

#new:
#5 is doublet plasma with t cell; 14 has cdc1 genes ; 11 is monocyte 
#12?
#1 is macrophage 
#pdc = 0,2,3,4,6,12
#plasma: 1,7,8,10,13
#cd3: 1,7,8,9,13

cluster.ids=c("pDC",'Macrophage','pDC','pDC','pDC','Plasma-T doublet','pDC',
              'Plasma-T doublet','Plasma-T doublet','Plasma-T doublet','Plasma',
              'Plasma-Monocyte','pDC','Plasma-T doublet','cDC1')

names(cluster.ids) <- levels(other_B)
other_B <- RenameIdents(other_B, cluster.ids)
DimPlot(other_B, reduction = "umap", label = TRUE, repel = TRUE)  
other_B$celltype = Idents(other_B)
other_B@meta.data$celltype %>% table


DotPlot(other_B, features = c( 'JCHAIN','IGHA1','IGHG1','IGHM', 
                             'MS4A1','BANK1','CD79A','CD79B',"CD19",
                             'MZB1','TOP2A','CDK1','MKI67', #PLASMABLAST
                             'TCL1A','FCER2'), cols = 'RdYlBu')
VlnPlot(other_B, features = c('CD3D','CD3E',
                              'CD40LG','CD8A','CD8B'))

VlnPlot(other_B, features = c("TOP2A",'MKI67','IGHE'))
VlnPlot(other_B, features=c('CLEC4C', 'CD1C','PTGDS',  'IL3RA', 'CLIC3'))

FeaturePlot(other_B, features=c('CLEC4C', 'CD1C','PTGDS',  'IL3RA', 'CLIC3'),  cols = c("gray",'red'), label = T, reduction ='umap')

FeaturePlot(other_B, Macro_mono_genes, cols = c("gray",'red'), label = T, reduction ='umap')
VlnPlot(other_B, features = c("CD14",'FCGR3A','CX3CR1','ITGAM','CD68','CD163'))

VlnPlot(other_B, features = c('JCHAIN','IGHA1','IGHG1','IGHM', 
                                     'MS4A1','BANK1','CD79A','CD79B',"CD19"))
VlnPlot(other_B, features = c("IL1B",'FCN1','CXCL8','S100A8','S100A9',
                              'CLEC10A','CLEC12A','CD1C'))
DotPlot(other_B, DC_genes, cols = 'RdYlBu')
DotPlot(other_B, features = c("C1QA",'C1QB','TNF','CCL3','CCL4','CCL3L1','CCL4L1',
                              'CD14','CD163','TREM2'), cols = 'RdYlBu')

#current version:
#plasma:2,14,18   
#pdc= 4,6,9
#20= cdc2a
#check 7,13; 
#SUBCLUSTER 19 AND 21 BECAUSE THERE SEEMS TO BE OTHER CELLS UNDER IG GENES, NOT NECEARRILY DOUBLET 
#monocyte?: 1, 8, 18
#still unknown: 18,8,7,13 

#old:
#3= cytocd8
#get rid of 8 and 16 
#5,15 = t cell (cytocd8); MAYBE 14 (MAYBE NK BUT VERY INFLAMMATORY); 5 MAYBE A DOUBLET 
#pdc: 4,6,7,8
#there are some myeloid/monocyte genes: 1 = mono,  16 HAS CCR7 AND LAMP3, 18 IS NOT JUNK 
#general T cells with cd3: 0,4,5,14,15.18,19,20? 

cluster11 = FindMarkers(other_B, ident.1 = 11, min.pct =0.25)


table(other_B$seurat_clusters)

#subcluster 5
#get rid of 1 , 16, 19
#second run: get rid of 19


#other_B= subset(other_B, idents = c(0,2,3:15,17,18,20,21))
#other_B= subset(other_B, idents = c(0:18, 20))
#other_B = subset(other_B, idents = c(0:7, 9:15,17,18))
#other_B=subset(other_B, idents = c(2,4,6,7,8,9,13,14,18,20))
other_B= subset(other_B, idents = c(0,1,2,3,5,6,7, 9,10,11,12,13,14))
other_B= subset(other_B, idents =c(0:7,9:14))
#subcluster 5 from other_b, it has both cytocd8 and b cells 

B_myeloid = subset(other_B, idents = c(0,12,19,21,1,11,3,15))
B_myeloid=  SeuratAnalysis1(B_myeloid)
ElbowPlot(B_myeloid, ndims = 50)
B_myeloid = SeuratAnalysis2(B_myeloid, harmony = F, n = 20)
DimPlot(B_myeloid,
        reduction = "umap",label = T,
        raster = FALSE ,label.size = 5)

table(B_myeloid$seurat_clusters)
VlnPlot(B_myeloid, Macro_mono_genes)
VlnPlot(B_myeloid, features = c('CD3D','CD3E','CD8A','CD8B','JCHAIN','CD19'))

DotPlot(B_myeloid, features = c('JCHAIN','IGHA1','IGHG1','IGHM', 
                            'MS4A1','BANK1','CD79A','CD79B',"CD19",
                            'CD3D','CD3E','CD8A','CD8B'
                            ), cols = 'RdYlBu')

#new:
#plasma: 8,10,12
#pdc = 11
#bcell = 3 
#myeloid : 4(mono), 7 (pvm)
cluster.ids = c("Macrophage",'Macrophage','Macrophage','BCell',
                'Monocyte','Macrophage','Macrophage','PVM',
                'Plasma','Macrophage','Plasma','pDC','Plasma','Macrophage')
names(cluster.ids) <- levels(B_myeloid)
B_myeloid <- RenameIdents(B_myeloid, cluster.ids)
DimPlot(B_myeloid, reduction = "umap", label = TRUE, repel = TRUE)  

B_myeloid$celltype = Idents(B_myeloid)
B_myeloid@meta.data$celltype %>% table


#plasma: 11,21,20,23; 9 is pdc 
#bcell: 17-19, 
#cytocd8 10 
#3,12,16 = monocyte; 4,7 and 18 is unclear myeloid 
#0 is unknown with some t cell 

#6,12,13,16 has cd14


B_myeloid=subset(B_myeloid,idents = c(0:10,12:25 ))
B_myeloid=subset(B_myeloid, idents = c(11,21,20,23,9,17,18,19,10,3,12,16,4,7,0))
B_myeloid = subset(B_myeloid, idents = c(0:9, 11,13,14))



DotPlot(B_myeloid, features=c('CLEC4C', 'CD1C','PTGDS',  'IL3RA', 'CLIC3'), cols='RdYlBu' ) 
FeaturePlot(B_myeloid, c(Macro_mono_genes,'TNF',
                         'JCHAIN','BANK1','CD79A','MS4A1'), cols = c("gray",'red'), label = F, reduction ='umap')
FeaturePlot(B_myeloid, features = c('JCHAIN','IGHA1','IGHG1','IGHM', 
                                    'MS4A1','BANK1','CD79A','CD79B',"CD19"), 
            cols = c("gray",'red'), label = T, label.size =3, reduction ='umap')

FeaturePlot(B_myeloid, features = c("TNF",'CCL3','CCL4','CCL3L1',
                                    'C1QC','IFITM1','ISG15','CXCL9'), 
            cols = c("gray",'red'), label = T, label.size =3, reduction ='umap')


DotPlot(B_myeloid, Macro_mono_genes, cols = 'RdYlBu')

VlnPlot(B_myeloid, features = c(Macro_mono_genes[1:7], "IL1B",'FCN1','CXCL8','S100A8','S100A9'
                              #'CLEC10A','CLEC12A','CD1C'
                              )) 

cluster9 = FindMarkers(B_myeloid, ident.1 = 9, min.pct=0.25)

#don't include all of these!
B_T = subset(other_B, idents = c(5,16,10,17))
B_T =  SeuratAnalysis1(B_T)
ElbowPlot(B_T, ndims = 50)
B_T = SeuratAnalysis2(B_T, harmony = T, 
                          harmony_by = 'orig.ident', n = 20)
DimPlot(B_T,
        reduction = "umap",label = T,
        raster = FALSE ,label.size = 5)
VlnPlot(B_T, features = c('JCHAIN','IGHA1','IGHG1','IGHM', 
                          'MS4A1','BANK1','CD79A','CD79B',"CD19",
                          'CD3D','CD3E','CD8A','CD8B'))
FeaturePlot(B_T, features = c('JCHAIN','IGHA1','IGHG1','IGHM', 
                          'MS4A1','BANK1','CD79A','CD79B',"CD19",
                          'CD3D','CD3E','CD8A','CD8B'), cols = c("gray",'red'),
            reduction = 'umap', label = T)
cluster5 = FindMarkers(B_T, ident.1 = 5, min.pct=0.25)

B_T= subset(B_T, idents = c(0,1,2,3,4,6,7,8))
saveRDS(B_T,file = 'BCell_TCell_doublet.rds')

#merge all b cell -----
BCell_final = merge(x= other_B, y = c(B_myeloid, plasma, sub_B))
BCell_final = SeuratAnalysis1(BCell_final)
ElbowPlot(BCell_final, ndims = 50)
BCell_final = SeuratAnalysis2(BCell_final, harmony = T, 
                              harmony_by = 'orig.ident', n = 15)
DimPlot(BCell_final,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 5, group.by = 'celltype')


table(BCell_final$celltype)
Idents(BCell_final)='celltype'
BCell_final = subset(BCell_final, idents = c('BCell','cDC1','CytoCD8','gdT',
                                             'IgE+ BCell','Monocyte','naive BCell',
                                             'nCD4','pDC','Plasma','Plasma-Monocyte',
                                             'Plasmablast','PVM'))
#the second time: 
BCell_final = subset(BCell_final, idents = c('BCell',
                                             'IgE+ BCell','naive BCell',
                                             'nCD4','pDC','Plasma',
                                             'Plasmablast','PVM'))


VlnPlot(BCell_final, features = c('JCHAIN','IGHA1','IGHG1','IGHM', 
                          'MS4A1','BANK1','CD79A','CD79B',"CD19",
                          'CLEC4C','PTGDS', 'IL3RA'
                         # 'CD3D','CD3E','CD8A','CD8B'
                          ), group.by = 'celltype')
VlnPlot(BCell_final, features = c('TOP2A','MKI67'), group.by = 'seurat_clusters')
DotPlot(BCell_final, features = c('JCHAIN','IGHA1','IGHG1','IGHM', 
                                  'MS4A1','BANK1','CD79A','CD79B',"CD19",
                                  'CLEC4C', 'CD1C','PTGDS',  'IL3RA',
                                  'CD3D','CD3E','CD40LG',
                                  'KLRD1','GNLY','TRDC','TRDV1','TRGC2',
                                  'CXCL8','IL1B','FCN1','ITGAX','ITGAM',
                                  'CD14','CD163','C1QA','SELENOP','LYVE1',
                                  'CLEC9A','XCR1'), cols = 'RdYlBu', group.by = 'celltype')+coord_flip()
DotPlot(BCell_final, features = c('JCHAIN','IGHA1','IGHG1','SDC1',
                                  'IGHE',
                                  'MZB1','TOP2A','CDK1','MKI67', #PLASMABLAST
                                  'IGHM','TCL1A','FCER2', # NAIVE 
                                  'IGHD','TNFRSF13B', #MEMORY 
                                  'MS4A1','BANK1','CD79A','CD79B',"CD19" ), cols = 'RdYlBu', group.by = 'celltype')+coord_flip()

VlnPlot(BCell_final, features = c('MS4A1','CD79A',"CD19",
                                  'TNFRSF13B', #MEMORY
                                  'CD22', #follicular
                                  'TCL1A','FCER2', #NAIVE 
                                  'IGLL1','RAG1','DNTT', #pro and pre-B
                                    'AICDA','BCL6','SUGCT', #germination center B cell
                                  'JCHAIN',  'XBP1','SDC1','CD38',  'CERCAM','ITGA8' #plasma and plasmablast, the last two are shown to be only in plasma 
), raster=F, group.by='seurat_clusters') 
#7 has aicda and sugct, maybe germination center b cell 

FeaturePlot(BCell_final, features = c('JCHAIN','IGHA1','IGHG1','IGHM', 
                              'MS4A1','BANK1','CD79A','CD79B',"CD19",
                              'CD3D','CD3E','CD8A','CD8B'), cols = c("gray",'red'),
            reduction = 'umap', label = T)

cluster7 = FindMarkers(BCell_final, ident.1 = 7, min.pct = 0.25)

saveRDS(BCell_final, file = 'BCell_final_0212.rds')


#merge all the lymphoid ------
lymphoid_0212= merge(BCell_final, T_NK_final_0128)
#lymphoid_0212= JoinLayers(lymphoid_0212)
table(lymphoid_0212$celltype)

lymphoid_0227= SeuratAnalysis1(lymphoid_0227)
ElbowPlot(lymphoid_0227, ndims = 50)
lymphoid_0227 = SeuratAnalysis2(lymphoid_0227, harmony = T, 
                        harmony_by = 'orig.ident', n = 20)

DimPlot(lymphoid_0227,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 3, group.by = 'celltype')

lymphoid_0227$cluster_celltype = paste0(lymphoid_0227$seurat_clusters, "_", lymphoid_0227$celltype)
table(lymphoid_0227$celltype)

#put the cd4+ tcm back
lymphoid_rest = subset(lymphoid_0212, idents = unique(lymphoid_0212$celltype)[c(1:18, 20:25)])
lymphoid_0212 = merge(lymphoid_rest, CD4_TCM)

saveRDS(lymphoid_0212, file = 'Lymphoid_final_v2_0213.rds')

#put the trnk back 
lymphoid_rest = subset(lymphoid_0212, idents = unique(lymphoid_0212$celltype)[c(1:18, 20:28)])
lymphoid_0227 = merge(lymphoid_rest, trNK)


lymphoid_0227 = JoinLayers(lymphoid_0227)
saveRDS(lymphoid_0227, file = 'Lymphoid_final_v2_0227.rds')


table(lymphoid_0212$seurat_clusters)
Idents(lymphoid_0212) ='celltype'
CD4_TCM_genes = FindMarkers(lymphoid_0212, ident.1 = 'CD4+ Tcm', min.pct = 0.25)

#reexamine CD4+ Tcm
CD4_TCM = subset(lymphoid_0212, idents = 'CD4+ Tcm')
CD4_TCM= SeuratAnalysis1(CD4_TCM)
ElbowPlot(CD4_TCM, ndims = 50)
CD4_TCM = SeuratAnalysis2(CD4_TCM, harmony = F,n=20)
DimPlot(CD4_TCM,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 5, group.by = 'seurat_clusters')

VlnPlot(CD4_TCM, features = c("CD3D",'CD3E','CD40LG','TCF7',
                              'KLRF1','PRF1',
                              'CCR7','RORA','RORC','CD8A','CD8B' ), group.by= 'seurat_clusters')
FeaturePlot(CD4_TCM, features = c( 'CD40LG','TCF7','KLRG1', 
                                   'ZNF683', 'KLRC2', 'TRDC', 'TRDV1', 'TRGC1', 'TRGV9','TRG-AS1','KIR2DL4',
                                  'CD3D','CD3E','CD8A','CD8B'), cols = c("gray",'red'),
            reduction = 'umap', label = T)


VlnPlot(CD4_TCM, features= c(#'AREG','SOX4','KIT','TNFRSF18','TBX21','ZNF683',
                             #'CXCR6',
                             'SOX4','KIT', 'TBX21', 'AQP3','RBPJ', 'HPGD', #ILC1 
                             'VEGFA','AFF3','GATA3', #ILC3
                             'PTGS2', #ILC2 
                             'SELL', 'TMPRSS11E', 'FXYD7', 'XCL1', 'TNFRSF4',#ILC PROGENITOR 
                             'CCR6','CD200R1','KIT','KLRB1')) 
VlnPlot(CD4_TCM, features = c('ISG15','IFITM1','IFITM2'))
VlnPlot(CD4_TCM, features = c('IL23R','IL1R1','KIT','CCR6','CXCR6',
                              'IL17R','CD4','TGFBR2','TCF7','RORA','RORC')) #ILC3 
VlnPlot(CD4_TCM, features= c( 
  'ZNF683','TBX21', 'AQP3','RBPJ', 'HPGD', #ILC1  
  'PTGS2', 'IL1RL1','GATA3', #ILC2 
  'SELL', 'IL1R1','TMPRSS11E', 'FXYD7', 'XCL1', 'TNFRSF4',#ILC PROGENITOR 
  'CCR6','CD200R1','KIT','KLRB1',#GENERAL 
  'NCAM1','EOMES')) #PROVE NOT NK 
DotPlot(CD4_TCM, features= c( 
  'ZNF683','IL32', 'CEBPB','RBPJ', 'GEM', #ILC1  
  'PTGS2', 'IL1RL1','HPGD', 'SLAMF1', #ILC2 
  'IL1R1','KIT','CCR6','AFF3', 'VEGFA',
  'SELL', 'TNFRSF4','TMPRSS11E', 'FXYD7', 'XCL1'  #ILC PROGENITOR  
  ), cols = 'RdYlBu', group.by = 'celltype') + coord_flip() #PROVE NOT NK 



cluster3 = FindMarkers(CD4_TCM, ident.1 = 3, min.pct = 0.25)
#1 IS GDT, 4 seems to be ifn primed but not sure t cells 
#0,5,6 MAY BE ILC 

#0,4 = ILC3; 2=ILC-PROGENITOR (HIGHER SELL) ; 1= ILC1 WITH HIGH ZNF683, 3 = ILC2 

cluster.ids = c('ILC3','ILC1','ILC-progenitor', 'ILC2','ILC3')
names(cluster.ids) <- levels(CD4_TCM)
CD4_TCM <- RenameIdents(CD4_TCM, cluster.ids)
DimPlot(CD4_TCM, reduction = "umap", label = TRUE, repel = TRUE)  

CD4_TCM$celltype = Idents(CD4_TCM)
CD4_TCM@meta.data$celltype %>% table

#CD4_TCM = subset(CD4_TCM, idents = c(0,1,2,3,4,6))
CD4_TCM = subset(CD4_TCM, idents = c(0,4,5,6,1))


#check trNK
#trNK = subset(lymphoid_0212, idents = 'trNK')
trNK= SeuratAnalysis1(trNK)
ElbowPlot(trNK, ndims = 50)
trNK = SeuratAnalysis2(trNK, harmony = F,n= 15)
DimPlot(trNK,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 3)

VlnPlot(trNK, features = c( 'KLRF1','PRF1',
                            'GZMB','CD160','XCL1','CD247', 'CXCR6','TOX2','TIGIT','IRF8',
                            'NCAM1','FCGR3A',  
                            'ADGRG1','FCER1G','LAIR2','SPON2',
                           'CCL3','TYROBP','CLIC3', 'ADGRG3','IL2RB','SPRY2','GZMK'  #FROM CELLTYPIST OF IMMATURE NK 
                           # 'CD3D','CD3E','CD8A','CD8B','CD40LG'
                            ))
VlnPlot(trNK, features = c('JCHAIN','IGHM','IGHD', 'IGHA1'))
FeaturePlot(trNK,  features = c( 'KLRF1','PRF1',
                                 'GZMB','CD160','XCL1','CD247'), reduction='umap',
            cols = c('gray','red'))
cluster3= FindMarkers(trNK, ident.1 = 3, min.pct = 0.25)
#3 AND 6 NOT TRNK 

table(trNK$seurat_clusters)
#trNK = subset(trNK, idents = c(0,1,2,4,5))
#trNK = subset(trNK, idents = c(0,3,4))

saveRDS(trNK, '~/jliang/Human tissue harmonization/Lymphoid/trNK_subcluster_v2.rds')


#the rest of the clusters in Lymphoid raw------
Unknown = subset(Lymphoid, idents = 
                   unique(Lymphoid$seurat_clusters) [!unique(Lymphoid$seurat_clusters) %in% c(0,30,19,34,18,9,4,5,12,32,33,2,3,15,1,8,11,17,39,6)])
Unknown= SeuratAnalysis1(Unknown)
ElbowPlot(Unknown, ndims = 50)
Unknown = SeuratAnalysis2(Unknown, harmony = T, 
                        harmony_by = 'orig.ident', n = 20)
DimPlot(Unknown,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 5, group.by = 'orig.ident')

DotPlot(Unknown, features = c( "CD19",'CD79A','BANK1','MS4A1','JCHAIN',
  'KLRF1','PRF1',
  'CD3D','CD3E',
  'CD14','CD163','IL1B','CXCL8', 'FCN1',
  'ITGAM','ITGAX'
  #'GZMB','CXCR6','CD160','XCL1','CD247', #MATURE 
  #'NCAM1','FCGR3A','ICAM1','KLRK1','NCR1', #IMMATURE 
 # 'KLRC2','CD2' #ADAPTIVE
), cols = 'RdYlBu')+
  coord_flip() + theme(axis.text.x = element_text(size = 8),
                       axis.text.y = element_text(size = 8))

#11 IS PLASMA + BCELL MIXED 
#T CELL: 0,4,9,10, 13-16, 27
#NK: 3 , 19, 25
#MYELOID: 1,2,6,7,8,12,17,18,24,26,20,21 
#UNKNOWN: 5, 12, 22,23  -- get rid of them 

cluster23 = FindMarkers(Unknown, ident.1 = 23, min.pct = 0.25)

table(Unknown$seurat_clusters)

#subcluster unknown plasma and b cell --------
Unknown_b= subset(Unknown, idents = 11)
Unknown_b = SeuratAnalysis1(Unknown_b)
ElbowPlot(Unknown_b, ndims = 50)
Unknown_b = SeuratAnalysis2(Unknown_b, harmony = T, 
                          harmony_by = 'orig.ident', n = 20)
DimPlot(Unknown_b,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 5, group.by = 'celltype')

VlnPlot(Unknown_b, features = c( 'JCHAIN','IGHA1','IGHG1','IGHM',
                                 'MS4A1','BANK1','CD79A','CD79B',"CD19",
                                 'CD3D','CD3E',
                                 'CD40LG','CD8A','CD8B',
                                 'KLRF1','PRF1','GZMB','CXCR6','CD160','XCL1','CD247'
))

DotPlot(Unknown_b, features = c( 'JCHAIN','IGHA1','IGHG1','IGHM',
                                     'MS4A1','BANK1','CD79A','CD79B',"CD19",
                                 'TOP2A','MKI67',
                                     'CD3D','CD3E',
                                     'CD40LG','CD8A','CD8B'
), cols = 'RdYlBu')


#v4:
#plasma: 0,4,9,10
#t cell: 2;
# 7 = TRNK  
#cdc1=6, cdc2 = 1
#3,5 = INFLAMMATORY MAC 
#myeloid: 2,5,7,8 --> THERE ARE MACROPHAGES 

cluster.ids = c("Plasma",'Monocyte-derived cDC2b','CytoCD8','Macrophage',
                'Plasma','Macrophage','cDC1','trNK','Monocyte',
                'Plasma','Plasma')
names(cluster.ids) <- levels(Unknown_b)
Unknown_b <- RenameIdents(Unknown_b, cluster.ids)
DimPlot(Unknown_b, reduction = "umap", label = TRUE, repel = TRUE)  

Unknown_b$celltype = Idents(Unknown_b)
Unknown_b@meta.data$celltype %>% table



#v3:
#plasma: 0,4,11,12
#t cell: 2,8 ; 10 = TRNK --> GET RID OF 8 
#cdc1=7, cdc2 = 1; pvm = 9
#myeloid: 1,2,3,5,6,9 --> THERE ARE MACROPHAGES 



#v2:
#plasma: 0,4,13,14 
#t cell: 1,5,11 --> get rid of 10 and 5,8 ; there may be nk cells (12)
#cdc1=7; cdc2: 2,8; pvm = 9


#1, 11, 12 = cytocd8 --> check 1 for doublet, get rid of it 
#1,2,5,8,9,13,14 --> myeloid 
#cdc1 = 7
#0,4,16 = plasma 
#10 = naive b ? --> get rid of 10 and 13 
#the other clusters may not be useful --> check generally or get rid of them 
#look into 10, more like b cell 


DotPlot(Unknown_b, features = c( 'JCHAIN','IGHA1','IGHG1','IGHM', 
                               'MS4A1','BANK1','CD79A','CD79B',"CD19",
                               'MZB1','TOP2A','CDK1','MKI67', #PLASMABLAST
                               'TCL1A','FCER2'), cols = 'RdYlBu')
VlnPlot(Unknown_b, Macro_mono_genes)
VlnPlot(Unknown_b, c(DC_genes,'CLEC12A','FCN1','CD14','LY6G6C','CD1A'))
VlnPlot(Unknown_b, c('LYVE1', 'SELENOP','LGMN', 'STAB1', 'MAN1A1',
                     #'LILRB5', 'MS4A6A',
                     'TNF','IFNG','CCL3','CCL4'
                     ))
VlnPlot(Unknown_b, c(#'KLRF1','PRF1','CD160','NCR1','NCAM1','XCL1',
                     'KLRF1','PRF1','GZMB','CXCR6','CD160','XCL1','CD247', #MATURE 
                     'NCAM1','FCGR3A','ICAM1','KLRK1','NCR1', #IMMATURE 
                     'KLRC2','CD2' )) #ADAPTIVE 

table(Unknown_b$seurat_clusters)

cluster9 = FindMarkers(Unknown_b, ident.1 =9, min.pct = 0.25)
#Unknown_b= subset(Unknown_b, idents = c(0,1,2,3,4,5,6,7,8,9,11,12,14:17))
#Unknown_b= subset(Unknown_b, idents = c(0,1,2,3,4,6,7,9,11,12,13,14))
#Unknown_b= subset(Unknown_b, idents = c(0,1,2,3,4,5,6,7,9,10,11,12))
saveRDS(Unknown_b, 'Unknown_bcell_subcluster.rds')


#T cell unknown ------------

Unknown_T = subset(Unknown, idents = c(0,4,9,10,13:16, 27))
Unknown_T = SeuratAnalysis1(Unknown_T)
ElbowPlot(Unknown_T, ndims = 50)
Unknown_T = SeuratAnalysis2(Unknown_T, harmony = T, 
                            harmony_by = 'orig.ident', n = 10)
DimPlot(Unknown_T,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 5)

FeaturePlot(Unknown_T, features = c( 'JCHAIN','IGHA1','IGHG1','IGHM','IGKC',
                                # 'MS4A1','BANK1','CD79A','CD79B',"CD19",
                                 'CD3D','CD3E',
                                 'CD40LG','CD8A','CD8B'
                                
), cols = c('gray',"red"), reduction = 'umap', label = T)
#contamination with ig genes and jchain : 6,9-12,16,18
table(Unknown_T$seurat_clusters)

#final:  

#cytocd8= 4
#0,3 =  naive cd8  
#cd4=1,2

cluster.ids = c("nCD8",'nCD4','nCD4','nCD8','CytoCD8')
names(cluster.ids) <- levels(Unknown_T)
Unknown_T <- RenameIdents(Unknown_T, cluster.ids)
DimPlot(Unknown_T, reduction = "umap", label = TRUE, repel = TRUE)  

Unknown_T$celltype = Idents(Unknown_T)
Unknown_T@meta.data$celltype %>% table

#6 = tc7+sell+ccr7+ -> tcm? 
#5= nk 

#messy: 4, 6,7 (7 is tcf7+sell+ naive double negative ish)


VlnPlot(Unknown_T,features = c('CD3D','CD3E','CD40LG',
                               'CD8A','CD8B','IFNG', 'ITGB1','GZMB','GZMK')) 
DotPlot(Unknown_T, features = c('CD3D','CD3E','CD40LG',
                                'CD8A','CD8B','IFNG', 'ITGB1','GZMB','GZMK',
                                'CCR7', 'SELL', 'TCF7',  #tcm 
                                'CXCR6','ITGAE','ITGA1', #trm 
                                "LAG3",'PDCD1','HAVCR2',
                                'ISG15','IFIT1','IFIT2', 'IFITM2',
                                'IFNGR1','SLC4A10','TRAV1-2','TRBV6-2',
                                'NCAM1','FCGR3A','NCR1','ICAM1','KLRK1','GNLY','PRF1','KLRF1'
), cols = 'RdYlBu')  + coord_flip()


VlnPlot(Unknown_T, features = c('ISG15','IFIT1','IFIT2', 'IFITM2'))
VlnPlot(Unknown_T, features = c( 'KLRF1','PRF1','GZMB','CXCR6','CD160','XCL1','CD247'))

cluster0 = FindMarkers(Unknown_T, ident.1 = 0, min.pct = 0.25)
#Unknown_T= subset(Unknown_T, idents = c(0:5,7,8,13,14,15,17,19))
#Unknown_T= subset(Unknown_T, idents = c(0,1,2,3,5,7,8,9,10,11,12,13,15))
#Unknown_T= subset(Unknown_T, idents = c(0,1,2,3,4,6,7,9,10,11,12,13,14))
#Unknown_T= subset(Unknown_T, idents = c(1,2,3,4,5,6,7,8,9,11,12,13))
#Unknown_T= subset(Unknown_T, idents = c(0,1,2,3,4,5,6,7,8,10,11))
#Unknown_T= subset(Unknown_T, idents = c(0:9))
#Unknown_T= subset(Unknown_T, idents = c(0,1,2,5,8,9))
#Unknown_T= subset(Unknown_T, idents = c(1:5))


Unknown_T_sub = subset(Unknown_T, idents = c(4,6,7,3))
Unknown_T_sub = SeuratAnalysis1(Unknown_T_sub)
ElbowPlot(Unknown_T_sub, ndims = 50)
Unknown_T_sub = SeuratAnalysis2(Unknown_T_sub, harmony = T, 
                            harmony_by = 'orig.ident', n = 10)
DimPlot(Unknown_T_sub,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 5)

#1,2,4 =cytocd8  
#0,3,5,6,7 = ifn-primed 
cluster.ids = c("IFN-primed T cell",'CytoCD8','CytoCD8',"IFN-primed T cell",
                'CytoCD8',"IFN-primed T cell","IFN-primed T cell","IFN-primed T cell")
names(cluster.ids) <- levels(Unknown_T_sub)
Unknown_T_sub <- RenameIdents(Unknown_T_sub, cluster.ids)
DimPlot(Unknown_T_sub, reduction = "umap", label = TRUE, repel = TRUE)  

Unknown_T_sub$celltype = Idents(Unknown_T_sub)
Unknown_T_sub@meta.data$celltype %>% table


FeaturePlot(Unknown_T_sub, features = c( 'JCHAIN','IGHA1','IGHG1','IGHM','IGKC',
                                     # 'MS4A1','BANK1','CD79A','CD79B',"CD19",
                                     'CD3D','CD3E',
                                     'CD40LG','CD8A','CD8B'
                                     
), cols = c('gray',"red"), reduction = 'umap', label = T)
VlnPlot(Unknown_T_sub,features = c('CD3D','CD3E','CD40LG',
                               'CD8A','CD8B','IFNG', 'ITGB1','GZMB','GZMK')) 
DotPlot(Unknown_T_sub, features = c('CD3D','CD3E','CD40LG',
                                'CD8A','CD8B','IFNG', 'ITGB1','GZMB','GZMK',
                                'CCR7', 'SELL', 'TCF7',  #tcm 
                                'CXCR6','ITGAE','ITGA1', #trm 
                                "LAG3",'PDCD1','HAVCR2',
                                'ISG15','IFIT1','IFIT2', 'IFITM2',
                                'IFNGR1','SLC4A10','TRAV1-2','TRBV6-2',
                                'NCAM1','FCGR3A','NCR1','ICAM1','KLRK1','GNLY','PRF1','KLRF1'
), cols = 'RdYlBu')  + coord_flip()


VlnPlot(Unknown_T_sub, features = c('ISG15','IFIT1','IFIT2', 'IFITM2'))
VlnPlot(Unknown_T_sub, features = c( 'KLRF1','PRF1','GZMB','CXCR6','CD160','XCL1','CD247'))

cluster6 = FindMarkers(Unknown_T_sub, ident.1 = 6, min.pct = 0.25)
Unknown_T_sub= subset(Unknown_T_sub, idents = c(1,2,7))

#NK unknown : --------
Unknown_NK = subset(Unknown, idents = c(3,19,25))
Unknown_NK = SeuratAnalysis1(Unknown_NK)
ElbowPlot(Unknown_NK, ndims = 50)
Unknown_NK = SeuratAnalysis2(Unknown_NK, harmony = T, 
                                harmony_by = 'orig.ident', n = 10)
DimPlot(Unknown_NK,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 5)

Unknown_NK[['percent.rp']] = PercentageFeatureSet(Unknown_NK, pattern = "^RP[SL]")
VlnPlot(Unknown_NK, features = 'percent.rp')

#final
#nk: 2,7 have some ncami and ncr1; 0,2,6,7,8 have high fcgr3a; so they are all mature 
#3 = trNK 
#cd8+: 1,4,5 ; 1 is not as inflammatory.cytotoxic
 
cluster.ids = c("Mature NK", "nCD8","Mature NK",'trNK',
                'CytoCD8','CytoCD8',"Mature NK","Mature NK","Mature NK")
names(cluster.ids) <- levels(Unknown_NK)
Unknown_NK <- RenameIdents(Unknown_NK, cluster.ids)
DimPlot(Unknown_NK, reduction = "umap", label = TRUE, repel = TRUE)  

Unknown_NK$celltype = Idents(Unknown_NK)
Unknown_NK@meta.data$celltype %>% table


#another subclustering find a cluster that is cd8+ T cells but with lots of RP genes 
#0 looks myeloid: thbs1, s100a8, s100a9 
#1 is nkt!! --> I put them in another subcluster 

#11 may be NKT--> get rd 
#cd8+: 0,2,4,6,11 
#nk: 3,5,7,9,10 (10 may be tissue resident)
 

VlnPlot(Unknown_NK,features = c('CD3D','CD3E','CD40LG',
                                   'CD8A','CD8B','IFNG', 'ITGB1','GZMB','GZMK','KLRB1','CD28'))
VlnPlot(Unknown_NK, c(#'KLRF1','PRF1','CD160','NCR1','NCAM1','XCL1',
  'KLRF1','PRF1','GZMB','CXCR6','CD160','XCL1','CD247', #MATURE 
  'NCAM1','FCGR3A','ICAM1','KLRK1','NCR1','AREG', 'CD44', #IMMATURE 
  'KLRC2','CD2'
  )) #ADAPTIVE 
DotPlot(Unknown_NK, c(#'KLRF1','PRF1','CD160','NCR1','NCAM1','XCL1',
  'KLRF1','PRF1','GZMB','CXCR6','CD160','XCL1','CD247', #MATURE 
  'NCAM1','FCGR3A','ICAM1','KLRK1','NCR1','AREG', 'CD44', #IMMATURE 
  'KLRC2','CD2'
), cols= 'RdYlBu') 
VlnPlot(Unknown_NK, c( 'ADGRG1','FCER1G','LAIR2','SPON2'))

cluster4 = FindMarkers(Unknown_NK, ident.1 = 4, min.pct = 0.25)
#Unknown_NK =subset(Unknown_NK, idents = c(0:5, 7:9))
#Unknown_NK = subset(Unknown_NK, idents = c(2:9))
#Unknown_NK = subset(Unknown_NK, idents = c(1:9))
Unknown_NK = subset(Unknown_NK, idents = c(0:8))


#the subcluster of NK does not look good ------

Unknown_NK_sub = subset(Unknown_NK, idents = c(0,1))
Unknown_NK_sub = SeuratAnalysis1(Unknown_NK_sub)
ElbowPlot(Unknown_NK_sub, ndims = 50)
Unknown_NK_sub = SeuratAnalysis2(Unknown_NK_sub, harmony = F, n = 10)
DimPlot(Unknown_NK_sub,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 5)

VlnPlot(Unknown_NK_sub,features = c('CD3D','CD3E','CD40LG',
                                'CD8A','CD8B','IFNG', 'ITGB1','GZMB','GZMK'))
#1,2,5= cytocd8  
#3 has ifng 

#3 = mature nk --> doublet actually 

cluster1 = FindMarkers(Unknown_NK_sub, ident.1 = 1, min.pct = 0.25)
VlnPlot(Unknown_NK_sub, c(#'KLRF1','PRF1','CD160','NCR1','NCAM1','XCL1',
  'KLRF1','PRF1','GZMB','CXCR6','CD160','XCL1','CD247', #MATURE 
  'NCAM1','FCGR3A','ICAM1','KLRK1','NCR1','AREG', 'CD44', #IMMATURE 
  'KLRC2','CD2'
)) #ADAPTIVE 
VlnPlot(Unknown_NK_sub, Macro_mono_genes) #some monocyte mixing 


#Unknown_NK_sub = subset(Unknown_NK_sub, idents = 0:4)
#Unknown_NK_sub = subset(Unknown_NK_sub, idents = c(0:2))


Unknown_NK_sub2 = subset(Unknown_NK, idents = 0)
Unknown_NK_sub2 = SeuratAnalysis1(Unknown_NK_sub2)
ElbowPlot(Unknown_NK_sub2, ndims = 50)
Unknown_NK_sub2 = SeuratAnalysis2(Unknown_NK_sub2, harmony = F, n = 50)
DimPlot(Unknown_NK_sub2,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 5)

VlnPlot(Unknown_NK_sub2, features = c('CD3D','CD3E','CD40LG',
                                      'CD8A','CD8B','IFNG', 'ITGB1','GZMB','GZMK'))
#all clusters seem to have high RP genes 
cluster1 = FindMarkers(Unknown_NK_sub2, ident.1 = 1, min.pct = 0.25)

Unknown_NK_sub2[['percent.rp']] = PercentageFeatureSet(Unknown_NK_sub2, pattern = "^RP[SL]")
VlnPlot(Unknown_NK_sub2, features = 'percent.rp')
Unknown_NK_sub2 = subset(Unknown_NK_sub2, percent.rp <=15) #but then not many cells left, so I think we can just get rid of the whole thing


#merge the unknown lymphoid -----
Unknown_lymphoid = merge(x = Unknown_b, y = c(Unknown_T, Unknown_T_sub, Unknown_NK))
Unknown_lymphoid = SeuratAnalysis1(Unknown_lymphoid)
ElbowPlot(Unknown_lymphoid, ndims = 50)
Unknown_lymphoid = SeuratAnalysis2(Unknown_lymphoid, harmony = T, 
                                   harmony_by='orig.ident', n = 20)
DimPlot(Unknown_lymphoid,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 5, group.by = 'celltype')

#4,5,10-13 are cycling 
VlnPlot(Unknown_lymphoid, features = c('TOP2A','MKI67','STMN1','CENPE'),group.by = 'seurat_clusters' )
Idents(Unknown_lymphoid) ='seurat_clusters'
Unknown_lymphoid = subset(Unknown_lymphoid, idents = c(0:3, 6:9))

VlnPlot(Unknown_lymphoid, features = c('JCHAIN','IGHA1','IGHG1','IGHM'),group.by = 'seurat_clusters' )
Unknown_lymphoid= subset(Unknown_lymphoid, idents = c(0:6, 8:10))
Unknown_lymphoid = merge(Unknown_lymphoid, sub_7)
Unknown_lymphoid = JoinLayers(Unknown_lymphoid)

saveRDS(Unknown_lymphoid, file ='Unknown_lymphoid_merged_0225.rds')
 
#7 is a little weird 
sub_7 = subset(Unknown_lymphoid, idents = 7)
sub_7 = SeuratAnalysis1(sub_7)
ElbowPlot(sub_7, ndims = 50)
sub_7 = SeuratAnalysis2(sub_7, harmony = T, 
                                   harmony_by='orig.ident', n = 20)
DimPlot(sub_7,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 5, group.by = 'orig.ident' )

VlnPlot(sub_7, features = c('CD3D','CD3E','CD40LG',
                            'CD8A','CD8B','IFNG','ITGB1'))
VlnPlot(sub_7, features = 'percent.mt') 

cluster3 = FindMarkers(sub_7, ident.1 = 3, min.pct = 0.25)
#only 3 is good 

cluster.ids = c('Unknown','Unknown','Unknown','CD8+ Tcm')
names(cluster.ids) <- levels(sub_7)
sub_7 <- RenameIdents(sub_7, cluster.ids)
DimPlot(sub_7, reduction = "umap", label = TRUE, repel = TRUE)  

sub_7$celltype = Idents(sub_7)
sub_7@meta.data$celltype %>% table

sub_7 = subset(sub_7, idents = 'CD8+ Tcm')
sub_7 = JoinLayers(sub_7)

#check the ifn-primed T cell, with Ig genes 
IFN_T= subset(Unknown_lymphoid, idents = 'IFN-primed T cell')



#myeloid unknown ---------
Unknown_myeloid = subset(Unknown, idents = c(1,2,6,7,8,12,17,18,24,26,20,21))
Unknown_myeloid = SeuratAnalysis1(Unknown_myeloid)
ElbowPlot(Unknown_myeloid, ndims = 50)
Unknown_myeloid = SeuratAnalysis2(Unknown_myeloid, harmony = T, 
                             harmony_by = 'orig.ident', n = 20)
DimPlot(Unknown_myeloid,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 5)

VlnPlot(Unknown_myeloid, features = Macro_mono_genes)
cluster1 = FindMarkers(Unknown_myeloid, ident.1 = 1, min.pct = 0.25)

#the single 10 cluster is the alveolar cell 


FeaturePlot(Unknown_myeloid, features = Macro_mono_genes, label = T, cols = c('gray',"red"), 
            reduction='umap')
#the right half is macrophage, left half is monocyte 

#check for contamination:
FeaturePlot(Unknown_myeloid, features = c( 'JCHAIN','IGHA1','IGHG1','IGHM','IGKC', 
                                         'CD3D','CD3E',
                                         'CD40LG','CD8A','CD8B'
), cols = c('gray',"red"), reduction = 'umap', label = T)

VlnPlot(Unknown_myeloid, features = c( 'JCHAIN','IGHA1','IGHG1','IGHM','IGKC', 
                                       'CD3D','CD3E',
                                       'CD40LG','CD8A','CD8B'))
#1,5,13 seems to be T cells ; 1 seems to be plasma 
#there are igkc background : 0,1, 9,14,15


#t cell: 
Unknown_myeloid_sub1 = subset(Unknown_myeloid, idents = c(5))
Unknown_myeloid_sub1 = SeuratAnalysis1(Unknown_myeloid_sub1)
ElbowPlot(Unknown_myeloid_sub1, ndims = 50)
Unknown_myeloid_sub1 = SeuratAnalysis2(Unknown_myeloid_sub1, harmony = F,n = 10)
DimPlot(Unknown_myeloid_sub1,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 5)
VlnPlot(Unknown_myeloid_sub1, features = c( #'JCHAIN','IGHA1','IGHG1','IGHM','IGKC', 
                                       'CD3D','CD3E',
                                       'CD40LG','CD8A','CD8B',
                                       'NCAM1','FCGR3A','KLRF1','PRF1'))
VlnPlot(Unknown_myeloid_sub1, Macro_mono_genes)
VlnPlot(Unknown_myeloid_sub1, features = c( "CDKN1C",'CX3CR1','FCGR3A', 'HES4',  #ncMo
                                            'LILRA1','SIGLEC10', 'LILRB1','LILRB2' ,
                                            'S100A12','CD300E','FCN1','EREG','CD14','IL1RN','TCF7L2' , #iMo  
                                            'FCER1A','FCGR2B'  )) #mo-1 
FeaturePlot(Unknown_myeloid_sub1, features = c( 'CD14', 'CCR2', 'CX3CR1','CD86','CD40',
                                                'HLA-DRA',
                                             'S100A12','VCAN','S100A8','S100A9','FCGR3A'
), cols = c('gray','red'),reduction = 'umap', label = T)  

#final: ?? doublet?? background?? not clear?? --> I decided not to include this 
#cd8 t = 1
#nk = 7 (immature)
#ncmo = 6


#v4:
#cd8 t: 3 
#ncmo : 5,8
#1 & 4 IS IMO, HAS CD40 AND CD86, LOWER S100A8/12, SOME FCGR3A 
#others are all monocytes 
 
VlnPlot(Unknown_myeloid_sub1, features= c('APOC1', 'CD9', 'TREM2', 'APOE', 'FABP5',#LAM 
                                          'LYVE1', 'SELENOP','LGMN', 'STAB1', 'MAN1A1', 'LILRB5', 'MS4A6A',  #PVM
                                          'TNF','CCL3L1','CCL4L2','CCL3'))

#v3
#cd8 t: 6,1 --> 1 is a doublet 
#c1q macrophage= 3,7 --> 7 is a doublet 
#3=ncmo 
# all others are iMo 


#v2:
#2,3 --> still weird 
#8 is macrophage with c1q ; 5,6= myeloid; 9 looks like ncmo with fcgr3a and hes4 
#0,5,6,8 have igkc
#junk: 0, 

#cd8 t cell = 2,4,8
#9 and 11 is plasma or junk 
cluster8 = FindMarkers(Unknown_myeloid_sub1, ident.1 = 8, min.pct = 0.25)


#Unknown_myeloid_sub1 = subset(Unknown_myeloid_sub1, idents = c(0:8, 10))
#Unknown_myeloid_sub1 = subset(Unknown_myeloid_sub1, idents = c(1:10))
#Unknown_myeloid_sub1 = subset(Unknown_myeloid_sub1, idents = c(0,2,3,4,5,6,8))
#Unknown_myeloid_sub1 = subset(Unknown_myeloid_sub1, idents = c(1:9))  #0 has EIF and cdc genes 
Unknown_myeloid_sub1 = subset(Unknown_myeloid_sub1, idents = c(0:7))


Unknown_myeloid_sub2 = subset(Unknown_myeloid, idents = c(0,1,20,13,17,16,12,10))
Unknown_myeloid_sub2 = SeuratAnalysis1(Unknown_myeloid_sub2)
ElbowPlot(Unknown_myeloid_sub2, ndims = 50)
Unknown_myeloid_sub2 = SeuratAnalysis2(Unknown_myeloid_sub2, harmony = F, n = 20)
DimPlot(Unknown_myeloid_sub2,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 5 )

cluster6 = FindMarkers(Unknown_myeloid_sub2, ident.1 = 6, min.pct = 0.25)

#final:
#3,5,7= monocyte; 7 has highest fcn1 and s100a8/9 
# 5 HAS cxcl9 and 10 and ccl3; 3 IS IMO (HIGHEST HLA-DRA, SIGLEC10, LOWER S100A8/9); 7 HAS LGALS2 --> 7 IS TIM 
#5 AND 9 HAVE CXCL10 AND 9 --> 9 IS INFLAMMATORY MAC 
#TREM2+ : 0,1,2,8 


#after subset out the unknown I check:
#10 is plasma, not TIM --> doublet  
#6 IS tim 
#chck 2,3,8
#4 5 HAVE BAD GENES IN IT 
#LAM AND INFLAMMATORY MAC ARE FINE 
#GET RID OF 3 

#Unknown_myeloid_sub2 = subset(Unknown_myeloid_sub2, idents = c(0,1,2,4:9, 11))
#Unknown_myeloid_sub2 = subset(Unknown_myeloid_sub2, idents = c(0:7))

#v6
#6,7,10, 11= MONO; 11 IS TIM 
#8,9 HAS FIBRO GENES BUT ALSO SOME MAYBE PERIVASCULAR  --> PVM  
#INFLAMMATORY: 12 
#trem2+: 1,2,3,6-9, 11-13
# 0,3,4,5 ?? --> a mixture of junk and monocyte/DC2 

cluster.ids = c( 'LAM','LAM','LAM','iMo','PVM','IFN-stimulated monocytes','PVM',
                 'Tumor-infiltrating monocytes','LAM',
                'Inflamamtory macrophage')
names(cluster.ids) <- levels(Unknown_myeloid_sub2)
Unknown_myeloid_sub2 <- RenameIdents(Unknown_myeloid_sub2, cluster.ids)
DimPlot(Unknown_myeloid_sub2, reduction = "umap", label = TRUE, repel = TRUE)  

Unknown_myeloid_sub2$celltype = Idents(Unknown_myeloid_sub2)
Unknown_myeloid_sub2@meta.data$celltype %>% table



#v5:
#5 and 10 are monocytes ; 6 is imo 
#inflamamtory mac/ifn-primed: 11 ,0(ifn) 
#trem2+: 0, 5, 7,8,11 
#t cell genes in 4 (klrb1, klrd1, gnly)--> cannot determine the type so get rid of it 
#3 is a mix of club cell and at2


#V4
#13= IMO, MAYBE 6,7 IS NCMO ; 6 is mono
#14 is inflammatory mac 
#trem2+ : 0,7,8,10,12,14 --> so 0, 8,10, 12, they have high expression of others 
#3 has t cells 

#v3:
#14 = iMo
#macrophages all others ; -- seems to be all trem2+ lams ?? not inflammatory no TNF
#8 and 11 are structural 
#0,6,7,8,12,15 have trem2 
#15 IS INFLAMMATORY MAC WITH CCL 

#8 and 10 are the more lung structural cells 
#17 seems to be nk? or in general a lymphoid 

#Unknown_myeloid_sub2= subset(Unknown_myeloid_sub2, idents = c(0:10,12,14,16,17))
#Unknown_myeloid_sub2= subset(Unknown_myeloid_sub2, idents = c(0,1,2,3,6:13))
#Unknown_myeloid_sub2= subset(Unknown_myeloid_sub2, idents = c(0:10,13:15))
#Unknown_myeloid_sub2= subset(Unknown_myeloid_sub2, idents = c(0:4, 6:14))  
#Unknown_myeloid_sub2= subset(Unknown_myeloid_sub2, idents = c(0,1,2,4:8,10,11)) 
#Unknown_myeloid_sub2= subset(Unknown_myeloid_sub2, idents = c(0,1,2,4:12)) 

VlnPlot(Unknown_myeloid_sub2, features = c('SCGB1A1','SCGB3A2', AT1_gene))
VlnPlot(Unknown_myeloid_sub2, features = c(# 'JCHAIN','IGHA1','IGHG1','IGHM','IGKC', 
  'CD3D','CD3E',
  'CD40LG','CD8A','CD8B','KLRD1' ,
  'ZNF683', 'TRDC', 'TRDV1', 'TRGC1','KIR2DL4'
  #'NCAM1','FCGR3A','KLRF1','PRF1'
  ))
FeaturePlot(Unknown_myeloid_sub2, Macro_mono_genes, cols = c('gray','red'),
            reduction = 'umap', label = T)
VlnPlot(Unknown_myeloid_sub2, features = c("ISG15",'IFITM1','IFITM2','IFI27','CD163','MARCO','C1QA','C1QB'))
VlnPlot(Unknown_myeloid_sub2, features = c('CLEC10A','CD1C','ITGAM', 'ITGAX','IRF4','ZBTB46'))

VlnPlot(Unknown_myeloid_sub2, features = c( 'CD14', 'CCR2', 'IL1B','FCN1','FOS',
                                           # 'S100A12','EREG','CD14','VCAN','S100A8','S100A9','FCGR3A',
                                           # 'CCL3','CCL4','CCR5','TNF','CXCL9','CXCL10','ISG15','STAT1','IRF1' #inflammatory signatures 
                                            'LGALS2','CEBPD','CD36','IFITM2','VCAN', 'LYZ','SELL',
                                           'EGR1','NFKB1','CCL4','CXCR4','AREG'   #TUMOR INFILTRATING MONOCYTE
                                            ))  

VlnPlot(Unknown_myeloid_sub2, features = c('LUM','DCN','PTPRC','ADGRE1','CD68', 'FCGR3A'))

VlnPlot(Unknown_myeloid_sub2, features = c( #"CDKN1C",'CX3CR1','FCGR3A', 'HES4',  #ncMo
                                            #'LILRA1','SIGLEC10', 'LILRB1','LILRB2' ,
                                            'S100A12','CD300E','FCN1','EREG', 
                                            'IL1RN','TCF7L2' , #iMo  
                                            'S100A8','S100A9','CCR2','CD14','FCGR3A','HLA-DRA')) #classical monocytes 


VlnPlot(Unknown_myeloid_sub2, features= c('APOC1', 'CD9', 'TREM2', 'APOE', 'FABP5',
                                          'TNF','CCL3L1','CCL4L2','CCL3','CXCL10','CXCL9',
                                          'LYVE1', 'SELENOP','LGMN', 'STAB1' ,'MAN1A1', 'LILRB5', 'MS4A6A'
                                          #'LUM','DCN'
                                          ))
VlnPlot(Unknown_myeloid_sub2, c(Neutrophil_gene, 'CXCR4','CXCR2','ICAM1','ELANE','MPO')) 



Unknown_myeloid_sub3 = subset(Unknown_myeloid, idents = c(7,4,2,9,18,8,19,11,14,15,3,21,2,6))
Unknown_myeloid_sub3 = SeuratAnalysis1(Unknown_myeloid_sub3)
ElbowPlot(Unknown_myeloid_sub3, ndims = 50)
Unknown_myeloid_sub3 = SeuratAnalysis2(Unknown_myeloid_sub3, harmony = T, 
                                       harmony_by = 'orig.ident', n = 20)
DimPlot(Unknown_myeloid_sub3,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 5, group.by = 'seurat_clusters')

cluster5 = FindMarkers(Unknown_myeloid_sub3, ident.1 = 5, min.pct = 0.25)


cluster.ids = c('ncMo','cMo','ncMo', 'cMo','ncMo','Inflammatory macrophage','cMo',
                'cMo','Inflammatory macrophage','Inflammatory macrophage', 'cMo', 'Inflammatory macrophage',
                'ncMo','ncMo','cMo','ncMo')
names(cluster.ids) <- levels(Unknown_myeloid_sub3)
Unknown_myeloid_sub3 <- RenameIdents(Unknown_myeloid_sub3, cluster.ids)
DimPlot(Unknown_myeloid_sub3, reduction = "umap", label = TRUE, repel = TRUE)  

Unknown_myeloid_sub3$celltype = Idents(Unknown_myeloid_sub3)
Unknown_myeloid_sub3@meta.data$celltype %>% table

#initial: 9 is a doublet, get rid of it first 

#final:
#ncmo: 0,12,13,15,2,4
#cmo: 6,7,10,14,1,3
#macrophages with ccl3 and ccl4 family: 5,9,11 8 

#v2
#5 maybe macrophage 
#ncmo: 0,13,12 
#imo: 6, 11
#classical monocyte: 10,8,7,1,5 (THE LEFT HALF)
#15,16 = JUNK WITH H3/RP GENES 


#v1:
#ncmo: 0,12,13,15,18,5,6
#imo: 2, 16 (the upper half)
#tim: 1,4,9
#GET RID OF 17, 19
#CLUSTER 14 HAS CLEC10A, CD14, SELL, EGR2, ITGAE 

VlnPlot(Unknown_myeloid_sub3, features = c( 'JCHAIN','IGHA1','IGHG1','IGHM','IGKC', 
  'CD3D','CD3E',
  #'CD40LG','CD8A','CD8B','IFNG'  ,
  Macro_mono_genes
))

DotPlot(Unknown_myeloid_sub3, c(Macro_mono_genes, 'VCAN','CCR2' ,'CD300E','FCN1'), cols = 'RdYlBu')
VlnPlot(Unknown_myeloid_sub3, features = c('TNF','CCL3L1','CCL4L2','CCL3','CXCL10','CXCL9',
                                           Macro_mono_genes), group.by = 'seurat_clusters')
FeaturePlot(Unknown_myeloid_sub3, features = c( "CDKN1C",'CX3CR1','FCGR3A', 'HES4',  #ncMo
                                           # 'LILRA1','SIGLEC10', 'LILRB1','LILRB2' ,
                                            'S100A12','CCR2','CD14','VCAN','S100A8','S100A9', #CLASSICAL MO HAS HIGH OF THESE
                                            # 'CD300E','FCN1',#GENERAL MONOCYTE 
                                            #'IL1RN','TCF7L2','CCL2','CXCL10','IFITM2' , #IMO GENES FROM THE PAPER
                                            'HLA-DRA','HLA-DPB1' ,'CD74','CD86','CD40'
                                            #'FCER1A','FCGR2B'  
                                            ), cols = c('gray','red'),
        reduction = 'umap', label = T)  
FeaturePlot(Unknown_myeloid_sub3, features = c('CCL3','CCL3L3','CCL4','CCL4L2','CD163','MARCO','TNF','ISG15'), cols = c('gray','red'),
        reduction = 'umap', label = T)
FeaturePlot(Unknown_myeloid_sub3, features = c('CD163','MARCO','C1QA', 'CD68',
                                               'CD14', 'CCR2', 'IL1B','FCN1'
                                            #'FOS','LGALS2','LYZ','SELL' 
                                            ), cols = c('gray','red'),
            reduction = 'umap', label = T)

FeaturePlot(Unknown_myeloid_sub3, features= c('CD14','CD163',#'MARCO','TREM2', 'APOE', 'FABP5',
   'TNF','CCL3L1','CCL4L2','CCL3','CXCL10','CXCL9'
), cols = c('gray','red'),reduction = 'umap', label = T) #not really expressed 

VlnPlot(Unknown_myeloid_sub3, features = DC_genes)

#Unknown_myeloid_sub3 = subset(Unknown_myeloid_sub3, idents = c(0:8, 10:19))
#Unknown_myeloid_sub3 = subset(Unknown_myeloid_sub3, idents = c(0:16, 18))
#Unknown_myeloid_sub3 = subset(Unknown_myeloid_sub3, idents = c(0:5, 8:11, 13:17))
#Unknown_myeloid_sub3 = subset(Unknown_myeloid_sub3, idents = c(0:14))


Unknown_structural= subset(Unknown_myeloid_sub2, idents = c(3,9))
Unknown_structural = SeuratAnalysis1(Unknown_structural)
ElbowPlot(Unknown_structural, ndims = 50)
Unknown_structural = SeuratAnalysis2(Unknown_structural, harmony = T, 
                                       harmony_by = 'orig.ident', n = 15)
DimPlot(Unknown_structural,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 5)
cluster10 = FindMarkers(Unknown_structural, ident.1 = 10, min.pct = 0.25)

FeaturePlot(Unknown_structural, Macro_mono_genes, cols = c('gray','red'),
            reduction = 'umap', label = T)
FeaturePlot(Unknown_structural, features = c('S100A12','CCR2','CD14','VCAN','S100A8','S100A9'), cols = c('gray','red'),
            reduction = 'umap', label = T)
VlnPlot(Unknown_structural, c("SCGB1A1",'SCGB3A1')) 
VlnPlot(Unknown_structural, AT2_gene) 
VlnPlot(Unknown_structural,c(Endo_genes,  'ACKR1', 'LYVE1','SELENOP')) 

#v2:
#1,2 = club, 4 = at2 ((maybe 1 and 2)); 5 = fibro, 6 = ackr1+ endo 

cluster.ids =c("Unknown",'Club cell','Club cell','Unknown', 'AT2','Fibro','Vein endothelial',
               'Unknown','Unknown')
names(cluster.ids) <- levels(Unknown_structural)
Unknown_structural <- RenameIdents(Unknown_structural, cluster.ids)
DimPlot(Unknown_structural, reduction = "umap", label = TRUE, repel = TRUE)  

Unknown_structural$celltype = Idents(Unknown_structural)
Unknown_structural@meta.data$celltype %>% table


#v1:
#0 = macrophage ; 9 is monocyte; 6 is also some kind of immune 
#1,2,7 = junk 
#5 endo but inflammatory ; 8 = fibro 
#4 and 10 = at2 ; 3 = CLUB 
Unknown_structural_markers = FindAllMarkers(Unknown_structural, min.pct = 0.25, only.pos = T)

Unknown_structural_markers_filtered = Unknown_structural_markers%>% 
  group_by(cluster) %>% 
  slice_max(order_by = avg_log2FC, n = 50, with_ties = T) 

#Unknown_structural = subset(Unknown_structural, idents = c(0,3,4,5,6,8,9,10))
#Unknown_structural = subset(Unknown_structural, idents = c(1:7))


Unknown_monocyte = subset(Unknown_myeloid_sub2, idents = 'Unknown')
Unknown_monocyte = SeuratAnalysis1(Unknown_monocyte)
ElbowPlot(Unknown_monocyte, ndims = 50)
Unknown_monocyte = SeuratAnalysis2(Unknown_monocyte, harmony = T, 
                                     harmony_by = 'orig.ident', n = 10)
DimPlot(Unknown_monocyte,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 5)

VlnPlot(Unknown_monocyte, features = c("ISG15",'IFITM1','IFITM2','IFI27','MARCO','C1QA','C1QB',
                                       'FCN1','VCAN','CD300E'),group.by = 'seurat_clusters')

VlnPlot(Unknown_monocyte, features = c('CXCR4','EGR1','CCR2','NFKB1','SELL','THBS1','MPEG1','LGALS2'))

FeaturePlot(Unknown_monocyte, features = c('S100A12','CCR2','CD14','VCAN','S100A8','S100A9'), cols = c('gray','red'),
            reduction = 'umap', label = T)
VlnPlot(Unknown_monocyte, features = c( 'CD14','CD163','MARCO',
                                        'TREM2', 'APOE', 'FABP5',
                                        'LYVE1','SELENOP','STAB1','LGMN','MS4A6A',
                                        'TNF','CCL3L1','CCL4L2','CCL3','CCL4',
                                        'CXCL10','CXCL9','STAT1','ISG15') )
FeaturePlot(Unknown_monocyte, features = c( "CDKN1C",'CX3CR1','FCGR3A', 'HES4',  #ncMo
                                                'S100A12','CCR2','CD14','VCAN','S100A8','S100A9', #CLASSICAL MO 
                                                'HLA-DRA','HLA-DPB1' ,'CD74','CD86','CD40'
), cols = c('gray','red'), reduction = 'umap', label = T)  
FeaturePlot(Unknown_monocyte, features = c('FABP4', 'MCEMP1', 'PPARG','C1QA','C1QB'), cols = c('gray','red'), reduction = 'umap', label = T)
VlnPlot(Unknown_monocyte, features = c("TPSAB1",'KIT','CPA3'))

markers = FindAllMarkers(Unknown_monocyte, min.pct = 0.25, only.pos = T)

Unknown_monocyte_markers_filtered = markers %>% 
  group_by(cluster) %>% 
  slice_max(order_by = avg_log2FC, n = 50, with_ties = T) 

cluster3 = FindMarkers(Unknown_monocyte, ident.1 = 3, min.pct = 0.25)
#V2: 6 IS MAST; 7 HAS CXCL9, cxcl10, stat1 
#0 IS MONOCYTE(INFLAMMATORY, ifitm genes and nfkb1a); others macrophage  

cluster.ids = c("IFN-stimulated monocyte", 'Macrophage','Macrophage','Macrophage','Macrophage',
                'Macrophage','Mast cell', 'IFN-primed Macrophage')
names(cluster.ids) <- levels(Unknown_monocyte)
Unknown_monocyte <- RenameIdents(Unknown_monocyte, cluster.ids)
DimPlot(Unknown_monocyte, reduction = "umap", label = TRUE, repel = TRUE)  

Unknown_monocyte$celltype = Idents(Unknown_monocyte)
Unknown_monocyte@meta.data$celltype %>% table

#v1: only 0 is confusing; 7 is mast cell 

#Unknown_monocyte = subset(Unknown_monocyte, idents = 1:7)


#merge the myeloid within unknown + merge with lymphid all together -----
Unknown_myeloid_sub2 = subset(Unknown_myeloid_sub2, idents = c('iMo','Tumor-infiltrating monocytes','IFN-stimulated monocytes',
                                                               'Inflamamtory macrophage','LAM'))
Unknown_myeloid_sub2 = JoinLayers(Unknown_myeloid_sub2)

Unknown_structural = subset(Unknown_structural, idents = c("Club cell",'AT2','Vein endothelial','Fibro'))
Unknown_structural = JoinLayers(Unknown_structural)

Unknown_myeloid_all = merge(x = Unknown_myeloid_sub2, 
                            y = c(Unknown_myeloid_sub3, Unknown_monocyte, Unknown_structural))
Unknown_myeloid_all = JoinLayers(Unknown_myeloid_all)
saveRDS(Unknown_myeloid_all, file = 'Unknown_myeloid_merge_0225.rds')

Idents(Unknown_lymphoid)='celltype' #get rid of the single monocyte from lymphoid as I double check it is not
Unknown_lymphoid = subset(Unknown_lymphoid, idents = c('cDC1','IFN-primed T cell','CytoCD8',
                                                       'Macrophage','Mature NK','Monocyte-derived cDC2b',
                                                       'nCD4','nCD8','Plasma','trNK'))
Unknown_lymphoid =  JoinLayers(Unknown_lymphoid)

Unknown_final= merge(Unknown_lymphoid, Unknown_myeloid_all)
Unknown_final = SeuratAnalysis1(Unknown_final)
ElbowPlot(Unknown_final, ndims = 50)
Unknown_final = SeuratAnalysis2(Unknown_final, harmony = T, 
                                   harmony_by = 'orig.ident', n =20)
DimPlot(Unknown_final,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'celltype')
table(Unknown_final$celltype)

Idents(Unknown_final)='seurat_clusters'
cluster8 = FindMarkers(Unknown_final, ident.1 = 8, min.pct = 0.25)


Idents(Unknown_final)='celltype'
Unknown_final= RenameIdents(Unknown_final, 'IFN-stimulated monocytes'= 'IFN-primed monocyte', 
                            'IFN-stimulated monocyte' ='IFN-primed monocyte',
                            'Inflamamtory macrophage' = 'Inflammatory macrophage') 
Unknown_final$celltype = Idents(Unknown_final)

Idents(Unknown_final)='celltype'
Unknown_final= RenameIdents(Unknown_final,  'IFN-primed monocyte' ='IFN-primed Macrophage' ) 
Unknown_final$celltype = Idents(Unknown_final)

Unknown_final = subset(Unknown_final, idents = c(1:15, 17:24))

sub_16 = subset(Unknown_final, idents=16) #rename it to alveolar macrophage
sub_16$celltype %>% unique()
sub_16$celltype = 'Alveolar macrophage'

Unknown_final= merge(Unknown_final, sub_16)

#Unknown_final = subset(Unknown_final, idents = unique(Unknown_final$celltype)[c(1:10, 12:20)] )
#Unknown_final = subset(Unknown_final, idents = unique(Unknown_final$celltype)[2:21])

#get rid of plasma  
#look into the iMo close to lam, and the general macrophage 
Idents(Unknown_final) ='seurat_clusters'
cluster16 = FindMarkers(Unknown_final, ident.1 =16, min.pct=0.25, only.pos = T)

markers = FindAllMarkers(Unknown_final, min.pct = 0.25, only.pos = T)

Unknown_final_markers = markers %>% 
  group_by(cluster) %>% 
  slice_max(order_by = avg_log2FC, n = 50, with_ties = T) 

saveRDS(Unknown_final, file = 'All_Lymphoid_remaining_clusters_v2_0225.rds')

VlnPlot(Unknown_final, features = c("LUM",'DCN'))
VlnPlot(Unknown_final, features = c("ISG15",'IFITM1','IFITM2','IFI27','STAT1','IRF1','CXCL9','CXCL10'
                                    #'MARCO','C1QA','C1QB', 'FABP4', 'MCEMP1', 'PPARG','CCL18'
                                    #'TNF','CXCL9','CXCL10','CCL4','CCL3','CCL3L1'
                                    ), group.by = 'celltype') #16 is not ifn-primed
#16 is alveolar macrophage, not ifn-primed 

#monocyte derived cdc2b & PLASMA HAVE top2a and cycling markers --> GET RID OF  

#ifn-primed monocyte population has trem2 and c1q 
#lam is not clear in trem2 ; pvm is not very clear 
#TIM has ig contamination 
#SO LOOK BACK AT MYELOID SUB 2

Unknown_dotplot = DotPlot(Unknown_final, features = c('CD163', 'MARCO','C1QA','C1QB', 'FABP4', 'MCEMP1', 'PPARG',
                                    'TREM2','APOE','APOC1','CD9',
                                    'TNF','CCL3L1','CCL4L2','CCL3','CCL4',
                                    'CXCL10','CXCL9','STAT1','ISG15','IFITM2','IFI27',
                                    'FCN1',
                                    "CDKN1C",'CX3CR1','FCGR3A', 'HES4',  #ncMo
                                   'S100A12','CCR2','CD14','VCAN','S100A8','S100A9', #CLASSICAL MO 
                                   'HLA-DRA','HLA-DPB1' , 
                                   'LGALS2','CEBPD','CD36','EGR1','NFKB1','AREG',
                                   'LUM','DCN','PECAM1','VWF','ACKR1',
                                   'SCGB3A1','SCGB1A1', AT2_gene,
                                   'JCHAIN','IGHA1','IGHG1','IGHM', 
                                   'CD3D','CD3E',
                                   'CD40LG','CD8A','CD8B',
                                   'CCR7','SELL',
                                   'KLRF1','PRF1','CD160','XCL1',
                                   'NCAM1', 'ADGRG1','FCER1G','LAIR2','SPON2',
                                   'TPSAB1','KIT'), group.by = 'celltype', cols = 'RdYlBu') + RotatedAxis() + coord_flip() + theme(axis.text.x = element_text(size = 8),
                                                                                                                                   axis.text.y = element_text(size = 8))

behind = Unknown_dotplot$data
behind = behind %>% filter(avg.exp.scaled > 1)
Unknown_dotplot$data= behind
Unknown_dotplot


#myeloid main clustering ------------
total$seurat_clusters %>% table()
Myeloid = subset(total, idents= c(4,28,36,12,18,35,29,5))
Myeloid = SeuratAnalysis1(Myeloid)
ElbowPlot(Myeloid, ndims = 50)
Myeloid = SeuratAnalysis2(Myeloid, harmony = T, 
                                harmony_by = 'orig.ident', n =15)
DimPlot(Myeloid,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'seurat_clusters')

saveRDS(Myeloid, file ='Myeloid_raw.rds')
#do a bunch of feature plots on these myeloid 

FeaturePlot(Myeloid, 
                 features = B_cell_genes, reduction='umap',
                 raster=F, label = T,  cols = c('grey','red'))
#not seen
p1 = FeaturePlot(Myeloid, 
                 features =T_NK_cell_genes, reduction='umap',
                 raster=F, label = T,  cols = c('grey','red'))

#11 and part 15 
ggsave('FeaturePlot_T_NK.pdf', plot = p1, width = 20, height = 20)


p2 = FeaturePlot(Myeloid, 
                 features =Macro_mono_genes,reduction='umap',
                 raster=F, label = T,  cols = c('grey','red'))

# 2,14,26,8,5,12,3,20,13,16
ggsave('FeaturePlot_macrophage_monocyte.pdf', plot =p2, width = 20, height = 20)

p3 = FeaturePlot(Myeloid, 
                 features =DC_genes, reduction='umap',
                 raster=F, label = T,  cols = c('grey','red'))

#9,21,24 = cdc1, 4 and 10 may be cdc2 

ggsave('FeaturePlot_DC.pdf', plot =p3, width = 20, height = 20)

FeaturePlot(Myeloid, 
                 features = 'EPCAM',  reduction='umap',
                 raster=F, label = T,  cols = c('grey','red'))

#6,7,17,23

FeaturePlot(Myeloid, 
                 features = Endo_genes,  reduction='umap',
                 raster=F, label = T,  cols = c('grey','red'))

#

FeaturePlot(Myeloid, 
                  features = fibro_genes,  reduction='umap',
                  raster=F, label = T,  cols = c('grey','red'))
#

FeaturePlot(Myeloid, 
                  features = SMC_genes,  reduction='umap',
                  raster=F, label = T,  cols = c('grey','red'))

# SMC: 8, 27, 37

#how I subcluster inside Myeloid: ------
#DC: 9,21,24 = cdc1, 4 =cdc2 
#macrophages and monocytes: 2,14,26,8,5,12,3,20,16, 19 , 22
#T and NK and B: 11 and part 15 , 1, 0, 25,18,
#EPCAM and structural: 6,7,17,20 
#unknown remaining: 10, 13, 23

#26 clusters in total 
cluster25 = FindMarkers(Myeloid, ident.1 = 25, min.pct = 0.25)

DotPlot(Myeloid, features = c(T_NK_cell_genes, Macro_mono_genes, DC_genes, 'EPCAM',
                              B_cell_genes), cols = 'RdYlBu')+ coord_flip() + RotatedAxis()

#cdc1 subcluster ------------
DC = subset(Myeloid, idents = c(9,21,24))
DC = SeuratAnalysis1(DC)
ElbowPlot(DC, ndims = 50)
DC = SeuratAnalysis2(DC, harmony = T, 
                          harmony_by = 'orig.ident', n = 15)
DimPlot(DC,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4)

VlnPlot(DC, features = c(DC_genes, 'SIRPA',
                         'MPEG1','WDFY4','BATF3','IRF8',
                         'CCR7','IL12B','FSCN1',
                         'CD274','PDCD1LG2', 
                         'LAMP3','CD40','CD80',
                         'RELB','NFKB1','STAT1',
                         'CTLA4','ICOSLG',
                         'IDO1')) #NO ACTIVATED, MAYBE 6 HAS IL12B
VlnPlot(DC, features = c("CD3D",'CD3E','CD8A','IFNG','GZMB','CD40LG','ITGB1')) 
VlnPlot(DC, features = c(Macro_mono_genes,'CTLA4','ICOSLG',
                         'IDO1', 'EBI3', 'PVR', 'DLL4','CCL22'))

cluster5 = FindMarkers(DC, ident.1 = 5, min.pct = 0.25)

FeaturePlot(DC, 
            features =DC_genes, reduction='umap',
            raster=F, label = T,  cols = c('grey','red'))

cluster.ids = c("cDC1",'cDC1','CD4+CD8+ T cell','cDC1','cDC1',
                'CytoCD8','cDC1','Activated cDC1','Activated cDC1'
                )

names(cluster.ids) <- levels(DC)
DC <- RenameIdents(DC, cluster.ids)
DimPlot(DC, reduction = "umap", label = TRUE, repel = TRUE)  

DC$celltype = Idents(DC)
DC@meta.data$celltype %>% table


#FINAL
#ACTIVATED CDC1: 7,8 WITH  LAMP3, NFKB, RELB, STAT1 
#CD8 T CELL =2,5 
#5 HAS LAMP3, NFKB, RELB, STAT1 
 
#GET RID OF 2,14,11 ; CHECK THE ISOLATED CLUSTERS (16, 14, 8, 11)
#general: a lot of cdc1 here 
#8 is cdc2/mono/macrophage? except for 1,2,14,15 
#there is a ifn-primed population 

table(DC$seurat_clusters)

#DC= subset(DC, idents = c(0,1,3:10, 12,13,15))
#DC= subset(DC, idents = c(0:7, 9:12))
#DC= subset(DC, idents = c(0:6, 8,9, 12))
#DC= subset(DC, idents = c(0, 2:10))
#DC= subset(DC, idents = c(0,1,2,3,4,5,9))
#DC= subset(DC, idents = c(0,2,4,5,6,7))

saveRDS(DC, 'DC1_subcluster.rds')

#cdc2 subcluster ------------
#DC2 = subset(Myeloid, idents = c(4))
DC2_all = SeuratAnalysis1(DC2_all)
ElbowPlot(DC2_all, ndims = 50)
DC2_all = SeuratAnalysis2(DC2_all, harmony = T, 
                     harmony_by = 'orig.ident', n = 25)
DimPlot(DC2_all,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'seurat_clusters')


#check if langerhans cells in hnscc--> seems to be a mix OF LUNG, NOT HNSCC 
#3 and 12 

DC2_all$cell_condition = paste0(DC2_all$celltype, "_", DC2_all$orig.ident)
table(DC2_all$cell_condition)


Idents(DC2_all)='celltype'
inflammatory_DC2_genes = FindMarkers(DC2_all, ident.1 = 'Inflammatory cDC2b', min.pct = 0.25)
DC2_Tem_genes = FindMarkers(DC2_all, ident.1 = 'Tem', min.pct = 0.25)



VlnPlot(DC2_all, features = c("FLT3",'CD14', 'CX3CR1','CCR2','ITGAX', #NOT DEPENDENT ON CX3CR1 OR CCR2 
                              'ZBTB46','LY6G6C','IRF4','ID2',
                              'FCGR3A','S100A8','S100A9')) 
#TRYING TO FIND MONO-DERIVED CDC2B: still express zbtb46  
#1 and 11 may be mono-derived 
VlnPlot(DC2_all, features = c("FLT3",'FCGR1A','C5AR1','MS4A3', #MORE MONOCYTE GENES 
                              'ITGAM','ITGAX'))
VlnPlot(DC2_all, features = c('SIRPA','DPP4','IRF4'), group.by = 'celltype') #inflammatory cdc2 markers from a paper
FeaturePlot(DC2_all, features = c(DC_genes, 'ZBTB46','SIRPA',
                         #'MPEG1','WDFY4','BATF3','IRF8',
                         #'CCR7','IL12B','FSCN1', 
                       # 'LAMP3','ESAM',
                        'CD207','CD1A',  'EPCAM','KRT19', 'ITGAE', 'HLA-DRA',#MAYBE WILL HAVE EPCAM 
                       'IL22RA2'
                       #'CD209','ISG15' #cd209 is dc-sign 
                      ), reduction ='umap', cols = c('gray', 'red'))  

#3 has no clec10a, 6 has high cd207 and cd1a 

VlnPlot(DC2_all, features = c( 'CCR7','FSCN1','LAMP3',
                           'CD274','PDCD1LG2', 'CD40','CD80', 
                            'RELB','NFKB1','STAT1',
                           'CTLA4','ICOSLG', 'IDO1',
                           'IL1B','MPEG1','SIRPA'),group.by = 'seurat_clusters')
VlnPlot(DC2_all, features = c( 'CCR7','FSCN1','LAMP3',
                           'CD274','PDCD1LG2', 'CD40','CD80', 
                           'RELB','NFKB1','STAT1'))#a paper describes the celltypes expressing il22bp have cd4 and sirpa
VlnPlot(DC2_all, features = c( 'CCR7','PDCD1LG2','FSCN1','LAMP3','CD274',
                               'CD86','CD40' ,'CD83',  'EBI3', 'IDO1','CCL22',
                               'IL12B','MPEG1','SIRPA'
                               #'CD163','MARCO','CD68','CD14',
                               #'C1QA','C1QB','TNF','CCL3','CCL4'
                               )) 

DotPlot(DC2_all,features = c(DC_genes, 'ZBTB46','SIRPA',
                            'CCR7','FSCN1', 'CD274','CD86','CD80','LAMP3', 'EBI3', 'IDO1','CCL22',
                             'CXCL8','CCL3L1', 'CCL4L2', 
                             'CD207','CD1A', 'IL22RA2',
                             'CD163','MARCO','CD68','CD14',
                             'C1QA','C1QB','TNF','CCL3','CCL4',
                            'TREM2', 'APOC1','APOE',
                            'CD3D','CD3E','GZMB','GZMM','KLRD1','RORA','CD8A','CD40LG'), group.by = 'celltype', cols = 'RdYlBu')+
  coord_flip() + RotatedAxis()


DC2_rest = subset(DC2_all, idents = c(0:8, 10, 12))
DC2_rename = subset(DC2_all, idents = c(9,11))
Idents(DC2_rename)='celltype'
DC2_rename = RenameIdents(DC2_rename, 'LAMP3+ DC' ='CCR7+ mature DC', 'cDC2b' ='Inflammatory cDC2b')
DC2_rename$celltype = Idents(DC2_rename)
DC2_all = merge(DC2_rest, DC2_rename)


Idents(DC2_all) ='celltype'
DC2_all = subset(DC2_all, idents=unique(DC2_all$celltype)[-9])
DC2_all = RenameIdents(DC2_all, 
                       'CCR7+ mature DC' ='Activated cDC2b')
DC2_all = RenameIdents(DC2_all, 
                       'Macrophage' ='Interstitial macrophages', 'Langerhans cells' = 'Langerin+ cDC2b')
DC2_all$celltype = Idents(DC2_all)
table(DC2_all$celltype)


saveRDS(DC2_all, file = 'cDC2_subcluster_v3.rds')
#ccr7: 5, 14, 16
#0 has high ctla4 and icoslg 
#sirpa level: 1,2 and 9, 10 and 13-15 not very high 

#v2: 
#6 IS A T CELL WITH BOTH CD8 AND CD4 
#4 has ccr7 
#11, 13 have lamp3, maybe cdc2a or langerhans 
#0 and 12 not cdc2--> 12 is langerhans!! 
#there are il22ra2; there are ifn-related genes 

#V3: 
#12,14= ISG15; 7,8,9 = IL22RA2 (have cd1a and cd207--> LANGERHANS )
#5,12,14 = CCR7 --> they seem all inflammatory (nkfb1, relb, cd40 and cd80 ) 
#4 = T CELL (cytocd8 with gzmb and ifng)
#14 IS JUNK, 10 HAS C1Q 
#MAYBE SUBCLUSTER 7,8,9 

#V4: 
#5-7 = IL22RA2 (have cd1a and cd207--> LANGERHANS ; 6 MAY NOT BE)
#11 has lamp3 (CDC2A)  
#4 HAS CD3, GZMB 
#C1Q MACROPHAGES: 13 (MONOCYTE?); 2 AND 12  

cluster.ids = c('cDC2b','cDC2b','Macrophage','cDC2b','Tem','IL22RA2+ Langerhans',
                'IL22RA2+ Langerhans','IL22RA2+ Langerhans','cDC2b', 'cDC2b',
                'cDC2b','LAMP3+ DC', 'Macrophage','Monocyte')
names(cluster.ids) <- levels(DC2)
DC2 <- RenameIdents(DC2, cluster.ids)
DimPlot(DC2, reduction = "umap", label = TRUE, repel = TRUE)  

DC2$celltype = Idents(DC2)
DC2@meta.data$celltype %>% table



cluster12 = FindMarkers(DC2, ident.1 = 12, min.pct = 0.25)
cluster4 = FindMarkers(DC2_all, ident.1 = 4, min.pct = 0.25) #11 is CCR7+ 
#DC2 = subset(DC2, idents = c(0:12, 14, 16))
#DC2 = subset(DC2, idents= 1:13)
#DC2 = subset(DC2, idents= 1:12)
#DC2 = subset(DC2, idents= C(0:13))
#DC2 = subset(DC2, idents= C(0:12)) #TRYING TO GET RID OF CYCLING CELLS 
#DC2 = subset(DC2, idents= c(0:12, 14))
#DC2 = subset(DC2, idents= c(0:6, 8:10))
#DC2 = subset(DC2, idents= c(0:7, 9:12))

FeaturePlot(DC2_all, features = c( 'LAMP3','CCR7', 
                                   'IL22RA2','CD207','CD1A',
                                   'CLEC10A','CD1C','CLEC12A'), reduction = 'umap',
            cols = c('gray','red'), label = T)
FeaturePlot(DC2, features = c( 'CLEC10A','CD1C','CLEC12A'), reduction = 'umap',
            cols = c('gray','red'), label = T)

VlnPlot(DC2, features = c("CD3D",'CD3E','CD8A','IFNG','GZMB','CD40LG','ITGB1',
                         'JCHAIN', 'MKI67', 'TOP2A')) 
VlnPlot(DC2, features = c(Macro_mono_genes, 'CXCL9','CXCL10'))

#look into the mono and macrophage within DC2
sub_DC2 = subset(DC2, idents = c("Macrophage",'Monocyte'))
sub_DC2 = SeuratAnalysis1(sub_DC2)
ElbowPlot(sub_DC2, ndims = 50)
sub_DC2 = SeuratAnalysis2(sub_DC2, harmony = T, 
                      harmony_by = 'orig.ident', n = 25)
DimPlot(sub_DC2,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4)

cluster5 = FindMarkers(sub_DC2, ident.1 = 5, min.pct = 0.25)
FeaturePlot(sub_DC2, features = c( 'LAMP3','CD207','CD1A',Macro_mono_genes), reduction = 'umap',
            cols = c('gray','red'), label = T) 
VlnPlot(sub_DC2, features = c(Macro_mono_genes, 'TNF','CCL3','CCL4'))
VlnPlot(sub_DC2, features = c('TREM2', 'APOE','APOC1','FABP4','CD9',
                              'LAMP3','CD207','CD1A','IL22RA2',
                              'CLEC10A','CD1C','CLEC12A'))

VlnPlot(sub_DC2, features = c( 'JCHAIN', 'MKI67', 'TOP2A','IGKC','IGHM'))

#6 is cMo; 0,2,3, macrophages (3 IS INFLAMMATORY; 2 HAS TREM2 ); 1 AND4 ARE CONTAIMINATION 
cluster.ids = c('Macrophage','Unknown','LAM', 'Inflammatory macrophages','Unknown',
                'Langerhans cells','cMo')
names(cluster.ids) <- levels(sub_DC2)
sub_DC2 <- RenameIdents(sub_DC2, cluster.ids)
DimPlot(sub_DC2, reduction = "umap", label = TRUE, repel = TRUE)  

sub_DC2$celltype = Idents(sub_DC2)
sub_DC2@meta.data$celltype %>% table


DC2 = subset(DC2, idents = c('cDC2b','Tem','IL22RA2+ Langerhans','LAMP3+ DC'))
sub_DC2= subset(sub_DC2, idents = c('Macrophage','LAM', 'Inflammatory macrophages',
                                    'Langerhans cells','cMo'))
DC2_all = merge(DC2, sub_DC2)
DC2_all= JoinLayers(DC2_all)
saveRDS(DC2_all, file = 'cDC2_subcluster.rds')


#macrophage and monocyte subcluster ----------
#under Macro, I further subset monocytes and macrophages separately to be Mono_final and Macro_final
Macro = subset(Myeloid, idents = c(2,14,26,8,5,12,3,20,16, 19,22))
Macro = SeuratAnalysis1(Macro)
ElbowPlot(Macro, ndims = 50)
Macro = SeuratAnalysis2(Macro, harmony = T, 
                          harmony_by = 'orig.ident', n = 25)
DimPlot(Macro,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4)

cluster0 = FindMarkers(Macro, ident.1 = 0, min.pct = 0.25)
table(Macro$seurat_clusters)
#Macro = subset(Macro, idents = c(0:11, 13:18, 20:22))
#Macro = subset(Macro, idents = c(0:9, 11:21))

#V2
#13 HAS IFN-PRIMED; 
#STRUCTURAL: 22,19
#MONO: 1,7,0,4, 17, 13, 14
#CDC2: 21


#V1
#2,1, 15 monocyte, 21 and 20 unknown ; 0 is a mix of mono and macro 
#20 = hepatocyte; 21= club and AT2 mixing 


FeaturePlot(Macro, features = c(Macro_mono_genes), reduction = 'umap',
            cols = c('gray','red'), label = T) 
VlnPlot(Macro, features = c(Macro_mono_genes, 'TNF','CCL3','CCL4'))
VlnPlot(Macro, features = c(#'TREM2', 'APOE','APOC1','FABP4','CD9',
                              'LAMP3','CD207','CD1A','IL22RA2',
                              'CLEC10A','CD1C','CLEC12A'))
VlnPlot(Macro, features = c('CXCL9','CXCL10','CCL3','CCL4','TNF',
                            'CCL4L2'))
#5, 15 HAVE HIGH CXCL9/10; TNF IN 12

VlnPlot(Macro, features = c(#AT2_gene, 'SCGB1A1','SCGB3A2', AT1_gene,
                            'ALDOB','HPR', 'F12',  'HGFAC', 'AGXT', 'SAA4','ARG1'))

VlnPlot(Macro, features = c( 'JCHAIN', 'MKI67', 'TOP2A','IGKC','IGHM'))


#monocytes---------- 
Mono = subset(Macro, idents = c( 1,7,0,4, 17, 13, 14))
Mono = SeuratAnalysis1(Mono)
ElbowPlot(Mono, ndims = 50)
Mono = SeuratAnalysis2(Mono, harmony = T, 
                        harmony_by = 'orig.ident', n = 25)
DimPlot(Mono,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4)

#5 is ifn-primed ; 3,6,7 NOT MACROPHAGES ; 3= CLASSICAL MONOCYTE  
#4 = INFLAMMATORY MAC 
FeaturePlot(Mono, features = c(Macro_mono_genes), reduction = 'umap',
            cols = c('gray','red'), label = T)


cluster7 = FindMarkers(sub_Mono, ident.1 = 7, min.pct = 0.25)
table(sub_Mono$seurat_clusters)

#CLUSTERS TOO BIG TO DIFFERENTIATE, FURTHER SUBCLUSTER 

sub_Mono = subset(Mono, idents =c(3,6,7))
sub_Mono = SeuratAnalysis1(sub_Mono)
ElbowPlot(sub_Mono, ndims = 50)
sub_Mono = SeuratAnalysis2(sub_Mono, harmony = T, 
                           harmony_by = 'orig.ident', n = 20)
DimPlot(sub_Mono,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4)

#v2:
#7 = ncmo; 0=iMo (high hla) ,3 = cMo( high s100a8/9)
#macrophages: 6 , 4 
#ifn: 1, 
#inflamamtory: 5 
# 0 and 7 have cxcl10  

cluster.ids= c('iMo', 'IFN-primed macrophages', 'Unknown','cMo',
               'Macrophages', 'Inflammatory macrophages','Macrophages',
               'ncMo','Unknown')
names(cluster.ids) <- levels(sub_Mono)
sub_Mono <- RenameIdents(sub_Mono, cluster.ids)
DimPlot(sub_Mono, reduction = "umap", label = TRUE, repel = TRUE)  

sub_Mono$celltype = Idents(sub_Mono)
sub_Mono@meta.data$celltype %>% table

Mono_unknown = subset(sub_Mono, idents = 'Unknown')
sub_Mono = subset(sub_Mono, idents = c('iMo', 'IFN-primed macrophages','cMo',
                                    'Macrophages', 'Inflammatory macrophages',  'ncMo'))
saveRDS(sub_Mono, file ='Monocyte subcluster.rds')
#v1: 
#4,9= C1Q macrophage
#10=ncmo ; 8 ifn-primed ; 1,2,5 = monocyte (high MHCII, CCR2+, high FCGR3A)
#3,6,0,8 

FeaturePlot(sub_Mono, features = c(Macro_mono_genes), reduction = 'umap',
            cols = c('gray','red'), label = T) 
VlnPlot(sub_Mono, features = c(Macro_mono_genes,
                               'TNF','CCL3','CCL4','KLF2','KLF6',
                               'ISG15','IFIT1','IFI27','CCL2','CCL8',
                               'CXCL10','CXCL9'))
VlnPlot(sub_Mono, features = c( 'CXCR4','IL1RN','SERPINB2','PPARG','THBS1'))
#VlnPlot(sub_Mono, features = c( 'JCHAIN','IGHA1','IGHG1','IGHM','IGKC', 
 # 'CD3D','CD3E',
 # 'CD40LG','CD8A','CD8B',
 # 'NCAM1','FCGR3A','KLRF1','PRF1')) #no doublet 
VlnPlot(sub_Mono, features = c( "CDKN1C",'CX3CR1','FCGR3A', 'HES4',  #ncMo
                                #   'LILRA1','SIGLEC10', 'LILRB1','LILRB2' ,
                                'S100A12','S100A8','S100A9','CCR2','EREG','CD14','IL1RN','TCF7L2' ,  
                                'FCER1A','FCGR2B'  )) #mo-1 
FeaturePlot(sub_Mono, features = c( 'CD14', 'CCR2','NFKB1','RELB', 'LGALS2',
                                    'CX3CR1','CD86','CD40',
                                    'HLA-DRA',
                                    'S100A12','VCAN','S100A8','S100A9','FCGR3A'), cols = c('gray','red'),reduction = 'umap', label = T)  

cluster8 = FindMarkers(sub_Mono, ident.1 = 8, min.pct = 0.25) 

#sub_Mono= subset(sub_Mono, idents = c(0:10)) #11 has MT genes 
#sub_Mono= subset(sub_Mono, idents = c(0:8))

#really bad the two unknown clusters  

#sub_Mono2= subset(Mono, idents = c(5,4,8,2,0,9,1,10))
sub_Mono2 = SeuratAnalysis1(sub_Mono2)
ElbowPlot(sub_Mono2, ndims = 50)
sub_Mono2 = SeuratAnalysis2(sub_Mono2, harmony = T, harmony_by = 'orig.ident', n = 20)
DimPlot(sub_Mono2,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4)

VlnPlot(sub_Mono2, Macro_mono_genes)
cluster6 = FindMarkers(sub_Mono2, ident.1 = 6, min.pct = 0.25) 

#sub_Mono2 = subset(sub_Mono2, idents= c(0,1,3:10))
#sub_Mono2 = subset(sub_Mono2, idents= c(0:11)) 
#sub_Mono2= subset(sub_Mono2, idents = c(0,2,3,4,5,6,7,10))
#sub_Mono2= subset(sub_Mono2, idents = c(0:8)) 
#sub_Mono2= subset(sub_Mono2, idents = c(0,1,3,4,6,7,8))
#sub_Mono2= subset(sub_Mono2, idents = c(0,1,2,3,5,6,7))
#sub_Mono2= subset(sub_Mono2, idents = 0:7) 
#sub_Mono2= subset(sub_Mono2, idents = c(0,1,2,4,6,8)) 
sub_Mono2= subset(sub_Mono2, idents = 0:5) 

#there is a treg cluster??
#the separated 3 is ifn-primed 
table(sub_Mono2$seurat_clusters)
VlnPlot(sub_Mono2, features = c(#Macro_mono_genes,
                               'TNF','CCL3','CCL4','KLF2','KLF6','GOS2','IL1RN','NEDD9',
                               'ISG15','IFIT1','IFI27','CCL2','CCL8',
                               'CXCL10','CXCL9',
                               'TREM2','CD9','APOE','FABP5'))
VlnPlot(sub_Mono2, features = c('TREM2','CD9','APOE','FABP5',
                                'STAB1','SELENOP','LYVE1'))
VlnPlot(sub_Mono2, features = c( 'JCHAIN','IGHA1','IGHG1','IGHM','IGKC', 
                                 'CD3D','CD3E',
                                 'CD40LG','CD8A','CD8B',
                                 'FOXP3','CTLA4','ICOS','PDCD1') ) 
FeaturePlot(sub_Mono2, features = c(Macro_mono_genes), reduction = 'umap',
            cols = c('gray','red'), label = T) 

VlnPlot(sub_Mono2, features = c( 'CLEC10A','CD1C','CLEC12A','CCL17','ZBTB46','SIRPA', 
                                 'ITGAM','STAT5A', 'C1QA','C1QB','MARCO',
                                'CCR2','FCN1'
                              #  'CX3CR1','CDKN1A','FCGR3A', 
                                #'S100A8','S100A9','HLA-DRA'
                                ))
 
#remaining 
#2,3 = ifn; 1,5  = trem2; 0,4 =inflammatory  

cluster.ids = c('Inflamamtory macrophages','LAM', 'IFN-primed macrophages',
                'IFN-primed macrophages', 'Inflammatory macrophages','LAM')
names(cluster.ids) <- levels(sub_Mono2)
sub_Mono2 <- RenameIdents(sub_Mono2, cluster.ids)
DimPlot(sub_Mono2, reduction = "umap", label = TRUE, repel = TRUE)  

sub_Mono2$celltype = Idents(sub_Mono2)
sub_Mono2@meta.data$celltype %>% table
saveRDS(sub_Mono2, file ='Macrophage_subcluster1.rds')


#v5:
#ifn prime: 4,6 ; also have cxcl10 and 9 
#trem2: 2,3,6 
#2,5 MORE LIKE ccl17+ CDC2 (have cd1e, icoslg, il1b, ccl20, il1r1, cxcl3,hla, s100a) 
#8 more like monocyte = cMo 

#v4:
#ifn prime: 5,6 ; inflammatory 8; 9 is a doublet seems like 
#trem2: 2,3,

#V2:
#IFN PRIMED: 5,6,12 ; 12 IS DOUBLET WITH iG 
#TREG: 8 ; 1,7 ARE T CELLS / NK 

#V3:
#T CELL: 1,8,9; TREG= 9; --> further subcluster (9 maybe a doublet because also have cxcl9 and 10)
#IFN PRIMED = 6,7 
#inflammatory: 10, 11 ? 0 is inflmmatory with nkfb 
#11 IS JUNK ; 2 IS DOUBLET ?? 
#3 AND 5 HAVE CLEC10A, CD1C (MAYBE THE DERIVED CDC2B?)



#sub_Mono_Tcell= subset(sub_Mono2, idents = c(1,8,9))
sub_Mono_Tcell = SeuratAnalysis1(sub_Mono_Tcell)
ElbowPlot(sub_Mono_Tcell, ndims = 50)
sub_Mono_Tcell = SeuratAnalysis2(sub_Mono_Tcell, harmony = T, harmony_by = 'orig.ident', n = 20)
DimPlot(sub_Mono_Tcell,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4)

VlnPlot(sub_Mono_Tcell, Macro_mono_genes)
cluster3 = FindMarkers(sub_Mono_Tcell, ident.1 = 3, min.pct = 0.25) 
cDC2b_genes = FindMarkers(sub_Mono_Tcell, ident.1 = 'cDC2b', min.pct = 0.25)

VlnPlot(sub_Mono_Tcell, DC_genes)
VlnPlot(sub_Mono_Tcell, features = c( 'JCHAIN','IGHA1','IGHG1','IGHM','IGKC', 
                                 'CD3D','CD3E',
                                 'CD40LG','CD8A','CD8B',
                                 'FOXP3','CTLA4','ICOS','PDCD1') ) 
VlnPlot(sub_Mono_Tcell, features = c('CLEC10A','CD1C', 'ZBTB46','SIRPA','CCL17',
                                     'KLRF1','PRF1','SPON2','LAIR2','FCGR3A','NCAM1'))
VlnPlot(sub_Mono_Tcell, features = c(
  'TNF','CCL3','CCL4','KLF2','KLF6','GOS2','IL1RN','NEDD9',
  'ISG15','IFIT1','IFI27','CCL2','CCL8',
  'CXCL10','CXCL9',
  'TREM2','CD9','APOE','FABP5'))

#V2:
#4-5 IS INFLAMMATORY , 6 IS MATURE NK , 2 IS CDC2B; 0 IS IFN-PRIMED;1 IS DOUBLET ; 3 IS C1Q MACROPHAGES
cluster.ids = c("IFN-primed macrophages", 'Unknown','cDC2b','Macrophages',
                'Inflammatory macrophages','Inflammatory macrophages','Mature NK')
names(cluster.ids) <- levels(sub_Mono_Tcell)
sub_Mono_Tcell <- RenameIdents(sub_Mono_Tcell, cluster.ids)
DimPlot(sub_Mono_Tcell, reduction = "umap", label = TRUE, repel = TRUE)  

sub_Mono_Tcell$celltype = Idents(sub_Mono_Tcell)
sub_Mono_Tcell@meta.data$celltype %>% table

sub_Mono_Tcell = subset(sub_Mono_Tcell, idents = c("IFN-primed macrophages",'cDC2b','Macrophages',
                                                   'Inflammatory macrophages','Mature NK'))
sub_Mono_Tcell = JoinLayers(sub_Mono_Tcell)
saveRDS(sub_Mono_Tcell, file = "Macrophage_subcluster2.rds")

#2,5 = cytocd8 (DOUBLET); 3 = treg; 1 contamination with Ig? 7 = nk 
#4 = CCL17+ CDC2 
#TREG DOUBLET WITH MACROPHAGES, 


#sub_Mono_Tcell= subset(sub_Mono_Tcell, idents = c(0:6, 8:10))
#sub_Mono_Tcell= subset(sub_Mono_Tcell, idents = c(0:6))
#sub_Mono_Tcell= subset(sub_Mono_Tcell, idents = c(0,1,4,6,7))
#sub_Mono_Tcell= subset(sub_Mono_Tcell, idents = c(0:6))


#sub_Mono3= subset(sub_Mono2, idents = c(2,5))
#SUB_Mono4 = subset(sub_Mono2, idents = c(3,5,7))

#sub_Mono3 = merge(sub_Mono3, SUB_Mono4)
sub_Mono3 = SeuratAnalysis1(sub_Mono3)
ElbowPlot(sub_Mono3, ndims = 50)
sub_Mono3 = SeuratAnalysis2(sub_Mono3, harmony = T, harmony_by = 'orig.ident', n = 20)
DimPlot(sub_Mono3,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'celltype')

VlnPlot(sub_Mono3, Macro_mono_genes)
cluster3 = FindMarkers(sub_Mono3, ident.1 = 3, min.pct = 0.25) 

VlnPlot(sub_Mono3, features = c( 'JCHAIN','IGHA1','IGHG1','IGHM','IGKC', 
                                      'CD3D','CD3E') ) #NO MORE DOUBLET 
VlnPlot(sub_Mono3, features = c('CCL17','CCR2','ICOSLG','IL1R1','CCL20','IL1B',
                                'CLEC10A','CD1C'), group.by = 'celltype')
VlnPlot(sub_Mono3, features = c('CLEC10A','CD1C', 'ZBTB46','SIRPA','CLEC9A')) 
#2,4,7 MAYBE CDC2B 
#cMo=6 ; lam = 2,8 ; 4=ifn primed ; 5 = PVM

cluster.ids = c("Unknown", 'Unknown', 'LAM', 'Unknown', 'IFN-primed macrophages',
                'PVM','cMo', 'Unknown','LAM')
names(cluster.ids) <- levels(sub_Mono3)
sub_Mono3 <- RenameIdents(sub_Mono3, cluster.ids)
DimPlot(sub_Mono3, reduction = "umap", label = TRUE, repel = TRUE)  

sub_Mono3$celltype = Idents(sub_Mono3)
sub_Mono3@meta.data$celltype %>% table

sub_Mono3 = subset(sub_Mono3, idents = c('LAM', 'IFN-primed macrophages',
                'PVM','cMo'))
sub_Mono3 = JoinLayers(sub_Mono3)
saveRDS(sub_Mono3, file = 'Macrophage_subcluster3.rds')


VlnPlot(sub_Mono3, Macro_mono_genes)
VlnPlot(sub_Mono3, features = c(
  #'TNF','CCL3','CCL4','KLF2','KLF6','NEDD9',
  'ISG15','IFIT1','IFIT2','CCL2','CCL8',
 # 'CXCL10','CXCL9',
  'TREM2','CD9','APOE','FABP5',
  'SELENOP','LYVE1'), group.by = 'celltype')

VlnPlot(sub_Mono3, features = c( #"CDKN1C",'CX3CR1','FCGR3A', 'HES4',  #ncMo
                                  # 'LILRA1','SIGLEC10', 'LILRB1','LILRB2' ,
                           'HLA-DRA','S100A8','FCN1','S100A9','CD14','FCGR3A', 'CCR2','CD86','CD40' #cMo 
                                 #'FCER1A','FCGR2B' 
                                            ), group.by = 'celltype')

VlnPlot(sub_Mono3, c(Neutrophil_gene, 'CXCR4','CXCR2','ICAM1','ELANE','MPO')) 


#sub_Mono3 = subset(sub_Mono3, idents = 0:7)
#sub_Mono3 = subset(sub_Mono3, idents = c(0:3, 5:7))
#sub_Mono3 = subset(sub_Mono3, idents = c(0:4))

#MERGE ALL MONOCYTE SUBCLUSTERS -------
Mono_final = merge(x = sub_Mono, y = c(sub_Mono2, sub_Mono3, sub_Mono_Tcell))
Mono_final = SeuratAnalysis1(Mono_final)
ElbowPlot(Mono_final, ndims = 50)
Mono_final = SeuratAnalysis2(Mono_final, harmony = T, harmony_by = 'orig.ident', n = 30)
DimPlot(Mono_final,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by='celltype')



Idents(Mono_final) = 'celltype'
Mono_final = RenameIdents(Mono_final, 'Inflamamtory macrophages' ='Inflammatory macrophages')
Mono_final$celltype = Idents(Mono_final)
Mono_final = RenameIdents(Mono_final, 'Macrophages' ='Interstitial macrophages')
Mono_final$celltype = Idents(Mono_final)
table(Mono_final$celltype)

Mono_nk = FindMarkers(Mono_final, ident.1 = "Mature NK", min.pct=0.25) #seems pretty confident 
#saveRDS(Mono_final, file = "Monocytes_v2.rds") 

VlnPlot(Mono_final, features = c(
  'C1QA','C1QB','MARCO','CD163',
  'ISG15','IFIT1','IFIT2','CCL2','CCL8',
  'CLEC10A','CD1C', 'ZBTB46','SIRPA',
  'KLRF1','PRF1','NKG7'
  ), group.by = 'celltype')

DotPlot(Mono_final, features = c(
  'C1QA','C1QB','MARCO','CD163',
  'TNF','CCL3','CCL4','KLF2','NEDD9',
  'ISG15','IFIT1','IFIT2','CCL2','CCL8',
  # 'CXCL10','CXCL9',
  'TREM2','CD9','APOE','FABP5',
  'SELENOP','LYVE1',
  'FCN1',"CDKN1C",'CX3CR1','FCGR3A','CD14','S100A8','S100A9','CCR2',
  'CLEC10A','CD1C', 'ZBTB46','SIRPA',
  'KLRF1','PRF1','NKG7'
), group.by = 'celltype', cols = 'RdYlBu') + coord_flip()+ RotatedAxis()

#using celltypist macrophage genes 
DotPlot(Mono_final, features = c('TNIP3','TNFAIP3','CCL2',
                                 'SPP1','VCAN','VEGFA',
                                 'FCN1'), cols = 'RdYlBu') #INTERMEDIATE MACROPHAGES 
#A LITTLE BIT SPP1, BUT NO OTHER ANGIOGENIC MARKERS 


#macrophages------------
#sub_Macro = subset(Macro, idents = c( 2,3,5,6,8,9,10,11,12,15,16,18,20))
sub_Macro = SeuratAnalysis1(sub_Macro)
ElbowPlot(sub_Macro, ndims = 50)
sub_Macro = SeuratAnalysis2(sub_Macro, harmony = T, 
                        harmony_by = 'orig.ident', n = 20)
DimPlot(sub_Macro,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'celltype')

cluster3 = FindMarkers(sub_Macro, ident.1 = 3, min.pct = 0.25)
table(sub_Macro$seurat_clusters)

VlnPlot(sub_Macro, features = c( 'JCHAIN','IGHA1','IGHG1','IGHM','IGKC', 
                                 'CD3D','CD3E','CD8A','IFNG','KLRD1') ) 
VlnPlot(sub_Macro, features = c('ISG15','IFIT1','IFI27','CCL2','CCL8',
                                'TREM2','CD9','APOE','FABP5',
                                'SELENOP','LYVE1'))
VlnPlot(sub_Macro, features = c('TNF','CCL3','CCL4','KLF2','NEDD9','CCL3L1',
                                'CCL4L2','KLF6'))
#NO CYCLING 
VlnPlot(sub_Macro, Macro_mono_genes)
VlnPlot(sub_Macro, c(fibro_genes))
VlnPlot(sub_Macro, AlveolarMac_gene)
DotPlot(sub_Macro, features = c(
  'C1QA','C1QB','MARCO','CD163',
  AlveolarMac_gene, 
  'TNF','CCL3','CCL4','KLF2', 
  'ISG15','IFIT1','IFIT2','CCL2','CCL8',
  # 'CXCL10','CXCL9',
  'TREM2','CD9','APOE','FABP5',
  'SELENOP','LYVE1',
 # 'FCN1',"CDKN1C",'CX3CR1','FCGR3A','CD14','S100A8','S100A9','CCR2',
 # 'CLEC10A','CD1C', 'ZBTB46','SIRPA',
  'CD3D','CD3E','CD8A','CD40LG', 
  'KLRF1','PRF1','NKG7',
  'LUM', 'DCN'
), group.by = 'seurat_clusters', cols = 'RdYlBu') + coord_flip()+ RotatedAxis()

DotPlot(sub_Macro, features = c(
  'C1QA','C1QB','MARCO','CD163',
  AlveolarMac_gene, 'GPNMB','BHLHE41', #THE LAST TWO FROM CELLTYPIST 
  'TNF','CCL3','CCL4','KLF2', 
  'ISG15','IFIT1','IFIT2','CCL2','CCL8', 
  'TREM2','CD9','APOE','FABP5',
  'SELENOP','LYVE1',
  #'CD5L','SPIC','SLC40A1', 'VCAM1','CXCL12', #FROM CELLTYPIST :DOES NOT REALLY HELP
  'CD3D','CD3E','CD8A','IFNG','GZMB'
), group.by = 'celltype', cols = 'RdYlBu') + coord_flip()+ RotatedAxis()



VlnPlot(sub_Macro, features = c("CDKN1C",'CX3CR1','FCGR3A', 'HES4',  #ncMo
                                'LILRA1','SIGLEC10', 'LILRB1','LILRB2' ,
                                'HLA-DRA','S100A8','FCN1','S100A9','CD14','FCGR3A', 'CCR2','CD86','CD40', #cMo 
                                'FCER1A','FCGR2B' ))
#final:
#1,11 = PVM; 6 = LAM, 8 = IFN-PRIMED, 4,7 = INFLAMMATORY 
#9= FIBRO, 3 IS CYTOCD8  
#0,5 are alveolar mac 
cluster.ids = c('Alveolar macrophage', 'PVM', 'CytoCD8', 'Unknown', 
                'Inflammatory macrophages',
                'Alveolar macrophage','LAM','Inflammatory macrophages',
                'IFN-primed macrophages',
                'Fibro', 'Unknown','PVM')  

names(cluster.ids) <- levels(sub_Macro)
sub_Macro <- RenameIdents(sub_Macro, cluster.ids)
DimPlot(sub_Macro, reduction = "umap", label = TRUE, repel = TRUE)  

sub_Macro$celltype = Idents(sub_Macro)
sub_Macro@meta.data$celltype %>% table

sub_Macro= subset(sub_Macro, idents = c('Alveolar macrophage', 'PVM', 'CytoCD8',  
                                        'Inflammatory macrophages', 'LAM', 
                                        'IFN-primed macrophages', 'Fibro'))
saveRDS(sub_Macro, file ="Macrophage subcluster 4.rds")
#v6:
#1, 12 = PVM; 6 = LAM, 10= IFN-PRIMED, 7,9 = INFLAMMATORY 
#2 AND 5 MAY BE MONOCYTE 
# 4 HAS CCR2, MAY ALSO BE MONOCYTE?? 
#11 IS FIBRO, 3 IS CYTOCD9  
#0 and 8 are alveolar mac 
#CHECK 4,2,5 

#v5
#5 = t cells (cytocd8) ; 12 is fibro; 1 may have endothelial but not strong 
#4  and 13 = pvm (13 may be inflammatory); 9 = lam; 1 = inflammatory mac (tnf) 
#11 is ifn-primed ; 3 may be ncmo 


#v3
#fibro =12 ; endo = 10, 2 = cytocd8, 15= plasma, 11 ALSO t CELLS
#monocyte maybe 0 and 8 
#16 very high c1q 

#V2
#14, 15 HAS JCHAIN (DOUBLET??), 2 IS CYTOCD8, 11 ALSO t CELLS --> GET RID OF 14 

#16 is fibro 
#8 has cd3 (cytocd8), 0 and 12.15 have ig genes (12 is endothelial)
#get rid of 15, 13 , 10, 17

#sub_Macro= subset(sub_Macro, idents = c(0:9,11,12,14,16))
#sub_Macro= subset(sub_Macro, idents = c(0:13, 15))
#sub_Macro= subset(sub_Macro, idents = c(0:15))
#sub_Macro= subset(sub_Macro, idents = c(1:10, 12))
#sub_Macro= subset(sub_Macro, idents = c(0,1,3:13))
#sub_Macro= subset(sub_Macro, idents = c(0,1,3:13))
#sub_Macro= subset(sub_Macro, idents = c(0,1,3,6,7,8,9,10,11,12))


#sub_Macro_unknown = subset(sub_Macro, idents ="Unknown")
sub_Macro_unknown = SeuratAnalysis1(sub_Macro_unknown)
ElbowPlot(sub_Macro_unknown, ndims = 50)
sub_Macro_unknown = SeuratAnalysis2(sub_Macro_unknown, harmony = T, 
                                 harmony_by = 'orig.ident', n = 10)
DimPlot(sub_Macro_unknown,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4)
table(sub_Macro_unknown$seurat_clusters)

cluster4 = FindMarkers(sub_Macro_unknown, ident.1 = 4, min.pct = 0.25)

VlnPlot(sub_Macro_unknown, features = c(Macro_mono_genes, AlveolarMac_gene, 'GPNMB','BHLHE41'))
VlnPlot(sub_Macro_unknown, features = c(
  #'C1QA','C1QB','MARCO','CD163',
  #AlveolarMac_gene, 'GPNMB','BHLHE41', #THE LAST TWO FROM CELLTYPIST 
  'TNF','CCL3','CCL4','KLF2', 
  'ISG15','IFIT1','IFIT2','CCL2','CCL8', 
  'TREM2','CD9','APOE','FABP5',
  'SELENOP','LYVE1'
  #'CD5L','SPIC','SLC40A1', 'VCAM1','CXCL12', #FROM CELLTYPIST :DOES NOT REALLY HELP
  ))

#only 0 and 3 are cd14+ 
#2 = trem2, 1= inflammatory,3 seems a mix 
 
cluster.ids = c("Unknown",'Unknown','Unknown', 'Inflammatory macrophages','Interstitial macrophages')
names(cluster.ids) <- levels(sub_Macro_unknown)
sub_Macro_unknown <- RenameIdents(sub_Macro_unknown, cluster.ids)
DimPlot(sub_Macro_unknown, reduction = "umap", label = TRUE, repel = TRUE)  

sub_Macro_unknown$celltype = Idents(sub_Macro_unknown)
sub_Macro_unknown@meta.data$celltype %>% table

sub_Macro_unknown = subset(sub_Macro_unknown, idents = c( 'Inflammatory macrophages','Interstitial macrophages'))


#sub_Macro_unknown= subset(sub_Macro_unknown, idents = c(0:4, 6))
#sub_Macro_unknown= subset(sub_Macro_unknown, idents = c(0:3, 5))
#sub_Macro_unknown= subset(sub_Macro_unknown, idents = c(0:3))


#sub_Macro_mono = subset(sub_Macro, idents = c(2,4,5))
sub_Macro_mono = SeuratAnalysis1(sub_Macro_mono)
ElbowPlot(sub_Macro_mono, ndims = 50)
sub_Macro_mono = SeuratAnalysis2(sub_Macro_mono, harmony = T, 
                            harmony_by = 'orig.ident', n = 20)
DimPlot(sub_Macro_mono,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4)

VlnPlot(sub_Macro_mono, features = c("CDKN1C",'CX3CR1','FCGR3A', 'HES4',  #ncMo
                                #'LILRA1','SIGLEC10', 'LILRB1','LILRB2' ,
                                'HLA-DRA','S100A8','FCN1','S100A9','CD14','FCGR3A', 'CCR2','CD86','CD40',
                                'CRISPLD2','CLU','MARC2', #cMo FROM CELLTYPIST
                                'FCER1A','FCGR2B' ))
FeaturePlot(sub_Macro_mono, c(Macro_mono_genes,'S100A8','CCR2','S100A9','CRISPLD2','CLU','MARC2'), reduction = 'umap', 
            cols = c('gray','red'), label = T)

VlnPlot(sub_Macro_mono, features = c("SPP1",'TNIP3','CCL2',
                                     #'VCAN','VEGFA','THBS1','FCN1','TIMP1', 'SERPINB2',
                                     #'ISG15','IFIT1','IFIT2','CCL2','CCL8',
                                    'TNF','CCL3','CCL4','KLF2'), group.by = 'seurat_clusters')
VlnPlot(sub_Macro_mono, features = c('TREM2','CD9','APOE','FABP5','APOC1'))
VlnPlot(sub_Macro_mono, features = c(Macro_mono_genes, AlveolarMac_gene, 'GPNMB','BHLHE41'))
DotPlot(sub_Macro_mono, features = c(
  'C1QA','C1QB','MARCO','CD163',
  AlveolarMac_gene, 'GPNMB','BHLHE41', #THE LAST TWO FROM CELLTYPIST 
  'TNF','CCL3','CCL4','KLF2', 
  'ISG15','IFIT1','IFIT2','CCL2','CCL8', 
  'TREM2','CD9','APOE','FABP5',
  'SELENOP','LYVE1',
  #'CD5L','SPIC','SLC40A1', 'VCAM1','CXCL12', #FROM CELLTYPIST :DOES NOT REALLY HELP
  'CD3D','CD3E','CD8A','IFNG','GZMB'
), cols = 'RdYlBu') + coord_flip()+ RotatedAxis()
VlnPlot(sub_Macro_mono, DC_genes)


cluster3 = FindMarkers(sub_Macro_mono, ident.1 = 3, min.pct = 0.25)

#v3:
#no ncmo BUT HIGH FCGR3A
#0 +4is alveolar, 1 has spp1/ccr2/cd274/cxcl10, 3=spp1+ 
#5 IS INFLAMMATORY, 2,6 IS IFN PRIMED 
 
cluster.ids = c("Alveolar macrophage", 'SPP1+ macrophages','IFN-primed macrophages',
                'SPP1+ macrophages',"Alveolar macrophage", 'Inflammatory macrophages',
                'IFN-primed macrophages' ) 

names(cluster.ids) <- levels(sub_Macro_mono)
sub_Macro_mono <- RenameIdents(sub_Macro_mono, cluster.ids)
DimPlot(sub_Macro_mono, reduction = "umap", label = TRUE, repel = TRUE)  

sub_Macro_mono$celltype = Idents(sub_Macro_mono)
sub_Macro_mono@meta.data$celltype %>% table


#sub_Macro_mono = subset(sub_Macro_mono, idents = 0:4)
#sub_Macro_mono = subset(sub_Macro_mono, idents = 1:4)

#merge macrophages 
#Macro_final = merge(x=sub_Macro, y=c(sub_Macro_mono, sub_Macro_unknown) ) 
Macro_final = merge(x = Macro_final, y = c(sub_DC, sub_structural))

bad= subset(Macro_final, celltype_condition =='Alveolar macrophage_Breast' | celltype_cluster == "10_IFN-primed macrophages" | celltype_cluster =="10_Inflammatory macrophages")
Macro_final$barcodes = rownames(Macro_final@meta.data)
Idents(Macro_final) ='barcodes'
Macro_final = subset(Macro_final, idents = Macro_final$barcodes[!Macro_final$barcodes %in% rownames(bad@meta.data)] )
Macro_final=JoinLayers(Macro_final)
Macro_final = merge(Macro_final, macro_clean)

Macro_final = SeuratAnalysis1(Macro_final)
ElbowPlot(Macro_final, ndims = 50)
Macro_final = SeuratAnalysis2(Macro_final, harmony = T,
                              harmony_by = 'orig.ident', n = 20)
DimPlot(Macro_final,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by='orig.ident')

saveRDS(Macro_final, file = 'Macrophages_merged_v2.rds')

table(Macro_final$seurat_clusters)
DotPlot(Macro_final, features = c(
  AT1_gene, AT2_gene,
  'CLEC10A','CLEC12A','CD1C','CD1E', #'CSF2RA', 'ZBTB46','SIRPA','ITGAM',
  'C1QA','C1QB',#'MARCO','CD163','CD68','MRC1', 
  'SPP1','CCR2','CXCL10','TGFB1','CD80', #BUT SPP1 DOES NOT HAVE VCAN, VEGFA, THBS1 FOR ANGIOGENIC 
  AlveolarMac_gene, 'GPNMB','BHLHE41', #THE LAST TWO FROM CELLTYPIST 
  'TNF','CCL3','CCL4','CCL3L1','CCL4L2',#'KLF2', 'NEDD9','G0S2',
  'ISG15','IFIT1','IFIT2','CCL2','CCL8', 
  'TREM2','CD9','APOE','FABP5','ACP5',
  'SELENOP','LYVE1','LGMN','STAB1',
  'CD3D','CD3E','CD8A','IFNG','GZMB',
   'LUM','DCN'
), cols = 'RdYlBu', group.by = 'celltype') + coord_flip()+ RotatedAxis()

DotPlot(Macro_final, features = c('TREM2', 'CD86', 'CCR2', 'MRC1','CD80'), cols = 'RdYlBu', group.by = 'celltype')

VlnPlot(Macro_final, features = c(
  'C1QA','C1QB','MARCO','CD163','CD68','MRC1'), group.by = 'celltype') 
  
cluster11 = FindMarkers(Macro_final, ident.1= 11, min.pct =0.25)
#10 has hepatocyte
Idents(Macro_final) ='celltype_condition'
Alveolar_genes  = FindMarkers(Macro_final, ident.1= 'Alveolar macrophage_Breast', min.pct =0.25)

Idents(Macro_final) ='celltype'
IFNmacro_genes = FindMarkers(Macro_final, ident.1='IFN-primed macrophages', min.pct=0.25)
Inflammatory_genes = FindMarkers(Macro_final, ident.1='Inflammatory macrophages', min.pct=0.25)


Macro_final$celltype_condition = paste0(Macro_final$celltype, "_", Macro_final$orig.ident)
table(Macro_final$celltype_condition)

#sub_structural = subset(Macro, idents = c(22, 19))
sub_structural = SeuratAnalysis1(sub_structural)
ElbowPlot(sub_structural, ndims = 50)
sub_structural = SeuratAnalysis2(sub_structural, harmony = F, n = 10)
DimPlot(sub_structural,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'celltype')

cluster1 = FindMarkers(sub_structural, ident.1 = 1, min.pct = 0.25)
table(sub_structural$seurat_clusters)

VlnPlot(sub_structural, features = c(#"MKI67",'TOP2A','LUM','DCN',
                                     'EPCAM',
                                     AT1_gene, AT2_gene, 'SCGB1A1','SCGB3A1'))
VlnPlot(sub_structural, features = c(
  'C1QA','C1QB','MARCO','CD163','CD68','MRC1', DC_genes)) 
VlnPlot(sub_structural, features = c('TNF','CCL3','CCL4','CCL3L1','CCL4L2','KLF2', 'NEDD9','G0S2','KLF6',
                                     'ISG15','IFIT1','IFIT2',#'CCL2','CCL8', 'LAMP3',
                                    # 'TREM2','CD9','APOE','FABP5',
                                    # 'SELENOP','LYVE1',
                                     'NKKB1','TNFAIP3'))
#V2: 2= CDC2; 5 IS AT1, 4 IS AT2; 0,1,3 ARE MACROPHAGES 

cluster.ids = c('IFN-primed macrophages','Inflammatory macrophages', 'cDC2b',
                'Inflammatory macrophages','Alveolar cell type 2','Alveolar cell type 1')
names(cluster.ids) <- levels(sub_structural)
sub_structural <- RenameIdents(sub_structural, cluster.ids)
DimPlot(sub_structural, reduction = "umap", label = TRUE, repel = TRUE)  

sub_structural$celltype = Idents(sub_structural)
sub_structural@meta.data$celltype %>% table


#4 HAS EPCAM (lung structural cells), 0,1,5 are macrophages; 3 = CDC2B
#6 is T cells  
#sub_structural= subset(sub_structural, idents = c(0,1,3,4,5))

#sub_DC = subset(Macro, idents = c(21)) 
sub_DC = SeuratAnalysis1(sub_DC)
ElbowPlot(sub_DC, ndims = 50)
sub_DC = SeuratAnalysis2(sub_DC, harmony = T, harmony_by = 'orig.ident',
                         n = 25)
DimPlot(sub_DC,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4)

cluster1 = FindMarkers(sub_DC, ident.1 = 1, min.pct = 0.25)
table(sub_DC$seurat_clusters)

VlnPlot(sub_DC, features = c(DC_genes, 'C1QA','C1QB','MARCO','CD163','CD68','MRC1'))
#1 is macrophages with high c1q and trem2 

VlnPlot(sub_DC, features = c(#'TREM2','CD9','APOE','FABP5',
                             #'ICAM1','CXCR4','CXCR2','S100A8','S100A9','CSF3R', #possible neutrophil
                             'S100A12','CD14','FCN1','S100A8','S100A9','FCGR3A')) #possible monocytes 
VlnPlot(sub_DC, features = c('CLEC10A','CD1C','ITGAM', 'ITGAX','IRF4','ZBTB46'))
#not really confident about the DC2 because has low zbtb46 and itgam, and the findmarkers reveals cycling genes 

#4 has IG genes 
#sub_DC= subset(sub_DC, idents = c(0,1,2,3))

cluster.ids = c('Unknown','LAM','Unknown')
names(cluster.ids) <- levels(sub_DC)
sub_DC <- RenameIdents(sub_DC, cluster.ids)
DimPlot(sub_DC, reduction = "umap", label = TRUE, repel = TRUE)  

sub_DC$celltype = Idents(sub_DC)
sub_DC@meta.data$celltype %>% table
sub_DC = subset(sub_DC, idents='LAM')

#clean up from Macro final
#macro_clean = subset(Macro_final, celltype_condition =='Alveolar macrophage_Breast' | celltype_cluster == "10_IFN-primed macrophages" | celltype_cluster =="10_Inflammatory macrophages")
macro_clean = SeuratAnalysis1(macro_clean)
ElbowPlot(macro_clean, ndims = 50)
macro_clean = SeuratAnalysis2(macro_clean, harmony = T, harmony_by = 'orig.ident',
                         n = 20)
DimPlot(macro_clean,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4)

cluster0 = FindMarkers(macro_clean, ident.1 = 0, min.pct = 0.25)

cluster.ids=c("Unknown", 'Inflammatory macrophages','IFN-primed macrophages')
names(cluster.ids) <- levels(macro_clean)
macro_clean <- RenameIdents(macro_clean, cluster.ids)
DimPlot(macro_clean, reduction = "umap", label = TRUE, repel = TRUE)  

macro_clean$celltype = Idents(macro_clean)
macro_clean@meta.data$celltype %>% table
macro_clean = subset(macro_clean, idents= c( 'Inflammatory macrophages','IFN-primed macrophages'))

#macro_clean= subset(macro_clean, idents = c(0,1))
#macro_clean= subset(macro_clean, idents = c(0,1,3)) 

VlnPlot(macro_clean, features = c(
  'C1QA','C1QB','MARCO','CD163','CD68','MRC1', DC_genes)) 
VlnPlot(macro_clean, features = c('TNF','CCL3','CCL4','CCL3L1','CCL4L2','KLF2', 'NEDD9','G0S2','KLF6',
                                     'ISG15','IFIT1','IFIT2'))



#T and NK and B in myeloid------ 
#T and NK and B: 11 and part 15 , 1, 0, 25,18,
#Myeloid_T = subset(Myeloid, idents = c(11, 15,1,0,25,18))
Myeloid_T = SeuratAnalysis1(Myeloid_T)
ElbowPlot(Myeloid_T, ndims = 50)
Myeloid_T = SeuratAnalysis2(Myeloid_T, harmony = T, 
                            harmony_by = 'orig.ident', n = 10)
DimPlot(Myeloid_T,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'celltype')

cluster6 = FindMarkers(Myeloid_T, ident.1 = 6, min.pct = 0.25)
table(Myeloid_T$seurat_clusters)
#all have very high RP or MT genes 

#Myeloid_T = subset(Myeloid_T, idents = c(0:17))
#Myeloid_T = subset(Myeloid_T, idents = c(0:8,10, 13:17,19,20))
#Myeloid_T = subset(Myeloid_T, idents = c(0, 2:19 ))
#Myeloid_T = subset(Myeloid_T, idents = c(0:7, 10:12,15:18))
#Myeloid_T= subset(Myeloid_T, idents =c(1,5,6,9,10,12))
#Myeloid_T= subset(Myeloid_T, idents =0:9)
#Myeloid_T= subset(Myeloid_T, idents =0:11)
#Myeloid_T= subset(Myeloid_T, idents =0:10)
#Myeloid_T= subset(Myeloid_T, idents = c(0:9))
#Myeloid_T= subset(Myeloid_T, idents = c(0:5,7,8,10))
#Myeloid_T= subset(Myeloid_T, idents = c(1,5,7,8,10))
#Myeloid_T=subset(Myeloid_T, idents=c(0,1,4:8))
#Myeloid_T=subset(Myeloid_T, idents=c(0:3, 5:8))

VlnPlot(Myeloid_T, features = c("IGKC",'IGHA1','JCHAIN',#'MKI67','TOP2A',
                                'CD3D','CD3E',#'CD40LG',
                                'CD8A','CD8B', 
                                'NCAM1','FCGR3A', 'KLRF1','PRF1',
                                'IFNG','IFNGR1', 'ITGB1','GZMB','GZMK'
))
VlnPlot(Myeloid_T, features = c(AT1_gene, AT2_gene,
                                'FABP4', 'MCEMP1', 'PPARG'))
FeaturePlot(Myeloid_T, features = c(Macro_mono_genes, 'CD3D','CD3E', 
                                    'CD8A','CD8B' ), reduction = 'umap',
            cols =c('gray','red'), label =T)
VlnPlot(Myeloid_T, Macro_mono_genes)
VlnPlot(Myeloid_T, features = c('TNF','CCL3','CCL4','CCL3L1','CCL4L2','KLF2', 'NEDD9','G0S2','KLF6',
                               'ISG15','IFIT1','IFIT2',#'CCL2','CCL8', 'LAMP3',
                               'TREM2','CD9','APOE','FABP5',
                                'SELENOP','LYVE1',
                               'SPP1','VCAN','VEGFA','THBS1'
                                     #'NKKB1','TNFAIP3'
                                     ))
#doublets of myeloid and T cells! 
#clean version:
#1,4,5,6,7 are t lymphoid;1,5,6 have cd8 ; 7 may be nk 
#the rest are c1q; 
#7 is nk-macrophage doublet ; 4 more like inflammatory mac

cluster.ids = c("Interstitial macrophages", 'CytoCD8-Interstitial macrophages',
                'Interstitial macrophages', 'Interstitial macrophages', 
                'Inflammatory macrophages',
                'SPP1+ macrophages', 'Unknown', 'Mature NK-macrophage')
names(cluster.ids) <- levels(Myeloid_T)
Myeloid_T <- RenameIdents(Myeloid_T, cluster.ids)
DimPlot(Myeloid_T, reduction = "umap", label = TRUE, repel = TRUE)  

Myeloid_T$celltype = Idents(Myeloid_T)
Myeloid_T@meta.data$celltype %>% table
Myeloid_T = subset(Myeloid_T, idents= c( 'Inflammatory macrophages',
                                         'CytoCD8-Interstitial macrophages',
                                         'Interstitial macrophages', 
                                         'SPP1+ macrophages', 'Mature NK-macrophage'))

saveRDS(Myeloid_T, 'Myeloid_T_subcluster1.rds')

#10 = nk; 1,7,8 have cd3; 8 has cd8 ; 5 HAS TRDC, KLRD1, GZM, a little cd8 
#except for 0,2,5,7,8 all have c1q  
#1, 7,8 is spp1+ tam; 6 AND 9 ARE MIXED OF C1Q AND SCGB LUNG GENES 
#3 HAS C1Q; 4 IS JUNK 

Myeloid_T_junk3 = subset(Myeloid_T, idents =c(0,2,3,4,6,9)) #either doublet or have very high RP genes 


#v5:
#9 = nk; 2 and 4 have cd3 (also 0,1,8,10,12); igkc in 1,5,9 
#2,4,1,5,7 are macrophages
#12 is doublet ; 2 and 4 are doublets 

#v4: 
#5,6,1 macrophages; 1 may have monocytes; 9 may be gdt; 12 is inflammatory 
#10 has high spp1 and vcan; 

#v3: 18 is macrophage (inflammatory maybe); 
#13-14 are t cells with some nk and jchain ; 19 is doublet; 8 and 9 also have cd3d 
#0-6 and 13 and 17-19 have c1q 

#cycling: 11 (maybe plasmablast)
#t cell: 0-6 ; 5 also has nk 
#V2:
#seems A MIXING OF B AND T CELLS 
#nk=8,6; t cells are 7,8,10,13-15; 1 and 20 have high ig genes ; 0-3 may be alveolar cells 
#14 have cd14; 8 as well --> doublet of myeloid and lymphoid 



#separate out the clusters more like T cell  
#Myeloid_T_sub= subset(Myeloid_T, idents = c(13,14,8,9))
Myeloid_T_sub = SeuratAnalysis1(Myeloid_T_sub)
ElbowPlot(Myeloid_T_sub, ndims = 50)
Myeloid_T_sub = SeuratAnalysis2(Myeloid_T_sub, harmony = T, 
                            harmony_by ='orig.ident', n = 15)
DimPlot(Myeloid_T_sub,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'celltype' )

cluster9 = FindMarkers(Myeloid_T_sub, ident.1 = 9, min.pct = 0.25)
table(Myeloid_T_sub$seurat_clusters)
#8 is mast cell; 7 is ncmo 
#2,3,5 are macrophages (c1q)

#v2: 1,2,7,9 macrophages; 4 and 8 are monocytes ; 
#7,8,9 doublet with T cell
#11 = MAST CELL, 8 = NCMO; 4 = CMO ; 7 = LAM; 2= SPP1+ MAC 
#3= CLUB CELL

cluster.ids = c('Unknown', 'Interstitial macrophages','SPP1+ macrophages',
                'Club cell','cMo', 'Unknown', 'Unknown', 'LAM','ncMo',
                'Interstitial macrophages', 'Unknown','Mast cell')
names(cluster.ids) <- levels(Myeloid_T_sub)
Myeloid_T_sub <- RenameIdents(Myeloid_T_sub, cluster.ids)
DimPlot(Myeloid_T_sub, reduction = "umap", label = TRUE, repel = TRUE)  

Myeloid_T_sub$celltype = Idents(Myeloid_T_sub)
Myeloid_T_sub@meta.data$celltype %>% table
Myeloid_T_sub = subset(Myeloid_T_sub, idents= c('Interstitial macrophages','SPP1+ macrophages',
                                            'Club cell','cMo', 'LAM','ncMo',
                                           'Mast cell'))

saveRDS(Myeloid_T_sub, 'Myeloid_T_subcluster2.rds')



#Myeloid_T_sub= subset(Myeloid_T_sub, idents=c(2,3,7,8,1,6))

VlnPlot(Myeloid_T_sub, features = c("IGKC",'IGHA1','JCHAIN','BANK1','MS4A1', 
                                    'CD14','CD163','CD68','MARCO','C1QA','C1QB',
                                'CD3D','CD3E','CD40LG',
                                'CD8A','CD8B', 'NCAM1','FCGR3A', 'KLRF1','PRF1'
                               
))
VlnPlot(Myeloid_T_sub, features = c('TNF','CCL3','CCL4','CCL3L1','CCL4L2','KLF2', 'NEDD9','G0S2','KLF6',
                                'ISG15','IFIT1','IFIT2',#'CCL2','CCL8', 'LAMP3',
                                'TREM2','CD9','APOE','FABP5',
                                'SELENOP','LYVE1',
                                'SPP1','VCAN','VEGFA','THBS1'
                                #'NKKB1','TNFAIP3'
))

VlnPlot(Myeloid_T_sub, Macro_mono_genes)
VlnPlot(Myeloid_T_sub, features = c( "CDKN1C",'CX3CR1','FCGR3A', 'HES4',  #ncMo
                                      'LILRA1','SIGLEC10', 'LILRB1','LILRB2' ,
                                      'S100A12','CD300E','FCN1','TCF7L2' ,'S100A8','S100A9',
                                     'TSPAB1', 'CPA3'
                                          ))
VlnPlot(Myeloid_T_sub, features = c('SCGB1A1','SCGB3A2', AT1_gene, AT2_gene))

Myeloid_T_junk = merge(x = Myeloid_T_junk, y = c(Myeloid_T_junk2, Myeloid_T_junk3))
#Myeloid_T_junk = subset(Myeloid_T, idents = c(11,18,12,9))
#Myeloid_T_junk2= subset(Myeloid_T, idents =c(0,14,4,11,15,13,7,8,2,3))
Myeloid_T_junk = SeuratAnalysis1(Myeloid_T_junk)
ElbowPlot(Myeloid_T_junk, ndims = 50)
Myeloid_T_junk = SeuratAnalysis2(Myeloid_T_junk, harmony = F, n = 10)
DimPlot(Myeloid_T_junk,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'seurat_clusters' )

cluster14 = FindMarkers(Myeloid_T_junk, ident.1 = 14, min.pct = 0.25)



VlnPlot(Myeloid_T_junk, features = c('SCGB1A1','SCGB3A1', AT2_gene,
                                'FABP4', 'MCEMP1', 'PPARG',  'GPNMB','BHLHE41'))
#check 4;
#1,5,6 macrophages; the rest are at2/club

VlnPlot(Myeloid_T_junk, features = c("IGKC",'IGHA1','JCHAIN','BANK1','MS4A1', 
                                    'CD14','CD163','CD68','MARCO','C1QA','C1QB',
                                    'CD3D','CD3E','CD40LG',
                                    'CD8A','CD8B', 'NCAM1','FCGR3A', 'KLRF1','PRF1'
                                    
)) #6 is macrophage + t cell doublet; 1 and 5 have c1q 
VlnPlot(Myeloid_T_junk, features = c('TNF','CCL3','CCL4','CCL3L1','CCL4L2','KLF2', 'NEDD9','G0S2','KLF6',
                                    'ISG15','IFIT1','IFIT2',#'CCL2','CCL8', 'LAMP3',
                                    'TREM2','CD9','APOE','FABP5',
                                    'SELENOP','LYVE1',
                                    'SPP1','VCAN','VEGFA','THBS1'
                                   
))

Myeloid_T_junk[['percent.rp']] = PercentageFeatureSet(Myeloid_T_junk, pattern = "^RP[SL]")
VlnPlot(Myeloid_T_junk, 'percent.rp' ) 
Myeloid_T_junk= subset(Myeloid_T_junk, percent.rp <= 20)

#12 and 16 are cycling 
#check 15,1,8,9,13,14 --> keep them

#Myeloid_T_junk= subset(Myeloid_T_junk, idents=c(15,1,8,9,13,14))

#structural in myeloid----------
#EPCAM and structural: 6,7,17,20 
#Myeloid_structural = subset(Myeloid, idents = c(6,7,17,20))
Myeloid_structural = SeuratAnalysis1(Myeloid_structural)
ElbowPlot(Myeloid_structural, ndims = 50)
Myeloid_structural = SeuratAnalysis2(Myeloid_structural, harmony = T, 
                            harmony_by = 'orig.ident', n = 20)
DimPlot(Myeloid_structural,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'celltype')

cluster7 = FindMarkers(Myeloid_structural, ident.1 = 7, min.pct = 0.25)
table(Myeloid_structural$seurat_clusters)

#V1: 20 IS CYCLING; 18 HAS AN INTERESTNG lnc-RNA about astrocyte and neurons; 19 may be hepatocytes; 11 has cd8+ 

#v2: no more cycling
#12 may be ciliated cells ; 

#t cells with cd3: 6,10,15 ; 15 has nk genes too ; 4 has nktr 
#6 and 10 have a lot of granzymes, ccl5, ccl4 (EFFECTOR T CELLS)
#10= CYTOCD8; 6= TEM (no obvious cd8 or cd4)

#17 like hepatocyte ; 5, 11, 8 14 like endothelial 

DotPlot(Myeloid_structural, features = c( 
                                     'CD3D','CD3E','CD40LG',
                                     'CD8A','CD8B',
                                     #'GZMB','GZMM','PRDM1','IFNG','CCL5',
                                     #'RORA','PTPRC','CCR7' #high ptprc but no ccr7 validates cluster 6 is terminal effector 
                                     'NCAM1','FCGR3A', 'KLRF1','PRF1'
                                     
), cols = 'RdYlBu')

VlnPlot(Myeloid_structural, features = c('MKI67','TOP2A','STMN1','EPCAM'))
VlnPlot(Myeloid_structural, features = c( 'GPBAR1','FPR2', 'CD300LF','LILRB2', 'SLC2A6', 'HK3', 'C19orf38', 'SLC11A1', 'FCN1')) #kupffer 
VlnPlot(Myeloid_structural, features = c(Macrophage_gene, 'SPP1','TNF','CCL3','CCL4'
                                         #'S100A12','CD300E','FCN1','EREG','CD14','IL1RN','TCF7L2'
                                          ))
DotPlot(Myeloid_structural, features = c(    'PDGFRB', 'PLN', 'RGS5', 'TAGLN', 'ACTA2', #SMC
                                           #  'IGFBP5', 'MYLK', 'RERGL', 'CASQ2', #vascular SMC
                                            # "CNN1",'MYH11', 'SYNM','ACTG2', #VASCULAR SMC 
                                             'LUM', 'DCN', 'OGN', 'POSTN', #fibro 
                                             "PECAM1",  "NOTCH4",'VWF', 'PLVAP', #general endo
                                            # 'SEMA3G', #artery 
                                            # 'LYVE1', 'PROX1',   'HSPG2', 'CLDN5',   'CLEC14A',  #lymphatic
                                           #  'CLEC4G', 'FCN2', 'CLEC4M', 'IFITM10',  'RELN', 'MRO', 'CLEC1B', 'OIT3', 'STAB2', 'FCN3',#sinusoidal
                                           #  'AIF1L', 'MGP', 'SLCO2A1',  'SPARCL1', 'CD34', 'EDNRB', 'EMCN', 'CLEC14A', #MACROVASCULAR 
                                             'PLPP2', 'CLDN4', 'FXYD2', 'FAM3B', 'SFRP5', 'SLC5A1','KRT19','B3GNT3','PRR15L','CLDN10', #cholangiocyte
                                             'HPR', 'F12',  'HGFAC', 'AGXT', 'SAA4','ARG1',#HEPATOCYTE
                                             'EPCAM' #epithelial
                                             ), cols = 'RdYlBu') + coord_flip()
VlnPlot(Myeloid_structural, features = c( "PECAM1",  "NOTCH4",'VWF', 'PLVAP', #general endo
                                          'LUM', 'DCN', 'OGN', 'POSTN' #general fibro 
                                          ))
VlnPlot(Myeloid_structural, features = c( 'ALDOB' , 'SLC34A1','PDZK1IP1','ANPEP',
                                          'CA9', 'NNMT','NDUFA4L2',
                                          "CADM1",  'NFASC','NRXN1','GRIK2' #schwann from lung
)) #PROXIMAL TUBULE CELLS; cancer cells FROM KIDNEY 

DotPlot(Myeloid_structural, features = c( "FOXJ1", 'RFX2', 'CFAP52', 'PIFO','DCDC2B','DTHD1','GSTA1',
                                          'AGER','RTKN2','CLIC5','HOPX', 'CAV1', #at1
                                          'SFTPA2','SFTPB', 'MUC1', 'SFTPC', #at2 
                                          'EPCAM', 'CDH1', 'KRT18','KRT19','KRT14', 'KRT15', 'KRT5', 'KRT17', #epithelial
                                          'SCGB1A1','SCGB3A2'
                                          ), cols = 'RdYlBu') + coord_flip()
#8 = club cells 
#0,1,3,6,9,12 = ciliated cells 
#13 = at2 ; 11 has high krt18 
#17 = hepatocyte 
#5 (c1q) and 14 = cmo (fcn1, vcan, s100a8) are macrophages/monocytes ; 2 is spp1+ macrophage 
#7 is junk ; 11 has the three genes from the kindey paper of kidney cancer cells 
#don't know what 4 and 16 are

cluster.ids = c('Ciliated cells','Ciliated cells', 'SPP1+ macrophages',
                'Ciliated cells', 'Unknown', 'Inflammatory macrophages',
                'Tem', 'Unknown', 'Club cells', 'Ciliated cells', 'CytoCD8',
                'Kidney cancer cells', 'Ciliated cells', 'Alveolar cell type 2',
                'cMo','NKT', 'Unknown', 'Hepatocyte')
names(cluster.ids) <- levels(Myeloid_structural)
Myeloid_structural <- RenameIdents(Myeloid_structural, cluster.ids)
DimPlot(Myeloid_structural, reduction = "umap", label = TRUE, repel = TRUE)  

Myeloid_structural$celltype = Idents(Myeloid_structural)
Myeloid_structural@meta.data$celltype %>% table

Myeloid_structural_sub2 = subset(Myeloid_structural, idents = c('Ciliated cells', 'SPP1+ macrophages',
                                                          'Inflammatory macrophages',
                                                           'Tem', 'Club cells', 'CytoCD8',
                                                           'Kidney cancer cells', 'Alveolar cell type 2',
                                                           'cMo','Hepatocyte'))

#Myeloid_structural_sub = subset(Myeloid_structural, idents = c(4, 16))
#Myeloid_structural = subset(Myeloid_structural, idents = 0:19)

#Myeloid_structural_sub = subset(Myeloid_structural, idents = c(4, 16))
Myeloid_structural_sub = SeuratAnalysis1(Myeloid_structural_sub)
ElbowPlot(Myeloid_structural_sub, ndims = 50)
Myeloid_structural_sub = SeuratAnalysis2(Myeloid_structural_sub, harmony = T, 
                                     harmony_by = 'orig.ident', n = 10)
DimPlot(Myeloid_structural_sub,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'seurat_clusters')

VlnPlot(Myeloid_structural_sub, features = c('MKI67','TOP2A','STMN1','EPCAM'))
DotPlot(Myeloid_structural_sub, features = c( "FOXJ1", 'RFX2', 'CFAP52', 'PIFO','DCDC2B','DTHD1','GSTA1',
                                          'AGER','RTKN2','CLIC5','HOPX', 'CAV1', #at1
                                          'SFTPA2','SFTPB', 'MUC1', 'SFTPC', #at2 
                                          'EPCAM', 'CDH1', 'KRT18','KRT19','KRT14', 'KRT15', 'KRT5', 'KRT17', #epithelial
                                          'SCGB1A1','SCGB3A2'
), cols = 'RdYlBu') + coord_flip()
#0,1,2,3 all ciliated 
#4 and 6 are epithelial in kidney; 5 is hepatocyte  
#reference the tissue specific lineage spreadsheet that 6 is collecting duct intercalated cells 

cluster.ids = c('Ciliated cells','Ciliated cells','Ciliated cells','Ciliated cells',
                'Renal epithelial cell','Hepatocyte', 'Collecting duct intercalated cell'
              )                
names(cluster.ids) <- levels(Myeloid_structural_sub)
Myeloid_structural_sub <- RenameIdents(Myeloid_structural_sub, cluster.ids)
DimPlot(Myeloid_structural_sub, reduction = "umap", label = TRUE, repel = TRUE)  

Myeloid_structural_sub$celltype = Idents(Myeloid_structural_sub)
Myeloid_structural_sub@meta.data$celltype %>% table


cluster6 = FindMarkers(Myeloid_structural_sub, ident.1 = 6, min.pct = 0.25)
 
Myeloid_structural_final = merge(Myeloid_structural_sub, Myeloid_structural_sub2)
saveRDS(Myeloid_structural_final, file = 'Myeloid_structural_subcluster_all.rds')

#myeloid unknown-----------
#unknown remaining: 10, 13, 23
#Myeloid_unknown = subset(Myeloid, idents = c(10, 13, 23))
Myeloid_unknown = SeuratAnalysis1(Myeloid_unknown)
ElbowPlot(Myeloid_unknown, ndims = 50)
Myeloid_unknown = SeuratAnalysis2(Myeloid_unknown, harmony = T, 
                                     harmony_by = 'orig.ident', n = 10)
DimPlot(Myeloid_unknown,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'orig.ident')
 

VlnPlot(Myeloid_unknown, features = c('MKI67','TOP2A','STMN1','EPCAM',
                                      'ALDOB','APOA1','ALB', 'ALDOA','AKR1B10'))
VlnPlot(Myeloid_unknown, features = c(Macrophage_gene, 
                                      'TNF','CCL3','CCL4','CCL3L1','CCL4L2',#'KLF2', 'NEDD9','G0S2','KLF6',
                                      'ISG15','IFIT1','IFIT2',#'CCL2','CCL8', 'LAMP3',
                                      #'TREM2','CD9','APOE','FABP5',
                                      #'SELENOP','LYVE1',
                                      'SPP1','VCAN','VEGFA','THBS1')) 

#11,12,17 may be macrophages ; 7 may be monocytes 
#3 = t cells 
#4 = at2, 9 and 10 = epi,
#14 and 16 may be hepatocyte or any liver structural cells  
  

#v2:
#7 and 9 have high alb, 5,9,7 all HCC ; 5 = cholangiocyte ; 7 IS HCC CANCER CELL WITH AKR1B10 
#0,4 has cd3, ccr7, gzmB, GZMM, low SELL/TCF7, low CX3CCR1, prdm1, ifng, ccl5--> CYTOCD8 
#1 = at2; 5 = epi, 2 = NCMO 
#3, 6, 10 = mac (10 has vegfa, vcan)

cluster.ids = c('CytoCD8', 'Alveolar cell type 2', 'ncMo', 'Macrophage','CytoCD8',
                'Cholangiocyte', 'Macrophage', 'Liver cancer cell', 'Unknown',
                'Hepatocyte', 'Angiogenic macrophage'
)
names(cluster.ids) <- levels(Myeloid_unknown)
Myeloid_unknown <- RenameIdents(Myeloid_unknown, cluster.ids)
DimPlot(Myeloid_unknown, reduction = "umap", label = TRUE, repel = TRUE)  

Myeloid_unknown$celltype = Idents(Myeloid_unknown)
Myeloid_unknown@meta.data$celltype %>% table

Myeloid_unknown = subset(Myeloid_unknown, idents =  c('CytoCD8', 'Alveolar cell type 2', 'ncMo', 
                                                      'Cholangiocyte', 'Macrophage', 'Liver cancer cell', 
                                                      'Hepatocyte', 'Angiogenic macrophage'
))
saveRDS(Myeloid_unknown, 'Myeloid_unknown_subcluster1.rds')


VlnPlot(Myeloid_unknown, features = c("FCN1",'S100A8','S100A9', 'FCGR3A','HES4','CX3CR1',"CDKN1C", 
                                      'CD14','CXCL8','IL1B',
                                      'ICAM1', 'MPO','CXCR4','NAMPT','CSF3R','CTSG','ELANE'
                                      ))

DotPlot(Myeloid_unknown, features = c( "FOXJ1", 'RFX2', 'CFAP52', 'PIFO','DCDC2B','DTHD1','GSTA1',
                                              'AGER','RTKN2','CLIC5','HOPX', 'CAV1', #at1
                                              'SFTPA2','SFTPB', 'MUC1', 'SFTPC', #at2 
                                              'EPCAM', 'CDH1', 'KRT18','KRT19','KRT14', 'KRT15', 'KRT5', 'KRT17', #epithelial
                                              'SCGB1A1','SCGB3A2',
                                              'HPR', 'F12',  'HGFAC', 'AGXT', 'SAA4','ARG1' 
), cols = 'RdYlBu') + coord_flip()
VlnPlot(Myeloid_unknown, features = c( 
  'CD3D','CD3E','CD40LG',
  'CD8A','CD8B',
  'GZMB','GZMK','ITGB1','PRDM1','IFNG','CCL5','CX3CR1',
  'RORA','PTPRC','CCR7', 'SELL','TCF7'  #high ptprc but no ccr7 validates cluster 6 is terminal effector 
 # 'NCAM1','FCGR3A', 'KLRF1','PRF1'
  ) )
DotPlot(Myeloid_unknown, features = c(    'PLPP2', 'CLDN4', 'FXYD2', 'FAM3B', 'SFRP5', 'SLC5A1','KRT19','B3GNT3','PRR15L','CLDN10', #cholangiocyte
                                             'HPR', 'F12',  'HGFAC', 'AGXT', 'SAA4','ARG1',#HEPATOCYTE
                                             'EPCAM' #epithelial
), cols = 'RdYlBu') + coord_flip()

 
VlnPlot(Myeloid_unknown, features = c( 'GPBAR1','FPR2', 'CD300LF','LILRB2', 'SLC2A6', 'HK3', 'C19orf38', 'SLC11A1', 'FCN1')) #kupffer \
VlnPlot(Myeloid_unknown, features = c( 'LYZ','CCR2','S100A8','S100A9', #VERY WEAK, BUT PRESENT 
                                       'FCER1A','CLEC10A','ITGAX'
                                       ))  
VlnPlot(Myeloid_unknown, features = c('FXYD3','SFN','TACSTD2','S100A14','KRT6A','ACTA1'
)) #EPI FROM HNSCC 

cluster9 = FindMarkers(Myeloid_unknown, ident.1 = 9, min.pct = 0.25)
table(Myeloid_unknown$seurat_clusters)
#Myeloid_unknown = subset(Myeloid_unknown, idents = c(0,1,2,4:8))




#Myeloid_unknown_sub = subset(Myeloid_unknown, idents = c(0,1,2,5,6,8,13,15))
Myeloid_unknown_sub = SeuratAnalysis1(Myeloid_unknown_sub)
ElbowPlot(Myeloid_unknown_sub, ndims = 50)
Myeloid_unknown_sub = SeuratAnalysis2(Myeloid_unknown_sub, harmony = F, n = 20)
DimPlot(Myeloid_unknown_sub,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'seurat_clusters')

#v2:
#8 is a single cluster in kidney, seems to be macrophage (lgmn, apoe, IFIT2, CCL3, some aldoa, icoslg )
#4= mono-cdc2 ; 2 = at2 or club; 6 = epi; 7 may be at1 (contaminate with ig)
#5 IS IL12B+WDFY4+ CDC1; 0 = mac 

cluster.ids = c("Macrophage", 'Unknown', 'Unknown', 'Unknown', 'Monocyte-DC2', 'IL12B+ activated cDC1',
                'Unknown', 'Unknown','Macrophage')
names(cluster.ids) <- levels(Myeloid_unknown_sub)
Myeloid_unknown_sub <- RenameIdents(Myeloid_unknown_sub, cluster.ids)
DimPlot(Myeloid_unknown_sub, reduction = "umap", label = TRUE, repel = TRUE)  

Myeloid_unknown_sub$celltype = Idents(Myeloid_unknown_sub)
Myeloid_unknown_sub@meta.data$celltype %>% table

Myeloid_unknown_sub = subset(Myeloid_unknown_sub, idents =  c('Macrophage', 'Monocyte-DC2', 'IL12B+ activated cDC1'))


VlnPlot(Myeloid_unknown_sub, features = c('MKI67','TOP2A','STMN1','EPCAM',
                                      'ALDOB','APOA1','ALB', 'ALDOA','AKR1B10'))
DotPlot(Myeloid_unknown_sub, features = unique(c(Macrophage_gene, 
                                     Monocyte_gene,'HLA-DRA')), cols = 'RdYlBu') + coord_flip() 


cluster6 = FindMarkers(Myeloid_unknown_sub, ident.1 = 6, min.pct = 0.25)
 
VlnPlot(Myeloid_unknown_sub, features = c('CLEC4F','TIMD4','VSIG4','CLEC1B', 'GPHBP1',
                                          'FCGR2A')) #OTHER KUPFFER CELL GENES FROM ONLINE 
VlnPlot(Myeloid_unknown_sub, features = c('ZBTB46','CLEC9A', 'XCR1', 'IL12B','MPEG1','WDFY4', 'CCR7', 'CD3D','CD3E'))
VlnPlot(Myeloid_unknown_sub, features = c('TNF','CCL3','CCL4','KLF2','KLF6',
                                          'ISG15','IFIT1','IFI27','CCL2','CCL8',
                                          'CXCL10','CXCL9'))

VlnPlot(Myeloid_unknown_sub, features = c( 'FCN1','LYZ','CCR2','S100A8','S100A9','CXCL8','CD14', #VERY WEAK, BUT PRESENT 
                                       'FCER1A','CLEC10A','ITGAX','CD1C', 'CD1B', 'CCL22','ISG15'
))  

DotPlot(Myeloid_unknown_sub, features = c( "FOXJ1", 'RFX2', 'CFAP52', 'PIFO','DCDC2B','DTHD1','GSTA1',
                                       'AGER','RTKN2','CLIC5','HOPX', 'CAV1', #at1
                                       'SFTPA2','SFTPB', 'MUC1', 'SFTPC', #at2 
                                       'EPCAM', 'CDH1', 'KRT18','KRT19','KRT14', 'KRT15', 'KRT5', 'KRT17', #epithelial
                                       'SCGB1A1','SCGB3A2',
                                       'HPR', 'F12',  'HGFAC', 'AGXT', 'SAA4','ARG1' 
), cols = 'RdYlBu') + coord_flip()

#Myeloid_unknown_sub = subset(Myeloid_unknown_sub, idents = c(0:3, 5:9))
  

#merge all myeloid and recheck lineages -----------
#Myeloid_final = merge(x = Macrophage, y = c(cDC1, cDC2, Monocyte, Myeloid_T1, Myeloid_T2, 
                                           # Myeloid_structural, Myeloid_unknown1, Myeloid_unknown_sub))
#Myeloid_final = readRDS("~/jliang/Human tissue harmonization/Myeloid/Myeloid_final_v7.rds")
Myeloid_final = JoinLayers(Myeloid_final)
Myeloid_final = SeuratAnalysis1(Myeloid_final)
ElbowPlot(Myeloid_final, ndims = 50)
Myeloid_final = SeuratAnalysis2(Myeloid_final, harmony = T, 
                                harmony_by = 'orig.ident', n = 20)
DimPlot(Myeloid_final,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'orig.ident')  

#rename orig ident
Myeloid_final$orig.ident= gsub("Breast",'BRCA', Myeloid_final$orig.ident)
Myeloid_final$orig.ident= gsub("hnscc",'HNSCC', Myeloid_final$orig.ident)
Myeloid_final$orig.ident= gsub("Kidney",'KIRC', Myeloid_final$orig.ident) 


#get the dataframe of umap behind and label the coordinates by celltype? 
umap_data <- FetchData(object = Myeloid_final, 
                       vars = c("umap_1", "umap_2", "orig.ident", "celltype")) #this is the key
head(umap_data)

label_centers <- umap_data %>% #create another dataframe telling where the label should be 
  group_by(celltype) %>%
  summarize(umap_1 = median(umap_1), umap_2 = median(umap_2))

ggplot(umap_data, aes(x = umap_1, y = umap_2)) +
  # Color the dots by your cancer condition
  geom_point(aes(color = orig.ident), size = 0.3, alpha = 1) + 
  # Add the manual labels at the calculated center points
  geom_text(data = label_centers, aes(label = celltype), 
            color = "black",  check_overlap = TRUE) +
  theme_bw() +
  labs(title = "All myeloid cells in 6 human cancers",
       x = "UMAP 1", y = "UMAP 2") +
  guides(color = guide_legend(override.aes = list(size = 3))) + 
  theme(plot.title = element_text(size = 15, hjust =0.5),
        panel.grid = element_blank(),
        panel.border = element_blank(), #get rid of the box 
        axis.line.x = element_line(linetype = "solid", colour = "black"),
        axis.line.y = element_line(linetype = "solid", colour = "black")) 


table(Myeloid_final$celltype)
Idents(Myeloid_final) ='celltype'
#Myeloid_final = RenameIdents(Myeloid_final,'SPP1+ macrophages' ='CXCL9/10/11+ LAM')
Myeloid_final = RenameIdents(Myeloid_final,'Inflammatory macrophages' ='Interstitial macrophages', 
                             'Inflammatory cDC2b' ='cDC2b','Activated cDC1'= 'Migratory DC1',
                             'MigDC' ='Migratory DC1')
Myeloid_final$celltype= Idents(Myeloid_final)


#Myeloid_final = subset(Myeloid_final, idents = unique(Myeloid_final$celltype)[c(1:24, 27:41)])
#Myeloid_final = subset(Myeloid_final, idents = unique(Myeloid_final$celltype)[c(1,2,4:10, 13,14,16:24,
                                                                               # 27, 29, 36, 37, 40, 41)])
#Myeloid_rest = subset(Myeloid_final, idents = unique(Myeloid_final$celltype)[c(1,2,3,5:11, 13:20, 22:24)])
#Myeloid_rest = subset(Myeloid_final, idents = unique(Myeloid_final$celltype)[c(1:19, 22:24)])
#Myeloid_final = merge(Myeloid_rest, Myeloid_mac)
#Myeloid_final = merge(Myeloid_rest, Myeloid_T)
#Myeloid_final = subset(Myeloid_final, idents = unique(Myeloid_final$celltype)[c(1:21, 23)])
#Myeloid_final = subset(Myeloid_final, idents = unique(Myeloid_final$celltype)[c(1:21)])
Myeloid_final$barcodes = rownames(Myeloid_final@meta.data)
Idents(Myeloid_final) ='barcodes'

#Myeloid_rest = subset(Myeloid_final, idents = Myeloid_final$barcodes [!Myeloid_final$barcodes %in% rownames(Myeloid_sub@meta.data)])
#Myeloid_final = merge(Myeloid_rest, Myeloid_sub_final)
#Myeloid_final = subset(Myeloid_final, idents = unique(Myeloid_final$celltype)[1:19])
#Myeloid_final = subset(Myeloid_final, idents = 0:19)
#Myeloid_rest = subset(Myeloid_final, idents = Myeloid_final$barcodes [!Myeloid_final$barcodes %in% rownames(Kupffer@meta.data)])
Myeloid_rest = subset(Myeloid_final, idents = unique(Myeloid_final$celltype)[c(2:21)])


#put back the neutrophils and kuppfer cells back 
Neutrophil_metadata = Neutrophil_final@meta.data[, 1:11]
Kupffer_metadata = Kupffer@meta.data[, 1:11]
neutrophil_export = CreateSeuratObject(counts = LayerData(Neutrophil_final, assay = "RNA", layer = "counts"), meta.data= Neutrophil_metadata)
kupffer_export = CreateSeuratObject(counts = LayerData(Kupffer, assay = "RNA", layer = "counts"), meta.data= Kupffer_metadata)
add = merge(neutrophil_export, kupffer_export)
head(add@meta.data)
Idents(add)='celltype'
add$celltype = Idents(add)

LAM_metadata = LAM@meta.data[, 1:11]
LAM_export = CreateSeuratObject(counts = LayerData(LAM, assay = "RNA", layer = "counts"), meta.data= LAM_metadata)
Myeloid_final =merge(Myeloid_rest, LAM_export)


#Myeloid_rest = subset(Myeloid_final, idents = Myeloid_final$barcodes [!Myeloid_final$barcodes %in% rownames(add@meta.data)])
Myeloid_rest = JoinLayers(Myeloid_rest)
#spp1_metadata= Spp1@meta.data[, 1:11]
#add =CreateSeuratObject(counts = LayerData(Spp1, assay = "RNA", layer = "counts"), meta.data= spp1_metadata)
Interstitial_metadata = Interstitial@meta.data[, 1:11]
Interstitial_export =CreateSeuratObject(counts = LayerData(Interstitial, assay = "RNA", layer = "counts"), meta.data= Interstitial_metadata)
DC_metadata = DC1_Mig@meta.data[, 1:11]
DC_export =CreateSeuratObject(counts = LayerData(DC1_Mig, assay = "RNA", layer = "counts"), meta.data= DC_metadata)

Myeloid_final = merge(x = Myeloid_rest, y = c(Interstitial_export, DC_export))

saveRDS(Myeloid_final, file = 'Myeloid_final_v7_rename.rds')

Myeloid_final$cluster_celltype = paste0(Myeloid_final$seurat_clusters,"_", Myeloid_final$celltype)
table(Myeloid_final$cluster_celltype)
table(Myeloid_final$seurat_clusters)

Myeloid_final$celltype_condition = paste0(Myeloid_final$celltype,"_", Myeloid_final$orig.ident)
table(Myeloid_final$celltype_condition)


cluster20 = FindMarkers(Myeloid_final, ident.1 = 20, min.pct = 0.25)

Myeloid_0401_all_genes = FindAllMarkers(Myeloid_final, group.by = 'celltype', min.pct = 0.25,
                                         only.pos = TRUE)  
Myeloid_0401_top50 = Myeloid_0401_all_genes %>% 
  group_by(cluster) %>% 
  slice_max(order_by = avg_log2FC, n = 50, with_ties = T)
write.csv(Myeloid_0401_all_genes, file = 'Myeloid_all_cell_all_findmarkers_v7_rename.csv')

Idents(Myeloid_final) ='celltype'
Mac_genes = FindMarkers(Myeloid_final, ident.1 = 'Macrophages', min.pct = 0.25)



Myeloid_final = RenameIdents(Myeloid_final, 'IL12B+ activated cDC1' ='Activated cDC1', 'Monocyte-DC2'='Migratory DC2',
                             'Activated cDC2b'='Migratory DC2')
Myeloid_final = RenameIdents(Myeloid_final, 'LAMP3+ DC' ='Migratory DC2')
Myeloid_final$celltype = Idents(Myeloid_final)

DotPlot(Myeloid_final, features = c( 
  'CD3D','CD3E','CD40LG',
  'CD8A','CD8B','ITGB1','IFNG',
  'CCR7', 'TCF7','SELL',
  'GZMB','GZMH','PRDM1','KLRG1','KLF2',  'RORA', 'CX3CR1', #TEM
  'ITGAE','ITGA1', 'TNIP3', #TRM 
  'NCAM1','FCGR3A', 'KLRF1','PRF1'
), cols = 'RdYlBu', group.by = 'celltype') + coord_flip() + RotatedAxis() #so no cx3cr1, low cd11b, 

VlnPlot(Myeloid_final, features = c())
 

#DC: seems to all makes sense; 
DotPlot(Myeloid_final, features = unique(c( 
  DC_gene,'ITGAM', 'ITGAX' ,'HLA-DRA','SIRPA',
  'MPEG1','WDFY4','BATF3','IRF8',
  'CCR7','IL12B','FSCN1',
  'CD274','PDCD1LG2', 'LAMP3','CD40','CD80','CD86', 'EBI3','CCL22', #more mig dc genes 
  'RELB','NFKB1','STAT1', #inflammatory 
  'CTLA4','ICOSLG',  'IDO1', #regulatory 
  'CD207','CD1A', 'IL22RA2' #langerin 
)), cols = 'RdYlBu', group.by = 'celltype') + coord_flip() + RotatedAxis()

DotPlot(Myeloid_final, features = unique(c( 
  'MPEG1','WDFY4','BATF3','IRF8',
  'CCR7','LAMP3','IL12B'
)), cols = 'RdYlBu', group.by = 'celltype') + coord_flip() + RotatedAxis()


VlnPlot(Myeloid_final, features = c(
  'ITGAM', 'ITGAX','SIRPA',
  'LAMP3', 'TBX21', 'ESAM','TCF4','CX3CR1' #trying to find cdc2a with esam and tbet but turns out they are more like migDC
), group.by = 'celltype') 

#check if LAMP3+ DC also has Tbet, tcf4, ESAM, low itgam, low cx3cr1--> i think it's more like a migdc2 instead of lamp3+ dc2 
#macrophages: just the "angiogenic" and "macrophage" do not make sense 

p1=DotPlot(Myeloid_final, c(Macrophage_gene[1:5],'C1QA','C1QB', 
                         'TNF','CCL3','CCL4','KLF2','KLF6','IL1RN',#'NEDD9',
                     'ISG15','IFIT1','IFI27','CCL8',
                     'CXCL10','CXCL9',
                     'TREM2','CD9','APOE','FABP5',
                     'LYVE1', 'SELENOP','LGMN', 'STAB1', 'MAN1A1', 'LILRB5', 'MS4A6A',  #PVM
                     AlveolarMac_gene, 'GPNMB','BHLHE41', #alveolar
                     'SPP1', 'VCAN','VEGFA','THBS1' , #ANGIOGENIC AND SPP1++ : VEGFA SHOULD BE HIGH IN ANGIOGENIC
                     'TIMD4','CD5L', 'VCAM1','HMOX1', 'CD36','CSF1R' #kupffer  
                     ), cols = 'RdYlBu', group.by = 'celltype')+ coord_flip() + RotatedAxis()

behind = p1$data
behind =behind %>% filter(id %in% c('Alveolar macrophage', 'IFN-primed macrophages', 
                                    'Inflammatory macrophages', 
                                    'Kupffer cells', 'LAM','PVM','SPP1+ LAM', "CXCL9/10/11+ LAM"))
p1$data = behind
p1

DotPlot(Myeloid_final, features = c('SPP1', 'VCAN','VEGFA','THBS1'), cols = 'RdYlBu', group.by = 'celltype')+ coord_flip() + RotatedAxis()

#MONOCYTE:
p2= DotPlot(Myeloid_final, features = c( "CDKN1C",'CX3CR1', 'HES4',  'MS4A7', # ncMo
                                     'CD14', 'FCGR3A',  #relative expression of these 2 
                                  'LILRA1','SIGLEC10', 'LILRB1','LILRB2' , , #ImO
                                  'S100A8','S100A9','CCR2','S100A12','CD300E' #'FCN1','EREG','IL1RN','TCF7L2' , #iMo  
                                    #'HLA-DRA', 'HLA-DPB1' 
                                  ), cols = 'RdYlBu', group.by = 'celltype')+ coord_flip() + RotatedAxis()
 
behind = p2$data
behind =behind %>% filter(id %in% c('cMo','iMo','ncMo'))
p2$data = behind
p2


#FROM DIFFERENT PAPERS, SEEM LIKE CD5L & TIMD4 most important for kupffer cells 
#kupffer seems to be within inflammatory mac 
VlnPlot(Myeloid_final, features = c('CLEC4F','TIMD4','CD5L', 'VCAM1','HMOX1', 'CD36','CSF1R', #'MRC1','MARCO', #'CETP','VSIG4','CLEC1B', 'GPHBP1', 'FCGR2A',
                                    'CCL2','TNIP3', #intermediate 
                                    'SEPP1' #KIDNEY-RESIDENT MAC? 
                                    ), group.by='seurat_clusters')
 #15,19,20 has kupffer!! 
VlnPlot(Myeloid_final, features = c(   'SEPP1', 'CD74','CD81','CSF2RA' #KIDNEY-RESIDENT MAC? 
), group.by='seurat_clusters')

#neutrophils:
DotPlot(Myeloid_final, features = c('SPARC','CSF3R')) #from andrew's il12b paper 
VlnPlot(Myeloid_final, features = c('LCN2', 'ORM1','MMP8','PAGE5','SYCP1','SYN3', 'CTAG2','DEFA3'
), group.by='seurat_clusters')
DotPlot(Myeloid_final, features = c('S100A8','S100A9','S100A6', 'CXCR4','CXCR2','CTSG','ELANE',
                                 'CD44','MPO','ZEB2', 'CST3','CSF3R','FCGR3B','CDA' ,'LY6D'
), group.by='celltype', cols = 'RdYlBu')+ coord_flip() + RotatedAxis()
VlnPlot(Myeloid_final, features = c('VEGFA','PLAUR', 'OLR1','CXCR2'
), group.by='celltype')
VlnPlot(Myeloid_final, features = c('CLC', 'CCR3', 
                                    'CXCR2','FCGR3B', 'OLR1','CD181' #eosinophil and neutrophil 
), group.by='seurat_clusters')


#still need to check spp1+ mac (have cd3 genes); differentiate cmo vs imo vs ncmo 


#structural_in_myeloid = subset(Myeloid_final, idents = unique(Myeloid_final$celltype)[c(3, 11, 12, 26, 28:33,
                                                                                      #  36,37)])
Idents(Myeloid_final)='celltype'
ciliated_in_myeloid = subset(Myeloid_final, idents = unique(Myeloid_final$celltype)[c(21,22)])
structural_in_myeloid= merge(structural_in_myeloid, ciliated_in_myeloid)
##hepatocyte in Spp1 subcluster!! 
saveRDS(structural_in_myeloid, file = 'Structural_cells_in_myeloid_all.rds')
saveRDS(ciliated_in_myeloid, file = 'Structural_cells_in_myeloid2.rds')

#clean up the T cells in merged myeloid ------
#Myeloid_T= subset(Myeloid_final, idents = c('CytoCD8-Interstitial macrophages','CytoCD8','Tem'))
Myeloid_T = SeuratAnalysis1(Myeloid_T)
ElbowPlot(Myeloid_T, ndims = 50)
Myeloid_T = SeuratAnalysis2(Myeloid_T, harmony = T, harmony_by = 'orig.ident', n = 10)
DimPlot(Myeloid_T,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'celltype')  


DotPlot(Myeloid_T, features = c( 
  'CD3D','CD3E','CD40LG',
  'CD8A','CD8B','ITGB1','IFNG',
  'CCR7', 'TCF7','SELL',
  'GZMB','GZMH','PRDM1','KLRG1','KLF2',  'RORA', 'CX3CR1', #TEM
  'ITGAE','ITGA1', 'TNIP3', #TRM 
  'NCAM1','FCGR3A', 'KLRF1','PRF1'
), cols = 'RdYlBu') + coord_flip()
VlnPlot(Myeloid_T, c(Macrophage_gene, 'LY6G6C', 'FCN1','S100A8'))
VlnPlot(Myeloid_T, c('TNF','CCL3','CCL4','KLF2','KLF6','IL1RN','NEDD9',
                     'ISG15','IFIT1','IFI27','CCL2','CCL8',
                     #'CXCL10','CXCL9',
                     'TREM2','CD9','APOE','FABP5'))
VlnPlot(Myeloid_T, c(DC_gene, 'IL12B','SIRPA','ITGAX','ZBTB46'))
VlnPlot(Myeloid_T, features = c( 'CCR7','FSCN1', 'CD274','CD86','CD80','LAMP3', 'EBI3', 'IDO1','CCL22'))
DotPlot(Myeloid_T, features = c( "FOXJ1", 'RFX2', 'CFAP52', 'PIFO','DCDC2B','DTHD1','GSTA1',
                                           'AGER','RTKN2','CLIC5','HOPX', 'CAV1', #at1
                                           'SFTPA2','SFTPB', 'MUC1', 'SFTPC', #at2 
                                           'EPCAM', 'CDH1', 'KRT18','KRT19','KRT14', 'KRT15', 'KRT5', 'KRT17', #epithelial
                                           'SCGB1A1','SCGB3A2',
                                           'HPR', 'F12',  'HGFAC', 'AGXT', 'SAA4','ARG1' 
), cols = 'RdYlBu') + coord_flip()

cluster9 = FindMarkers(Myeloid_T, ident.1 = 9, min.pct = 0.25)

#4,5,6 all ciliated, 9 can be at2 mixed with at1 (doublet)
#8 = cdc1; 1= cdc2; 2 and 11 with ccr7 and lamp3 (MIG DC)
#0, 3, 10 = cytocd8 (doublet with mac);  7 unknown 
#0 = LAM 

cluster.ids = c('LAM', 'cDC2', 'MigDC', 'LAM', 'Ciliated cells',
                'Ciliated cells','Ciliated cells',
                'Unknown', 'cDC1', 'Unknown', 'Unknown', 'MigDC')
names(cluster.ids) <- levels(Myeloid_T)
Myeloid_T <- RenameIdents(Myeloid_T, cluster.ids)
DimPlot(Myeloid_T, reduction = "umap", label = TRUE, repel = TRUE)  

Myeloid_T$celltype = Idents(Myeloid_T)
Myeloid_T@meta.data$celltype %>% table

Myeloid_T = subset(Myeloid_T, idents =  c('LAM', 'cDC2', 'MigDC',  'Ciliated cells',
                                           'cDC1'))
saveRDS(Myeloid_T, 'Myeloid_T_subcluster_final.rds')

#Myeloid_T = subset(Myeloid_T, idents = c(0:11))
#Myeloid_T = subset(Myeloid_T, idents = c(0, 2:10))


#clean up the macrophages in merged macrophages ------ 
Myeloid_mac = subset(Myeloid_final, idents = c("Macrophage", 'Angiogenic macrophage'))
Myeloid_mac = SeuratAnalysis1(Myeloid_mac)
ElbowPlot(Myeloid_mac, ndims = 50)
Myeloid_mac = SeuratAnalysis2(Myeloid_mac, harmony = F, n = 25)
DimPlot(Myeloid_mac,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'seurat_clusters')  

VlnPlot(Myeloid_mac, c(Macrophage_gene,'SPP1','VEGFA','VCAN','THBS1',
                       'CCL2','TNIP3')) #NOTHING FOR THE INTERMEDIATE GENE 

cluster5 = FindMarkers(Myeloid_mac, ident.1 = 5, min.pct = 0.25)
cluster.ids = c("Unknown","Unknown","Unknown","Unknown",'Macrophage','Macrophage')
names(cluster.ids) <- levels(Myeloid_mac)
Myeloid_mac <- RenameIdents(Myeloid_mac, cluster.ids)
DimPlot(Myeloid_mac, reduction = "umap", label = TRUE, repel = TRUE)  

Myeloid_mac$celltype = Idents(Myeloid_mac)
Myeloid_mac@meta.data$celltype %>% table

Myeloid_mac = subset(Myeloid_mac, idents =  'Macrophage')
 
#clean up cluster 7 and cDC2 from myeloid merged ---- 
Myeloid_sub = subset(Myeloid_final, seurat_clusters == '7' | celltype == 'cDC2')
Myeloid_sub = SeuratAnalysis1(Myeloid_sub)
ElbowPlot(Myeloid_sub, ndims = 50)
Myeloid_sub = SeuratAnalysis2(Myeloid_sub, harmony = F, n = 15)
DimPlot(Myeloid_sub,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'celltype') 


VlnPlot(Myeloid_sub, c(Macrophage_gene,'SPP1','VEGFA','VCAN','THBS1',
                       'TNF','CCL3','CCL4','CCL3L1','CCL4L1'
                       #'CD3D','CD3E',
                       #'SIRPA','CLEC10A','CD1C'
                       
                       ))  
VlnPlot(Myeloid_sub, c('ITGAM', 'ITGAX', 'ZBTB46', DC_gene,
                       'CD207', 'CD1A'))
DotPlot(Myeloid_sub, features = c( "FOXJ1", 'RFX2', 'CFAP52', 'PIFO','DCDC2B','DTHD1','GSTA1',
                                 'AGER','RTKN2','CLIC5','HOPX', 'CAV1', #at1
                                 'SFTPA2','SFTPB', 'MUC1', 'SFTPC', #at2 
                                 'EPCAM', 'CDH1', 'KRT18','KRT19','KRT14', 'KRT15', 'KRT5', 'KRT17', #epithelial
                                 'SCGB1A1','SCGB3A2',
                                 'HPR', 'F12',  'HGFAC', 'AGXT', 'SAA4','ARG1' 
), cols = 'RdYlBu') + coord_flip()



#5,6,7 all doublet junk 
cluster3 = FindMarkers(Myeloid_sub, ident.1 = 3, min.pct = 0.25)

#4 is a cdc2 ; 8 has il22ra2 + CD1C BUT NO CLEC10A  (IL22RA2+ Langerhans)   
#1 macrophages with c1q ; 0 is ciliated ; 2 = at2 

cluster.ids = c("Ciliated cells", 'Macrophages', 'Alveolar cells type 2',
                'Unknown', 'cDC2b','Unknown', 'Unknown', 'Unknown', 
                'IL22RA2+ Langerhans', 'Unknown')
names(cluster.ids) <- levels(Myeloid_sub)
Myeloid_sub <- RenameIdents(Myeloid_sub, cluster.ids)
DimPlot(Myeloid_sub, reduction = "umap", label = TRUE, repel = TRUE)  

Myeloid_sub$celltype = Idents(Myeloid_sub)
Myeloid_sub@meta.data$celltype %>% table

Myeloid_sub_final = subset(Myeloid_sub, idents =  c("Ciliated cells", 'Macrophages', 'Alveolar cells type 2',
                                              'cDC2b',
                                              'IL22RA2+ Langerhans'))

rownames(Myeloid_sub@meta.data) %in% rownames(Myeloid_final@meta.data) %>% sum()

#kupffer and neutrophil----
Kupffer = subset(Myeloid_final, idents = c(15,19))
Kupffer = SeuratAnalysis1(Kupffer)
ElbowPlot(Kupffer, ndims = 50)
Kupffer = SeuratAnalysis2(Kupffer, harmony = F, n = 25)
DimPlot(Kupffer,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'seurat_clusters') 

#0,1,2,3 all from HCC 

VlnPlot(Kupffer, features = c('CLEC4F','TIMD4','CD5L', 'VCAM1','HMOX1', 'CD36','CSF1R',
                              Macrophage_gene, 'LAMP1','ADGRE1','CX3CR1'
), group.by='seurat_clusters')

VlnPlot(Kupffer, features = c(#'CLEC4F','TIMD4','CD5L', 'VCAM1','HMOX1', 'CD36','CSF1R',
                              'CETP','TNF','CCL3','CCL4','KLF2','KLF6','IL1RN','NEDD9'
), group.by='seurat_clusters')


cluster0 = FindMarkers(Kupffer, ident.1 = 0, min.pct = 0.25)

Kupffer = subset(Kupffer, idents = c(0,3))
Kupffer$celltype ='Kupffer cells'

Neutrophil = subset(Myeloid_final, idents =  c(9)) #i also try to subcluster both 8 and 9 but turns out not able to find more but make it more blurry
Neutrophil = SeuratAnalysis1(Neutrophil)
ElbowPlot(Neutrophil, ndims = 50)
Neutrophil = SeuratAnalysis2(Neutrophil, harmony = F, n = 30)
DimPlot(Neutrophil,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'celltype') 

VlnPlot(Neutrophil, features = c('LCN2', 'ORM1','MMP8','PAGE5','SYCP1','SYN3', 'CTAG2','DEFA3'
), group.by='seurat_clusters')

VlnPlot(Neutrophil, features = c('S100A8','S100A9','S100A6', 'CXCR4','CXCR2','CTSG','ELANE',
                                 'CD44','MPO','ZEB2', 'CST3','CSF3R','FCGR3B','CDA' ,'LY6D'
), group.by='seurat_clusters')
#7 IS LIKELY BECAUSE IT HAS LY6G, HLA-G, S100A8/9/16 (LY6D SIGNALS NEUTROPHIL RECRUITMENT)
VlnPlot(Neutrophil, features = Neutrophil_gene)


cluster3 = FindMarkers(Neutrophil, ident.1 = 3, min.pct = 0.25)
#Neutrophil= subset(Neutrophil, idents = 1:7)
Neutrophil_final = subset(Neutrophil, idents = 5)
Neutrophil_final$celltype ='Neutrophil'

#clena up spp1+ mac ------ 
#Spp1 = subset(Myeloid_final, idents =  "SPP1+ macrophages") #i also try to subcluster both 8 and 9 but turns out not able to find more but make it more blurry
Spp1 = SeuratAnalysis1(Spp1)
ElbowPlot(Spp1, ndims = 50)
Spp1 = SeuratAnalysis2(Spp1, harmony = F, n = 20)
DimPlot(Spp1,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'seurat_clusters') 

VlnPlot(Spp1, features = c( 
  'CD3D','CD3E','CD40LG',
  'CD8A','CD8B','ITGB1','IFNG',
  'CCR7', 'TCF7','SELL',
  'GZMB','GZMH','PRDM1','KLRG1','KLF2',  'RORA', 'CX3CR1', #TEM
  'ITGAE','ITGA1', 'TNIP3', #TRM 
  'NCAM1','FCGR3A', 'KLRF1','PRF1'
))

VlnPlot(Spp1, c(Macrophage_gene,  'SPP1','VEGFA','VCAN','THBS1'))
#v2:
#1,2,3 have cd3; especially 3; 2,6,7 are not mac 

VlnPlot(Spp1, c('TNF','CCL3','CCL4','KLF2','KLF6','IL1RN','NEDD9',
                     'ISG15','IFIT1','IFI27','CCL2','CCL8',
                     'CXCL10','CXCL9',
                     'TREM2','CD9','APOE','FABP5'))

cluster3 = FindMarkers(Spp1, ident.1 = 3, min.pct = 0.25)
#1 are doublet with cytocd8; look more into 4, 7; 8 is hepatocyte 
#Spp1 = subset(Spp1, idents = c(0, 2:8))
#Spp1 = subset(Spp1, idents = c(0,1,4:7))

cluster.ids = c("SPP1+ macrophages", 'SPP1+ macrophages','Unknown', 'SPP1+ macrophages',
                'SPP1+ macrophages','Hepatocyte')
names(cluster.ids) <- levels(Spp1)
Spp1 <- RenameIdents(Spp1, cluster.ids)
DimPlot(Spp1, reduction = "umap", label = TRUE, repel = TRUE)  

Spp1$celltype = Idents(Spp1)
Spp1@meta.data$celltype %>% table

Spp1 = subset(Spp1, idents = c("SPP1+ macrophages", 'Hepatocyte'))
saveRDS(Spp1, file = 'Spp1 mac subcluster.rds')


#check interstitial mac-----
#Interstitial = subset(Myeloid_final, idents =  "Interstitial macrophages") #i also try to subcluster both 8 and 9 but turns out not able to find more but make it more blurry
Interstitial = SeuratAnalysis1(Interstitial)
ElbowPlot(Interstitial, ndims = 50)
Interstitial = SeuratAnalysis2(Interstitial, harmony = F, harmony_by='orig.ident',
                               n = 30)
DimPlot(Interstitial,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'seurat_clusters') 

VlnPlot(Interstitial, c(Macrophage_gene
                       # 'CD3D','CD3E','CD40LG',  'CD8A','CD8B' 
                        ))
FeaturePlot(Interstitial, Macrophage_gene, cols = c("gray",'red'), label =T)
# 0 is pure spp1; 2,4 = spp1+ lam; 1= lam; 3 = cdc2b 

cluster2 = FindMarkers(Interstitial, ident.1 = 2, min.pct = 0.25)

VlnPlot(Interstitial, c('TNF','CCL3','CCL4','KLF2','KLF6','IL1RN','NEDD9',
                'ISG15','IFIT1','IFI27','CCL2','CCL8',
                'CXCL10','CXCL9',
                'TREM2','CD9','APOE','FABP5',
                'CCL2','TNIP3',
                'SPP1' , 'VCAN','VEGFA','THBS1')) 
VlnPlot(Interstitial, features = c(DC_gene, 'CD1A','CD207'))
VlnPlot(Interstitial, features = c('S100A8','S100A9','S100A6', 'CXCR4','CXCR2','CTSG','ELANE',
                                 'CD44','MPO','ZEB2', 'CST3','CSF3R','FCGR3B','CDA' ,'LY6D'
), group.by='seurat_clusters')


table(Interstitial$celltype)

cluster.ids = c( 'LAM',  'SPP1+ LAM',
                'IL22RA2+ Langerhans')
names(cluster.ids) <- levels(Interstitial)
Interstitial <- RenameIdents(Interstitial, cluster.ids)
DimPlot(Interstitial, reduction = "umap", label = TRUE, repel = TRUE)  


Interstitial$celltype = Idents(Interstitial)
Interstitial@meta.data$celltype %>% table 
saveRDS(Interstitial, 'Interstitial_macrophage_subcluster.rds')


#Interstitial = subset(Interstitial, idents = c(0,1,3,4,6:10))
Interstitial = subset(Interstitial, idents = c(3:8))
Interstitial = subset(Interstitial, idents = c(0,2:5))
Interstitial = subset(Interstitial, idents = 1:4)

#clean up lam (may have T cell)-----  
#LAM = subset(Myeloid_final, idents = 'LAM')
LAM = SeuratAnalysis1(LAM)
ElbowPlot(LAM, ndims = 50)
LAM = SeuratAnalysis2(LAM, harmony =F, n = 20)
DimPlot(LAM,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'seurat_clusters') 

cluster6 = FindMarkers(LAM, ident.1 = 6, min.pct = 0.25)
#5,4,7 are separated

VlnPlot(LAM, features =c(Macrophage_gene,
                         #'CD3D','CD3E','CD40LG',  'CD8A','CD8B',
                         'TREM2','CD9','APOE','FABP5'))
VlnPlot(LAM, features = c('S100A8','S100A9','S100A6', 'CXCR4','CXCR2','CTSG','ELANE',
                          'CD44','MPO','ZEB2', 'CST3','CSF3R','FCGR3B','CDA' ,'LY6D'
) )

DotPlot(LAM, c(#'TNF','CCL3','CCL4','KLF2','KLF6','IL1RN','NEDD9',
                       # 'ISG15','IFIT1','IFI27','CCL2','CCL8',
                       # 'CXCL10','CXCL9',
                        'TREM2','CD9','APOE','FABP5',
                        DC_gene, 'ITGAX','ITGAM','ZBTB46'
                       # 'CCL2','TNIP3',
                       # 'SPP1' , 'VCAN','VEGFA','THBS1'
                       ), cols = 'RdYlBu') 
#5 is ifn-stimulated; 6 IS CDC2
cluster.ids = c("LAM", "LAM","LAM","LAM","LAM",'IFN-primed macrophages',
                "LAM","LAM")
names(cluster.ids) <- levels(LAM)
LAM <- RenameIdents(LAM, cluster.ids)
DimPlot(LAM, reduction = "umap", label = TRUE, repel = TRUE)  


LAM$celltype = Idents(LAM)
LAM@meta.data$celltype %>% table 
saveRDS(LAM, 'LAM_macrophage_subcluster.rds')


#first get rid of the T cells 
#LAM = subset(LAM, idents = c(1:3, 5:10))
#LAM= subset(LAM, idents = c(0:3, 6,8,9))


#structural cells----------
structural = subset(total, idents = c(0, 8,11, 15,16, 17, 19, 23, 24,25,26,27,30:34, 37:40))
structural = SeuratAnalysis1(structural)
ElbowPlot(structural, ndims = 50)
structural = SeuratAnalysis2(structural, harmony = F, n = 30)
DimPlot(structural,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'seurat_clusters')  

table(structural$seurat_clusters)
saveRDS(structural, file='~/jliang/Human tissue harmonization/structural/structural_raw.rds')
#first get rid of the cycling cells 
VlnPlot(structural, features = c('MKI67','TOP2A'), raster=F) #13, 24, 25 

#subcluster by cancer / tissue type
#hcc structural ----------
HCC = subset(structural, idents=c(42,39,50,24, 37, 41,44,40,19,34,51))
HCC = SeuratAnalysis1(HCC)
ElbowPlot(HCC, ndims = 50)
HCC = SeuratAnalysis2(HCC, harmony = T, harmony_by = 'Patient', n = 20)
DimPlot(HCC,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'seurat_clusters')  

table(HCC$seurat_clusters)
#there are a few hnscc 
#v1: 6 = cycling? 2 and 14 = CHOLANGIOCYTE 
#0 and 10 and 13 hepatocytes; 4 = SINUSOIDAL ; 3,5,9,12,15 = macrovascular
#10 is the hepatocyte version that is more "malignant" (SAA1... and MT genes)

#v2
#0, 12 cholangiocyte; 4,6,8,11 = macrovascular ; 5= sinusoidal ; 1,2,3,7,10,13 = hepatocyte 
cluster.ids =c('Cholangiocyte', 'Hepatocyte','Hepatocyte','Hepatocyte','Macrovascular endothelial',
               'Sinusoidal endothelial', 'Macrovascular endothelial', 'Hepatocyte', 
               'Macrovascular endothelial', 'Hepatocyte', 'Hepatocyte', 'Fibro',
               'Cholangiocyte', 'Hepatocyte')
names(cluster.ids) <- levels(HCC)
HCC <- RenameIdents(HCC, cluster.ids)
DimPlot(HCC, reduction = "umap", label = TRUE, repel = TRUE)  

HCC$celltype = Idents(HCC)
HCC@meta.data$celltype %>% table
saveRDS(HCC, 'HCC_structural.rds')
#the single cluster 12 is still cholangiocyte because the top 1 gene CHST4 is exclusive to cholangiocyte; no weird and cycling gene 


VlnPlot(HCC, features = c('MKI67','TOP2A','EPCAM'), raster=F)  
VlnPlot(HCC, features = c('ALDOB','APOA1','ALB', 'ALDOA','AKR1B10'))
VlnPlot(HCC, features = c('EPCAM','CD24', 'KRT19',
                                'PLPP2', 'CLDN4', 'FXYD2', 'FAM3B', 
                                'SFRP5', 'SLC5A1','KRT19','B3GNT3','PRR15L','CLDN10')) #CHOLANGIOCYTE
VlnPlot(HCC, features = c('CLEC4G', 'FCN2', 'CLEC4M', 'IFITM10', 
                                  'RELN', 'MRO', 'CLEC1B', 'OIT3', 'STAB2', 'FCN3')) #sinusoidal 
VlnPlot(HCC, features = c('PLVAP', 'AIF1L', 'MGP', 'SLCO2A1', 'VWF', 
                          'SPARCL1', 'CD34', 'EDNRB', 'EMCN', 'CLEC14A',
                          fibro_gene)) #macrovascular
VlnPlot(HCC, features = c( 'HPR', 'F12',  'HGFAC', 'AGXT', 'SAA4','ARG1' ) ) 


DotPlot(HCC, features = c( 'EPCAM', 'KRT19', 'PLPP2', 'CLDN4', 'FXYD2', 'FAM3B', 'SLC5A1', 'CLDN10',
                           'CLEC4G', 'FCN2', 'CLEC4M', 'IFITM10', 'STAB2', 'FCN3',
                           'PLVAP', 'AIF1L', 'MGP', 'SLCO2A1', 'VWF', 'EMCN', 'CLEC14A',
                           'HPR', 'F12',  'ALB', 'ALDOB', 'AGXT', 'SAA4','ARG1',
                           fibro_gene
), group.by = 'celltype', cols = 'RdYlBu') + coord_flip() + RotatedAxis()


cluster12 = FindMarkers(HCC, ident.1 = 12, min.pct =0.25)
#HCC = subset(HCC, idents ='HCC')
#HCC= subset(HCC, idents = c(0:5, 7:14))



#hnscc structural---------- 
Idents(structural) ='orig.ident'
#hnscc= subset(structural, idents='hnscc')
hnscc = SeuratAnalysis1(hnscc)
ElbowPlot(hnscc, ndims = 50)
hnscc = SeuratAnalysis2(hnscc, harmony = T, harmony_by = 'Patient', n = 25)
DimPlot(hnscc,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'celltype' )  

table(hnscc$seurat_clusters)
nrow(hnscc@meta.data)

VlnPlot(hnscc, features = c('CD3D','CD3E','CD40LG',
                          'CD8A','CD8B','CD14', 'CD163', 'CD68', 'MARCO', 'MRC1',
                          'C1QA', 'C1QB',
                          Cycling_genes))
VlnPlot(hnscc, c('MYL9','CDH5', #SMC 
                 'MYH1', 'MYOG', 'TNNC1', 'TNNI1', #SKELETAL 
                 'PAX7','MYF5')) #PROGENITOR 

#13 may be cycling; 12 may be myeloid; 16 and 3 and 5 may be T cells (16 seems to be a mix of cytocd8 and treg?)
#3 is a mix of myeloid and epithelial (krt); 5 is also a mix

#v2: 3,4,15 can be immune mix 
#3,10, 17, 18 = epithelial
#7,8,13= endo
#13 lymphatic; 8 with sema3g (artery); 7 with ackr1 (venous)
#fibro: 0-2, 4-6, 9, 15-16 fibro (the upper blob and 4 and 15)
#14= mast; 3 is epi 
#11 = MUSCLE STEM; 9 = SMC; 12 = SKELETAL MUSCLE 
 
cluster.ids = c('Fibro','Fibro','Fibro','Epithelial', 'Fibro','Fibro','Fibro',
                'Vein endothelial','Artery endothelial' ,'Smooth muscle cell',  'Epithelial',
                'Muscle stem cell','Skeletal muscle cell','Lymphatic endothelial',
                'Mast cell', 'Immune', 'Fibro','Epithelial',  "Epithelial"
)
names(cluster.ids) <- levels(hnscc)
hnscc <- RenameIdents(hnscc, cluster.ids)
DimPlot(hnscc, reduction = "umap", label = TRUE, repel = TRUE)  

hnscc$celltype = Idents(hnscc)
hnscc@meta.data$celltype %>% table
hnscc= subset(hnscc, idents = c('Fibro',  'Vein endothelial','Artery endothelial' ,'Smooth muscle cell',  'Epithelial',
                                'Muscle stem cell','Mast cell', 'Skeletal muscle cell','Lymphatic endothelial'))
saveRDS(hnscc, 'hnscc_structural.rds')

VlnPlot(hnscc,features =c('EPCAM', 'CDH1', 'KRT18','KRT19','KRT14', 'KRT15', 'KRT5', 'KRT17',
                          endo_gene) )
FeaturePlot(hnscc, features = c(fibro_gene, 'LUM'), label =T, 
            reduction ='umap', cols=c('gray','red'))
DotPlot(hnscc, features = c('EPCAM', 'CDH1', 'KRT18','KRT19','KRT14', 'KRT15', 'KRT5', 'KRT17',
                            endo_gene, 
                            fibro_gene, 'LUM',  
                            SMC_gene,'MYH1', 'MYOG', 'TNNC1', 'TNNI1', #SKELETAL 
                            'PAX7','MYF5',
                            'TPSAB1','CPA3','KIT'
                            ), group.by = 'celltype', cols = 'RdYlBu') + coord_flip() + RotatedAxis()

cluster18 = FindMarkers(hnscc, ident.1 = 18, min.pct =0.25)

#hnscc= subset(structural, idents='hnscc')
#hnscc= subset(hnscc, idents = c(0:13, 15:19))
#hnscc= subset(hnscc, idents = c(0:5, 7:17))
#hnscc= subset(hnscc, idents = c(0:17, 19))

#kidney structural---- 
kidney = subset(structural, idents= 'Kidney')
kidney_final = SeuratAnalysis1(kidney_final)
ElbowPlot(kidney_final, ndims = 50)
kidney_final = SeuratAnalysis2(kidney_final, harmony = T, harmony_by = 'Patient', n = 20)
DimPlot(kidney_final,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'celltype')  

saveRDS(kidney_final, file = 'Kidney_structural.rds')
table(kidney$seurat_clusters)

cluster14 = FindMarkers(kidney, ident.1 = 14, min.pct =0.25)
#v1: 18 = mast; 15= cd8+ T cell; 16 with FOXI1 may be renal epi 
#1 AND 11 = SMC, 17 is fibro
#5 is vein, 6=  artery; 13= RENAL GLOMERULAR CELL; 2 also endo (may be lymphatic endo); 19 = cap endo (with rgcc) 
#7,10,16 = epi
#check 0,3,8,9, 14 

#16 = COLLECTING DUCT INTERCALATED; 10 = HENLE THICK ASCENDING  ; 
#7 has saa1, saa2, fga, fgb, fgg -> fibrosis?inflammatory?
#4 and 12 = PROXIMAL TUBULE  

cluster.ids = c('Kidney cancer cell', 'Smooth muscle cell', 'Unknown', 'Kidney cancer cell', 
                'Proximal tubule cell', 'Vein endothelial', 'Artery endothelial', 'Kidney cancer cell', 
                'Kidney cancer cell', 'Kidney cancer cell', 'Loop of henle thick ascending cell',
                'Smooth muscle cell', 'Proximal tubule cell', 'Renal glomerular cell', 'Kidney cancer cell',
                'CD8+ T cell', 'Collecting duct intercalated cell', 'Fibro', 'Mast cell', 'Capillary endothelial')

names(cluster.ids) <- levels(kidney)
kidney <- RenameIdents(kidney, cluster.ids)
DimPlot(kidney, reduction = "umap", label = TRUE, repel = TRUE)  

kidney$celltype = Idents(kidney)
kidney@meta.data$celltype %>% table


kidney= subset(kidney, idents = unique(kidney$celltype)[-3])
kidney_final = merge(kidney, kidney_sub_endo)

DotPlot(kidney, features = c('EPCAM', 'CDH1', 'KRT18','KRT19','KRT14', 'KRT15', 'KRT5', 'KRT17',
                            endo_gene, 
                            fibro_gene, 'LUM',  
                            SMC_gene,'MYH1', 'MYOG', 'TNNC1', 'TNNI1', #SKELETAL 
                            'PAX7','MYF5',
                            'TPSAB1','CPA3','KIT'
), group.by = 'seurat_clusters', cols = 'RdYlBu') + coord_flip() + RotatedAxis()
DotPlot(kidney, c('TAGLN','ACTA2', 'MYL9','CDH5',
                  'RGS5','PDGFRB', 'PLN',
                  'ANPEP','CSPG4','MCAM',
                  fibro_gene), group.by = 'seurat_clusters', cols = 'RdYlBu') + coord_flip() + RotatedAxis() 
 
VlnPlot(kidney, c("EMCN", 'SOST', 'PLAT', endo_gene))
VlnPlot(kidney, features = c("CLCNKB",'EPCAM', 'SLC8A1', 'STC1','CALB1','ADM', #COLLECTING TUBE EPI 
                             'HS3ST5','SLC14A2', 'KSR2', 'CNGA1', 'ELF5','KCNK3','CAPN8', #COLLECTING DUCT PRINCIPAL 
                             'NXPH2', 'CLNK','DMRT2', 'SLC4A9' #COLLECTING DUCT INTERCALATED 
                             )) 
VlnPlot(kidney, features = c('PLG', 'AGXT2','ACSM2B','ALDOB','ACSM2A', 'FUT6', 'SLC22A6', 'SLC13A3', 'SLC36A2', 'SLC22A8'  #proximal tubule 
                             #'KCNK2' ,'MOXD1', 'CLMP', #INTERSTITIAL FIBRO 
                             #'SLC12A3','TRIM50', 'TRPM6', 'SALL3' #DISTAL CONVOLUTED TUBULE EPI 
                             )) 
VlnPlot(kidney, features=c('CA9', 'NNMT', 'NDUFA4L2' )) #CANCER CELL MARKERS FROM ORIGINAL PAPERS -- 0,3,7,8,9,14 
VlnPlot(kidney, features= c('UMOD', 'GP2', 'SLC12A1', 'ACP3', #HENLE THICK 
                           # 'SH3GL3', 'ABCA4', 'EPHA6', 'COLEC10', 'CCNP', #LOOP OF  HENLE THIN ASCNEDING
                            #'SLC12A3', 'TRIM50', 'TRPM6','SALL3', #DISTAL CONVOLUTED TUBE 
                           # 'SLCO1A2', 'ABCA13', 'KLRG2', 'SLC44A5', #HENLE THINN DESCENDING
                            'PTPRQ', 'PTPRO', 'NPHS2', 'ADAMTS19', 'NPHS1'
                            ))

#final dotplot:
DotPlot(kidney_final, features =unique( c('PLG', 'AGXT2','ACSM2B','ALDOB','ACSM2A', 'FUT6', 'SLC22A6',
                             'NXPH2', 'CLNK','DMRT2', 'SLC4A9' , #collecting duct 
                             'UMOD', 'GP2', 'SLC12A1', 'ACP3',#loop of henle 
                             "EMCN", 'SOST', 'PLAT', endo_gene, #renal glomerular 
                             'TAGLN','ACTA2', 'MYL9', 'RGS5','PDGFRB', 'PLN',
                             fibro_gene,'LUM', 
                             'CA9', 'NNMT', 'NDUFA4L2' ,'EPCAM', 'KRT19',
                             'CD3D','CD3E','CD8A','GZMB', 'IFNG', 
                             'TPSAB1','CPA3','KIT'
                             )),
        group.by = 'celltype', cols = 'RdYlBu') + coord_flip()+ RotatedAxis()

kidney_final$celltype = factor(kidney_final$celltype, levels = c(unique(kidney_final$celltype)[c(11, 12, 7, 9, 4,6,10,13, 3, 8, 1, 2, 5)] ))
#subcluster 2 

#kidney_sub_endo = subset(kidney, idents = 2)
kidney_sub_endo = SeuratAnalysis1(kidney_sub_endo)
ElbowPlot(kidney_sub_endo, ndims = 50)
kidney_sub_endo = SeuratAnalysis2(kidney_sub_endo, harmony = T, harmony_by = 'Patient', n = 20)
DimPlot(kidney_sub_endo,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'seurat_clusters')  

VlnPlot(kidney_sub_endo, c(endo_gene, 'TAGLN','ACTA2', 'MYL9','CDH5',
                           'RGS5','PDGFRB', 'EPCAM'))
VlnPlot(kidney_sub_endo, c("EMCN", 'SOST', 'PLAT'))
table(kidney_sub_endo$seurat_clusters)

cluster4 = FindMarkers(kidney_sub_endo, ident.1 = 4, min.pct =0.25)

#kidney_sub_endo= subset(kidney_sub_endo, idents= 0:7)
cluster.ids = rep("Endothelial cell", 8)
names(cluster.ids) <- levels(kidney_sub_endo)
kidney_sub_endo <- RenameIdents(kidney_sub_endo, cluster.ids)
DimPlot(kidney_sub_endo, reduction = "umap", label = TRUE, repel = TRUE)  

kidney_sub_endo$celltype = Idents(kidney_sub_endo)
kidney_sub_endo@meta.data$celltype %>% table


#breast structural---- 
breast = subset(structural, idents= 'Breast')
breast = SeuratAnalysis1(breast)
ElbowPlot(breast, ndims = 50)
breast = SeuratAnalysis2(breast, harmony = T, harmony_by = 'Patient', n = 20)
DimPlot(breast,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'seurat_clusters')  

table(breast$seurat_clusters)
cluster23 = FindMarkers(breast, ident.1 = 23, min.pct =0.25)
DotPlot(breast, features = c('CD3D','CD3E','CD40LG',
                            'CD8A','CD8B','CD14', 'CD163', 'CD68', 'MARCO', 'MRC1',
                            'C1QA', 'C1QB',
                            Cycling_genes), cols='RdYlBu') + coord_flip()
#immune: 5,20 

breast[['percent.rp']] = PercentageFeatureSet(breast, pattern = "^RP[SL]")
VlnPlot(breast, 'percent.rp', raster=F, group.by = 'orig.ident')
breast= subset(breast, percent.rp <= 20)

breast_immune =  subset(breast, idents =c(5, 20)) #immune
breast_immune = SeuratAnalysis1(breast_immune)
ElbowPlot(breast_immune, ndims = 50)
breast_immune = SeuratAnalysis2(breast_immune, harmony = F, n = 10)
DimPlot(breast_immune,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'seurat_clusters')  

VlnPlot(breast_immune, features = c('CD3D','CD3E','CD40LG',
                             'CD8A','CD8B','CD14', 'CD163', 'CD68', 'MARCO', 'MRC1',
                             'C1QA', 'C1QB'), raster=F)  
table(breast_immune$seurat_clusters)

cluster4 = FindMarkers(breast_immune, ident.1 = 4, min.pct =0.25)
breast_immune= subset(breast_immune, idents =3)
Idents(breast_immune)  ='celltype'
breast_immune$celltype= 'Fibro'
Idents(breast_immune) = breast_immune$celltype

#3= fibro matrix 


#breast_lum = subset(breast, idents = c(0,11,12,16,17))
breast_lum = SeuratAnalysis1(breast_lum)
ElbowPlot(breast_lum, ndims = 50)
breast_lum = SeuratAnalysis2(breast_lum, harmony = F, n = 10)
DimPlot(breast_lum,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'Patient')  

DotPlot(breast_lum, features = c('KRT14','KRT17', 'DST','KRT5','SAA1','SFN',  'ACTA2','TP63','EPCAM','TAGLN','TPM2', #BASAL
                             'AREG','MUCL1','DBI','AZGP1','PIP','KRT18','AGR2',  'ESR1','AR', 'PGR', #LUM-hr
                             'SCGB2A2','SLPI','WFDC2','LTF','KRT15','MMP7','SCGB3A1','MYLK' #LUM-SEC
), cols='RdYlBu') + coord_flip()
VlnPlot(breast_lum, features = c('KRT14','KRT17','EPCAM', 
                                 'AREG','MUCL1','DBI','AZGP1','PIP'))
DotPlot(breast_lum, features = c("CFO",'FOS','GPX3','PLA2G2A,','WSP2','MGST1','MFAP5',#prematrix
                                    'SFRP4','CLU','GOS2','MGP','OGN','IGFBP5',#sfrp4
                                    'MMP3','CXCL1','GEM','CEBPB','TNFAIP6',  #fibro major 
                                    'COL3A1','POSTN','COL1A1','IGF1','IGFBP2','ADAM12','TNC' #fibro matrix 
), cols='RdYlBu') + coord_flip()

cluster9 = FindMarkers(breast_lum, ident.1 = 9, min.pct =0.25)

#8 = basal with tpm2 and tagln; 6,10 = lum-sec ;11 =fibro 
cluster.ids = c("Luminal epithelial-hr","Luminal epithelial-hr","Luminal epithelial-hr","Luminal epithelial-hr",
                "Luminal epithelial-hr","Luminal epithelial-hr",'Luminal epithelial-secretory',
                "Luminal epithelial-hr", 'Basal', "Luminal epithelial-hr", 'Luminal epithelial-secretory',
                "Fibro")
names(cluster.ids) <- levels(breast_lum)
breast_lum <- RenameIdents(breast_lum, cluster.ids)
DimPlot(breast_lum, reduction = "umap", label = TRUE, repel = TRUE)  

breast_lum$celltype = Idents(breast_lum)
breast_lum@meta.data$celltype %>% table


#breast_lum= subset(breast_lum, idents = c(0, 3:10))
#breast_lum= subset(breast_lum, idents = c(0:5,7:13))

breast_fibro = subset(breast, idents = c(1,6,13,18,19)) #also has smc 
breast_fibro = SeuratAnalysis1(breast_fibro)
ElbowPlot(breast_fibro, ndims = 50)
breast_fibro = SeuratAnalysis2(breast_fibro, harmony = F, n = 20)
DimPlot(breast_fibro,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'seurat_clusters')  

DotPlot(breast_fibro, features = c("CFO",'FOS','GPX3','PLA2G2A,','WSP2','MGST1','MFAP5',#prematrix
                                 'SFRP4','CLU','GOS2','MGP','OGN','IGFBP5',#sfrp4
                                 'MMP3','CXCL1','GEM','CEBPB','TNFAIP6',  #fibro major 
                                 'COL3A1','POSTN','COL1A1','IGF1','IGFBP2','ADAM12','TNC' #fibro matrix 
), cols='RdYlBu') + coord_flip()
VlnPlot(breast_fibro, features = SMC_gene)
VlnPlot(breast_fibro, features = c('GYPA', 'HBB', 'HBA1','FTL', 'ALAS2',
                                   'HBM','HBQ1'))
VlnPlot(breast_fibro, features = c('COL1A1', 'COL1A2', 'DCN',
                                   'ITGBL1','ASPN','LMOD1','NALCN','ITGBL1','WISP2','F2R','WIF1','CCDC68'))

cluster10 = FindMarkers(breast_fibro, ident.1 = 10, min.pct =0.25)
#9 = sema3g (artery); 12 =ackr1 (venous); 3? 
#11 = red blood cell ? 6 = smc with rgs5; 7 with mmp3 exclusively 

cluster.ids = c("Fibro", "Fibro","Fibro","Fibro","Fibro","Fibro",
                'SMC',"Fibro","Myofibroblast",'Artery endothelial', "Fibro",'Erythrocyte',
                'Vein endothelial')
names(cluster.ids) <- levels(breast_fibro)
breast_fibro <- RenameIdents(breast_fibro, cluster.ids)
DimPlot(breast_fibro, reduction = "umap", label = TRUE, repel = TRUE)  

breast_fibro$celltype = Idents(breast_fibro)
breast_fibro@meta.data$celltype %>% table


#breast_endo = subset(breast, idents = c(9,10,14))
breast_endo = SeuratAnalysis1(breast_endo)
ElbowPlot(breast_endo, ndims = 50)
breast_endo = SeuratAnalysis2(breast_endo, harmony = F, n = 20)
DimPlot(breast_endo,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'seurat_clusters')  

cluster4 = FindMarkers(breast_endo, ident.1= 4, min.pct= 0.25) 
#7= smc 
VlnPlot(breast_endo, features = c(SMC_gene))
VlnPlot(breast_endo, endo_gene)
#0,2,3 = lymphatic; 1,5,6 = artery (sema3g) ; 4 with rgcc (capillary )

cluster.ids = c("Lymphatic endothelial", 'Artery endothelial',"Lymphatic endothelial", "Lymphatic endothelial", 
               "Capillary endothelial",'Artery endothelial', 'Artery endothelial', 'SMC' )
names(cluster.ids) <- levels(breast_endo)
breast_endo <- RenameIdents(breast_endo, cluster.ids)
DimPlot(breast_endo, reduction = "umap", label = TRUE, repel = TRUE)  

breast_endo$celltype = Idents(breast_endo)
breast_endo@meta.data$celltype %>% table



#breast_endo = subset(breast_endo, idents = c(1:10))
#breast_endo = subset(breast_endo, idents = c(0,1,2,4:7, 9))


breast_unknown = subset(breast, idents = c(2,3,4,7,8,15))
breast_unknown = SeuratAnalysis1(breast_unknown)
ElbowPlot(breast_unknown, ndims = 50)
breast_unknown = SeuratAnalysis2(breast_unknown, harmony = F, n = 10)
DimPlot(breast_unknown,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'seurat_clusters')  

cluster4 = FindMarkers(breast_unknown, ident.1= 4, min.pct= 0.25) 
 
#0= lum-hr; 3 = lymphatic endo; 1,2,6 = fibro  ; 5,7,9 = erythrocyte ; 4,8=adipocyte
cluster.ids = c("Luminal epithelial-hr", 'Fibro', 'Fibro', 'Lymphatic endothelial',
                'Adipocyte', 'Erythrocyte', 'Fibro', 'Erythrocyte', 'Adipocyte',
                'Erythrocyte')
names(cluster.ids) <- levels(breast_unknown)
breast_unknown <- RenameIdents(breast_unknown, cluster.ids)
DimPlot(breast_unknown, reduction = "umap", label = TRUE, repel = TRUE)  

breast_unknown$celltype = Idents(breast_unknown)
breast_unknown@meta.data$celltype %>% table


VlnPlot(breast_unknown, c(SMC_gene, endo_gene, fibro_gene))
VlnPlot(breast_unknown, 'percent.rp')
VlnPlot(breast_unknown, features = c('GYPA', 'HBB', 'HBA1','FTL', 'ALAS2',
                                   'HBM','HBQ1'))

VlnPlot(breast_unknown, features = c("PLIN4",'ADIPOQ','FABP4','LEP','LPL')) #ADIPOCYTES -- NO CLUSTERS EXPRESS ALL THESE GENES 

#breast_unknown= subset(breast_unknown, idents = c(1,10,5,7,9,11))
#breast_unknown= subset(breast_unknown, idents = c(0,1,3,5,6,8,7))


breast_all = merge(x = breast_unknown, y = c(breast_endo, breast_lum,
                                             breast_fibro, breast_immune))
breast_all = SeuratAnalysis1(breast_all)
ElbowPlot(breast_all, ndims = 50)
breast_all = SeuratAnalysis2(breast_all, harmony = T, harmony_by = 'Patient', n = 20)
DimPlot(breast_all,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'celltype')  

Idents(breast_all) ='celltype'
breast_all=RenameIdents(breast_all, 'Myofibroblast' = 'Fibro')
breast_all$celltype = Idents(breast_all)


table(breast_all$celltype)
cluster17 = FindMarkers(breast_all, ident.1= 17, min.pct= 0.25)
#breast_all= subset(breast_all, idents= c(0:16,18,19)) #17 is a doublet of t cell and macrophage
saveRDS(breast_all, file ='Breast_structural.rds')
 
DotPlot(breast_all, features = c( "PLIN4",'ADIPOQ','FABP4','LEP','LPL', #adipocyte
                                  endo_gene, 
                                   'KRT14','KRT17', 'DST','KRT5','SAA1','SFN',  'ACTA2','TP63','EPCAM', #BASAL
                                  'HBB', 'HBA1', 'ALAS2',
                                  fibro_gene, 'LUM',
                                  'ITGBL1','ASPN','LMOD1','NALCN','WISP2', #myofibro
                                  'AREG','MUCL1','DBI','AZGP1','PIP','KRT18','AGR2',  'ESR1','AR', 'PGR', #LUM-hr
                             'SCGB2A2','SLPI','WFDC2','LTF','KRT15','MMP7','SCGB3A1','MYLK', #LUM-SEC 
                             SMC_gene
                            ) %>% unique(), cols='RdYlBu', group.by = 'celltype') + coord_flip() + RotatedAxis()
#0, 11, 12, 16 = lum-hr; 17= lum-sec + basal;  18 has muscle gene 
unique(breast_all$celltype)
breast_all$celltype = factor(breast_all$celltype, levels = c(unique(breast_all$celltype)[c(2, 7,6,4,11,9, 1, 
                                                                                         3,5,10,8)]))
DotPlot(breast_all, features = c('PI16','CD34','MFAP5', #UNIVERSAL
                            'CCL19','CXCL12','CH25H','APOE','EFEMP1','APOC1','C7','PPARG','MYOC', #SUPERFICIAL, PERIVASCULAR, FIBROBLASTIC RETICULAR CELL, CLOSE TO IMMUNE CELLS
                            'ASPN','COL11A1', 'MKX', 'TNMD', #FOLLICLE ASSOCIATED 
                            'SCN7A','FMO2','FGFBP2','OLFML2A','NGFR','RAMP1','CGRP', #SCHWANN, NEURAL 
                            'WNT2','COL10A1','LAMP5','NIRG1','OGN','ACTA2', 'TAGLN', 'RUNX2','SULF1','ADAM12','COL18A1','CDH2','LRRC17','PIEZO2','ASPN','COMP','COL11A1','SFRPA','CTHRC1'
) %>% unique(), group.by = 'celltype', cols = 'RdYlBu')+ RotatedAxis() + coord_flip()
                                                                             


DotPlot(breast_all, features = c("CFO",'FOS','GPX3','PLA2G2A,','WSP2','MGST1','MFAP5',#prematrix
                                 'SFRP4','CLU','GOS2','MGP','OGN','IGFBP5',#sfrp4
                                 'MMP3','CXCL1','GEM','CEBPB','TNFAIP6',  #fibro major 
                                 'COL3A1','POSTN','COL1A1','IGF1','IGFBP2','ADAM12','TNC' #fibro matrix 
), cols='RdYlBu') + coord_flip()
#1, 6, 13 = fibro; 18,19 POSTN+ 

DotPlot(breast_all, features = c(  'ITGBL1','ASPN','LMOD1','NALCN','ITGBL1','WISP2','F2R','WIF1','CCDC68' #myofibro
) %>% unique(),cols = 'RdYlBu') + coord_flip() 

#9, 10, 14, 19 = endothelial ; 13, 17,18, 1, 6 = smc

VlnPlot(breast_all, features = c("CNN1",'MYH11', 'SYNM','ACTG2')) #VASCULAR SMC 


#lung------
#Lung = subset(structural, idents=c("LUAD", 'LUSC'))
Lung = SeuratAnalysis1(Lung)
ElbowPlot(Lung, ndims = 50)
Lung = SeuratAnalysis2(Lung, harmony = T, harmony_by = 'Patient', n = 20)
DimPlot(Lung,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'seurat_clusters')  


Lung[['percent.rp']] = PercentageFeatureSet(Lung, pattern = "^RP[SL]")
VlnPlot(Lung, 'percent.rp', raster=F, group.by = 'orig.ident')
#Lung= subset(Lung, percent.rp <= 20)

#28 may be schwann (neural ish); 24 is cycling 
table(Lung$seurat_clusters)
#Lung = subset(Lung, idents=c(0:23, 25:29))

#7,11 = at2; 15= at1; 9 and 12 is a mix of ciliated and at 
#immune: 0,1,21, 3,4
#endo: 6, 10 
#fibro: 16,17 ; 16 may be myofibro 
#12 may be goblet; epi: 2, 23, 21,22, 18, 14, 26
#5,13,24 = mast cells 

cluster25 = FindMarkers(Lung, ident.1 = 25, min.pct=0.25)

DotPlot(Lung, features =c('CD3D','CD3E','CD8A','GZMB', 'IFNG', Macrophage_gene), cols = 'RdYlBu')+ coord_flip()+ RotatedAxis()
VlnPlot(Lung, features = c('ASPN', 'CTHRC1', 'COL1A1',
                           'ITGBL1','COL3A1', 'CCN5' , fibro_gene)) #myofibroblast
DotPlot(Lung, features = c('MUC5AC', 'SERPINB11', 'CEACAM5', 'FER1L6',  'CYP2F1','MUC16', 'MSMB', #goblet 
                           'TP63','KRT14', 'KRT15', 'KRT5', 'KRT17','EPCAM'), cols = 'RdYlBu') + coord_flip()+ RotatedAxis() #epi   
VlnPlot(Lung, Cycling_genes)  
VlnPlot(Lung, endo_gene) 
VlnPlot(Lung, Mast_genes)
DotPlot(Lung, c(AT2_gene, AT1_gene, 'SCGB1A1', 'SCGB3A1'), cols = 'RdYlBu') + coord_flip()+ RotatedAxis()


#Lung_epi = subset(Lung, idents =c( 2, 23, 21,22, 18, 14, 26))
Lung_epi = SeuratAnalysis1(Lung_epi)
ElbowPlot(Lung_epi, ndims = 50)
Lung_epi = SeuratAnalysis2(Lung_epi, harmony = T, harmony_by = 'orig.ident', n = 20)
DimPlot(Lung_epi,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'seurat_clusters')  

cluster3 = FindMarkers(Lung_epi, ident.1= 3, min.pct= 0.25)
VlnPlot(Lung_epi, c('MUC5AC', 'SERPINB11', 'CEACAM5', 'FER1L6',  'CYP2F1','MUC16', 'MSMB', #goblet 
                    'TP63','KRT14', 'KRT15', 'KRT5', 'KRT17','EPCAM'
                   # endo_gene
                    )) 
VlnPlot(Lung_epi, c("TFF1", 'TFF2','CLDN18','S100P',
                    'KRT6A', 'KRT6B'))
VlnPlot(Lung_epi, c(Macrophage_gene, "SPP1",
                    'CD3D','CD3E', 'CD40LG','CD8A'))

#after getting rid of the immune cells:
#epi: except for 0,5 have epcam; 6,9 = GOBLET; 0,1,5,7 HAVE KRT6  
#3 HAS TFF1, tff2, cldn18  
 
cluster.ids = c('Squamous epithelial', 'Squamous epithelial','Basal',
                'TFF1+S100P+ tumor cell','Epithelial','Squamous epithelial', 'Goblet cell', 
                "Squamous epithelial", 'Epithelial','Goblet cell','Basal')
names(cluster.ids) <- levels(Lung_epi)
Lung_epi <- RenameIdents(Lung_epi, cluster.ids)
DimPlot(Lung_epi, reduction = "umap", label = TRUE, repel = TRUE)  

Lung_epi$celltype = Idents(Lung_epi)
Lung_epi@meta.data$celltype %>% table


table(Lung_epi$seurat_clusters)
 
#Lung_epi = subset(Lung_epi, idents= c(3,5:16))
#Lung_epi = subset(Lung_epi, idents= c(0:11))
#Lung_epi = subset(Lung_epi, idents= c(0, 2:12))
#Lung_epi = subset(Lung_epi, idents= c(1,2,3,5:12))
#Lung_epi = subset(Lung_epi, idents= c(0:3, 5:9))
#Lung_epi = subset(Lung_epi, idents= c(0, 2:10))
#Lung_epi = subset(Lung_epi, idents= c(0, 1,3:9))
#Lung_epi = subset(Lung_epi, idents= c(0, 2,3,4,5,7,8,9,10))


#Lung_epi_sub= subset(Lung_epi, idents=c(0,1,2,4))
Lung_epi_sub = SeuratAnalysis1(Lung_epi_sub)
ElbowPlot(Lung_epi_sub, ndims = 50)
Lung_epi_sub = SeuratAnalysis2(Lung_epi_sub, harmony = T, 
                               harmony_by = 'orig.ident', n = 20)
DimPlot(Lung_epi_sub,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'seurat_clusters')  


VlnPlot(Lung_epi_sub, c(
                    'CD3D','CD3E', 'CD40LG','CD8A',
                   'CXCL10','CXCL9','CXCL2','CXCL8','G0S2', 'NFKBIA', 'CCL2') )

VlnPlot(Lung_epi_sub, features = c( 'TP63','KRT14',
                                    'KRT15', 'KRT5', 'KRT17','EPCAM' ,
                                    'KRT6A', 'KRT6B'))

VlnPlot(Lung_epi_sub, c(AT1_gene, AT2_gene, fibro_gene))

#check NEUTROPHIL
VlnPlot(Lung_epi_sub, features = c('LCN2', 'ORM1','MMP8','PAGE5','SYCP1','SYN3', 'CTAG2','DEFA3'
), group.by='seurat_clusters')
DotPlot(Lung_epi_sub, features = c('S100A8','S100A9','S100A6', 'CXCR4','CXCR2','CTSG','ELANE',
                                    'CD44','MPO','ZEB2', 'CST3','CSF3R','FCGR3B','CDA' ,'LY6D'
), cols = 'RdYlBu')+ coord_flip() + RotatedAxis()
VlnPlot(Lung_epi_sub, features = c('VEGFA','PLAUR', 'OLR1','CXCR2'))

cluster0 = FindMarkers(Lung_epi_sub, ident.1= 0, min.pct= 0.25)

#Lung_epi_sub = subset(Lung_epi_sub, idents = c(0,2,3,4,6,7,9))
#Lung_epi_sub = subset(Lung_epi_sub, idents = c(1:9))
#Lung_epi_sub = subset(Lung_epi_sub, idents = c(0,1,2,4,5,6,7))

 
#4,5 = cytocd8; 4 is epi (krt6 and 17 high); 1,3 is basal 
#epcam all over the place  ; 0 AND 6  ARE INFLAMMATORY 
cluster.ids = c('Inflamamtory epithelial', 'Basal', 'Epithelial', 'Basal', 
                'CytoCD8', 'CytoCD8', 'Epithelial')
names(cluster.ids) <- levels(Lung_epi_sub)
Lung_epi_sub <- RenameIdents(Lung_epi_sub, cluster.ids)
DimPlot(Lung_epi_sub, reduction = "umap", label = TRUE, repel = TRUE)  

Lung_epi_sub$celltype = Idents(Lung_epi_sub)
Lung_epi_sub@meta.data$celltype %>% table

table(Lung_epi_sub$seurat_clusters)


#Lung_at = subset(Lung, idents = c(9, 12,7,11,15))
Lung_at = SeuratAnalysis1(Lung_at)
ElbowPlot(Lung_at, ndims = 50)
Lung_at = SeuratAnalysis2(Lung_at, harmony = T, 
                               harmony_by = 'orig.ident', n = 20)
DimPlot(Lung_at,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'seurat_clusters')  


VlnPlot(Lung_at, c(
  'CD3D','CD3E', 'CD40LG','CD8A','PRDM1', 'GZMB', 'RORA', 'KLRD1',
  Macrophage_gene))

VlnPlot(Lung_at, features = c( 'TP63','KRT14',
                              'KRT15', 'KRT5', 'KRT17','EPCAM' ,
                             'KRT6A', 'KRT6B'))
VlnPlot(Lung_at, c(
                   'MUC5AC', 'SERPINB11', 'CEACAM5', 'FER1L6',  'CYP2F1','MUC16', 'MSMB' #goblet 
                   ))

#6,13 = at1; 0,1,3,8,10,12 = at 2
#7 = goblet 
#11 = t cell (Tem, no obvious cd8 or cd4); 9 = mac 
#2, 4,5,14= CILIATED; 

cluster.ids = c('Alveolar cell type 2', 'Alveolar cell type 2', "Ciliated cell", 
               'Alveolar cell type 2', 'Ciliated cell','Ciliated cell',
               'Alveolar cell type 1', 'Goblet cell',
                'Alveolar cell type 2', 'Interstitial macrophage', 'Alveolar cell type 2',
                'Tem','Alveolar cell type 2','Alveolar cell type 1','Ciliated cell')
names(cluster.ids) <- levels(Lung_at)
Lung_at <- RenameIdents(Lung_at, cluster.ids)
DimPlot(Lung_at, reduction = "umap", label = TRUE, repel = TRUE)  

Lung_at$celltype = Idents(Lung_at)
Lung_at@meta.data$celltype %>% table


VlnPlot(Lung_at, c(AT1_gene, AT2_gene, 'SCGB3A1','SCGB1A1'))

VlnPlot(Lung_at, endo_gene)

#Lung_at= subset(Lung_at, idents = 0:14)
Lung_at= subset(Lung_at, idents = c(0:3, 5, 7:14))

cluster8 = FindMarkers(Lung_at, ident.1= 8, min.pct= 0.25)




#Lung_fibro_endo = subset(Lung, idents=c(6,10,16,17))
Lung_fibro_endo = SeuratAnalysis1(Lung_fibro_endo)
ElbowPlot(Lung_fibro_endo, ndims = 50)
Lung_fibro_endo = SeuratAnalysis2(Lung_fibro_endo, harmony = T, 
                          harmony_by = 'orig.ident', n = 20)
DimPlot(Lung_fibro_endo,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'seurat_clusters')  

#1,2 = lymphatic; 14 has sema3g ; 4,8 also endo (ackr1)
#fibro= 3,5,6,10,15 ; 5,10,15 myofibro 
#0 = at2; 9,16 = club cell 
#7, 12 =interstitial mac;  11= cytocd8 


cluster.ids = c("Alveolar cell type 2",  'Lymphatic endothelial' , 'Lymphatic endothelial' ,
                'Fibro','Vein endothelial','Myofibroblast','Fibro', 'Alveolar macrophage',
                'Vein endothelial',  'Club cell','Myofibroblast', 'CytoCD8', 
                'Interstitial macrophage','Unknown', 'Artery endothelial', 'Myofibroblast',
                'Ciliated cell')
names(cluster.ids) <- levels(Lung_fibro_endo)
Lung_fibro_endo <- RenameIdents(Lung_fibro_endo, cluster.ids)
DimPlot(Lung_fibro_endo, reduction = "umap", label = TRUE, repel = TRUE)  

Lung_fibro_endo$celltype = Idents(Lung_fibro_endo)
Lung_fibro_endo@meta.data$celltype %>% table
Lung_fibro_endo = subset(Lung_fibro_endo, idents = unique(Lung_fibro_endo$celltype)[-6])

VlnPlot(Lung_fibro_endo, c(fibro_gene, endo_gene,'EPCAM'))
VlnPlot(Lung_fibro_endo, c(AT1_gene, AT2_gene, 'SCGB3A1','SCGB1A1'))

VlnPlot(Lung_fibro_endo, c(
  'CD3D','CD3E', 'CD40LG','CD8A','PRDM1', 'GZMB', 'RORA', 'KLRD1',
  Macrophage_gene))
VlnPlot(Lung_fibro_endo, features = c('FABP4', 'MCEMP1', 'PPARG',  'GPNMB','BHLHE41'))
VlnPlot(Lung_fibro_endo, SMC_gene)
VlnPlot(Lung_fibro_endo, features = c("FOXJ1", 'RFX2', 'CFAP52',
                                  'PIFO','DCDC2B','DTHD1','GSTA1'))


cluster15 = FindMarkers(Lung_fibro_endo, ident.1= 15, min.pct=0.25)
#Lung_fibro_endo = subset(Lung_fibro_endo, idents =0:18)
#Lung_fibro_endo = subset(Lung_fibro_endo, idents=c(0, 2:18))
#Lung_fibro_endo = subset(Lung_fibro_endo, idents=c(0:15,17))



#Lung_immune = subset(Lung, idents = c(0,1,3,4, 5,13,24))
Lung_immune = SeuratAnalysis1(Lung_immune)
ElbowPlot(Lung_immune, ndims = 50)
Lung_immune = SeuratAnalysis2(Lung_immune, harmony = T, 
                                  harmony_by = 'orig.ident', n = 10)
DimPlot(Lung_immune,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'celltype')  

VlnPlot(Lung_immune, features = c('TPSAB1', 'CPA3', 'KIT'))
VlnPlot(Lung_immune, c(
  'CD3D','CD3E', 'CD40LG','CD8A','PRDM1', 'GZMB', 'RORA', 'KLRD1'))
VlnPlot(Lung_immune, features =c(Macrophage_gene, 'SPP1','VEGFA','VCAN' ))
VlnPlot(Lung_immune, c('ISG15' ,'IFIT1','IFIT2',
                       'TNF','CCL3','CCL4'))

 
#1 , 10= cytocd8 (a mix with neutrophil? ); everything other than 1 and 9= mast cell;  
#8 = spp1+ mac (with s100a8/9, csf3r, vegfa, cebpb) 

cluster.ids = c('Mast cell','Mast cell','Mast cell','Mast cell','Mast cell',
                'Mast cell','Mast cell','Mast cell','SPP1+ macrophages','Mast cell','Mast cell')
names(cluster.ids) <- levels(Lung_immune)
Lung_immune <- RenameIdents(Lung_immune, cluster.ids)
DimPlot(Lung_immune, reduction = "umap", label = TRUE, repel = TRUE)  

Lung_immune$celltype = Idents(Lung_immune)
Lung_immune@meta.data$celltype %>% table




VlnPlot(Lung_immune, features = c('S100A8','S100A9','S100A6', 'CXCR4','CXCR2','CTSG','ELANE',
                                    'CD44','MPO','ZEB2', 'CST3','CSF3R','FCGR3B','CDA' ,'LY6D'
)) 
VlnPlot(Lung_immune, features = c( #"CD19",'CD79A','BANK1','MS4A1','JCHAIN',
  'KLRF1','PRF1','GZMB','CXCR6','CD160','XCL1','CD247', #MATURE 
  'NCAM1','FCGR3A','ICAM1','KLRK1','NCR1', #IMMATURE 
  'KLRC2','CD2' #ADAPTIVE
))

cluster9 = FindMarkers(Lung_immune, ident.1= 9, min.pct=0.25)
#Lung_immune = subset(Lung_immune, idents= c(0:10,12))
#Lung_immune = subset(Lung_immune, idents= c(0:4, 6:10, 12,13))
#Lung_immune = subset(Lung_immune, idents= c(0:8, 10:11))
#Lung_immune = subset(Lung_immune, idents= c(1:3, 5:10)) 
#Lung_immune = subset(Lung_immune, idents= c(3:11)) 
Lung_immune = subset(Lung_immune, idents= c(0,1,3,4,5,6,8,9)) 
Lung_immune = subset(Lung_immune, idents= c(0,2:9, 11)) 



#Lung_immune_sub = subset(Lung_immune, idents = 0:2)
#Lung_immune_sub2= subset(Lung_immune, idents=c(1,10))
#Lung_immune_sub = merge(Lung_immune_sub, Lung_immune_sub2)
Lung_immune_sub = SeuratAnalysis1(Lung_immune_sub)
ElbowPlot(Lung_immune_sub, ndims = 50)
Lung_immune_sub = SeuratAnalysis2(Lung_immune_sub, harmony = T, 
                              harmony_by = 'orig.ident', n = 20)
DimPlot(Lung_immune_sub,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'seurat_clusters')  
table(Lung_immune_sub$seurat_clusters)

VlnPlot(Lung_immune_sub, c(
  'CD3D','CD3E', 'CD40LG','CD8A','PRDM1', 'GZMB', 'RORA', 'KLRD1'))
VlnPlot(Lung_immune_sub, features =c(Macrophage_gene, 'SPP1','VEGFA','VCAN' ))
VlnPlot(Lung_immune_sub, c('ISG15' ,'IFIT1','IFIT2',
                       'TNF','CCL3','CCL4'))
VlnPlot(Lung_immune_sub, features = c('LCN2', 'ORM1','MMP8','PAGE5','SYCP1','SYN3', 'CTAG2','DEFA3'
), group.by='seurat_clusters')
VlnPlot(Lung_immune_sub, features = c('S100A8','S100A9','S100A6', 'CXCR4','CXCR2','CTSG','ELANE',
                                  'CD44','MPO','ZEB2', 'CST3','CSF3R','FCGR3B','CDA' ,'LY6D'
)) 


Lung_immune_sub[['percent.rp']] = PercentageFeatureSet(Lung_immune_sub, pattern = "^RP[SL]")
VlnPlot(Lung_immune_sub, features = 'percent.rp', group.by = 'orig.ident')
#Unknown_NK_sub2 = subset(Unknown_NK_sub2, percent.rp <=15) #but then not many cells left, so I think we can just get rid of the whole thing


cluster1 = FindMarkers(Lung_immune_sub, ident.1= 1, min.pct=0.25)
# 7 = cytocd8; ALL OTHERS SEEM TO BE JUNK WITH HIGH RP GENES 
#0 is junk with ig genes; 1 has a lot of mt genes 

Lung_immune_sub= subset(Lung_immune_sub, idents=c(7))
Lung_immune_sub= subset(Lung_immune_sub, idents=c(2:4))

#so the remaining is just cytocd8 

cluster.ids = rep("CytoCD8", 4)
names(cluster.ids) <- levels(Lung_immune_sub)
Lung_immune_sub <- RenameIdents(Lung_immune_sub, cluster.ids)
DimPlot(Lung_immune_sub, reduction = "umap", label = TRUE, repel = TRUE)  

Lung_immune_sub$celltype = Idents(Lung_immune_sub)
Lung_immune_sub@meta.data$celltype %>% table


#Lung_others =subset(Lung, idents = c(8,19,20,25 ))
Lung_others = SeuratAnalysis1(Lung_others)
ElbowPlot(Lung_others, ndims = 50)
Lung_others = SeuratAnalysis2(Lung_others, harmony = T, 
                                  harmony_by = 'orig.ident', n = 10)
DimPlot(Lung_others,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'seurat_clusters')  


VlnPlot(Lung_others, c(
  'CD3D','CD3E', 'CD40LG','CD8A','PRDM1', 'GZMB', 'RORA', 'KLRD1',
  Macrophage_gene)) 
#SCHWANN 
#4 is at2 
#ALL CLUB CELLS EXCEPT FOR 2 AND 7 

VlnPlot(Lung_others, c(SMC_gene, AT2_gene, 'SCGB3A1', 'SCGB1A1'))
VlnPlot(Lung_others, features = c("FOXJ1", 'RFX2', 'CFAP52', 
                                 'PIFO','DCDC2B','DTHD1','GSTA1'))
VlnPlot(Lung_others, features =c(Macrophage_gene, 'SPP1','VEGFA','VCAN' ))
VlnPlot(Lung_others, c("SCN7A", 'NCAM1', 'CHD2', 'L1CAM', 'NGFR')) #9 may be schwann?? 
#schwann: also have CHGA/CHGB family, SCGN family that points to neural identity 
table(Lung_others$seurat_clusters)

cluster6 = FindMarkers(Lung_others, ident.1= 6, min.pct= 0.25)
#Lung_others= subset(Lung_others, idents=c(0:6, 8))

cluster.ids = c("Club cell", 'SPP1+ macrophages', 'Club cell',
                'Club cell','Alveolar cell type 2','Ciliated cell',
                'Club cell','Nonmyelinating schwann cell')
names(cluster.ids) <- levels(Lung_others)
Lung_others <- RenameIdents(Lung_others, cluster.ids)
DimPlot(Lung_others, reduction = "umap", label = TRUE, repel = TRUE)  

Lung_others$celltype = Idents(Lung_others)
Lung_others@meta.data$celltype %>% table



#merging lung structural ----- 
Lung_final = merge(x= Lung_epi, y=c(Lung_epi_sub, Lung_fibro_endo, Lung_at,
                                    Lung_immune, Lung_immune_sub, Lung_others))

Lung_final = SeuratAnalysis1(Lung_final)
ElbowPlot(Lung_final, ndims = 50)
Lung_final = SeuratAnalysis2(Lung_final, harmony = T, 
                              harmony_by = 'orig.ident', n = 25)
DimPlot(Lung_final,
        reduction = "umap",label = T,repel = T,
        raster = FALSE ,label.size = 4, group.by = 'celltype')  

cluster19 = FindMarkers(Lung_final, ident.1= 19, min.pct= 0.25)

Lung_final$cell_cluster= paste0(Lung_final$seurat_clusters, '_', Lung_final$celltype)
table(Lung_final$cell_cluster)
table(Lung_final$seurat_clusters)

Lung_structural_genes = FindAllMarkers(Lung_final, group.by = 'celltype', min.pct = 0.25,
                                         only.pos = TRUE)  
Lung_structural_genes_filtered = Lung_structural_genes %>% 
  group_by(cluster) %>% 
  slice_max(order_by = avg_log2FC, n = 20, with_ties = T)

Lung_final$celltype= factor(Lung_final$celltype, levels = c(
  names(table(Lung_final$celltype))[c(1,2,3,7,10, 16,
                                      18, 5, 11,19,8, 13, 4,
                                      9, 15, 6, 17, 12, 14)]
))

DotPlot(Lung_final, features = c( 
  AT1_gene[C(1:5)],
  AT2_gene, 
  'TP63','KRT14', 'KRT15', 'KRT5', 'KRT17','EPCAM' ,   
  'CXCL10','CXCL9',
  'KRT6A', 'KRT6B', 
  'TFF1', 'TFF2','S100P', 
  'SCGB1A1','SCGB3A2', #club
  "PECAM1", "NOTCH4", "TEK",'VWF', #endo 
  'LYVE1', 'PROX1','CCL21','FLT4', #LYMPHATIC
  'ACKR1', 'HDAC9','PKHD1L1','PTGS1','ZNF385D','OLFM1','TPD52L1', #VEIN  
   "COL1A1", 'COL1A2', 'SPARC', 'DCN','LUM', #fibro 
  'ITGBL1','ASPN','LMOD1','NALCN','WISP2', 'ACTA2', 'TAGLN',  #MYOFIBRO 
  "FOXJ1", 'RFX2', 'CFAP52', 'PIFO','DCDC2B','DTHD1','GSTA1', #CILIATED
  #SMC_gene, 
  'MUC5AC','TSPAN8', 'CEACAM5','VSIG2', #goblet  
  "CADM1",  'NFASC','NRXN1','GRIK2', #schwann 
  'CD3D','CD3E','CD8A','PRDM1', 'GZMB', 'RORA', 'KLRD1', 
  "TPSAB1",'KIT','CPA3',
  'ICAM1','CXCR4','CXCR2','S100A8','S100A9','CSF3R'
) %>% unique(), cols='RdYlBu', group.by = 'celltype') + coord_flip() + RotatedAxis()+ 
  theme(axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 8))

saveRDS(Lung_final, file ='Lung_structural_v2.rds')


table(Lung_final$celltype)
VlnPlot(Lung_final, features= c('CD3D','CD3E','CD8A','PRDM1', 'GZMB', 'RORA', 'KLRD1'), group.by = 'celltype')

#clean up 
Idents(Lung_final) ='celltype'
#Lung_mac = subset(Lung_final, idents= c('Alveolar macrophage', 'SPP1+ macrophages', 
                                       # 'Interstitial macrophage', 'Fibro'))
Lung_mac = SeuratAnalysis1(Lung_mac)
ElbowPlot(Lung_mac, ndims = 50)
Lung_mac = SeuratAnalysis2(Lung_mac, harmony = T, harmony_by = 'Patient', n = 15)
DimPlot(Lung_mac,
        reduction = "umap",label = T, repel = T,
        raster = FALSE ,label.size = 4, group.by = 'seurat_clusters')  

table(Lung_mac$seurat_clusters)
VlnPlot(Lung_mac, features=c( 'FABP4', 'MCEMP1', 'PPARG',  'GPNMB','BHLHE41',
                              'C1QA','C1QB', 'CD163', 'MRC1',
                              'SPP1','VCAN', 'VEGFA', 
                              fibro_gene))
VlnPlot(Lung_mac, endo_gene)
VlnPlot(Lung_mac, c(AT1_gene, AT2_gene, 'SCGB1A1','SCGB3A1'))
#0,3,6= fibro ; 1 = lymphatic; 4 = vein endothelial ; 5=normal endo  ; 2,7 = AT2
cluster.ids = c("Fibro", 'Lymphatic endothelial', 'Alveolar cell type 2',
                'Fibro', 'Vein endothelial', 'Endothelial', 'Fibro',  'Alveolar cell type 2')
names(cluster.ids) <- levels(Lung_mac)
Lung_mac <- RenameIdents(Lung_mac, cluster.ids)
DimPlot(Lung_mac, reduction = "umap", label = TRUE, repel = TRUE)  

Lung_mac$celltype = Idents(Lung_mac)
Lung_mac@meta.data$celltype %>% table

cluster7 = FindMarkers(Lung_mac, ident.1 = 7, min.pct=0.25)

#Lung_mac = subset(Lung_mac, idents= c(0,1,2,3,4,5))
#Lung_mac = subset(Lung_mac, idents= 1:7)
Lung_mac = subset(Lung_mac, idents= c(0,2:8)) 
Lung_mac = subset(Lung_mac, idents= c(0:6, 8))
Lung_mac = subset(Lung_mac, idents= c(0:4, 6))


#Lung_spp1 = subset(Lung_mac, idents = 0)
Lung_spp1 = SeuratAnalysis1(Lung_spp1)
ElbowPlot(Lung_spp1, ndims = 50)
Lung_spp1 = SeuratAnalysis2(Lung_spp1, harmony = T, harmony_by = 'Patient', n = 10)
DimPlot(Lung_spp1,
        reduction = "umap",label = T, repel = T,
        raster = FALSE ,label.size = 4, group.by = 'seurat_clusters')  

table(Lung_spp1$seurat_clusters)
VlnPlot(Lung_spp1, features=c( #'FABP4', 'MCEMP1', 'PPARG',  'GPNMB','BHLHE41',
                              'C1QA','C1QB', 'CD163', 'MRC1',
                              'SPP1','VCAN', 'VEGFA', 
                              'EPCAM','KRT19', 'PECAM1','VWF', 
                              Cycling_genes
                              ))
VlnPlot(Lung_spp1, features = c('ICAM1','CXCR4','CXCR2','S100A8','S100A9','CSF3R'))

cluster0 = FindMarkers(Lung_spp1, ident.1 = 0, min.pct=0.25)
#2,3 = neutrophil ; 0 and 1 hs epcam 
 
cluster.ids = c('Epithelial', 'Epithelial',  'Neutrophil' ,'Neutrophil')
names(cluster.ids) <- levels(Lung_spp1)
Lung_spp1 <- RenameIdents(Lung_spp1, cluster.ids)
DimPlot(Lung_spp1, reduction = "umap", label = TRUE, repel = TRUE)  

Lung_spp1$celltype = Idents(Lung_spp1)
Lung_spp1@meta.data$celltype %>% table

#Lung_spp1 = subset(Lung_spp1, idents = c(0,1,3,5))

#merge back to lung 
Idents(Lung_final) ='celltype'
Lung_rest = subset(Lung_final, idents = unique(Lung_final$celltype)[c(1:9, 12,13,15:19,21,22)])
Lung_final = merge(Lung_rest, c(Lung_spp1, Lung_mac))

#still a lot of background immune genes in basal / epithelial; artery endo has mixture of smc genes 
#club cells have mix  
Lung_final = subset(Lung_final, idents = unique(Lung_final$celltype)[c(1:12, 14:20)])


#merging all structural --------
Structural_0501 = merge(x = Lung_final, y = c(kidney_final, HCC, hnscc, breast_all))
Structural_0501 =  JoinLayers(Structural_0501)
Idents(Structural_0501) ='celltype'
unique(Structural_0501$celltype)
Structural_0501= RenameIdents(Structural_0501, 'CD8+ T cell' ='CytoCD8', 'SMC' = 'Smooth muscle cell'
                              )
Structural_0501$celltype = Idents(Structural_0501)

#Structural_0504 = merge(Structural_rest, check)
Structural_0506 = SeuratAnalysis1(Structural_0506)
ElbowPlot(Structural_0506, ndims = 50)
Structural_0506 = SeuratAnalysis2(Structural_0506, harmony = F, n = 30)
DimPlot(Structural_0506,
        reduction = "umap",label = T, repel = T,
        raster = FALSE ,label.size = 4, group.by = 'celltype')  

#saveRDS(Structural_0501, file ='All_structural_v1.rds')
#saveRDS(Structural_0504, file ='All_structural_v2.rds')
#saveRDS(Structural_0506, file = 'All_structural_v3.rds')

structural_genes = FindAllMarkers(Structural_0506, group.by = 'celltype', min.pct = 0.25,
                                       only.pos = TRUE)  
structural_genes_filtered = structural_genes %>% 
  group_by(cluster) %>% 
  slice_max(order_by = avg_log2FC, n = 20, with_ties = T)
write.csv(structural_genes_filtered, file = 'Structural_cell_top20_markers.csv')

#endothelial cell?? cholangiocyte? 
#skeletal muscle is valid, just with a few RP in top 20  

DotPlot(Structural_0506, features= c('AREG','MUCL1','DBI','AZGP1','PIP','KRT18','AGR2',  'ESR1','AR', 'PGR', #LUM-hr
                                    'SCGB2A2','SLPI','WFDC2','LTF','KRT15','MMP7','SCGB3A1','MYLK'),
        group.by = 'celltype', cols = 'RdYlBu') + coord_flip() + RotatedAxis()

DotPlot(Structural_0506, features = c( 'EPCAM', 'KRT19', 'PLPP2', 'CLDN4', 'FXYD2', 'FAM3B',  #cholangiocyte 
                           'CLEC4G', 'FCN2', 'CLEC4M', 'IFITM10', 'STAB2', 'FCN3', #sinusoidal 
                           'AIF1L', 'MGP', 'SLCO2A1', 'CLEC14A', #macrovascular 
                           'HPR', 'F12',  'ALB', 'ALDOB', 'AGXT', 'SAA4','ARG1' 
), group.by = 'celltype', cols = 'RdYlBu') + coord_flip() + RotatedAxis()

DotPlot(Structural_0506, c('MYL9','CDH5', SMC_gene, #SMC 
                 'MYH1', 'MYOG', 'TNNC1', 'TNNI1', #SKELETAL 
                 'PAX7','MYF5'), 
        group.by = 'celltype', cols = 'RdYlBu') + coord_flip() + RotatedAxis()  

VlnPlot(Structural_0504, features = c("SCN7A", 'NCAM1', 'CHD2', 'L1CAM', 'NGFR'), 
        group.by = 'celltype') #non-myelinating schwann seem not too confident 

#Structural_rest = subset(Structural_0501, 
                       #  idents = unique(Structural_0501$celltype)[-29])
Structural_rest = subset(Structural_0504, cell_condition !='Fibro_hnscc' & celltype != 'Tem')
Fibro_metadata = Fibro_hnscc@meta.data[, 1:12]
Fibro_export = CreateSeuratObject(counts = LayerData(Fibro_hnscc, assay = "RNA", layer = "counts"), 
                                    meta.data= Fibro_metadata) 

Structural_0506= merge(Structural_rest, Fibro_export)

#look into fibro of hnscc ------
Fibro_hnscc = subset(Structural_0506, cell_condition =='Fibro_hnscc')
head(Fibro_hnscc@meta.data)
Fibro_hnscc@meta.data = Fibro_hnscc@meta.data %>% mutate(
  Fibro_type = case_when(
    Fibro1 >=0.5 ~ "hi",
    TRUE ~ 'low'
  )
)

#Fibro_hnscc_hi = subset(Fibro_hnscc, Fibro1 >= median(Fibro1))
Fibro_hnscc = SeuratAnalysis1(Fibro_hnscc)
ElbowPlot(Fibro_hnscc, ndims = 50)
Fibro_hnscc = SeuratAnalysis2(Fibro_hnscc, harmony = T, harmony_by = 'Patient', n= 20)
DimPlot(Fibro_hnscc,
        reduction = "umap",label = T, repel = T,
        raster = FALSE ,label.size = 4, group.by = 'seurat_clusters')

VlnPlot_scCustom(Fibro_hnscc, features = 'Fibro1',
                 plot_median = T, group.by='Fibro_type' , raster=F) +NoLegend()# +ylim(-0.5, 1.5)

cluster8 = FindMarkers(Fibro_hnscc, ident.1 = 8, min.pct=0.25)
FeaturePlot(Fibro_hnscc, c(fibro_gene, 'LUM'), cols = c('gray','red'), label = T)
VlnPlot(Fibro_hnscc, c(fibro_gene, SMC_gene), group.by = 'seurat_clusters')
VlnPlot(Fibro_hnscc, features = c('PLA2G2A','CFD'), group.by = 'seurat_clusters')

#Fibro_hnscc = subset(Fibro_hnscc, idents = 0:11)
Fibro_hnscc = subset(Fibro_hnscc, idents = 0:12)
saveRDS(Fibro_hnscc, 'Fibroblast_hnscc.rds')

Idents(Fibro_hnscc) ='Fibro_type'
Fibro_lo_gene = FindMarkers(Fibro_hnscc, ident.1 = 'low', min.pct=0.25) %>% 
  filter(avg_log2FC > 0.5)
write.csv(Fibro_lo_gene, file = 'Fibro_hnscc_low_gene.csv')

Fibro_hi_gene = FindMarkers(Fibro_hnscc, ident.1 = 'hi', min.pct=0.25) %>% 
  filter(avg_log2FC > 0.5)
write.csv(Fibro_hi_gene, file = 'Fibro_hnscc_hi_gene.csv')


#look into artery vs vein endo 
artery = structural_genes_filtered %>% filter(cluster=='Artery endothelial')
vein = structural_genes_filtered %>% filter(cluster=='Vein endothelial')
artery$gene[artery$gene %in% vein$gene]

Tem = structural_genes_filtered %>% filter(cluster=='Tem')
CytoCD8 = structural_genes_filtered %>% filter(cluster=='CytoCD8')
CytoCD8$gene[!CytoCD8$gene %in% Tem$gene]


#merge everything-----------
Lymphoid= readRDS("~/jliang/Human tissue harmonization/Lymphoid/Lymphoid_final_v3_0312.rds")
colnames(Lymphoid@meta.data)
Lymphoid_metadata = Lymphoid@meta.data[, 1:13]
Lymphoid_export = CreateSeuratObject(counts = LayerData(Lymphoid, assay = "RNA", layer = "counts"), 
                                     meta.data= Lymphoid_metadata)


Myeloid = readRDS("~/jliang/Human tissue harmonization/Myeloid/Myeloid_final_v7_rename.rds") 
Myeloid_metadata = Myeloid@meta.data[, 1:11]
Myeloid_export = CreateSeuratObject(counts = LayerData(Myeloid, assay = "RNA", layer = "counts"), 
                                    meta.data= Myeloid_metadata)


Structural = readRDS("~/jliang/Human tissue harmonization/Structural/All_structural_v3.rds")
Structural = JoinLayers(Structural)
Structural_metadata = Structural@meta.data[, 1:12]
Structural_export = CreateSeuratObject(counts = LayerData(Structural, assay = "RNA", layer = "counts"), 
                                       meta.data= Structural_metadata)

#Total = subset(Total,idents = c(0:48))
#Total_rest = subset(Total, ! celltype %in% c('trNK', 'Mature NK'))


Total_rest = subset(Total, ! celltype %in% c('Inflammatory macrophage','Th2',
                                             'Exhausted T cell'))
Total_metadata = Total_rest@meta.data
Total_rest_export =  CreateSeuratObject(counts = LayerData(Total_rest, assay = "RNA", layer = "counts"), 
                                        meta.data= Total_metadata)
#NK = readRDS("~/jliang/Human tissue harmonization/check NK/Mature_trNK_final.rds")
InflamMac = readRDS("~/jliang/Human tissue harmonization/Myeloid/Inflammatory_macrophage_clean.rds")
Exhausted = readRDS("~/jliang/Human tissue harmonization/Exhausted_Tcell_clean.rds")
#Total = merge(Total_rest_export, NK)
Total = merge(Total_rest_export, c(InflamMac, Exhausted))

Total = merge(x = Lymphoid_export, y = c(Lymphoid2_export, Myeloid_export, Structural_export))
table(Total$celltype)
table(Total$orig.ident)
table(Total$Patient)

Idents(Total) ='orig.ident'
Total = RenameIdents(Total, 'hnscc' = 'HNSCC', 'Kidney' ='KIRC', 'Breast' ='BRCA')
Total$orig.ident = Idents(Total)
Idents(Total) ='celltype'
Total = RenameIdents(Total, 'IFN-primed Macrophage' = 'IFN-primed macrophages', 'AT2' ='Alveolar cell type 2',
                     'Macrophage' ='Interstitial macrophages',
                     'IgE+ BCell' ='BCell','Plasmablast' ='Plasma')
Total = RenameIdents(Total, 'Interstitial macrophages' ='Interstitial macrophage')
Total$celltype = Idents(Total)

Total = RenameIdents(Total, 'ILC1' = 'Immature NK')
Total$celltype = Idents(Total)

saveRDS(Total, "~/jliang/Human tissue harmonization/Merged_object_v3_UMAP.rds") 


Total$celltype_condition = paste0(Total$celltype, '_', Total$orig.ident)
table(Total$celltype_condition)


Total = SeuratAnalysis1(Total)
ElbowPlot(Total, ndims = 50)
Total = SeuratAnalysis2(Total, harmony = T,
                        harmony_by = 'orig.ident', n = 30)
p1= DimPlot(Total,
             reduction = "umap",label = T,repel = T,
             raster = FALSE ,label.size = 4, group.by = 'celltype')  
p1 #celltype
ggsave(p1, file = 'Total_celltype_umap_legend_neighbor20_v3.png', width = 30, height= 15)

p4 = DimPlot(Total,
             reduction = "umap",label = F,repel = T,
             raster = FALSE ,label.size = 4, group.by = 'orig.ident') + 
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 4)))    

ggsave(p4, file = 'Total_cancer_type_umap_neighbor20_v3.png', width = 25, height= 15)


p5 = DimPlot(Total,
             reduction = "umap",label = F,repel = T,
             raster = FALSE ,label.size = 4, group.by = 'condition') + 
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 4)))    

ggsave(p5, file = 'Total_paired_condition_umap_neighbor20_v3.png', width = 25, height= 15)


p3 = DimPlot(Total,
            reduction = "umap",label = T,repel = T,
            raster = FALSE ,label.size = 4, group.by = 'seurat_clusters')  + NoLegend() 
ggsave(p3, file = 'Total_seurat_clusters_umap_neighbor20_v3.png', width = 20, height= 15)

table(Total$celltype)



#check the merged object-----
 
DotPlot(Total, features = c(    'PDGFRB', 'PLN', 'RGS5', 'TAGLN', 'ACTA2', #SMC
                               'IGFBP5', 'MYLK', 'RERGL', 'CASQ2', "CNN1",'MYH11', 'SYNM','ACTG2'),
        group.by = 'seurat_clusters', cols ='RdYlBu')
                          
DotPlot(Total, features = c( 'APOC1', 'CD9', 'TREM2', 'APOE', 'FABP5', 
                             'C1QA', 'C1QB', 'CD163','CD68','MRC1',
                            # 'GPNMB','BHLHE41', 'PPARG', 'MCEPM1',
                             'CXCL9','CXCL10','CXCL11',
                             'SPP1','VCAN','VEGFA'),
        group.by = 'celltype', cols ='RdYlBu') + coord_flip() + RotatedAxis()

DotPlot(Total, features = c( 'TBX21','IFNG','TNF','CD40LG','FASLG','ANXA3','RUNX3',
                             'IL2','LTA'),
        group.by = 'celltype', cols ='RdYlBu') + coord_flip() + RotatedAxis()

DotPlot(Total, features = c('EGR1','NFKB1','CCL4','CXCR4','AREG'  ),
        group.by = 'celltype', cols ='RdYlBu') + coord_flip() + RotatedAxis()


table(Total$seurat_clusters)                                      
#check 23, 43, 49 
cluster47 = FindMarkers(Total, ident.1 = 47, min.pct=0.25)

Total$celltype_cluster = paste0(Total$seurat_clusters, '_', Total$celltype)
table(Total$celltype_cluster)


#first output lineage genes -----------
All_lineages = FindAllMarkers(Total, group.by = 'celltype', min.pct = 0.25,
                              only.pos = TRUE) 
All_lineages_filtered = All_lineages %>% 
  group_by(cluster) %>% 
  slice_max(order_by = avg_log2FC, n = 50, with_ties = T)
write.csv(All_lineages, file = 'Total_findallmarkers_v3.csv')

gene_list <- split(All_lineages_filtered$gene, All_lineages_filtered$cluster)
max_len <- 50
gene_list_padded <- lapply(gene_list, function(x) {
  length(x) <- max_len
  return(x)
})
final_df <- as.data.frame(gene_list_padded)
write.csv(final_df, file = '~/jliang/Human tissue harmonization/Total_lineage_top50_v3.csv')

#sanity check gene modules and lineages 
Myeloid_genes = read.csv("~/jliang/Human tissue harmonization/Myeloid/Myeloid_all_cell_all_findmarkers_v7_rename.csv",
                         row.names = 1)
Myeloid_genes = Myeloid_genes %>% 
  group_by(cluster) %>% 
  slice_max(order_by = avg_log2FC, n = 50, with_ties = T)

interstitial_gene = Myeloid_genes %>% filter(cluster == 'Interstitial macrophages')
LAM_gene = Myeloid_genes %>% filter(cluster == 'LAM')
Alveolar_gene = Myeloid_genes %>% filter(cluster == 'Alveolar macrophage')

Inflammatory_gene = read.csv("Total_findallmarkers.csv", row.names = 1) %>% 
  filter(cluster == 'Inflammatory macrophage')%>% 
  slice_max(order_by = avg_log2FC, n = 50, with_ties = T)

#compare them 
interstitial_gene$gene[interstitial_gene$gene %in% Inflammatory_gene$gene ]

Total = AddModuleScore(object = Total,
                                 features =list(c(Inflammatory_gene$gene,'TNF', 'CCL3','CCL3L1','CCL4L1','NEDD9','KLF6')),
                                 ctrl = 5,
                                 name = 'Inflammatory_macrophage')
  

VlnPlot_scCustom(Total, features = 'Inflammatory_macrophage1',
                 plot_median = T, group.by='celltype' , raster=F) +NoLegend()# +ylim(-0.5, 1.5)


gene_list <- split(Myeloid_genes$gene, Myeloid_genes$cluster)
max_len <- 50
gene_list_padded <- lapply(gene_list, function(x) {
  length(x) <- max_len
  return(x)
})
myeloid_df <- as.data.frame(gene_list_padded)
write.csv(myeloid_df, file = '~/jliang/Human tissue harmonization/Myeloid_lineage_top50.csv')


#module score----------
module_score = read.csv('~/jliang/Human tissue harmonization/gene_modules_final.csv', row.names = 1,
                        check.names = F) #check.names can prevent the trailing dots

#and then modify colnames and rownames 
unique(Total$celltype) [! unique(Total$celltype) %in% colnames(module_score) ] 
colnames(module_score) [ ! colnames(module_score) %in% unique(Total$celltype)]

#write.csv(module_score, file = 'gene_modules_final.csv')

#rename after examining first round of checking the module scores 
Total$celltype2 = Total$celltype
Idents(Total) ='celltype2'

Total = RenameIdents(Total, 'CD8+ Tcm' ='Activated CD8+ T cells', 
                     'CytoCD8' = 'Activated CD8+ T cells',
                    'TRTex' ='CD8+ Trm',
                    'CD4+ Tem' ='gdT' )
Total$celltype2 = Idents(Total)
table(Total$celltype2)
Total$celltype = Total$celltype2

CytoCD8_genes = gsub('^ ', '', module_score$CytoCD8) 
Total= AddModuleScore(object = Total,
                      features =list(CytoCD8_genes),
                      ctrl = 5,
                      name = 'Activated CD8+ T cells')


head(Total@meta.data)
colnames(Total@meta.data)
Total@meta.data= Total@meta.data[, c(1:92)]

for (i in 1:length(unique(Total$celltype))){ 
  genes = gsub('^ ', '', module_score[, i]) 
  cell_name = colnames(module_score)[i]
  Total= AddModuleScore(object = Total,
                                features =list(genes),
                                ctrl = 5,
                                name = cell_name)
}


VlnPlot_scCustom(Total, features = c('gdT_score_v2' ),
                 plot_median = T, group.by='celltype' , raster=F) +NoLegend()


VlnPlot_scCustom(Total, features = c(  'Activated CD8+ T cells1' ), 
                 plot_median = T, group.by='celltype2' , raster=F) +NoLegend()


#make a heatmap 
hmat = Total@meta.data

library(tibble) 
colnames(Total@meta.data)[93:168] = gsub('1$', '_score_v2', colnames(Total@meta.data)[93:168] ) 

heat_df <- Total@meta.data %>%
  group_by(celltype) %>%
  summarise(across(ends_with('_score_v2'), median, na.rm = TRUE)) #calculate the median score for each module score per cell type

heat_mat <- heat_df %>%
  column_to_rownames("celltype") %>%
  as.matrix()

heat_mat_scaled <- t(scale(heat_mat)) #for one score, what cell type is it represented the most
names = gsub('_score_v2', '',rownames(heat_mat_scaled)) 
names = gsub('_', ' ', names) 
sum(names %in% colnames(heat_mat_scaled))
heat_mat_scaled <- heat_mat_scaled[, names] #order 

library(pheatmap)
pheatmap(
  heat_mat_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("lightblue", "white", "red"))(100),
  border_color = NA,
  main = 'Median module scores across all celltypes'
) 
 

FeaturePlot(Total, features = c(colnames(Total@meta.data)[14:20]), cols = c('gray','red'),
            label = F, raster = F, reduction = 'umap')



VlnPlot_scCustom(Total, features = c("Alveolar macrophage_score", 'LAM_score',"SPP1+ LAM_score", 'CXCL9/10/11+ LAM_score'
                                     #'Interstitial macrophage_score', 'IFN-primed macrophages_score',
                                     #'Kupffer cells_score', 'PVM_score'
                                     ),
                 plot_median = T, group.by='celltype' , raster=F) +NoLegend() +  
  theme(axis.text.x = element_text(size = 5)) 

Stacked_VlnPlot(Total, features = c("Alveolar macrophage_score", 'LAM_score',"SPP1+ LAM_score", 'CXCL9/10/11+ LAM_score'),
                 group.by='celltype' , raster=F)+NoLegend() +  
  theme(axis.text.x = element_text(size = 8, angle = 45)) 



