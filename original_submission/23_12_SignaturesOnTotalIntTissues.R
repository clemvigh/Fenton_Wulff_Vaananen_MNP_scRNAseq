###### Libraries ######
library(gplots)
library(dplyr)
library(stringr)
library(scales)
library(ggplot2)
library(viridis)
library(scales)
library(Seurat)
library(ggdendro)
library(plotly)
library(gprofiler2)
library(ggrastr)

### Clear and read in data ###
rm(list=ls())

### variables
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)
# colour string for imputation and overlays
mycols_b <- c("#bdbdbd","#d9d9d9","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F")
mycols <- rev(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))
# project
project <- "MNP_LP"
#trajectory subset
round="CD14CD1C_prol"
#int res 0.1, cl 3 = cDC1
F2_col_v3 <- c("cDC1"="#F564E3","cDC2"="#9590FF","cDC3"="#D89000","amb"="lightgrey","HLA low"="darkgrey","mono/mac"="#6BB100")


#### data read in and variables ####
MNP <- readRDS("/Volumes/LWulffExD/Projects/FentonWulffData/R6/Rdata/R2_MNPs_MNP_LP_wCITEdsbnorm.rds")
MNP <- UpdateSeuratObject(MNP)
#meta
CD14CD1Cmeta <- readRDS("/Volumes/LWulffExD/Projects/FentonWulffData/R6/Rdata_v2/R3_CD14CD1Cprol_MNPLP_METAONLY_AUG22.rds")
monomac <- readRDS("/Volumes/LWulffExD/Projects/FentonWulffData/R6/Rdata_v2/monomac_traj_versionApril2022.rds")
monomac <- monomac@meta.data

MNP@meta.data$F2F5 <- NA
MNP@meta.data[rownames(CD14CD1Cmeta),]$F2F5 <- CD14CD1Cmeta$semisupDC_v5
MNP@meta.data[MNP@meta.data$integrated_snn_res.0.1==3,]$F2F5 <- "cDC1"
MNP@meta.data[MNP@meta.data$integrated_snn_res.0.1==6,]$F2F5 <- "CD16+ mono"
MNP@meta.data[MNP@meta.data$integrated_snn_res.0.1==5,]$F2F5 <- "pDC"
MNP@meta.data[rownames(monomac),]$F2F5 <- as.character(monomac$tSP_clustering_F4)

DimPlot(MNP, group.by = "F2F5", label = T)

Idents(MNP) <- "F2F5"
DefaultAssay(MNP) <- 'RNA'
DEGs <- FindAllMarkers(MNP, logfc.threshold = 0.5, only.pos = T)
DEGs$ENSG <- gconvert(DEGs$gene,organism="hsapiens",target="ENSG",filter_na = F,mthreshold = 1)$target

extras <- c("HIST1H4C"="ENSG00000197061","KIAA0101"="ENSG00000166803","H2AFZ"="ENSG00000164032",
            "LINC01272"="ENSG00000224397", "SEPP1"="ENSG00000250722")
DEGs[DEGs$gene %in% names(extras),]$ENSG <- extras[DEGs[DEGs$gene %in% names(extras),]$gene]

top50 <- DEGs %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

#### Read in external data set and add each DEG set as module ####
Elm <- readRDS("/Users/linewulff/Documents/work/projects/FentonWulff_LP_MNP/Lit_data/ElmentaiteNature2021.rds")
# Elm <- subset(Elm, cells = rownames(Elm@meta.data[Elm@meta.data$Diagnosis=="Healthy adult",]))
head(Elm@meta.data)
DimPlot(Elm, group.by = "category", label = T)+NoLegend()

for (clus in unique(top50$cluster)){
  print(clus)
  clus_top50 <- top50[top50$cluster==clus,]$ENSG
  if (length(clus_top50)>10){
    nbin = length(clus_top50)} 
  else {nbin = 10}
  print(length(clus_top50))
  Elm <- AddModuleScore(Elm, features = list(clus_top50), name = clus, nbin = nbin)
}

HLAs <- rownames(MNP@assays$RNA@data)[startsWith(rownames(MNP@assays$RNA@data), "HLA")]
HLAs <- gconvert(HLAs,organism="hsapiens",target="ENSG",filter_na = F)$target
Elm <- AddModuleScore(Elm, features = list(HLAs), name = "HLAs", nbin = nbin)

spec_genes <- c("CD1C","CD207","PLAC8","PKIB","PPA1","SLC38A1","FCER1A","CD1E","LSP1","LTB")
spec_genes <- gconvert(spec_genes,organism="hsapiens",target="ENSG",filter_na = F)$target
Elm <- AddModuleScore(Elm, features = list(spec_genes), name = "cDC2sup", nbin = 10)

FeaturePlot(Elm, features = top50[top50$cluster=="HLA low",]$ENSG)

Elm <- cbind(Elm@meta.data, Elm@reductions$umap@cell.embeddings)
colnames(Elm)[35:55] <- str_sub(colnames(Elm)[35:55], start = 1, end = -2)
colnames(Elm)[35:55] <- str_replace(colnames(Elm)[35:55], " ","")
colnames(Elm)[35:55] <- str_replace(colnames(Elm)[35:55], "/","")

for (clus in colnames(Elm)[35:54]){
  moduleplot <- ggplot(Elm, aes(x = UMAP_1, y = UMAP_2, colour = Elm[,clus]))+
  geom_point_rast()+
  scale_color_gradientn(colours = mycols, name = clus)+
  theme_classic()
  
  pdf(paste(dato,"Elmentaiteetal_DEGsbasedModulescores_",clus,".pdf",sep = ""),height = 5, width = 5)
  print(moduleplot)
  dev.off()
}
for (clus in colnames(Elm)[35:54]){
  moduleplot <- ggplot(Elm[Elm$Diagnosis=="Healthy adult",], aes(x = UMAP_1, y = UMAP_2, colour = Elm[Elm$Diagnosis=="Healthy adult",clus]))+
    geom_point_rast()+
    scale_color_gradientn(colours = mycols, name = clus)+
    theme_classic()
  
  pdf(paste(dato,"Elmentaiteetal_HEALTHYADULT_DEGsbasedModulescores_",clus,".pdf",sep = ""),height = 5, width = 5)
  print(moduleplot)
  dev.off()
}

#### Subcluster on myeloid and repeat ####
Elm <- readRDS("/Users/linewulff/Documents/work/projects/FentonWulff_LP_MNP/Lit_data/ElmentaiteNature2021.rds")
# Elm <- subset(Elm, cells = rownames(Elm@meta.data[Elm@meta.data$Diagnosis=="Healthy adult",]))
Elm <- subset(Elm, cells = rownames(Elm@meta.data[Elm@meta.data$category=="Myeloid",]))

for (clus in unique(top50$cluster)){
  print(clus)
  clus_top50 <- top50[top50$cluster==clus,]$ENSG
  if (length(clus_top50)>10){
    nbin = length(clus_top50)} 
  else {nbin = 10}
  print(length(clus_top50))
  Elm <- AddModuleScore(Elm, features = list(clus_top50), name = clus, nbin = nbin)
}

HLAs <- rownames(MNP@assays$RNA@data)[startsWith(rownames(MNP@assays$RNA@data), "HLA")]
HLAs <- gconvert(HLAs,organism="hsapiens",target="ENSG",filter_na = F)$target
Elm <- AddModuleScore(Elm, features = list(HLAs), name = "HLAs", nbin = nbin)

spec_genes <- c("CD1C","CD207","PLAC8","PKIB","PPA1","SLC38A1","FCER1A","CD1E","LSP1","LTB")
spec_genes <- gconvert(spec_genes,organism="hsapiens",target="ENSG",filter_na = F)$target
Elm <- AddModuleScore(Elm, features = list(spec_genes), name = "cDC2sup", nbin = 10)


## Reintegrate data
# split the dataset into a list of two seurat objects (stim and CTRL)
Elm.list <- SplitObject(Elm, split.by = "Age_group")

# normalize and identify variable features for each dataset independently
Elm.list <- lapply(X = Elm.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = Elm.list)
Elm.list <- lapply(X = Elm.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)})

immune.anchors <- FindIntegrationAnchors(object.list = Elm.list, anchor.features = features, reduction = "rpca")
# this command creates an 'integrated' data assay
Elm <- IntegrateData(anchorset = immune.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(Elm) <- "integrated"

# Run the standard workflow for visualization 
Elm <- FindVariableFeatures(Elm, nfeatures = 2000)
Elm <- ScaleData(Elm)
Elm <- RunPCA(Elm)
ElbowPlot(Elm, ndims = 30)

Elm <- RunUMAP(Elm,reduction = "pca", dims = 1:20)
DimPlot(Elm, group.by = c("author_cell_type","Age_group"))

## Some cDC2 specific genes
VlnPlot(Elm, #subset(Elm_orig, cells = rownames(Elm_orig@meta.data[Elm_orig@meta.data$Diagnosis=="Healthy adult"])), 
        group.by = "author_cell_type", features = "cDC2sup1",
        ncol = 1,
        pt.size = 0)+NoLegend()

## MPO progenitor genes
MPOs <- gconvert(c("MPO","AZU1","ELANE"),organism="hsapiens",target="ENSG",filter_na = F)$target
FeaturePlot(subset(Elm, cells = rownames(Elm@meta.data[Elm@meta.data$Diagnosis=="Healthy adult",])), features = Bcell[4])+
  scale_colour_gradientn(colours=mycols)


Elm_orig <- Elm

Elm <- cbind(Elm_orig@meta.data, Elm_orig@reductions$umap@cell.embeddings)
colnames(Elm)[35:55] <- str_sub(colnames(Elm)[35:55], start = 1, end = -2)
colnames(Elm)[35:55] <- str_replace(colnames(Elm)[35:55], " ","")
colnames(Elm)[35:55] <- str_replace(colnames(Elm)[35:55], "/","")

for (clus in colnames(Elm)[35:55]){
  moduleplot <- ggplot(Elm, aes(x = UMAP_1, y = UMAP_2, colour = Elm[,clus]))+
    geom_point_rast()+
    scale_color_gradientn(colours = mycols, name = clus)+
    theme_classic()
  
  pdf(paste(dato,"ReclustMyeloid_Elmentaiteetal_DEGsbasedModulescores_",clus,".pdf",sep = ""),height = 5, width = 5)
  print(moduleplot)
  dev.off()
}
for (clus in colnames(Elm)[35:55]){
  moduleplot <- ggplot(Elm[Elm$Diagnosis=="Healthy adult",], aes(x = UMAP_1, y = UMAP_2, colour = Elm[Elm$Diagnosis=="Healthy adult",clus]))+
    geom_point_rast()+
    scale_color_gradientn(colours = mycols, name = clus)+
    theme_classic()
  
  pdf(paste(dato,"ReclustMyeloid_Elmentaiteetal_HEALTHYADULT_DEGsbasedModulescores_",clus,".pdf",sep = ""),height = 5, width = 5)
  print(moduleplot)
  dev.off()
}
for (clus in colnames(Elm_orig@meta.data)[35:55]){
  Vlnplot <- VlnPlot(Elm_orig, #subset(Elm_orig, cells = rownames(Elm_orig@meta.data[Elm_orig@meta.data$Diagnosis=="Healthy adult"])), 
                     group.by = "author_cell_type", features = clus, 
                     pt.size = 0)+NoLegend()

  pdf(paste(dato,"Vln_ReclustMyeloid_Elmentaiteetal_HEALTHYADULT_DEGsbasedModulescores_",clus,".pdf",sep = ""),height = 5, width = 5)
  print(Vlnplot)
  dev.off()
}


pdf(paste(dato,"Vln_ReclustMyeloid_Elmentaiteetal_semisupcDC2genes.pdf",sep = ""),height = 16, width = 12)
VlnPlot(Elm_orig, #subset(Elm_orig, cells = rownames(Elm_orig@meta.data[Elm_orig@meta.data$Diagnosis=="Healthy adult"])), 
        group.by = "author_cell_type", features = spec_genes,
        ncol = 4,
        pt.size = 0)+NoLegend()
dev.off()

Idents(Elm_orig) <- "author_cell_type"
DEGs_ELm <- FindAllMarkers(Elm_orig, only.pos = T, logfc.threshold = 0.5)
DEGs_Elm$gene_ID <- gconvert(DEGs_Elm$gene,organism="hsapiens",target="ENTREZGENE",filter_na = F,mthreshold = 1)$target



ggplot(Elm[Elm$Age_group=="Adult",], aes(x = UMAP_1, y = UMAP_2, colour = author_cell_type))+
  geom_point_rast()+
  theme_classic()
ggplot(Elm[Elm$Age_group=="Adult",], aes(x = UMAP_1, y = UMAP_2, colour = cDC2sup))+
  geom_point_rast()+
  scale_color_gradientn(colours = mycols, name = "cDC2sup")+
  theme_classic()

#### #####
Elm <- readRDS("/Users/linewulff/Documents/work/projects/FentonWulff_LP_MNP/Lit_data/ElmentaiteNature2021.rds")
# Elm <- subset(Elm, cells = rownames(Elm@meta.data[Elm@meta.data$Diagnosis=="Healthy adult",]))
head(Elm@meta.data)
DimPlot(Elm, group.by = "category", label = T)+NoLegend()

Bcell <- gconvert(c("CD19","CD79A","CD79B","IGHA1"),organism="hsapiens",target="ENSG",filter_na = F)$target
FeaturePlot(subset(Elm, cells = rownames(Elm@meta.data[Elm@meta.data$Diagnosis=="Healthy adult",])), features = Bcell[4])+
  scale_colour_gradientn(colours=mycols)
DimPlot(subset(Elm, cells = rownames(Elm@meta.data[Elm@meta.data$Diagnosis=="Healthy adult",])), 
        group.by = "category", label = T)+NoLegend()



for (clus in unique(top50$cluster)){
  print(clus)
  clus_top50 <- top50[top50$cluster==clus,]$ENSG
  if (length(clus_top50)>10){
    nbin = length(clus_top50)} 
  else {nbin = 10}
  print(length(clus_top50))
  Elm <- AddModuleScore(Elm, features = list(clus_top50), name = clus, nbin = nbin)
}

HLAs <- rownames(MNP@assays$RNA@data)[startsWith(rownames(MNP@assays$RNA@data), "HLA")]
HLAs <- gconvert(HLAs,organism="hsapiens",target="ENSG",filter_na = F)$target
Elm <- AddModuleScore(Elm, features = list(HLAs), name = "HLAs", nbin = nbin)

spec_genes <- c("CD1C","CD207","PLAC8","PKIB","PPA1","SLC38A1","FCER1A","CD1E","LSP1","LTB")
spec_genes <- gconvert(spec_genes,organism="hsapiens",target="ENSG",filter_na = F)$target
Elm <- AddModuleScore(Elm, features = list(spec_genes), name = "cDC2sup", nbin = 10)

spec_sub <- c("ILC3","Mast cell","NK T cell","CLC+ Mast cell","ILC2","CLP","MPO+ mono-neutrophil")     
Idents(Elm) <- 'author_cell_type'
Elm <- subset(Elm, idents = spec_sub)

for (clus in colnames(Elm@meta.data)[35:55]){
  Vlnplot <- VlnPlot(subset(Elm, cells = rownames(Elm@meta.data[Elm@meta.data$Diagnosis=="Healthy adult",])), 
                     group.by = "author_cell_type", features = clus, 
                     pt.size = 0)+NoLegend()
  if (str_detect(clus,'/')){clus <- str_replace(clus,"/",".")}
  pdf(paste(dato,"Vln_ReclustPosCont_Elmentaiteetal_HEALTHYADULT_DEGsbasedModulescores_",clus,".pdf",sep = ""),height = 5, width = 5)
  print(Vlnplot)
  dev.off()
}

Elm <- cbind(Elm@meta.data, Elm@reductions$umap@cell.embeddings)
colnames(Elm)[35:55] <- str_sub(colnames(Elm)[35:55], start = 1, end = -2)
colnames(Elm)[35:55] <- str_replace(colnames(Elm)[35:55], " ","")
colnames(Elm)[35:55] <- str_replace(colnames(Elm)[35:55], "/","")

ggplot(Elm[Elm$Diagnosis=="Healthy adult",], aes(x=UMAP_1, y=UMAP_2))+
  geom_point_rast(colour="lightgrey")+
  geom_point_rast(data=Elm[Elm$Diagnosis=="Healthy adult" & Elm$author_cell_type=="Mast cell",], 
                  aes(x=UMAP_1, y=UMAP_2), colour="red")+
  theme_classic()


#### Finalizing the gene sets ####
sign_list <- list()
for (clus in unique(top50$cluster)){
  sign_list[[clus]] <- top50[top50$cluster==clus,]$gene
}
sign_list[["cDC2"]] <- c("CD1C","CD207","PLAC8","PKIB","PPA1","SLC38A1","FCER1A","CD1E","LSP1","LTB")
sign_list <- vectorlist(sign_list)
write.csv(sign_list,"2312_FentonWulff_MNPsignatures.csv")
