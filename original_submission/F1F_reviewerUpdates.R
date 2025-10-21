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

### Clear and read in data ###
rm(list=ls())

#### data read in and variables ####
MNP <- readRDS("/Volumes/LWulffExD/Projects/FentonWulffData/R6/Rdata_v2/R3_CD14CD1Cprol_MNPLP.rds")
MNP <- UpdateSeuratObject(MNP)

### variables
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)
# colour string for imputation and overlays
mycols_b <- c("#bdbdbd","#d9d9d9","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F")
# project
project <- "MNP_LP"
#trajectory subset
round="CD14CD1C_prol"

ord <- c(12,24,22,4,33,16,19,18,35,17,28,0,10,15,2,5,29,21,3,13,37,34,36,38,26,7,8,30,32,27,11,20,9,1,14,6,31,25,23)
ord <- as.character(ord)
names(ord) <- c(1:length(ord))

MNP@meta.data$integrated_snn_res.2.8 <- as.character(MNP@meta.data$integrated_snn_res.2.8)
MNP@meta.data$res2.8_ord <- NA
for (clus in unique(MNP@meta.data$integrated_snn_res.2.8)){
  MNP@meta.data[MNP@meta.data$integrated_snn_res.2.8==clus,]$res2.8_ord <- names(ord)[which(ord==clus)]
}
Idents(MNP) <- 'res2.8_ord'

### 
MNP_av <- AverageExpression(MNP, assays = "integrated", return.seurat = T, features = VariableFeatures(MNP))
MNP_av@meta.data

hc <- hclust(dist(t(MNP_av@assays$integrated@scale.data)))
p <- ggdendrogram(hc, rotate = FALSE, size = 2)

ggplotly(p)


## By Segura gene signatures
segura_sign <- read.csv("/Users/linewulff/Documents/work/projects/FentonWulff_LP_MNP/Lit_data/Nat_Huang_Segura_2018_sheet1.csv")
colnames(segura_sign)
segura_sign <- segura_sign[2:431,c(1,3,9,11)]
bloodCD14 <- as.character(segura_sign$blood.CD14..monocytes)[as.character(segura_sign$blood.CD14..monocytes)!=""]
bloodDC2 <- as.character(segura_sign$blood.cDC2)[as.character(segura_sign$blood.cDC2)!=""]
invitroMoMac <- as.character(segura_sign$in.vitro.mo.Mac)[as.character(segura_sign$in.vitro.mo.Mac)!=""]
tissueDC2 <- as.character(segura_sign$tissue.cDC2)[as.character(segura_sign$tissue.cDC2)!=""]

segura_sign <- c(bloodCD14, bloodDC2, invitroMoMac, tissueDC2)
segura_sign[!segura_sign %in% rownames(MNP@assays$integrated@data)]
segura_sign <- unique(segura_sign[segura_sign %in% rownames(MNP@assays$integrated@data)])

cur_lin_genes <- segura_sign

MNP_av <- AverageExpression(MNP, assays = "integrated", return.seurat = T, features = cur_lin_genes)
MNP_av@meta.data

hc <- hclust(dist(t(MNP_av@assays$integrated@scale.data)))
p <- ggdendrogram(hc, rotate = FALSE, size = 2)

ggplotly(p)

## Mulder et al
Mulder_sign <- read.csv("/Users/linewulff/Documents/work/projects/FentonWulff_LP_MNP/Lit_data/Mulder2021Imm_MNPSign.csv")
colnames(Mulder_sign)
Macro <- Mulder_sign$Macro
cMo <- Mulder_sign$cMo
DC2.3 <- Mulder_sign$DC2.DC3
mReg <- Mulder_sign$mregDC

cur_lin_genes <- c(Macro,cMo,DC2.3,mReg)

MNP_av <- AverageExpression(MNP, assays = "integrated", return.seurat = T, features = cur_lin_genes)
MNP_av@meta.data

hc <- hclust(dist(t(MNP_av@assays$integrated@scale.data)))
p <- ggdendrogram(hc, rotate = FALSE, size = 2)

pdf(paste(dato,"Dendrogram_Pseudobulk_res2.8_MulderGenes.pdf",sep = "_"), height = 5, width = 7)
p
dev.off()

ggplotly(p)

Idents(MNP) <- factor(MNP@meta.data$res2.8_ord, levels = seq(1,39))
VlnPlot(subset(MNP, idents = c("1","4","5","6","7","36","37","38","39")), features = c("CD1C","FCER1A","CD14","C5AR1","C1QA"), pt.size = 0, ncol=2)

#### Mulder et al correlation ####
obj2 <- readRDS("/Users/linewulff/Documents/work/projects/FentonWulff_LP_MNP/Lit_data/2021_MNP_Verse.RDS")
Idents(obj2) <- "Clusters"
DefaultAssay(obj2) <- "integrated"
DimPlot(obj2, label = T)

MNPcor <- cluster_corr(MNP, "res2.8_ord", "integrated",MNPverse,"Clusters","integrated")

png(paste(dato,project,))
heatmap.2(as.matrix(corr_mat), scale = "none", col =  bluered(150), trace="none", density ="none",
          dendrogram = "column",
          margins = c(8,5))
dev.off()


