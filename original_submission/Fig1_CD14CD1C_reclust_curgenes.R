###### Libraries ######
library(gplots)
library(ggplot2)
library(dplyr)
library(stringr)
library(scales)
library(Rmagic)
library(ggplot2)
library(viridis)
library(clustree)
library(ccRemover)
library(rgl)
library(Seurat , lib.loc = "/Library/Frameworks/R.framework/Versions/3.5/Resources/other_libs")
library(scales)
library(rgl)
library(ggbiplot) #R version 3.5
library(chemometrics)
library(ChemometricsWithR)

### Clear and read in data ###
rm(list=ls())

########### After supercomputer part #####
#### data read in and variables ####
MNP_LP <-readRDS("/Volumes/Mucosal-Immunology/WA group/Tom and Line data/cLP_SILP_merged/R6/Rdata/R3_CD14CD1Cprol_MNP_LP.rds")
MNP_LP_imp <- readRDS("/Volumes/Mucosal-Immunology/WA group/Tom and Line data/cLP_SILP_merged/R6/Rdata/R3_CD14CD1Cprol_MNP_LP_imputed.rds")
MNP_LP_imp@meta.data$integrated_snn_res.2.8 <- MNP_LP@meta.data$integrated_snn_res.2.8

MNP_LP_sc_pca <- MNP_LP@reductions$pca@feature.loadings
print(MNP_LP[["pca"]], dims = 1:10, nfeatures = 30)
pdf("PCA_loadings.pdf", height = 80, width = 30)
VizDimLoadings(MNP_LP, dims = 1:20, reduction = "pca", nfeatures = 100)
dev.off()

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

MNP_LP@meta.data$integrated_snn_res.2.8 <- as.character(MNP_LP@meta.data$integrated_snn_res.2.8)
MNP_LP@meta.data$integrated_snn_res.2.8 <- factor(MNP_LP@meta.data$integrated_snn_res.2.8, levels=ord)
Idents(MNP_LP) <- as.character(MNP_LP@meta.data$integrated_snn_res.2.8)

col_pal <- hue_pal()(39)
names(col_pal) <- ord

MNP_LP@meta.data$lineage <- NA
for (cl in levels(MNP_LP@meta.data$integrated_snn_res.2.8)){
  if (cl %in% levels(MNP_LP@meta.data$integrated_snn_res.2.8)[1:3]){
    MNP_LP@meta.data[MNP_LP@meta.data$integrated_snn_res.2.8==cl,]$lineage <- "mono"}
  else if (cl %in% levels(MNP_LP@meta.data$integrated_snn_res.2.8)[4:20]){
    MNP_LP@meta.data[MNP_LP@meta.data$integrated_snn_res.2.8==cl,]$lineage <- "mac"}
  else {MNP_LP@meta.data[MNP_LP@meta.data$integrated_snn_res.2.8==cl,]$lineage <- "cDC2-like"}
}
##### Curated_lineage genes #####
#cur_lin_genes <- c('HEY1','KLF7','MAFB','BHLHE41','NR1H3','MAF','HIC1','IRF4','ZEB1','PLD3','CTSD','MERTK','ABCA1','AP1S3','BTLA','ADAM19','SPINT2','P2RY10','FLT3','DPP4','KIT','SEPT6','PLCB1','AQP9','SLC2A6','DNAAF1','TNFAIP6','MS4A7','PSAP','FTL','CTSB','FOLR2','CSF1R','SLC40A1','MAF','DAB2','SLCO2B1','PDK4','STAB1','GADD45G','CTSD','FUCA1','APOE','TMEM176B','LGALS3','PLA2G7','TMEM176A','PLD3','CD207','CST7','DSG2','PKIB','PPP1R14A','PPA1','CCL22','SLC38A1','LTB','DUSP5','FAM110A','HIC1','GNG7','TRAF5','G3BP2','CD52','ENHO','AIM1','SPIB','TUBA1A','NRARP','TCTN3','CST7','UPF2','ADAM8','SDPR')
  
# unique(c('HEY1','MYC','KLF7','MAFB','BHLHE41','NR1H3','MAF','CEBPA','HIC1','IRF4','ZEB1','PLD3','CTSD','MERTK','TCN2','SEPP1',
#                           'ABCA1','AP1S3','BTLA','ADAM19','SPINT2','P2RY10','FLT3','DPP4','KIT','SEPT6',#'S100A9','S100A8','FCN1','VCAN','SERPINA1','APOBEC3A','SLC11A1','PLCB1','AQP9','SLC2A6','S100A12','DNAAF1','THBS1','TNFAIP6','CLEC4E','SOD2','IER3','MS4A7','PSAP','FTL',
#                           'CTSB','FOLR2','CSF1R','SEPP1','SLC40A1','MAF','DAB2','SLCO2B1','PDK4','STAB1','GADD45G','CTSD','SEPP1','FUCA1',
#                           'APOE','TMEM176B','LGALS3','PLA2G7','PRDX1','CTSC','TMEM176A','PLD3','SERPINF1','CD207','CST7','DSG2','PKIB','PPP1R14A',
#                           'PPA1','CCL22','SLC38A1','LTB','DUSP5','FAM110A','HIC1','RALA','GNG7','TRAF5','NAP1L1','G3BP2','CD52','ENHO','AIM1',
#                           'SPIB','TUBA1A','NRARP','TCTN3','CST7','UPF2','ADAM8','SDPR'))
cur_lin_genes <-unique(c('S100A9','S100A8','FCN1','VCAN','SERPINA1','APOBEC3A','SLC11A1','S100A6','S100A4','FPR1','PLCB1','AQP9','SLC2A6','S100A12','DNAAF1','THBS1','TNFAIP6','CLEC4E','SOD2','CXCL8','IL1B','PSAP','FTL','MS4A7','IER3','FOLR2','CTSB','KLF2','SEPP1','C1QA','C1QB','SLC40A1','MAF','DAB2','SLCO2B1','PDK4','C1QC','STAB1','GADD45G','CTSD','C1QB','SEPP1','FUCA1','APOE','C1QC','C1QA','TMEM176B','LGALS3','PLA2G7','PRDX1','CTSC','TMEM176A','PLD3','SERPINF1','DUSP1','CD207','CST7','DSG2','PKIB','PPP1R14A','PPA1','CCL22','SLC38A1','LTB','DUSP5','FAM110A','HIC1','RALA','GNG7','HLA-DQB1','PLAC8','TRAF5','NAP1L1','G3BP2','CD52','ENHO','AIM1','ICAM3','CREM','SPIB','TUBA1A','NRARP','TCTN3','CST7','ARL4C','EZR','UPF2','ADAM8','IL1R2','CD1C','VDR','IFITM2','TRMT6','SDPR','C1orf162','XCR1','CLEC9A','S100B','SOX4','CCR7','LAMP3','CD40', 'ZBTB16','HEY1','MYC','KLF7','NFKBIA','NFKBIZ','MAFB','BHLHE41','NR1H3','MAF','CEBPA','ID3','NR2F2','BATF3','GTF2B','HIC1','IRF4','ZEB1','ZNF165', 'FLT3'))
length(cur_lin_genes)
length(cur_lin_genes[cur_lin_genes %in% rownames(MNP_LP@assays$integrated)]) # all there

#length(cur_lin_genes[cur_lin_genes %in% rownames(MNP_LP_imp@assays$imputated)]) # all there
#### Segura gene signatures for curated list of genes ####
segura_sign <- read.csv("/Users/lwc/Google Drive/PhD/Literature/Immunology/MNP/Nat_Huang_Segura_2018_sheet1.csv")
colnames(segura_sign)
segura_sign <- segura_sign[2:431,c(1,3,9,11)]
bloodCD14 <- as.character(segura_sign$blood.CD14..monocytes)[as.character(segura_sign$blood.CD14..monocytes)!=""]
bloodDC2 <- as.character(segura_sign$blood.cDC2)[as.character(segura_sign$blood.cDC2)!=""]
invitroMoMac <- as.character(segura_sign$in.vitro.mo.Mac)[as.character(segura_sign$in.vitro.mo.Mac)!=""]
tissueDC2 <- as.character(segura_sign$tissue.cDC2)[as.character(segura_sign$tissue.cDC2)!=""]

segura_sign <- c(bloodCD14, bloodDC2, invitroMoMac)#, tissueDC2)
segura_sign[!segura_sign %in% rownames(MNP_LP@assays$integrated@data)]
segura_sign <- unique(segura_sign[segura_sign %in% rownames(MNP_LP@assays$integrated@data)])
  
cur_lin_genes <- segura_sign

#### Av exp pca ####
Idents(MNP_LP) <- MNP_LP@meta.data$integrated_snn_res.2.8
MNP_LP_av <- AverageExpression(MNP_LP, features = cur_lin_genes, assays = "integrated", return.seurat = T)
MNP_LP_av_pc <- t(MNP_LP_av$integrated@data)

leaveouts <- c()
for (cluster in rownames(MNP_LP_av_pc)){
  navs <- colnames(MNP_LP_av_pc)[which(is.na(MNP_LP_av_pc[cluster,]))]
  if (length(navs)>0){
    print(paste(cluster,navs))
    leaveouts <- append(leaveouts,navs)
  }}
saves <- colnames(MNP_LP_av_pc)[!colnames(MNP_LP_av_pc) %in% leaveouts]
MNP_LP_av_pc <- MNP_LP_av_pc[,saves]
dim(MNP_LP_av_pc) #162 genes

MNP_pca <- prcomp(MNP_LP_av_pc, center = TRUE,scale = TRUE)
MNP_df <- as.data.frame(MNP_pca$x)
MNP_PCA_2 <- PCA(scale(MNP_LP_av_pc))

MNP_df$lineage <- NA
MNP_df[rownames(MNP_df) %in% c('12','24','22'),]$lineage <- "monocytes"
MNP_df[rownames(MNP_df) %in% c('4','33','16','19'),]$lineage <- "monocytes" #"mono-ints"
MNP_df[rownames(MNP_df) %in% c('18','35','17','28','0','10','15','2','5','29','21','3','13'),]$lineage <- "macrophages"
MNP_df[rownames(MNP_df) %in% c('37','34','36'),]$lineage <- "DC2-like"#"DC-precursors"
MNP_df[rownames(MNP_df) %in% c('38'),]$lineage <- "DC2-like"#"mig-DC2"
MNP_df[rownames(MNP_df) %in% c('26','7','8','30','32','27'),]$lineage <- "DC2-like"#"DC2A"
MNP_df[rownames(MNP_df) %in% c('11','20','9','1','14','6','31','25','23'),]$lineage <- "DC2-like"#"DC2B/3"

dev.off()
pdf(paste(dato,project,"Elbow_CurGenes.pdf",sep="_"),height = 5,width = 5)
plot(1:length(MNP_PCA_2$var), MNP_PCA_2$var/MNP_PCA_2$totalvar, cex=2, ylab = "% Variance explained", xlab = "PCs")
lines(1:length(MNP_PCA_2$var), MNP_PCA_2$var/MNP_PCA_2$totalvar)
dev.off()

# ggplot(MNP_df, aes(x=PC1,y=PC2,colour=rownames(MNP_df)))+geom_point(alpha=0)+scale_color_manual(values=col_pal)+
#   geom_text(label=rownames(MNP_df),size=6)+theme_classic()+theme(legend.position = "none")
# dev.off()
# 
# ggbiplot(MNP_pca, obs.scale = 1, var.scale = 1, alpha=0 ,choices = c(1,2), ellipse = F, circle = FALSE, var.axes = F)+theme_classic()+
#   geom_text(label=rownames(MNP_df), aes(color=rownames(MNP_df)))+
#   theme(legend.position = "none")+scale_color_manual(values = col_pal)
# 
# ggbiplot(MNP_pca, obs.scale = 1, var.scale = 1, alpha=0 ,choices = c(1,2), ellipse = T, circle = FALSE, var.axes = F, groups = MNP_df$lineage)+theme_classic()+
#   geom_text(label=rownames(MNP_df))+
#   theme(legend.position = "none")+scale_color_manual(values = col_pal)
# 
# col_pal_ex <- append(col_pal, c("lightslateblue","forestgreen","darkorange2"))
# names(col_pal_ex) <- append(names(col_pal),c("DC2-like","macrophages","monocytes"))
# ggbiplot(MNP_pca, obs.scale = 1, var.scale = 1, groups = MNP_df$lineage, alpha=0 ,choices = c(1,2), ellipse = T, circle = FALSE, var.axes = F)+
#   geom_text(label=rownames(MNP_df), aes(color=rownames(MNP_df)))+
#   theme_classic()+scale_color_manual(values=col_pal_ex)+
#   guides(color=FALSE)

col_pal_ex <- append(col_pal, c("lightslateblue","forestgreen","darkorange2"))
names(col_pal_ex) <- append(names(col_pal),c("DC2-like","macrophages","monocytes"))
pdf(paste(dato,project,"pca_pseudobulk_curated_gene_list.pdf",sep="_"),height = 6, width = 6)
plot1 <- ggbiplot(MNP_pca, obs.scale = 1, var.scale = 1, groups = MNP_df$lineage,alpha=0 ,choices = c(1,2), ellipse = T, circle = FALSE, var.axes = F)+
  geom_text(label=rownames(MNP_df), aes(color=rownames(MNP_df)), show.legend  = F)+
  theme_classic()+theme(legend.title = element_text(face="bold"))
plot1+ scale_colour_manual(values = col_pal_ex, limits = c("DC2-like","macrophages","monocytes"), name="Cell types")+
  guides(colour = guide_legend(override.aes = list(shape=16)))#+geom_vline(xintercept = -0.3)
dev.off()


#############################################################################
## Av exp pca RANDOM LIST OF GENES
Idents(MNP_LP) <- MNP_LP@meta.data$integrated_snn_res.2.8

non_zeros <- rownames(MNP_LP@assays$integrated@data[complete.cases(as.matrix(MNP_LP_imp@assays$integrated@data)),])
set.seed(8)
rand_genes <- sample(non_zeros,120)
MNP_LP_av <- AverageExpression(MNP_LP_imp, features = rand_genes, assays = "integrated", return.seurat = T)
MNP_LP_av_pc <- t(as.matrix(t(MNP_LP_av$integrated@data)))
MNP_LP_av_pc <- t(na.omit(MNP_LP_av_pc))
dim(MNP_LP_av_pc)
rand_genes_used <- colnames(MNP_LP_av_pc)
write.csv(rand_genes_used, file="rand_genes.csv")

MNP_pca <- prcomp(MNP_LP_av_pc, center = TRUE,scale = TRUE)
MNP_df <- as.data.frame(MNP_pca$x)
MNP_PCA_2 <- PCA(scale(MNP_LP_av_pc))

MNP_df$lineage <- NA
MNP_df[rownames(MNP_df) %in% c('12','24','22'),]$lineage <- "monocytes"
MNP_df[rownames(MNP_df) %in% c('4','33','16','19'),]$lineage <- "monocytes" #"mono-ints"
MNP_df[rownames(MNP_df) %in% c('18','35','17','28','0','10','15','2','5','29','21','3','13'),]$lineage <- "macrophages"
MNP_df[rownames(MNP_df) %in% c('37','34','36'),]$lineage <- "DC2-like"#"DC-precursors"
MNP_df[rownames(MNP_df) %in% c('38'),]$lineage <- "DC2-like"#"mig-DC2"
MNP_df[rownames(MNP_df) %in% c('26','7','8','30','32','27'),]$lineage <- "DC2-like"#"DC2A"
MNP_df[rownames(MNP_df) %in% c('11','20','9','1','14','6','31','25','23'),]$lineage <- "DC2-like"#"DC2B/3"

dev.off()
pdf(paste(dato,project,"Elbow_RandGenes.pdf",sep="_"),height = 5,width = 5)
plot(1:length(MNP_PCA_2$var), MNP_PCA_2$var/MNP_PCA_2$totalvar, cex=2, ylab = "% Variance explained", xlab = "PCs")
lines(1:length(MNP_PCA_2$var), MNP_PCA_2$var/MNP_PCA_2$totalvar)
dev.off()

col_pal_ex <- append(col_pal, c("lightslateblue","forestgreen","darkorange2"))
names(col_pal_ex) <- append(names(col_pal),c("DC2-like","macrophages","monocytes"))
pdf(paste(dato,project,"pca_pseudobulk_RANDOM_gene_list.pdf",sep="_"),height = 6, width = 6)
plot1 <- ggbiplot(MNP_pca, obs.scale = 1, var.scale = 1, groups = MNP_df$lineage,alpha=0 ,choices = c(1,2), ellipse = T, circle = FALSE, var.axes = F)+
  geom_text(label=rownames(MNP_df), aes(color=rownames(MNP_df)), show.legend  = F)+
  theme_classic()+theme(legend.title = element_text(face="bold"))
plot1+ scale_colour_manual(values = col_pal_ex, limits = c("DC2-like","macrophages","monocytes"), name="Cell types")+
  guides(colour = guide_legend(override.aes = list(shape=16)))
dev.off()

## ALl variable genes (non NA values) - NON-imputated
Idents(MNP_LP) <- MNP_LP@meta.data$integrated_snn_res.2.8

non_zeros <- rownames(MNP_LP@assays$integrated@data[complete.cases(as.matrix(MNP_LP_imp@assays$integrated@data)),])
MNP_LP_av <- AverageExpression(MNP_LP_imp, assays = "integrated", return.seurat = T)
MNP_LP_av_pc <- t(as.matrix(t(MNP_LP_av$integrated@data)))
MNP_LP_av_pc <- t(na.omit(MNP_LP_av_pc))
dim(MNP_LP_av_pc)

MNP_pca <- prcomp(MNP_LP_av_pc, center = TRUE,scale = TRUE)
MNP_df <- as.data.frame(MNP_pca$x)
MNP_PCA_2 <- PCA(scale(MNP_LP_av_pc))

MNP_df$lineage <- NA
MNP_df[rownames(MNP_df) %in% c('12','24','22'),]$lineage <- "monocytes"
MNP_df[rownames(MNP_df) %in% c('4','33','16','19'),]$lineage <- "monocytes" #"mono-ints"
MNP_df[rownames(MNP_df) %in% c('18','35','17','28','0','10','15','2','5','29','21','3','13'),]$lineage <- "macrophages"
MNP_df[rownames(MNP_df) %in% c('37','34','36'),]$lineage <- "DC2-like"#"DC-precursors"
MNP_df[rownames(MNP_df) %in% c('38'),]$lineage <- "DC2-like"#"mig-DC2"
MNP_df[rownames(MNP_df) %in% c('26','7','8','30','32','27'),]$lineage <- "DC2-like"#"DC2A"
MNP_df[rownames(MNP_df) %in% c('11','20','9','1','14','6','31','25','23'),]$lineage <- "DC2-like"#"DC2B/3"

dev.off()
pdf(paste(dato,project,"Elbow_AllvarGenes.pdf",sep="_"),height = 5,width = 5)
plot(1:length(MNP_PCA_2$var), MNP_PCA_2$var/MNP_PCA_2$totalvar, cex=2, ylab = "% Variance explained", xlab = "PCs")
lines(1:length(MNP_PCA_2$var), MNP_PCA_2$var/MNP_PCA_2$totalvar)
dev.off()

col_pal_ex <- append(col_pal, c("lightslateblue","forestgreen","darkorange2"))
names(col_pal_ex) <- append(names(col_pal),c("DC2-like","macrophages","monocytes"))
pdf(paste(dato,project,"pca_pseudobulk_AllVarGenes_list.pdf",sep="_"),height = 6, width = 6)
plot1 <- ggbiplot(MNP_pca, obs.scale = 1, var.scale = 1, groups = MNP_df$lineage,alpha=0 ,choices = c(1,2), ellipse = T, circle = FALSE, var.axes = F)+
  geom_text(label=rownames(MNP_df), aes(color=rownames(MNP_df)), show.legend  = F)+
  theme_classic()+theme(legend.title = element_text(face="bold"))
plot1+ scale_colour_manual(values = col_pal_ex, limits = c("DC2-like","macrophages","monocytes"), name="Cell types")+
  guides(colour = guide_legend(override.aes = list(shape=16)))
dev.off()

#### MNP LP UMAP to make better sense with different metrics for umap calculation ####
library(umap) #R 3.5

umap.conf <- umap.defaults
umap.conf$metric <- 'pearson2'
set.seed(42)
umap.ts <- umap::umap(MNP_LP@reductions$pca@cell.embeddings[,c(1:3,5,7:15)], config = umap.conf) #avoid PC4 and 6
colnames(umap.ts$layout) <- c('UMAP_1', 'UMAP_2')

MNP_LP@reductions$umap@cell.embeddings <- umap.ts$layout
DimPlot(MNP_LP, reduction = "umap",group.by = "integrated_snn_res.2.8", label = T, label.size = 6)

saveRDS(MNP_LP, file="/Volumes/Mucosal-Immunology/WA group/Tom and Line data/cLP_SILP_merged/R6/Rdata/R3_CD14CD1Cprol_MNP_LP.rds")
FeaturePlot(MNP_LP, features = "cite_CD14", reduction = "umap", cols = mycols_b, cells = rownames(MNP_LP@meta.data[MNP_LP@meta.data$cite=="YES",]))

#####################
## Rerun PCA and UMAP with input genes cur_lin_genes
# MNP_LP_cp <- RunPCA(MNP_LP_cp, npcs = 30, verbose = FALSE, features=cur_lin_genes) #
# 
# ElbowPlot(MNP_LP_cp, ndims = 30)
# dims=10
# 
# MNP_LP_imp <- RunUMAP(MNP_LP_imp, dims = c(1:dims), reduction = "pca",)
# 
# DimPlot(MNP_LP_imp, group.by = "integrated_snn_res.2.8", reduction = "umap", label = T, pt.size = 0.3, label.size = 6,cols = col_pal)
# FeatureScatter(MNP_LP_cp, feature1 = "PC_1", feature2 = "PC_2", cols = col_pal, pt.size =0.3)
# FeatureScatter(MNP_LP, feature1 = "PC_1", feature2 = "PC_5", cols = col_pal, pt.size =0.3)#, cells = save_cells)
# 
# 
# MNP_df <- MNP_LP@reductions$umap@cell.embeddings
# MNP_df <- cbind(MNP_df, MNP_LP@meta.data)
# MNP_df <- cbind(MNP_df, MNP_LP@reductions$pca@cell.embeddings)
# write.csv(MNP_df,file="/Volumes/Mucosal-immunology/WA group/Tom and Line data/cLP_SILP_merged/R6/DimReductions/MNP_LP_6D.csv")
# 
# save_cells <- CellSelector(FeatureScatter(MNP_LP_imp, feature1 = "PC_1", feature2 = "PC_2", cols = col_pal, pt.size =0.3))
# other_cells <- colnames(MNP_LP)[!colnames(MNP_LP) %in% save_cells]
# 
# DimPlot(MNP_LP, group.by = "integrated_snn_res.2.8", reduction = "umap", label = T, pt.size = 0.3, label.size = 6,cols = col_pal, cells = save_cells)
# FeatureScatter(MNP_LP_imp, feature1 = "PC_5", feature2 = "PC_4", cols = col_pal, pt.size =0.3, cells=save_cells)
# FeatureScatter(MNP_LP_imp, feature1 = "PC_1", feature2 = "PC_2", cols = col_pal, pt.size =0.3)

#### F1 update to fit HM (renaming clusters ...) 22.07.04 ####
ord <- c(12,24,22,4,33,16,19,18,35,17,28,0,10,15,2,5,29,21,3,13,37,34,36,38,26,7,8,30,32,27,11,20,9,1,14,6,31,25,23)
ord <- as.character(ord)

MNP_LP@meta.data$integrated_snn_res.2.8 <- as.character(MNP_LP@meta.data$integrated_snn_res.2.8)
MNP_LP@meta.data$integrated_snn_res.2.8 <- factor(MNP_LP@meta.data$integrated_snn_res.2.8, levels=ord)
Idents(MNP_LP) <- as.character(MNP_LP@meta.data$integrated_snn_res.2.8)

col_pal <- hue_pal()(39)
names(col_pal) <- ord

MNP_LP@meta.data$lineage <- NA
for (cl in levels(MNP_LP@meta.data$integrated_snn_res.2.8)){
  if (cl %in% levels(MNP_LP@meta.data$integrated_snn_res.2.8)[1:3]){
    MNP_LP@meta.data[MNP_LP@meta.data$integrated_snn_res.2.8==cl,]$lineage <- "mono"}
  else if (cl %in% levels(MNP_LP@meta.data$integrated_snn_res.2.8)[4:20]){
    MNP_LP@meta.data[MNP_LP@meta.data$integrated_snn_res.2.8==cl,]$lineage <- "mac"}
  else {MNP_LP@meta.data[MNP_LP@meta.data$integrated_snn_res.2.8==cl,]$lineage <- "cDC2-like"}
}

## gene lists 
segura_sign <- read.csv("/Users/lwc/Google Drive/PhD/Literature/Immunology/MNP/Nat_Huang_Segura_2018_sheet1.csv")
colnames(segura_sign)
segura_sign <- segura_sign[2:431,c(1,3,9,11)]
bloodCD14 <- as.character(segura_sign$blood.CD14..monocytes)[as.character(segura_sign$blood.CD14..monocytes)!=""]
bloodDC2 <- as.character(segura_sign$blood.cDC2)[as.character(segura_sign$blood.cDC2)!=""]
invitroMoMac <- as.character(segura_sign$in.vitro.mo.Mac)[as.character(segura_sign$in.vitro.mo.Mac)!=""]
tissueDC2 <- as.character(segura_sign$tissue.cDC2)[as.character(segura_sign$tissue.cDC2)!=""]

segura_sign <- c(bloodCD14, bloodDC2, invitroMoMac)#, tissueDC2)
segura_sign[!segura_sign %in% rownames(MNP_LP@assays$integrated@data)]
segura_sign <- unique(segura_sign[segura_sign %in% rownames(MNP_LP@assays$integrated@data)])

cur_lin_genes <- segura_sign

## av pca
Idents(MNP_LP) <- MNP_LP@meta.data$integrated_snn_res.2.8
MNP_LP_av <- AverageExpression(MNP_LP, features = cur_lin_genes, assays = "integrated", return.seurat = T)
MNP_LP_av_pc <- t(MNP_LP_av$integrated@data)

leaveouts <- c()
for (cluster in rownames(MNP_LP_av_pc)){
  navs <- colnames(MNP_LP_av_pc)[which(is.na(MNP_LP_av_pc[cluster,]))]
  if (length(navs)>0){
    print(paste(cluster,navs))
    leaveouts <- append(leaveouts,navs)
  }}
saves <- colnames(MNP_LP_av_pc)[!colnames(MNP_LP_av_pc) %in% leaveouts]
MNP_LP_av_pc <- MNP_LP_av_pc[,saves]
dim(MNP_LP_av_pc) #162 genes

MNP_pca <- prcomp(MNP_LP_av_pc, center = TRUE,scale = TRUE)
MNP_df <- as.data.frame(MNP_pca$x)
MNP_PCA_2 <- PCA(scale(MNP_LP_av_pc))

MNP_df$lineage <- NA
MNP_df[rownames(MNP_df) %in% c('12','24','22'),]$lineage <- "monocytes"
MNP_df[rownames(MNP_df) %in% c('4','33','16','19'),]$lineage <- "monocytes" #"mono-ints"
MNP_df[rownames(MNP_df) %in% c('18','35','17','28','0','10','15','2','5','29','21','3','13'),]$lineage <- "macrophages"
MNP_df[rownames(MNP_df) %in% c('37','34','36'),]$lineage <- "DC2-like"#"DC-precursors"
MNP_df[rownames(MNP_df) %in% c('38'),]$lineage <- "DC2-like"#"mig-DC2"
MNP_df[rownames(MNP_df) %in% c('26','7','8','30','32','27'),]$lineage <- "DC2-like"#"DC2A"
MNP_df[rownames(MNP_df) %in% c('11','20','9','1','14','6','31','25','23'),]$lineage <- "DC2-like"#"DC2B/3"

dev.off()
pdf(paste(dato,project,"Elbow_CurGenes.pdf",sep="_"),height = 5,width = 5)
plot(1:length(MNP_PCA_2$var), MNP_PCA_2$var/MNP_PCA_2$totalvar, cex=2, ylab = "% Variance explained", xlab = "PCs")
lines(1:length(MNP_PCA_2$var), MNP_PCA_2$var/MNP_PCA_2$totalvar)
dev.off()

col_pal_ex <- append(col_pal, c("lightslateblue","forestgreen","darkorange2"))
names(col_pal_ex) <- append(names(col_pal),c("DC2-like","macrophages","monocytes"))
# pdf(paste(dato,project,"pca_pseudobulk_curated_gene_list.pdf",sep="_"),height = 6, width = 6)
# plot1 <- ggbiplot(MNP_pca, obs.scale = 1, var.scale = 1, groups = MNP_df$lineage,alpha=0 ,choices = c(1,2), ellipse = T, circle = FALSE, var.axes = F)+
#   geom_text(label=rownames(MNP_df), aes(color=rownames(MNP_df)), show.legend  = F)+
#   theme_classic()+theme(legend.title = element_text(face="bold"))
# plot1+ scale_colour_manual(values = col_pal_ex, limits = c("DC2-like","macrophages","monocytes"), name="Cell types")+
#   guides(colour = guide_legend(override.aes = list(shape=16)))#+geom_vline(xintercept = -0.3)
# dev.off()

rownames(MNP_pca$x) <- seq(1,39,1)
col_pal_ex_new <- col_pal_ex
colnames(col_pal_ex_new)[1:39] <- seq(1,39,1)

pdf(paste(dato,project,"pca_pseudobulk_curated_gene_list_reord.pdf",sep="_"),height = 6, width = 6)
plot1 <- ggbiplot(MNP_pca, obs.scale = 1, var.scale = 1, groups = MNP_df$lineage,alpha=0, choices = c(1,2), ellipse = T, circle = FALSE, var.axes = F)+
  geom_text(label=rownames(MNP_pca$x), aes(color=rownames(MNP_df)), show.legend  = F)+
  theme_classic()+theme(legend.title = element_text(face="bold"))
plot1+ scale_colour_manual(values = col_pal_ex_new, limits = c("DC2-like","macrophages","monocytes"), name="Cell types")+
  guides(colour = guide_legend(override.aes = list(shape=16)))#+geom_vline(xintercept = -0.3)
dev.off()

