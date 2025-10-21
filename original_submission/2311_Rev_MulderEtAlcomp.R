



library(gplots)
library(ggplot2)
library(dplyr)
library(stringr)
library(scales)
library(viridis)
library(Seurat)

#### Read in data ####
GSE_count <- read.table("/Users/linewulff/Documents/work/projects/FentonWulff_LP_MNP/Lit_data/GSE178209_SMARTSeq2_TPM.txt", header = T, row.names = 1)
GSE_count[1:5,1:5]
ncol(GSE_count)

GSE_met <- read.csv("/Users/linewulff/Documents/work/projects/FentonWulff_LP_MNP/Lit_data/GSE178209_Metadata_SMARTseq2.csv", fill = T, header = T, row.names = 1)
GSE_met[1:5,]

# check matches
nrow(GSE_met)==ncol(GSE_count)
colnames(GSE_count)[!colnames(GSE_count) %in% rownames(GSE_met)]

#### Seurat object ####
GSE <- CreateSeuratObject(counts = GSE_count, 
                          meta.data = GSE_met,
                          project = "Mulder et al.")

GSE <- subset(GSE, cells = rownames(GSE@meta.data[GSE@meta.data$FACS_Indexed_AllTissuesBatch1.2_tSNE1>(-25000),]))

FeatureScatter(GSE, feature1 = "RNA35PC_UMAP1",feature2 = "RNA35PC_UMAP2", group.by = "Seurat_clusters_35PCsRes3.25")
FeatureScatter(GSE, feature1 = "FACS_Indexed_AllTissuesBatch1.2_tSNE1",feature2 = "FACS_Indexed_AllTissuesBatch1.2_tSNE2", group.by = "Seurat_clusters_35PCsRes3.25")
FeatureScatter(GSE, feature1 = "X_RNA_tSNE1p30",feature2 = "X_RNA_tSNE2p30")

ggplot(GSE_met[Cells(GSE),], aes(x = RNA35PC_UMAP1, y = RNA35PC_UMAP2, colour = CD141.or.CADM1.APC.A))+
  geom_point()+
  scale_colour_viridis_c()+
  theme_classic()

GSE@meta.data$S2B <- NA
GSE@meta.data[GSE@meta.data$Seurat_clusters_35PCsRes3.25 %in% c(19),]$S2B <- "pre-DC"
GSE@meta.data[GSE@meta.data$Seurat_clusters_35PCsRes3.25 %in% c(8),]$S2B <- "cDC1"
GSE@meta.data[GSE@meta.data$Seurat_clusters_35PCsRes3.25 %in% c(0,17),]$S2B <- "pDC"
GSE@meta.data[GSE@meta.data$Seurat_clusters_35PCsRes3.25 %in% c(18),]$S2B <- "CD1C+CD1A+"
GSE@meta.data[GSE@meta.data$Seurat_clusters_35PCsRes3.25 %in% c(1,2,3,4,5,6,15),]$S2B <- "cDC2"
GSE@meta.data[GSE@meta.data$Seurat_clusters_35PCsRes3.25 %in% c(11),]$S2B <- "CD14hi mono"
GSE@meta.data[GSE@meta.data$Seurat_clusters_35PCsRes3.25 %in% c(16,9,14,10,12),]$S2B <- "CD16lo Mono"
GSE@meta.data[GSE@meta.data$Seurat_clusters_35PCsRes3.25 %in% c(13),]$S2B <- "T/NK cells"
GSE@meta.data[GSE@meta.data$Seurat_clusters_35PCsRes3.25 %in% c(7),]$S2B <- "B cells"
unique(GSE@meta.data$S2B)
FeatureScatter(GSE, feature1 = "RNA35PC_UMAP1",feature2 = "RNA35PC_UMAP2", group.by = "S2B")

pdf(paste(dato,"MulderEtAl_VlnProteinExp.pdf",sep="_"),height = 15,width = 10)
VlnPlot(GSE, group.by = "S2B",
        features = c("CD141.or.CADM1.APC.A",
                     "CD16.APC.Cy7.A",
                     "CD1a.Alexa.Fluor.700.A",
                     "CD123.BUV395.A",
                     "CD11b.BUV737.A",
                     "HLA.DR.BV.785.A",
                     "CD1c.BV421.A",
                     "CD34.or.CD163.BV605.A",
                     "CD3.14.19.20.BV650.A",
                     "CD5.BV711.A",
                     "L.D.DAPI.A",
                     "CD45RA.FITC.A",
                     "CD169.PE.A",
                     "CD206.PE.CF594.A",
                     "AUTOFLUO.PE.Cy5.A",
                     "CD88.PE.Cy7.A",
                     "FceRIa.PerCP.A",
                     "CD45.V500.A"),
        ncol = 4, pt.size = 0)
dev.off()


MNPverse <- readRDS("/Users/linewulff/Documents/work/projects/FentonWulff_LP_MNP/Lit_data/2021_MNP_Verse.RDS")
Idents(MNPverse) <- "Clusters"
DimPlot(MNPverse, label = T)
