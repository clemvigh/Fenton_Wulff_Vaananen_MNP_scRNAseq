library(Seurat)
library(openxlsx)
library(ggplot2)
library(scales)
library(ggbiplot)
library(pheatmap)
library(ComplexHeatmap)
library(colorRamp2)
library(RColorBrewer)


MNP_LP <- readRDS("R3_CD14CD1Cprol_MNPLP.rds")
MNP_LP@images <- list()


Idents(MNP_LP) <- "integrated_snn_res.2.8"
MNP_LP@active.ident <- factor(MNP_LP@active.ident, levels=c('24','12','22','16','4','33','19','17','2','15','35','3','18','13','21','5','28','0','10','29','37','34','36','38','26','7','8','30','32','27','11','20','9','1','14','6','31','25','23'))

# heatmap
genes.to.plot <- c("ZBTB16", "C5AR1", "SOD2", "CLEC4E", "CD52", "ITGAL", "TLR4", "ABCA1", "CD14", "S100A9", 
                   "MAFB", "HES1", "CTSD", "PLD3", "SEPP1", "MERTK", "MAF", "NR1H3", "ID3", "CD163", "C1QA", 
                   "CSF1R", "CCR2", "AP1S3", "FLT3", "SEPT6", "IRF4", "CST7", "SLC38A1", "DUSP4", "CD86", 
                   "CSF2RA", "CSF2RB", "CD207", "CD1C", "CCR7", "MKI67")

MNP.averages <- AverageExpression(MNP_LP, assay = "RNA", return.seurat = TRUE)
MNP.averages <- ScaleData(MNP.averages)


## re-numbering for dimplot
# Extract current cluster identities
current_clusters <- Idents(MNP_LP)

new_order <- c('24','12','22','16','4','33','19','17','2','15','35','3','18','13','21','5','28','0',
               '10','29','37','34','36','38','26','7','8','30','32','27','11','20','9','1','14','6','31','25','23')

# Create a mapping of old to new cluster numbers
cluster_mapping <- setNames(seq_along(new_order), new_order)

# Map the current cluster identities to the new order
new_clusters <- as.character(cluster_mapping[as.character(current_clusters)])

# Update the identities in the Seurat object
MNP_LP <- SetIdent(MNP_LP, value = new_clusters)
MNP_LP$clusters_reorder <- MNP_LP@active.ident

Idents(MNP_LP) <- "clusters_reorder"
MNP.averages@active.ident <- factor(MNP.averages@active.ident, levels = seq(1, 39, 1))

pdf("fig1G_ordered.pdf")
DimPlot(MNP_LP, label = T, label.size = 2, 
        raster = TRUE, raster.dpi = c(512, 512), pt.size = 2) + coord_fixed(ratio = 1)
dev.off()


heatmap_cols <- rev(c('#D7191C','#FDAE61','#FFFFBF','#ABD9E9','#2C7BB6'))
p2 <- DoHeatmap(MNP.averages, features = genes.to.plot, draw.lines = FALSE, size = 3.5, group.by = "orig.ident") + NoLegend() + scale_fill_gradientn(colours = heatmap_cols)

pdf("fig1E_heatmap_extended_genelist_ordered.pdf", height = 5, width = 7)
print(p2)
dev.off()


