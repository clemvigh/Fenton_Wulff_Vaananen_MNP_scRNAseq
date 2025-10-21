### Plots for monomac figure 2 + 3 + S2 ###
library(Seurat)
library(tidyverse)
library(VennDiagram)
library(cowplot)
library(ggplot2)
library(scCustomize)
library(RColorBrewer)
library(DESeq2)
library(patchwork)
library(UCell)
library(ComplexHeatmap)
library(circlize)
library(SingleR)

## Loading data
monomac <- readRDS("~/Downloads/monomac_traj_versionApril2022.rds")
monomac@images <- list()

Idents(monomac) <- "tSP_clustering_F4"
macs <- subset(monomac, idents = c("M6", "M7", "M8"))

### DGE between colon and ileum for M6, M7 and M8
# first subset to only include patients that have both segments
Idents(macs) <- "orig.ident"
macs.paired <- subset(macs, idents = c("cLP_pat5", "cLP_pat6"), invert = T)

run_deseq2_paired <- function(cluster) {
  
  sub <- subset(macs.paired, idents = as.character(cluster))
  counts <- AggregateExpression(sub, return.seurat = F, verbose = T, assays = "RNA",
                                group.by = "orig.ident")
  counts <- counts$RNA
  
  anno_df <- data.frame(IDs = colnames(counts), 
                              segment = factor(str_extract(colnames(counts),  "^[^-]+")),
                              patient = factor(str_extract(colnames(counts),  "[^-]+$")))
  
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = anno_df,
                                design = ~ segment + patient)
  
  dds <- DESeq(dds)
  res <- results(dds)
  
  return(list(res = res, dds = dds))
}

Idents(macs.paired) <- "tSP_clustering_F4"

M6.dge <- run_deseq2_paired("M6")
M7.dge <- run_deseq2_paired("M7")
M8.dge <- run_deseq2_paired("M8")


M6_silp_clp <- data.frame(results(M6.dge$dds, contrast = c("segment", "SILP", "cLP")))
M7_silp_clp <- data.frame(results(M7.dge$dds, contrast = c("segment", "SILP", "cLP")))
M8_silp_clp <- data.frame(results(M8.dge$dds, contrast = c("segment", "SILP", "cLP")))

write.xlsx(M6_silp_clp, file = "M6_silp_clp_dge.xlsx", sep = "\t", rowNames = T, quote = F)
write.xlsx(M7_silp_clp, file = "M7_silp_clp_dge.xlsx", sep = "\t", rowNames = T, quote = F)
write.xlsx(M8_silp_clp, file = "M8_silp_clp_dge.xlsx", sep = "\t", rowNames = T, quote = F)

# saving up in each segment
M6_ileum <- M6_silp_clp |> mutate(gene = rownames(M6_silp_clp)) |> 
  filter(log2FoldChange > 0, padj < 0.05) |> arrange(desc(log2FoldChange)) |> select(gene, log2FoldChange, padj)
M6_colon <- M6_silp_clp |> mutate(gene = rownames(M6_silp_clp)) |> 
  filter(log2FoldChange < 0, padj < 0.05) |> arrange(log2FoldChange) |> select(gene, log2FoldChange, padj)

M7_ileum <- M7_silp_clp |> mutate(gene = rownames(M7_silp_clp)) |> 
  filter(log2FoldChange > 0, padj < 0.05) |> arrange(desc(log2FoldChange)) |> select(gene, log2FoldChange, padj)
M7_colon <- M7_silp_clp |> mutate(gene = rownames(M7_silp_clp)) |> 
  filter(log2FoldChange < 0, padj < 0.05) |> arrange(log2FoldChange) |> select(gene, log2FoldChange, padj)

M8_ileum <- M8_silp_clp |> mutate(gene = rownames(M8_silp_clp)) |> 
  filter(log2FoldChange > 0, padj < 0.05) |> arrange(desc(log2FoldChange)) |> select(gene, log2FoldChange, padj)
M8_colon <- M8_silp_clp |> mutate(gene = rownames(M8_silp_clp)) |> 
  filter(log2FoldChange < 0, padj < 0.05) |> arrange(log2FoldChange) |> select(gene, log2FoldChange, padj)


## volcano plots -- Figure 2G
fc_cutoff <- 0.58
p_cutoff <- 0.05

res <- M6_silp_clp # done for each of the three subsets M6, M7 and M8


color_scale <- c("ileum" = "#F3766E", "colon" = "#1CBDC2", "ns" = "grey")
keyvals <- ifelse(
  res$log2FoldChange < -0.58 & res$padj < 0.05, '#1CBDC2',
  ifelse(res$log2FoldChange > 0.58 & res$padj < 0.05, '#F3766E',
         'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == '#1CBDC2'] <- 'colon'
names(keyvals)[keyvals == 'grey'] <- 'ns'
names(keyvals)[keyvals == '#F3766E'] <- 'ileum'

pdf("M6_cLP_SILP_volcano_new.pdf", width = 4, height = 4)
ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = names(keyvals))) +
  geom_point(size = 2) +
  scale_color_manual(values = color_scale) +
  # Add horizontal and vertical threshold lines
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed") +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
  #geom_text(data = top_genes, aes(label = gene), vjust = -1, hjust = 0.5, color = "black", size = 3) + #here can adjust top_genes gened or sig_genes
  # Labels and titles
  labs(
    title = "M6",
    x = expression(Log[2]~"fold change"),
    y = expression(-Log[10]~italic(P))
  ) +
  coord_cartesian(xlim = c(-7.5, 7.5), ylim = c(0,25)) + 
  # Customize theme for a clean look
  theme_classic() 
dev.off()



## venn diagrams -- Figure 2H
pdf("cLP_M6-8_VennDiagram.pdf")
# Create the Venn diagram and store it as an object
venn <- venn.diagram(
  x = list(M6_colon$gene, M7_colon$gene, M8_colon$gene),
  category.names = c("M6", "M7", "M8"),
  fill = c("#00c1ab", "#00bbda", "#00acfc"),
  alpha = 0.5,
  lwd = c(0, 0, 0),
  #lty = "blank",
  filename = NULL,
  output = TRUE
)
grid.draw(venn)
dev.off()

pdf("SILP_M6-8_VennDiagram.pdf")
# Create the Venn diagram and store it as an object
venn <- venn.diagram(
  x = list(M6_ileum$gene, M7_ileum$gene, M8_ileum$gene),
  category.names = c("M6", "M7", "M8"),
  fill = c("#59C1AE", "#54B8DD", "#4BA4F1"),
  alpha = 0.5,
  lwd = c(0, 0, 0),
  lty = "blank",
  filename = NULL,
  output = TRUE
)
grid.draw(venn)
dev.off()


## heatmap of intersecting genes -- Figure 2I
intersect.genes <- c(Reduce(intersect, list(M6_colon$gene, M7_colon$gene, M8_colon$gene)),
                     Reduce(intersect, list(M6_ileum$gene, M7_ileum$gene, M8_ileum$gene)))

data.to.plot <- AverageExpression(macs, return.seurat = T, features = intersect.genes,
                                    group.by = c("tSP_clustering_F4", "segment"))

mycols <- rev(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))

pdf("2I_cLPvsSILP_M6-8_overlap_heatmap.pdf", height = 16)
DoHeatmap(data.to.plot, features = intersect.genes, 
          draw.lines = F) + scale_fill_gradientn(colours = mycols, na.value = "white")
dev.off()


## violin plots for M6, M7 and M8 -- Figure 2D
genes.to.plot <- c("LYVE1", "SIGLEC1", "CD4", "ADAMDEC1", "LIPA", "PLD3", "FOLR2", 
                   "MAF", "C2", "IL4I1", "LGALS3", "CD68", "MARCO", "COLEC12")

pdf("macs_vlnplot.pdf", height = 16)
VlnPlot(sub, features = genes.to.plot, pt.size = 0, split.by = "segment", cols = c("#55B8C7", "#EA7D6C"), ncol = 2)
dev.off()

## more genes to plot
genes_enriched_M6 <- c("CLEC12A", "CCR2", "TNFSF10", "LYST", "S100Z", "TYK2", "TLR10", "AIM2", "CLEC7A")
genes_enriched_M7 <- c("CD163", "MRC1", "MS4A7", "VEGFA", "PDGFB", "TGFB1", "MRC1", "CD93",
                       "MERTK", "PTGS1", "PTGS2", "ARL5B")
genes_enriched_M8 <- c("FABP4", "FABP5", "CD63", "S100A11", "APOC1", "MTG1", "MT1E", "MMP12")

pdf("macs_M6_genes_vlnplot.pdf", height = 16)
VlnPlot(sub, features = genes_enriched_M6, pt.size = 0, split.by = "segment", cols = c("#55B8C7", "#EA7D6C"), ncol = 2)
dev.off()

pdf("macs_M7_genes_vlnplot.pdf", height = 16)
VlnPlot(sub, features = genes_enriched_M7, pt.size = 0, split.by = "segment", cols = c("#55B8C7", "#EA7D6C"), ncol = 2)
dev.off()

pdf("macs_M8_genes_vlnplot.pdf", height = 16)
VlnPlot(sub, features = genes_enriched_M8, pt.size = 0, split.by = "segment", cols = c("#55B8C7", "#EA7D6C"), ncol = 2)
dev.off()

pdf("macs_LGALS3+9_vlnplot.pdf", height = 16)
VlnPlot(sub, features = c("LGALS3", "LGALS9"), pt.size = 0, split.by = "segment", cols = c("#55B8C7", "#EA7D6C"), ncol = 2)
dev.off()


genes.to.plot <- c("MT1G", "MT1X", "MT2A", "FOLR2", "LGALS9", "FABP5")
pdf("macs_vlnplot_final.pdf", height = 24)
VlnPlot(macs, features = genes.to.plot, pt.size = 0, split.by = "segment", cols = c("#55B8C7", "#EA7D6C"), ncol = 1)
dev.off()


genes.to.plot <- c("IL4I1", "ACP5", "MT2A", "FOLR2", "LGALS9", "FABP5")
pdf("macs_vlnplot_IL4I1+ACP5.pdf", height = 24)
VlnPlot(macs, features = genes.to.plot, pt.size = 0, split.by = "segment", cols = c("#55B8C7", "#EA7D6C"), ncol = 1)
dev.off()



## Comparing to external datasets, frode et al. -- Figure S2D 
lpm <- readRDS("paper1-singleCell/LpM.rds")
lpm <- UpdateSeuratObject(lpm)
lpm@images <- list()

# getting the markers 
M1 <- markers |> filter(cluster == "M1") |> filter(gene %in% rownames(lpm)) |> arrange(p_val_adj, desc(avg_log2FC)) |> pull(gene) |> head(50)
M2 <- markers |> filter(cluster == "M2") |> filter(gene %in% rownames(lpm)) |> arrange(p_val_adj, desc(avg_log2FC)) |> pull(gene) |> head(50)
M3 <- markers |> filter(cluster == "M3") |> filter(gene %in% rownames(lpm)) |> arrange(p_val_adj, desc(avg_log2FC)) |> pull(gene) |> head(50)
M4 <- markers |> filter(cluster == "M4") |> filter(gene %in% rownames(lpm)) |> arrange(p_val_adj, desc(avg_log2FC)) |> pull(gene) |> head(50)
M5 <- markers |> filter(cluster == "M5") |> filter(gene %in% rownames(lpm)) |> arrange(p_val_adj, desc(avg_log2FC)) |> pull(gene) |> head(50)
M6 <- markers |> filter(cluster == "M6") |> filter(gene %in% rownames(lpm)) |> arrange(p_val_adj, desc(avg_log2FC)) |> pull(gene) |> head(50)
M7 <- markers |> filter(cluster == "M7") |> filter(gene %in% rownames(lpm)) |> arrange(p_val_adj, desc(avg_log2FC)) |> pull(gene) |> head(50)
M8 <- markers |> filter(cluster == "M8") |> filter(gene %in% rownames(lpm)) |> arrange(p_val_adj, desc(avg_log2FC)) |> pull(gene) |> head(50)
M9 <- markers |> filter(cluster == "M9") |> filter(gene %in% rownames(lpm)) |> arrange(p_val_adj, desc(avg_log2FC)) |> pull(gene) |> head(50)
M10 <- markers |> filter(cluster == "M10") |> filter(gene %in% rownames(lpm)) |> arrange(p_val_adj, desc(avg_log2FC)) |> pull(gene) |> head(50)
M11 <- markers |> filter(cluster == "M11") |> filter(gene %in% rownames(lpm)) |> arrange(p_val_adj, desc(avg_log2FC)) |> pull(gene) |> head(50)

lpm <- AddModuleScore(lpm, features = list(M1), name = "M1_signature")
lpm <- AddModuleScore(lpm, features = list(M2), name = "M2_signature")
lpm <- AddModuleScore(lpm, features = list(M3), name = "M3_signature")
lpm <- AddModuleScore(lpm, features = list(M4), name = "M4_signature")
lpm <- AddModuleScore(lpm, features = list(M5), name = "M5_signature")
lpm <- AddModuleScore(lpm, features = list(M6), name = "M6_signature")
lpm <- AddModuleScore(lpm, features = list(M7), name = "M7_signature")
lpm <- AddModuleScore(lpm, features = list(M8), name = "M8_signature")
lpm <- AddModuleScore(lpm, features = list(M9), name = "M9_signature")
lpm <- AddModuleScore(lpm, features = list(M10), name = "M10_signature")
lpm <- AddModuleScore(lpm, features = list(M11), name = "M11_signature")


Idents(lpm) <- "seurat_clusters"
p1 <- DimPlot(lpm) + coord_fixed(ratio = 0.75)
p2 <- FeaturePlot(lpm, features = "M6_signature1", order = T) + 
  coord_fixed(ratio = 0.75) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p3 <- FeaturePlot(lpm, features = "M7_signature1", order = T) + 
  coord_fixed(ratio = 0.75) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
p4 <- FeaturePlot(lpm, features = "M8_signature1", order = T) + 
  coord_fixed(ratio = 0.75) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

pdf("frode_comparisons_top50.pdf", width = 12)
p1 + p2 + p3 + p4 + plot_layout(ncol = 4)
dev.off()