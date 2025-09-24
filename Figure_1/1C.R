library(Seurat)
library(ggplot2)

# loading data
obj <- readRDS("R3_CD14CD1Cprol_MNPLP.rds")
obj@images <- list()
DefaultAssay(obj) <- "RNA"


## featureplot for figure 1C
genes.to.plot <- c("S100A8", "CD14", "C1QA", "CD1C", "CSF1R", "FLT3")
mycols_blue <-  c("#bdbdbd","#d9d9d9","skyblue3","steelblue","steelblue4","#1B3346")

pdf("fig1D_additional_genes.pdf")
patchwork::wrap_plots(FeaturePlot(obj, features = genes.to.plot, combine = FALSE, coord.fixed = 1,
                                  raster = TRUE, raster.dpi = c(512, 512), pt.size = 2)) & 
  scale_color_gradientn(colours=mycols_blue, na.value = "#bdbdbd")

dev.off()
