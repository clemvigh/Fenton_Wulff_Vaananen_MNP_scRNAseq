## 3D umap plot of monomac
library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(htmlwidgets)
library(plotly)

## monocytes and macrohages 
monomac <- readRDS("~/Downloads/monomac_traj_versionApril2022.rds")
monomac@images <- list()

monomac <- UpdateSeuratObject(monomac)
Idents(monomac) <- "tSP_clustering_F4"


### make flat umap based on fixed rotation
scalar1 <- function(x) {x / sqrt(sum(x^2))}
original_data <- as.data.frame(monomac@reductions$umap_3d@cell.embeddings, row.names = F)

rot <- c(-99.701588586417, -0.607350399017409, -19.8461393107943)

rot <- rev(rot) * (-pi/180)
x_axis <- c(1,0,0)
y_axis <- c(0,1,0)
ref.plane <- c(0,0,1)
rot.mat <- matrix(c(cos(rot[1])*cos(rot[2]),
                    sin(rot[1])*cos(rot[2]),
                    -sin(rot[2]),
                    cos(rot[1])*sin(rot[2])*sin(rot[3])-sin(rot[1])*cos(rot[3]),
                    sin(rot[1])*sin(rot[2])*sin(rot[3])+cos(rot[1])*cos(rot[3]),
                    cos(rot[2])*sin(rot[3]),
                    cos(rot[1])*sin(rot[2])*cos(rot[3])+sin(rot[1])*sin(rot[3]),
                    sin(rot[1])*sin(rot[2])*cos(rot[3])-cos(rot[1])*sin(rot[3]),
                    cos(rot[2])*cos(rot[3])), nrow = 3)

doRotation <- function(x) c(rot.mat %*% x)
x_axis.rot <- scalar1(doRotation(x_axis))
y_axis.rot <- scalar1(doRotation(y_axis))
ref.plane.rot <- scalar1(doRotation(ref.plane))

projected_data <- t(apply(original_data, 1, function(x) {
  c(x_axis.rot %*% (x - ref.plane.rot), y_axis.rot %*% (x - ref.plane.rot))
}))

dat.toplot <- as.data.frame(projected_data); colnames(dat.toplot) <- c("umap3D_1","umap3D_2");rownames(dat.toplot) <- colnames(monomac)

# adding to object
monomac[["umap3D.tpc"]] <- CreateDimReducObject(embeddings = as.matrix(dat.toplot), 
                                                key = "umap3D_", 
                                                assay = DefaultAssay(monomac))
# plotting clusters and pseudotime
pdf("monomac_2D_rotation_clusters.pdf")
DimPlot(monomac, reduction = "umap3D.tpc", group.by = 'tSP_clustering_F4', raster = T, pt.size = 2, raster.dpi = c(1024, 1024)) + xlim (c(-5, 2.5)) + ylim (c(-8.5, 7.5)) + NoLegend() + coord_fixed(ratio = 0.5)
dev.off()

pdf("monomac_2D_rotation_pseudotime_wlegend.pdf")
FeaturePlot(monomac, reduction = "umap3D.tpc", feature = "pseudotime", raster = T, pt.size = 2, raster.dpi = c(1024, 1024)) + xlim (c(-5, 2.5)) + ylim (c(-8.5, 7.5)) + coord_fixed(ratio = 0.5) + scale_color_gradientn(colours = c('#4D8E8A', 'lightgray', '#6B669B'))
dev.off()
