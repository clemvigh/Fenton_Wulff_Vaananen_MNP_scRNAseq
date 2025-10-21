######### Rotating 3D UMAP: Save as Interactive HTML
library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(htmlwidgets)
library(plotly)

## load and update the seurat object 
obj <- readRDS("~/Downloads/monomac_seurat.rds")
obj <- UpdateSeuratObject(obj)

## extract 3D UMAP and cluster labels
plot.data <- FetchData(object = obj, 
                       layer = "umap_3d",
                       vars = c("UMAP_1", "UMAP_2", "UMAP_3", "tSP_clustering_F4"))
plot.data$label <- paste(rownames(plot.data))

## color palette for clusters
cluster_colors <- c(
  "#E87D72", "#D58E31", "#B99D33", "#58B434", "#55BB77",
  "#56BEAC", "#54B8D6", "#4CA9F5", "#888FF1", "#C87AF7", "#E76CD7"
)

## create interactive 3D scatter plot
fig <- plot_ly(data = plot.data, 
               x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
               color = ~tSP_clustering_F4,
               jitter = T,
               colors = cluster_colors,
               type = "scatter3d", 
               mode = "markers", 
               marker = list(size = 4, width = 1), 
               text = ~label, 
               hoverinfo = "text") %>%
  layout(scene = list(camera = list(
    eye = list(
      x = 1.25,
      y = 1.25,
      z = 1.25
    ),
    center = list(x = 0,
                  y = 0,
                  z = 0)
  ))) %>%
  onRender("
      function(el, x){
  var id = el.getAttribute('id');
  var gd = document.getElementById(id);
  Plotly.update(id).then(attach);
  function attach() {
    var cnt = 0;
    
    function run() {
      rotate('scene', Math.PI / 1000);
      requestAnimationFrame(run);
    } 
    run();
    
    function rotate(id, angle) {
      var eye0 = gd.layout[id].camera.eye
      var rtz = xyz2rtz(eye0);
      rtz.t += angle;
      
      var eye1 = rtz2xyz(rtz);
      Plotly.relayout(gd, id + '.camera.eye', eye1)
    }
    
    function xyz2rtz(xyz) {
      return {
        r: Math.sqrt(xyz.x * xyz.x + xyz.y * xyz.y),
        t: Math.atan2(xyz.y, xyz.x),
        z: xyz.z
      };
    }
    
    function rtz2xyz(rtz) {
      return {
        x: rtz.r * Math.cos(rtz.t),
        y: rtz.r * Math.sin(rtz.t),
        z: rtz.z
      };
    }
  };
}
    ") %>% layout(legend = list(itemsizing = "constant")) # this one sizes the legend markers according to markers in plot

## save as standalone HTML
saveWidget(partial_bundle(fig), file = "~/Desktop/monomac_tspace_3d.HTML", selfcontained = TRUE)
