######### Rotating 3D UMAP: Save as Interactive HTML
library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(htmlwidgets)
library(plotly)
library(scales)

## load and update the seurat object 
obj <- readRDS("~/Downloads/cDC_seurat.rds")
obj <- UpdateSeuratObject(obj)

## extract 3D UMAP and cluster labels
plot.data <- FetchData(object = obj, vars = c("umaptpc_1",
                                              "umaptpc_2", 
                                              "umaptpc_3", 
                                              "combined_clusters"))
plot.data$label <- paste(rownames(plot.data))

## color palette for clusters
cluster_colors <- c(
  "cDC1" = "#F8766D",
  "cDC2" = "#B79F00",
  "Ambiguous" = "#00BA38",
  "cDC3" = "#00BFC4",
  "CCR7+ DC cluster 1" = "#A23E84",
  "CCR7+ DC cluster 2" = "#E5A6C8"
)

## create interactive 3D scatter plot
fig <- plot_ly(data = plot.data, 
               x = ~umaptpc_1, y = ~umaptpc_2, z = ~umaptpc_3,
               color = ~combined_clusters,
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
saveWidget(partial_bundle(fig), file = "~/Desktop/cDC_tspace_3d.HTML", selfcontained = TRUE)
