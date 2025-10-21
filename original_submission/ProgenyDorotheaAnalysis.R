### Progeny and Dorothea analysis
## same done for mature DC clusters and monomac clusters
## Author: Thomas Fenton

######## ----- Progeny ----- #######
# Libraries
library('progeny')
library('limma')
library('tidyr')
library('pheatmap')
library('tibble')

## At this point, the data should be split into colonic LP and ileal LP for each cluster##
## We create a data frame with the specification of the cells that belong to 
## each cluster to match with the Progeny scores
CellsClusters <- data.frame(Cell = names(Idents(mommac)), 
                            CellType = as.character(Idents(mommac)),
                            stringsAsFactors = FALSE)

## Compute the Progeny activity scores
mommac <- progeny(mommac, scale=FALSE, organism="Human", top=500, perm=1, 
                  return_assay = TRUE)

## Scale the pathway activity scores
mommac <- Seurat::ScaleData(mommac, assay = "progeny") 

## Transform Progeny scores into a data frame 
progeny_scores_df <- 
  as.data.frame(t(GetAssayData(mommac, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 

## Match Progeny scores with the cell clusters.
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

## Summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity),  std = sd(Activity))

## Prepare the data for the plot
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 

#reorder
head(summarized_progeny_scores)
head(x=summarized_progeny_scores[[]])
colnames(summarized_progeny_scores_df)
x<- c('7_SILP','7_cLP', '5_SILP','5_cLP', '3_SILP','3_cLP', '0_SILP','0_cLP', '1_SILP','1_cLP', '6_SILP','6_cLP', '2_SILP','2_cLP', '4_SILP','4_cLP', '9_SILP','9_cLP', '8_SILP','8_cLP')
y <- c( 'Trail','PI3K','Hypoxia','Estrogen', 'NFkB' ,'TNFa','EGFR', 'MAPK',   'VEGF','JAK-STAT', 'WNT','TGFb', 'p53' )

summarized_progeny_scores_df <- summarized_progeny_scores_df[x,]
summarized_progeny_scores_df <- summarized_progeny_scores_df[,y]

paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out=floor(paletteLength/2)))
progeny_hmap = pheatmap(t(summarized_progeny_scores_df[,]),fontsize=14, 
                        fontsize_row = 10, 
                        color=myColor, breaks = progenyBreaks, 
                        main = "PROGENy (500)", angle_col = 45,
                        treeheight_col = 0,  border_color = NA, cluster_cols =  FALSE, cluster_rows = F)




######## DOROTHEA #######
library(dorothea)
library(bcellViper)
library(dplyr)
library(viper)
library(pheatmap)

browseVignettes("dorothea")

## At this point, the data should be split into colonic LP and ileal LP for each cluster##
mommac <- mommac

################DOROTHEA#############################
# accessing expression data from bcellViper
data(bcellViper, package = "bcellViper")
# acessing (mouse) dorothea regulons
dorothea_regulon_mouse <- get(data("dorothea_hs", package = "dorothea"))

#filter to only include high quality regulons (could also be just A and B)
regulon <- dorothea_regulon_mouse %>%
  dplyr::filter(confidence %in% c("A","B"))

## We compute Viper Scores
mommac <- run_viper(mommac, regulon,
                  options = list(method = "scale", minsize = 4,
                                 eset.filter = FALSE, cores = 1,
                                 verbose = FALSE))

#scale data
DefaultAssay(mommac)<- 'dorothea'
mommac <- ScaleData(mommac)


## We transform Viper scores, scaled by seurat, into a data frame to better
## handling the results
viper_scores_df <- GetAssayData(mommac, slot = "scale.data",
                                assay = "dorothea") %>%
  data.frame(check.names = F) %>%
  t()

#to get original seurat clusters
DefaultAssay(mommac)<-'integrated'

## We create a data frame containing the cells and their clusters
CellsClusters <- data.frame(cell = names(Seurat::Idents(mommac)),
                            cell_type = as.character(Seurat::Idents(mommac)),
                            check.names = F)


## We create a data frame with the Viper score per cell and its clusters
viper_scores_clusters <- viper_scores_df  %>%
  data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  tidyr::gather(tf, activity, -cell) %>%
  inner_join(CellsClusters)

## We summarize the Viper scores by cellpopulation
summarized_viper_scores <- viper_scores_clusters %>%
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))

# We select the 20 most variable TFs. (20*4 populations = 80)
highly_variable_tfs <- summarized_viper_scores %>%
  group_by(tf) %>%
  mutate(var = var(avg))  %>%
  ungroup() %>%
  top_n(390, var) %>%
  distinct(tf)

## We prepare the data for the plot
summarized_viper_scores_df <- summarized_viper_scores %>%
  semi_join(highly_variable_tfs, by = "tf") %>%
  dplyr::select(-std) %>%  
  tidyr::spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE)

palette_length = 100
my_color = colorRampPalette(c("darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(min(summarized_viper_scores_df), 0,
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(summarized_viper_scores_df)/palette_length,
                   max(summarized_viper_scores_df),
                   length.out=floor(palette_length/2)))

#reorder
x<- c('7_SILP','7_cLP', '5_SILP','5_cLP', '3_SILP','3_cLP', '0_SILP','0_cLP', '1_SILP','1_cLP', '6_SILP','6_cLP', '2_SILP','2_cLP', '4_SILP','4_cLP', '9_SILP','9_cLP', '8_SILP','8_cLP')

summarized_viper_scores_df <- summarized_viper_scores_df[x,]
y <- c('WT1','FOSL2','CEBPB','HIF1A','FOS','NFKB1','REL','RELA','FOXL2','ATF4','ELK1','ATF2','CREB1','SP1','HSF1','REST','RFX5','FOSL1','STAT1', 'NR5A1')
summarized_viper_scores_df <- summarized_viper_scores_df[,y]

#heatmap
viper_hmap <- pheatmap(t(summarized_viper_scores_df),fontsize=14,
                       fontsize_row = 10,
                       color=my_color, breaks = my_breaks,
                       main = "DoRothEA (AB)", angle_col = 45,
                       treeheight_col = 0,  border_color = NA, cluster_cols =  FALSE, cluster_rows = F)


