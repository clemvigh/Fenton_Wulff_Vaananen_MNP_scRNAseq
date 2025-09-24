# module load tools openssl/1.1.1b intel/perflibs/2020_update4 gcc/8.2.0 R/3.6.1
library(ggradar)
library(openxlsx)

df <- read.csv("cDC2_ileum_colon_GOterms_RADARPLOT.csv", header = T, row.names = 1)
df <- read.xlsx("GOterms_cDCsubsets_RADARPLOT_final.xlsx", colNames = T)
df <- read.xlsx("GOterms_M6-8_subsets_radarplot_new.xlsx", colNames = T)

df[, 2:ncol(df)] <- sqrt(-log10(df[, 2:ncol(df)]))
colnames(df) <- c("group", paste("Var", 1:7, sep = " "))
colnames(df) <- c("group", paste("Var", 1:12, sep = " "))
colnames(df) <- c("group", paste("Var", 1:10, sep = " "))

# Color for the lines
lcols <- c("#1CBDC2", "#F3766E")
lcols <- c("#DD6DAA", "#8D8AC2", "#BD992F")
lcols <- c("#59C1AE", "#54B8DD", "#4BA4F1")

g <- ggradar(df,
        background.circle.colour = "white",
        grid.min = 0,
        grid.max = max(df[,2:ncol(df)]),
        gridline.min.linetype = 1,
        values.radar = c(sqrt(-log10(0.05))),
        gridline.mid.colour = "black",
        group.colours = lcols, 
        label.gridline.min = F,
        label.gridline.mid = F,
        label.gridline.max = F,
        group.point.size = 0)


ggsave("cDC2_ileum_colon_GOterms_radarplot_sqrt.pdf", plot = g, device = cairo_pdf)
ggsave("cDCsubsets_GOterms_radarplot_sqrt_final.pdf", plot = g, device = cairo_pdf)
ggsave("M6-8_subsets_GOterms_radarplot_sqrt_final.2.pdf", plot = g, device = cairo_pdf)

