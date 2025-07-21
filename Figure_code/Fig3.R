############# Fig3A - UMAP plot of B cell subclusters

library(Seurat)
library(ggplot2)

# Load B cell Seurat object
sce.Bcell <- readRDS("scRNA-seq/GC_AF/data/sce.Bcell.rds")

# Define colors for clusters
my_cols <- c('#87cdcb','#d5e9e8','#abaad4','#e7e7f3','#f4b3cb','#96a8d6',
             '#FDB28F','#B49FDC','#FE9392','#3367A2','#FAE28C','#84a3c7','#FDE6D9')

# Plot UMAP colored by CellType
DimPlot(sce.Bcell, group.by = "CellType", cols = my_cols, pt.size = 1) +
  labs(x = "UMAP1", y = "UMAP2", title = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 8))) +
  annotate("text", x = -Inf, y = -Inf, label = "n = 4237 cells", 
           hjust = -0.1, vjust = -0.5, size = 5, color = "black", fontface = "italic")

