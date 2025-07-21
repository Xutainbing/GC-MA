############# Fig3A - UMAP plot of T/NK subclusters

library(Seurat)
library(ggplot2)

sce.Bcell <- readRDS("scRNA-seq/GC_AF/data/sce.Bcell.rds")

my_cols = c('#87cdcb','#d5e9e8','#abaad4','#e7e7f3','#f4b3cb','#96a8d6','#FDB28F','#B49FDC','#FE9392','#3367A2','#FAE28C','#84a3c7','#FDE6D9')

DimPlot(sce.Bcell, group.by = c("CellType"), reduction = "fastMNNUMAP2D", cols = my_cols, pt.size = 1) +
  labs(x = "UMAP1", y = "UMAP2",title= NULL) +
  theme(
    panel.border = element_rect(fill=NA, color="black", size=2, linetype="solid"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    text = element_text(size = 20)  # 设置注释字体大小
  ) +
    guides(colour = guide_legend(override.aes = list(size = 8)))+
    annotate(
        "text", 
        x = -Inf, y = -Inf,  # 指定左下角位置
        label = "n = 4237 cells", 
        hjust = -0.1, vjust = -0.5,  # 调整位置
        size = 5, color = "black", fontface = "italic"
    )
