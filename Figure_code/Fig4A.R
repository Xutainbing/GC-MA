############# Fig4A - UMAP plot of Myeloid cell subclusters

library(Seurat)
library(ggplot2)

# Load B cell Seurat object
sce.Myeloid <- readRDS("scRNA-seq/GC_AF/data/Myeloid.rds")

# Define colors for clusters
my_cols <- c('#87cdcb','#d5e9e8','#abaad4','#e7e7f3','#f4b3cb','#96a8d6',
             '#FDB28F','#B49FDC','#FE9392','#3367A2','#FAE28C','#84a3c7','#FDE6D9')

my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', 
             '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', 
             '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', 
              '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')

# Plot UMAP colored by CellType
DimPlot(sce.Myeloid, group.by = "CellType", cols = my36colors, pt.size = 1) +
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
  annotate("text", x = -Inf, y = -Inf, label = "n = 85655 cells", 
           hjust = -0.1, vjust = -0.5, size = 5, color = "black", fontface = "italic")


############# Fig4B - Expression of marker genes across Myeloid subclusters
library(Seurat)
library(ggplot2)

# Set identity to CellType2 and specify order
sce.Myeloid$CellType2 <- factor(sce.Myeloid$CellType2, levels = c('M01','M02','M03','M04','M05','M06',
                                                  'M07','M08','M09','M10','F01'))
Idents(sce.Myeloid) <- sce.Myeloid$CellType2

# Define markers
markers <- c('CD14','CD163','LILRA4','PTCRA','CLEC9A','XCR1','CD1C','CLEC10A','LAMP3',
             'CCL19','CX3CR1','CDKN1C','FCN1','S100A9','CD3E','CD247','IL1B','CD83',
             'MRC1','C1QA','C1QC','APOE','PLTP','STAB1','FOLR2','OLR1','EREG','DCN','COL1A2')

# Plot violin plots
VlnPlot(sce.Myeloid, features = markers, group.by = "CellType2", pt.size = 0,
        stacked = TRUE, direction = "horizontal", cols = my36colors) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.ticks.y = element_blank())

