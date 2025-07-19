############# Fig1B - UMAP plot

library(Seurat)
library(ggplot2)

# Load data
sce.all <- readRDS("scRNA-seq/GC_AF/data/celltype_20240415.rds")

# Define color palette
my_cols <- c('#FE9392','#B49FDC','#87cdcb','#d5e9e8','#abaad4','#f4b3cb',
             '#FDB28F','#FAE28C','#96a8d6','#3367A2','#84a3c7','#FDE6D9','#e7e7f3')

my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', 
             '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', 
             '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', 
              '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')

# Set cell type order
sce.all$celltype <- factor(sce.all$celltype, 
                           levels = c('T cells','T_NK cells','B cells','Plasma',
                                      'Myeloid','Epithelial','Fibroblast','Platelets'))

# Plot UMAP
DimPlot(sce.all, group.by = "celltype", reduction = "umap", cols = my_cols, raster = FALSE) +
    labs(x = "UMAP1", y = "UMAP2", title = NULL) +
    theme_minimal(base_size = 12) +
    theme(
        panel.border = element_rect(fill = NA, color = "black", size = 1),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank()
    ) +
    annotate(
        "text",
        x = -Inf, y = -Inf,
        label = "n = 190392 cells",
        hjust = -0.1, vjust = -0.5,
        size = 5, color = "black", fontface = "italic"
    )


############# Fig1C - Marker gene expression (violin plot)

library(MySeuratWrappers)

# Define marker genes by lineage
genes_to_marker <- c(
    'PTPRC','CD3D','CD3E',        # T cells
    'GNLY','GZMA',                # NK cells
    'MS4A1','CD79A',              # B cells
    'IGHG1','MZB1',               # Plasma cells
    'CD163','CD14',               # Myeloid cells
    'EPCAM','KRT8',               # Epithelial
    'DCN','COL1A2',               # Fibroblasts
    'TUBB1','GP9'                 # Platelets
)

# Set identity
Idents(sce.all) <- sce.all$celltype

# Plot stacked violin plot
VlnPlot(sce.all,
        features = genes_to_marker,
        stacked = TRUE,
        pt.size = 0,
        cols = my_cols,
        direction = "horizontal",
        x.lab = "", y.lab = "") +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    )


############# Fig1F - Ro/e enrichment heatmap

library(ComplexHeatmap)
library(circlize)
library(grid)

# Load Ro/e distribution function
source("distribution_Roe.R")

# Extract metadata
meta <- sce.all@meta.data

# Calculate Ro/e enrichment
roe <- calTissueDist(
    dat.tb = meta,
    byPatient = FALSE,
    colname.cluster = "celltype",
    colname.patient = "patients",
    colname.tissue = "group1",
    method = "chisq",
    min.rowSum = 0
)

# Convert to data frame
roe_df <- as.data.frame.matrix(roe)

# Generate annotation symbols
symbols <- roe_df
symbols[symbols > 1] <- "+++"
symbols[symbols > 0.8 & symbols <= 1] <- "++"
symbols[symbols >= 0.2 & symbols <= 0.8] <- "+"
symbols[symbols > 0 & symbols < 0.2] <- "+/−"
symbols[symbols == 0] <- "−"

# Define color scale
col_fun2 <- colorRamp2(
    c(0, 0.2, 0.8, 1, max(roe_df)),
    c("#3367A2", "#7FD0D8", "#FAE28C", "#FDB28F", "#FE9392")
)

# Plot heatmap
Heatmap(
    data,
    name = "Ro/e",
    col = col_fun2,
    # Legend styling
    heatmap_legend_param = list(
        legend_height = unit(3, "cm"),
        grid_width = unit(0.4, "cm"),
        labels_gp = gpar(
            col = "gray20",
            fontsize = 8
        )
    ),
    # Cell annotation styling
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(
            p[i, j], x = x, y = y, vjust = 0.7,
            gp = gpar(fontsize = 13, col = "black")
        )
    }
)


############# Fig1E - Cell composition per sample and total barplot

library(ggplot2)
library(tidyverse)
library(patchwork)
library(dplyr)

# Extract relevant metadata
df <- sce.all@meta.data[, c("celltype", "patients", "sampleid", "group1")]

# Calculate cell type proportions per sample
df_ratio <- df %>%
    count(sampleid, celltype) %>%
    group_by(sampleid) %>%
    mutate(percentage = n / sum(n) * 100) %>%
    ungroup() %>%
    mutate(ID = str_split_fixed(sampleid, "_", 2)[,1],
           group = str_split_fixed(sampleid, "_", 2)[,2],
           celltype = factor(celltype, levels = c('T cells', 'T_NK cells', 'B cells', 'Plasma', 
                                                  'Myeloid', 'Epithelial', 'Fibroblast', 'Platelets')))

# Stacked barplot: percentage per sample
p <- ggplot(df_ratio, aes(x = sampleid, y = percentage, fill = celltype)) +
    geom_bar(stat = "identity", position = "fill", width = 0.8, color = "black", size = 0) +
    scale_fill_manual(values = my_cols) +
    labs(x = "Sample", y = "Relative Frequency") +
    theme_classic(base_size = 13) +
    theme(
        axis.text.x = element_text(angle = 60, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks.length = unit(0.1, "cm"),
        text = element_text(color = "black")
    )

# Total cell number per cell type
total_cells_by_type <- df_ratio %>%
    group_by(celltype) %>%
    summarise(total = sum(n)) %>%
    mutate(celltype = factor(celltype, levels = rev(levels(celltype))))

side_p <- ggplot(total_cells_by_type, aes(x = total, y = celltype)) +
    geom_bar(stat = "identity", fill = my_cols[2]) +
    labs(x = "Cell Number", y = NULL) +
    theme_classic(base_size = 13) +
    theme(
        axis.text.x = element_text(angle = 60, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none",
        text = element_text(color = "black")
    )

# Combine main plot and side bar
combined_plot <- p + side_p + 
    plot_layout(nrow = 1, widths = c(7, 3))


############# Fig1F
library(ComplexHeatmap)
library(circlize)

source("distribution_Roe.R")
meta = sce.all@meta.data
roe= calTissueDist(
  dat.tb = meta,
  byPatient = F,
  colname.cluster = "celltype",
  colname.patient = "patients",
  colname.tissue = "group1",
  method = "chisq",
  min.rowSum = 0
)

roe_df <- as.data.frame.matrix(roe)

p=roe_df
p[p > 1] <- "+++"
p[p> 0.8 & p <= 1] <- "++"
p[p>=0.2 & p <= 0.8] <- "+"
p[p> 0 & p < 0.2] <- "+/−"
p[p==0] <- "−"

col_fun2 = colorRamp2(c(0, 0.2,0.8,1, 2.023587), c("#3367A2", "#7FD0D8", "#FAE28C", "#FDB28F", "#FE9392"))

cellwidth = 1
cellheight = 1
cn = dim(roe_df)[2]
rn = dim(roe_df)[1]
w=cellwidth*cn
h=cellheight*rn

Heatmap(roe_df,name ="Ro/e", col = col_fun2,
        heatmap_legend_param = list(legend_height = unit(3, "cm"),
                                    grid_width = unit(0.4, "cm"),
                                    labels_gp = gpar(col = "gray20",
                                                     fontsize = 8)),
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(p[i, j], x, y, vjust = 0.7,
                      gp = gpar(fontsize = 13,col="black"))
        })
		



############# Fig1G - Relative frequencies comparison (T_NK cells & Plasma)

library(ggplot2)
library(ggpubr)
library(ggsignif)
library(patchwork)

# Define sample order and point shapes
sample_names <- unique(df$ID)
custom_shapes <- c(16, 17)

# Function to draw boxplot with statistics
plot_group_comparison <- function(data, celltype_name, y_limit) {
    df_sub <- data %>% filter(celltype == celltype_name)
    
    ggplot(df_sub, aes(x = group, y = percentage)) +
        geom_boxplot(fill = NA, width = 0.5, position = position_dodge(0.75)) +
        geom_jitter(
            aes(color = ID, shape = group),
            size = 3, alpha = 0.5,
            position = position_jitter(width = 0.2)
        ) +
        labs(
            x = NULL,
            y = "Relative Frequency (%)",
            title = celltype_name
        ) +
        theme_minimal(base_size = 13) +
        theme(
            panel.grid = element_blank(),
            axis.line = element_line(color = "black"),
            axis.ticks = element_line(color = "black"),
            plot.title = element_text(hjust = 0.5),
            legend.position = "none",
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
        ) +
        ylim(0, y_limit) +
        scale_color_manual(values = my_cols, limits = sample_names) +
        scale_shape_manual(values = custom_shapes) +
        stat_compare_means(
            comparisons = list(c("MA", "PBMC")),
            method = "wilcox.test",
            label = "p.format",
            tip.length = 0.01,
            size = 3
        )
}

# Draw plots
p1 <- plot_group_comparison(df_ratio, "T_NK cells", y_limit = 80)
p2 <- plot_group_comparison(df_ratio, "Plasma", y_limit = 8)

# Combine
p1 | p2
