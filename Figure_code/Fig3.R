############# Fig3A - UMAP plot of B cell subclusters

library(Seurat)
library(ggplot2)

# Load B cell Seurat object
sce.Bcell <- readRDS("scRNA-seq/GC_AF/data/sce.Bcell.rds")

# Define colors for clusters
my_cols <- c('#87cdcb','#d5e9e8','#abaad4','#e7e7f3','#f4b3cb','#96a8d6',
             '#FDB28F','#B49FDC','#FE9392','#3367A2','#FAE28C','#84a3c7','#FDE6D9')

my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', 
             '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', 
             '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', 
              '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')

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


############# Fig3B - DotPlot of B cell subtype markers

library(Seurat)
library(ggplot2)
library(dplyr)

# Define marker genes for B cell subtypes
Bcells_markers <- list(
  Naive = c('MS4A1', 'CCR7', 'TCL1A', 'IL4R'),
  Mem = c('AIM2', 'TNFRSF13B'),
  Plasma = c('CD38', 'MZB1', 'XBP1'),
  Cycling = c('MKI67', 'TOP2A'),
  Atm = c('FCRL5', 'ITGAX', 'ZEB2', 'FGR')
)

# Generate DotPlot and get data
p_all_markers <- DotPlot(sce.Bcell, features = Bcells_markers, assay = "RNA", group.by = "CellType") + 
  coord_flip()

data <- p_all_markers$data

# Rename columns for clarity
colnames(data) <- c("AvgExpUnscaled", "PctExpressed", "Feature", "CellType", "AvgExpression", "Group")

# Set factor levels for ordering
data$CellType <- factor(data$CellType, levels = rev(c("B naive", "B mem", "Plasma", "Plasma cycling", "AtM B")))

# Maintain feature order within groups
data <- data %>%
  group_by(Group) %>%
  mutate(Feature = factor(Feature, levels = unique(Feature))) %>%
  ungroup()

# Plot
ggplot(data, aes(x = Feature, y = CellType)) +
  geom_point(aes(size = PctExpressed, fill = AvgExpression), shape = 21, stroke = 0.2) +
  scale_size_continuous(range = c(1, 10)) +
  scale_fill_gradient2(low = "#F8F8FF", mid = "grey", high = "#E54924", midpoint = 0) +
  facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 14),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_blank(),
    panel.spacing.x = unit(1, "lines")
  )


############# Fig3C - Relative frequencies comparison (B mem & AtM B cells)

library(ggplot2)
library(ggpubr)
library(ggsignif)
library(patchwork)
library(tidyverse)
library(dplyr)

# Extract relevant metadata
df <- sce.Bcell@meta.data[, c("CellType", "patients", "sampleid", "group1")]

# Calculate cell type proportions per sample
df_ratio <- df %>%
    count(sampleid, CellType) %>%
    group_by(sampleid) %>%
    mutate(percentage = n / sum(n) * 100) %>%
    ungroup() %>%
    mutate(ID = str_split_fixed(sampleid, "_", 2)[,1],
           group = str_split_fixed(sampleid, "_", 2)[,2],
           celltype = factor(CellType, levels = c("B naive", "B mem", "Plasma", "Plasma cycling", "AtM B")))


# Define sample order and point shapes
sample_names <- unique(df$sampleid)
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
        scale_color_manual(values = my36colors, limits = sample_names) +
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
p1 <- plot_group_comparison(df_ratio, "B mem", y_limit = 70)
p2 <- plot_group_comparison(df_ratio, "AtM B", y_limit = 20)

# Combine
p1 | p2


############# Fig3D










