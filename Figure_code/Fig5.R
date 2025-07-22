############# Fig5A - UMAP plot of Myeloid cell subclusters

library(Seurat)
library(ggplot2)

# Load B cell Seurat object
epi <- readRDS("scRNA-seq/GC_AF/data/sce.epi.rds")

# Define colors for clusters
my_cols=  c('#AB3282', '#8C549C', '#5F3D69')

my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', 
             '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', 
             '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', 
              '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')

# Plot UMAP colored by CellType
DimPlot(epi, group.by = "CellType", cols = my_cols, pt.size = 1) +
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
  annotate("text", x = -Inf, y = -Inf, label = "n = 36282 cells", 
           hjust = -0.1, vjust = -0.5, size = 5, color = "black", fontface = "italic")


############# Fig5B - Expression of marker genes in epithelial cells (UMAP)

library(scRNAtoolVis)

# Plot selected marker genes on UMAP with consistent formatting
featureCornerAxes(
  object = epi,                        
  reduction = 'umap',                  
  groupFacet = NULL,                   
  relLength = 0.5,                     
  relDist = 0.2,                       
  features = c("EPCAM", "KRT8", "KRT18", "PTPRC", "CD14", "C1QB"),  
  aspect.ratio = 0.5,                 
  themebg = 'bwCorner',               
  pSize = 0.05,                       
  minExp = 0,                         
  maxExp = 4,                         
  nLayout = 3                      
)


############# Fig5C - Proportional distribution of epithelial cell subtypes by patient

library(ggplot2)
library(dplyr)
library(stringr)
library(patchwork)

# Extract sample information
epi@meta.data <- epi@meta.data %>%
  mutate(
    patients = str_extract(orig.ident, "^[^_]+"),
    groups = str_extract(orig.ident, "(?<=_).*"),
    samples = orig.ident
  )

# Count and normalize cell subtype frequency by patient
df_ratio <- epi@meta.data %>%
  select(CellType, patients) %>%
  group_by(patients, CellType) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(CellType) %>%
  mutate(relative_freq = n / sum(n)) %>%
  ungroup()

# Set factor levels for consistent ordering
df_ratio$CellType <- factor(df_ratio$CellType, levels = unique(df_ratio$CellType))
df_ratio$patients <- factor(df_ratio$patients, levels = unique(df_ratio$patients))

# Bar plot: raw count proportions
p1 <- ggplot(df_ratio, aes(x = CellType, y = n, fill = patients)) +
  geom_bar(stat = "identity", position = "fill", width = 0.95) +
  scale_fill_manual(values = my36colors) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11)
  )

# Bar plot: normalized proportions
p2 <- ggplot(df_ratio, aes(x = CellType, y = relative_freq, fill = patients)) +
  geom_bar(stat = "identity", position = "fill", width = 0.95) +
  scale_fill_manual(values = my36colors) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    legend.position = "none"
  )

# Combine plots using patchwork
p1 + p2 + plot_layout(ncol = 2)


############# Fig5D - Volcano-like rank plot of DEGs

library(Seurat)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(circlize)
library(dplyr)

# Set identity for DEG analysis
Idents(epi) <- epi$CellType

# Differential expression between epithelial subtypes
epi.markers <- FindAllMarkers(epi, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
epi.markers$pct.diff <- epi.markers$pct.1 - epi.markers$pct.2
epi.markers <- epi.markers %>% filter(p_val_adj < 0.05)

# Use external DEG result (deg) for plotting
deg <- epi.markers %>% 
  arrange(desc(avg_log2FC)) %>%
  mutate(
    rank = row_number(),
    p_val_adj = ifelse(p_val_adj == 0, min(p_val_adj[p_val_adj != 0]), p_val_adj),
    adj = -log10(p_val_adj)
  )

# Select top 10 upregulated and downregulated genes
top <- bind_rows(
  head(deg, 10),
  tail(deg, 10)
) %>% mutate(Gene = rownames(.))

# Color gradient function
colors <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100)
col_fun <- colorRamp2(seq(0, 303, length.out = 100), colors)

# Plot
p <- ggplot(deg, aes(x = rank, y = avg_log2FC, color = adj)) +
  geom_point(size = 1) +
  geom_text_repel(
    data = top, aes(label = Gene),
    size = 4, color = "#5f70b3",
    force = 20, point.padding = 0.5, segment.size = 0.3,
    segment.color = "grey20", segment.alpha = 0.8, nudge_y = -0.1
  ) +
  scale_color_gradientn(colors = col_fun(seq(0, 303, length.out = 100))) +
  scale_x_continuous(breaks = seq(0, 3000, by = 1000)) +
  labs(
    x = "DEGs", y = "avg_log2FC",
    title = "Epi vs Epi-Immu", color = "Adjusted\nP-value"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.border = element_rect(color = "black", size = 1.2, fill = NA),
    axis.text = element_text(face = "bold"),
    axis.title = element_text(size = 14)
  )

p


############# Fig5E - Correlation heatmap between cell types

library(Seurat)
library(pheatmap)
library(corrplot)
library(RColorBrewer)

# Set RNA as default assay
DefaultAssay(mac_epi) <- "RNA"

# Merge macrophage subclusters into a unified 'Mac' label
mac_epi$CellType3 <- ifelse(mac_epi$CellType %in% c('M08_Mac-IL1B', 'M09_Mac-C1QC', 'M10_Mac-PLTP', 'M11_Mac-SPP1'),
                            "Mac", mac_epi$CellType)

# Set new cell type identities
Idents(mac_epi) <- mac_epi$CellType3

# Compute average expression per cell type
avg_exp <- AverageExpression(mac_epi, return.seurat = FALSE)$RNA

# Select top 2000 most variable genes
top_genes <- names(sort(apply(avg_exp, 1, sd), decreasing = TRUE))[1:2000]
avg_exp_filtered <- avg_exp[top_genes, ]

# Compute Spearman correlation matrix between cell types
cor_mat <- cor(avg_exp_filtered, method = "spearman")

# Draw combined correlation heatmap (circle + square)
corrplot.mixed(corr = cor_mat,
               lower = "circle",
               upper = "square",
               order = "hclust",
               lower.col = rev(colorRampPalette(brewer.pal(11, "Spectral"))(200)),
               upper.col = rev(colorRampPalette(brewer.pal(11, "Spectral"))(200)),
               tl.col = "black",
               tl.cex = 1,
               cl.cex = 1,
               pch.cex = 1.5,
               addCoef.col = "black",
               outline = TRUE)


############# Fig5G - CNV scoring and violin plot

# Load inferCNV object and annotation
load("scRNA-seq/GC_AF/data/inferCNV/infercnv_obj.Rdata")
Anno <- read.table("scRNA-seq/GC_AF/data/inferCNV/annoFile.txt", header = FALSE, sep = '\t')
rownames(Anno) <- Anno$V1

# Extract CNV expression matrix
expr <- infercnv_obj@expr.data

# Define CNV reference baseline using B cells and Plasma cells
ref_cells <- cbind(
  expr[, infercnv_obj@reference_grouped_cell_indices$`B cells`],
  expr[, infercnv_obj@reference_grouped_cell_indices$`Plasma`]
)

# Set CNV thresholds
mean_expr <- rowMeans(ref_cells)
sd_expr <- apply(ref_cells, 1, sd)
down <- mean(mean_expr) - 2 * mean(sd_expr)
up   <- mean(mean_expr) + 2 * mean(sd_expr)
oneCopy <- up - down
a1 <- down - 2 * oneCopy
a2 <- down - oneCopy
a3 <- up + oneCopy
a4 <- up + 2 * oneCopy

# Initialize scoring matrix
cnv_score_table <- expr
cnv_score_mat <- as.matrix(expr)

# Assign categorical CNV scores
cnv_score_table[cnv_score_mat > 0  & cnv_score_mat < a2] <- "A"  # complete loss
cnv_score_table[cnv_score_mat >= a2 & cnv_score_mat < down] <- "B"  # one copy loss
cnv_score_table[cnv_score_mat >= down & cnv_score_mat < up] <- "C"  # neutral
cnv_score_table[cnv_score_mat >= up & cnv_score_mat <= a3] <- "D"  # one copy gain
cnv_score_table[cnv_score_mat > a3 & cnv_score_mat <= a4] <- "E"  # two copy gain
cnv_score_table[cnv_score_mat > a4] <- "F"  # more than two copies

# Convert to numeric scores
cnv_score_pts <- cnv_score_mat
cnv_score_pts[cnv_score_table == "A"] <- 2
cnv_score_pts[cnv_score_table == "B"] <- 1
cnv_score_pts[cnv_score_table == "C"] <- 0
cnv_score_pts[cnv_score_table == "D"] <- 1
cnv_score_pts[cnv_score_table == "E"] <- 2
cnv_score_pts[cnv_score_table == "F"] <- 2

# Calculate per-cell CNV scores
cell_scores_CNV <- data.frame(cnv_score = colSums(cnv_score_pts))
Anno$cnv_score <- cell_scores_CNV$cnv_score

# Add cell type annotation for violin plot
Anno$CellType <- factor(
  Anno$CellType,
  levels = c('B cells', 'Plasma', 'M02_cDC1-CLEC9A', 'Epi-Immu', 'Epi')
)

# Define color palette
my_cols <- c('#AB3282', '#8C549C', '#5F3D69', '#68A180', '#3A6963', '#53A85F')

# Plot violin and boxplot of CNV scores
ggplot(Anno, aes(x = CellType, y = cnv_score)) +
  geom_violin(aes(fill = CellType), cex = 1.2) +
  geom_boxplot(width = 0.1, cex = 1.2) +
  scale_fill_manual(values = my_cols) +
  theme_classic(base_size = 20) +
  theme(
    axis.text = element_text(color = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = 'none',
    panel.border = element_rect(fill = NA, color = "black", size = 2),
    text = element_text(size = 20),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 12, face = "italic")
  ) +
  labs(x = NULL) +
  geom_signif(
    comparisons = list(c("Epi", "Epi-Immu"), c("Epi-Immu", "M02")),
    map_signif_level = TRUE,
    textsize = 6
  )





