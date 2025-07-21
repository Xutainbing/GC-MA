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


############# Fig3D - Ro/e enrichment heatmap

library(ComplexHeatmap)
library(circlize)
library(grid)

# Load Ro/e distribution function
source("distribution_Roe.R")

# Extract metadata
meta <- sce.Bcell@meta.data

# Calculate Ro/e enrichment
roe <- calTissueDist(
    dat.tb = meta,
    byPatient = FALSE,
    colname.cluster = "CellType",
    colname.patient = "patients",
    colname.tissue = "group1",
    method = "chisq",
    min.rowSum = 0
)

# Convert to data frame
roe_df <- as.data.frame.matrix(roe)

# Generate annotation symbols
p=roe_df
p[p > 1] <- "+++"
p[p> 0.8 & p <= 1] <- "++"
p[p>=0.2 & p <= 0.8] <- "+"
p[p> 0 & p < 0.2] <- "+/−"
p[p==0] <- "−"

# Define color scale
col_fun2 = colorRamp2(c(0, 0.2,0.8,1, 1.7), c("#3367A2", "#7FD0D8", "#FAE28C", "#FDB28F", "#FE9392"))

cellwidth = 1
cellheight = 1
cn = dim(roe_df)[2]
rn = dim(roe_df)[1]
w=cellwidth*cn
h=cellheight*rn

# Plot heatmap
Heatmap(roe_df,name ="Ro/e", col = col_fun2,
        heatmap_legend_param = list(legend_height = unit(3, "cm"),
                                    grid_width = unit(0.4, "cm"),
                                    labels_gp = gpar(col = "gray20",
                                                     fontsize = 8)),
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(p[i, j], x, y, vjust = 0.7,
                      gp = gpar(fontsize = 13,col="black"))
        })


############# Fig3E - KEGG enrichment and GOcircle plot for B mem cells

library(Seurat)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

# Subset B memory cells and set identity
sce.Bcell_bm <- subset(sce.Bcell, CellType == "B mem")
Idents(sce.Bcell_bm) <- "group1"

# Differential gene expression between MA and PBMC groups
deg_bm <- FindMarkers(sce.Bcell_bm, ident.1 = "MA", ident.2 = "PBMC", 
                      min.pct = 0.25, logfc.threshold = 0.25)
deg_bm <- deg_bm[order(deg_bm$avg_log2FC, decreasing = TRUE), ]
deg_genes <- rownames(deg_bm)

# Map gene symbols to ENTREZ IDs
entrez_ids <- mapIds(org.Hs.eg.db, 
                     keys = deg_genes,
                     keytype = "SYMBOL", 
                     column = "ENTREZID")

# KEGG enrichment
kegg_res <- enrichKEGG(gene = entrez_ids,
                       organism = "hsa",
                       keyType = "kegg",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)

# Select specific KEGG terms by row index (can be adjusted)
df <- kegg_res@result[c(1, 5, 7, 13, 23, 27, 28, 45, 90, 92), ]

# Combine DEGs and KEGG terms
enrich_df <- deg_bm %>%
  mutate(SYMBOL = rownames(.)) %>%
  inner_join(diff.df, by = "SYMBOL") %>%  # Make sure diff.df is defined beforehand
  dplyr::select(SYMBOL, ENTREZID, avg_log2FC)

# Format for GOplot
kegg_plot_df <- data.frame(
  Category = "KEGG",
  ID = df$ID,
  Term = df$Description,
  Genes = gsub("/", ", ", df$geneID),
  adj_pval = df$pvalue
)

genelist <- data.frame(
  ID = enrich_df$ENTREZID,
  logFC = enrich_df$avg_log2FC
)

# Generate GOcircle input and plot
circ_data <- circle_dat(kegg_plot_df, genelist)

GOBubble(circ_data,
         colour = c('skyblue', 'pink', 'red'),
         display = 'single',
         labels = 1,
         table.legend = FALSE) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    aspect.ratio = 1,
    panel.spacing = unit(0.2, "lines")
  ) +
  coord_cartesian(xlim = c(-3, 3), ylim = c(0, 12), expand = FALSE)


############# Fig3G - Pseudotime trajectory analysis for B cells using monocle3
## The code is the same as Fig. 2H.


############# Fig3H - UMAP distribution of clonal expansion in B cells
## The code is the same as Fig. 2D.


############# Fig3K - Isotype distribution in memory B cells

library(dplyr)
library(ggplot2)
library(stringr)
library(ggpubr)

# Merge clonotype and B cell metadata
df <- do.call(rbind, lapply(combined, as.data.frame))
rownames(df) <- df$barcode
meta <- sce.Bcell@meta.data

df <- df[rownames(df) %in% rownames(meta), ]
meta <- meta[rownames(df), ]
stopifnot(identical(rownames(df), rownames(meta)))

df3 <- cbind(df[, c("IGH", "Tissue")], meta[, c("CTgene", "sampleid", "group1", "CellType")])
df3 <- df3[!is.na(df3$IGH), ]

# Assign isotype based on IGH sequence
df3 <- df3 %>%
  mutate(Isotype = case_when(
    str_detect(IGH, "IGHM") ~ "IgM",
    str_detect(IGH, "IGHA") ~ "IgA",
    str_detect(IGH, "IGHD") ~ "IgD",
    str_detect(IGH, "IGHG") ~ "IgG",
    TRUE ~ "Other"
  ))

# Calculate frequency
stats <- df3 %>%
  group_by(sampleid, CellType, Isotype) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sampleid, CellType) %>%
  mutate(frequency = count / sum(count)) %>%
  ungroup()

# Plot for memory B cells
ggplot(stats %>% filter(CellType == "B mem"), 
       aes(x = Isotype, y = frequency)) +
  geom_boxplot(fill = my_cols[5], width = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 2, color = "black") +
  stat_compare_means(method = "kruskal.test", label = "p.format", label.y = 0.2) +
  theme_minimal(base_size = 14) +
  labs(title = "Isotype distribution in memory B cells",
       x = NULL, y = "Relative Frequency") +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )


############# Fig3L - Hallmark pathway enrichment (GSVA) in B cell subtypes

library(Seurat)
library(GSVA)
library(GSEABase)
library(msigdbr)
library(readxl)
library(dplyr)
library(stringr)
library(pheatmap)

# Load hallmark gene sets
hallmark_sets <- msigdbr(species = "human", category = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)

# Subset B cells in malignant ascites (MA), excluding AtM B
sce_bcell_ma <- subset(sce.Bcell, subset = group1 == "MA" & CellType != "AtM B")

# Get normalized expression matrix
expr_mat <- as.matrix(sce_bcell_ma@assays$RNA@data)

# Run GSVA
gsva_scores <- gsva(expr_mat, 
                    gset.idx.list = hallmark_sets, 
                    method = "gsva", 
                    kcdf = "Gaussian", 
                    mx.diff = TRUE)

# Format GSVA result
gsva_scores <- t(gsva_scores)
colnames(gsva_scores) <- substring(colnames(gsva_scores), 10)
gsva_scores <- gsva_scores[rownames(sce_bcell_ma@meta.data), , drop = FALSE]

# Add cell type annotation
gsva_scores <- data.frame(gsva_scores)
gsva_scores$CellType <- sce_bcell_ma$CellType

# Calculate average enrichment per cell type
gsva_avg <- gsva_scores %>%
  group_by(CellType) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames("CellType") %>%
  t() %>%
  as.data.frame()

# Format row names
rownames(gsva_avg) <- gsub("_", " ", str_to_title(rownames(gsva_avg)))

# Load hallmark category annotations
anno <- read_excel("scRNA-seq/GC_AF/data/hallmark.xlsx") %>%
  column_to_rownames("Hallmark.Name") %>%
  arrange(Process.Category)
rownames(anno) <- rownames(anno)
anno <- anno["Process.Category", drop = FALSE]

# Match annotation order to matrix
gsva_avg <- gsva_avg[rownames(anno), ]

# Plot heatmap
pheatmap(gsva_avg,
         scale = "row",
         clustering_method = "ward.D",
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         annotation_row = anno,
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_names_row = FALSE,
         gaps_row = c(3,9,12,19,26,31,37,50),
         add_lines_row = c(3,9,12,19,26,31,37,50),
         color = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(100)),
         border = FALSE,
         fontsize = 14,
         legend = TRUE,
         legend_args = list(legend = "right"))


############# Fig3M-N - 
## The code is the same as Fig. 2L.









