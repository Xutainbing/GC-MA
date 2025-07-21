############# Fig2A - UMAP plot of T/NK subclusters

library(Seurat)
library(ggplot2)

# Load T/NK cell Seurat object
sce.Tcell <- readRDS("scRNA-seq/GC_AF/data/sce.Tcell.rds")

# Define custom color palette for subclusters
pair_my_cols <- c(
    'T01_CD4T' = '#87cdcb',
    'T02_CD4T' = '#d5e9e8',
    'T03_CD4T' = '#abaad4',
    'T04_CD4T' = '#e7e7f3',
    'T05_CD4T' = '#f4b3cb',
    'T01_CD8T' = '#96a8d6',
    'T02_CD8T' = '#FDB28F',
    'T03_CD8T' = '#B49FDC',
    'T04_CD8T' = '#FE9392',
    'γδT'      = '#FAE28C',
    'NK'       = '#84a3c7'
)

# Plot UMAP
DimPlot(sce.Tcell, group.by = "CellType2", label = TRUE, pt.size = 0.5, cols = pair_my_cols) +
    labs(x = "UMAP1", y = "UMAP2", title = NULL) +
    scale_color_manual(
        values = pair_my_cols,
        breaks = names(pair_my_cols),
        labels = c(
            'T01_CD4_Tn_CCR7', 'T02_CD4_Tcm_ANXA1', 'T03_CD4_Tfh_TOX2', 'T04_CD4_Treg_FOXP3', 'T05_CD4_Temra_GNLY',
            'T01_CD8_Tn_CCR7', 'T02_CD8_Tem/Trm_ITGA1', 'T03_CD8_Temra_CX3CR1', 'T04_CD8_MAIT_SLC4A10',
            'γδT', 'NK'
        )
    ) +
    theme_minimal(base_size = 12) +
    theme(
        panel.border     = element_rect(fill = NA, color = "black", size = 1),
        axis.text        = element_blank(),
        axis.ticks       = element_blank(),
        axis.line        = element_blank(),
        panel.grid       = element_blank()
    ) +
    guides(colour = guide_legend(override.aes = list(size = 8))) +
    annotate(
        "text",
        x = -Inf, y = -Inf,
        label = "n = 61976 cells",
        hjust = -0.1, vjust = -0.5,
        size = 5, color = "black", fontface = "italic"
    )


############# Fig2B - Marker gene expression heatmap for T/NK subclusters

library(SCP)

# Set subcluster factor levels
sce.Tcell$CellType2 <- factor(
    sce.Tcell$CellType2,
    levels = c(
        'T01_CD4T', 'T02_CD4T', 'T03_CD4T', 'T04_CD4T', 'T05_CD4T',
        'T01_CD8T', 'T02_CD8T', 'T03_CD8T', 'T04_CD8T',
        'γδT', 'NK'
    )
)

# Define key marker genes
all_Tcell_markers <- c(
    'CD4', 'CD8A', 'CCR7', 'LEF1', 'TCF7',            # Tn
    'CD69', 'IL7R', 'ANXA1',                          # Tcm
    'PDCD1', 'CTLA4', 'TOX2', 'ICOS',                 # Th1
    'FOXP3', 'IL2RA', 'TIGIT',                        # Treg
    'GNLY', 'NKG7', 'PRF1',                           # Temra (CD4)
    'ITGA1', 'ITGAE', 'CXCR4',                        # Trm
    'CX3CR1', 'FGFBP2', 'KLRG1',                      # Temra (CD8)
    'SLC4A10', 'TRAV1-2', 'RORC',                     # MAIT
    'TRGV9', 'TRDV1', 'TRDC1',                        # γδT
    'KLRF1', 'NCR1', 'KLRC1'                          # NK
)

# Generate heatmap
ht <- GroupHeatmap(
    srt = sce.Tcell,
    group.by = "CellType2",
    features = all_Tcell_markers,
    height = 8,
    width = 4,
    group_palette = "simspec",
    feature_split_palette = "simspec",
    group_palcolor = NULL,
    feature_split_palcolor = NULL,
    nlabel = 0,
    show_row_names = TRUE,
    row_names_side = "left"
)

ht$plot


############# Fig2C - Ro/e distribution for T/NK subclusters between AF and PBMC

library(ggplot2)

# Load Ro/e distribution calculation function
source("distribution_Roe.R")

# Extract metadata
meta <- sce.Tcell@meta.data

# Compute Ro/e enrichment score
roe <- calTissueDist(
    dat.tb = meta,
    byPatient = FALSE,
    colname.cluster = "CellType2",
    colname.patient = "patients",
    colname.tissue = "group1",
    method = "chisq",
    min.rowSum = 0
)

# Format result as data frame
roe_df <- as.data.frame.matrix(roe)
roe_df$CellType <- rownames(roe_df)

# Reshape for plotting
long_data <- roe_df %>%
    pivot_longer(
        cols = c("MA", "PBMC"),
        names_to = "Group",
        values_to = "Expression"
    )

# Add enrichment/depletion annotation
long_data$Type <- ifelse(long_data$Expression > 1, "Enrichment", "Depletion")

# Plot bubble chart
ggplot(long_data, aes(x = Group, y = CellType, size = Expression, fill = Type)) +
    geom_point(shape = 21, color = "black", alpha = 0.85) +
    scale_size_continuous(
        range = c(3, 12),
        breaks = c(0.2, 0.8, 1, 1.5, 2, 3),
        name = "Ro/e"
    ) +
    scale_fill_manual(
        values = c("Enrichment" = "#f4b3cb", "Depletion" = "#abaad4"),
        name = NULL
    ) +
    labs(
        x = NULL,
        y = NULL
    ) +
    theme_minimal(base_family = "STHeiti") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
    )


############# Fig2D - UMAP distribution of clonal expansion in T/NK cells

library(scRepertoire)
library(ggplot2)
library(Seurat)

# Load all TCR contig files (19 samples from MA and PBMC)
file_paths <- list.files(
    "scRNA-seq/GC_AF/raw_data/TCR",
    pattern = "_filtered_contig_annotations.csv$",
    full.names = TRUE,
    recursive = TRUE
)
contig_list <- lapply(file_paths, read.csv)

# Define patient IDs and tissue origins (must match order of contig_list)
samples <- c("RM001", "ZQQ003", "ZQQ003", "LFQ004", "S005ZJM",
             "S006SWQ", "S006SWQ", "S007LXX", "S007LXX", "S008ZJF",
             "S008ZJF", "S009ZHX", "S009ZHX", "S010YGF", "S010YGF",
             "S012WRJ", "S012WRJ", "RM001", "S005ZJM")

tissues <- c("PBMC", "MA", "PBMC", "PBMC", "PBMC",
             "MA", "PBMC", "MA", "PBMC", "MA",
             "PBMC", "MA", "PBMC", "MA", "PBMC",
             "MA", "PBMC", "MA", "MA")

# Combine all TCR data
combined <- combineTCR(
    contig_list,
    samples = samples,
    ID = tissues
)

# Add tissue type as a metadata variable
combined <- addVariable(
    combined,
    variable.name = "Tissue",
    variables = tissues
)

# Assign clonal expansion status to Seurat object
sce.Tcell <- combineExpression(
    combined,
    sce.Tcell,
    cloneCall = "strict",
    group.by = "sample",
    proportion = FALSE,
    cloneSize = c(Single = 1,
                  Small = 5,
                  Medium = 20,
                  Large = 100,
                  Hyperexpanded = 500)
)

# Define custom color palette for cloneSize
clone_colors <- rev(hcl.colors(n = 7, palette = "inferno", fixup = TRUE)[c(1, 3, 4, 5, 7)])

# Plot clone size distribution on UMAP
DimPlot(sce.Tcell, group.by = "cloneSize") +
    scale_color_manual(values = clone_colors) +
    labs(
        x = "UMAP1",
        y = "UMAP2",
        title = NULL,
        color = "Clone Size"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        panel.border     = element_rect(fill = NA, color = "black", size = 1),
        axis.text        = element_blank(),
        axis.ticks       = element_blank(),
        axis.line        = element_blank(),
        panel.grid       = element_blank()
    )  +
    guides(colour = guide_legend(override.aes = list(size = 4)))


############# Fig2E - TCR diversity (Gini.simpson) in T/NK subclusters between AF and PBMC

library(ggpubr)
library(ggplot2)
library(gridExtra)

# Define color palette and shape settings
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')
custom_shapes <- c(16, 17)

# Function to calculate and plot TCR diversity
plot_diversity <- function(celltype = NULL, title = "All T cells") {
    meta <- sce.Tcell@meta.data
    if (!is.null(celltype)) {
        meta <- meta[meta$CellType2 == celltype, ]
    }

    combined_filtered <- lapply(combined, function(df) {
        df[df$barcode %in% rownames(meta), ]
    })

    df <- clonalDiversity(
        combined_filtered,
        cloneCall = "strict",
        x.axis = "Tissue"
    )[["data"]]
    df <- df[df$variable == "gini.simpson", ]
    sample_names <- df$Group

    ggplot(df, aes(x = Tissue, y = value)) +
        geom_boxplot(fill = NA, width = 0.5, position = position_dodge(0.75)) +
        geom_jitter(aes(color = Group, shape = Tissue), size = 3, alpha = 0.5,
                    position = position_jitter(width = 0.2)) +
        labs(x = NULL, y = "Gini.simpson", title = paste("TCR diversity of", title)) +
        theme_minimal() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(color = "black"),
            axis.ticks = element_line(color = "black"),
            plot.title = element_text(hjust = 0.5),
            legend.position = "none",
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
        ) +
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

# Generate plots for all T cells and selected subclusters
p1 <- plot_diversity(NULL, "all T cells")
p2 <- plot_diversity("T01_CD4T", "T01_CD4T")
p3 <- plot_diversity("T02_CD8T", "T02_CD8T")
p4 <- plot_diversity("T03_CD8T", "T03_CD8T")

# Arrange plots in a grid
grid.arrange(p1, p2, p3, p4, ncol = 2)


############# Fig2F - STARTRAC analysis

library(Startrac)
library(tidyverse)
library(Seurat)
library(data.table)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

# Subset cell types
sce.Tcell2 <- subset(sce.Tcell, CellType %in% c(
  "T01_CD4_Tn_CCR7", "T01_CD8_Tn_CCR7", "T02_CD4_Tcm_ANXA1", "T02_CD8_Tem/Trm_ITGA1",
  "T03_CD4_Th1_TOX2", "T04_CD4_Treg_FOXP3", "T03_CD8_Temra_CX3CR1", "T05_CD4_Temra_GNLY",
  "T04_CD8_MAIT_SLC4A10"))

# Load and combine all VDJ data
files <- list.files("scRNA-seq/GC_AF/raw_data/TCR/", recursive = TRUE, full.names = TRUE, pattern = "filtered_contig_annotations.csv")
labels <- gsub("_filtered_contig_annotations.csv", "", basename(files))

vdj <- map2_dfr(files, labels, ~ {
  df <- read.csv(.x)
  df$file <- .y
  df
}) %>%
  filter(high_confidence == "True", productive == "True", chain %in% c("TRA", "TRB")) %>%
  mutate(Cell_name = paste0(file, "_", barcode))

# Process best TRA and TRB
get_best_chain <- function(vdj_data, chain_type) {
  vdj_chain <- vdj_data %>% filter(chain == chain_type) %>% arrange(desc(umis), desc(reads))
  best <- vdj_chain %>% group_by(Cell_name) %>% summarise(reads = max(reads), umis = max(umis))
  inner_join(vdj_chain, best)
}

vdj_a <- get_best_chain(vdj, "TRA")
vdj_b <- get_best_chain(vdj, "TRB")

final_vdj <- full_join(vdj_a, vdj_b, by = "Cell_name", suffix = c(".TRA", ".TRB"))

# Merge with metadata
sce.Tcell2@meta.data <- sce.Tcell2@meta.data %>%
  mutate(Cell_name = rownames(.)) %>%
  inner_join(final_vdj, by = "Cell_name") %>%
  filter(productive.TRA == "True", productive.TRB == "True") %>%
  mutate(Clone_AA = paste(cdr3.TRA, cdr3.TRB, sep = "_")) %>%
  arrange(Clone_AA)

# Generate Clone_ID and Clone_NUM
tmp <- sce.Tcell2@meta.data %>%
  group_by(Clone_AA) %>%
  summarise(Clone_NUM = n()) %>%
  mutate(Clone_ID = paste0("Clone_", row_number()))

sce.Tcell2@meta.data <- left_join(sce.Tcell2@meta.data, tmp)

# Prepare Startrac metadata
sce.Tcell2.meta <- sce.Tcell2@meta.data %>%
  transmute(
    Cell_Name = Cell_name,
    clone.id = Clone_ID,
    clone.status = ifelse(Clone_NUM > 1, "Clonal", "NoClonal"),
    orig.ident,
    patient = patients,
    stype = ifelse(CellType2 %in% c("T01_CD4T","T02_CD4T","T03_CD4T","T04_CD4T","T05_CD4"), "CD4", "CD8"),
    majorCluster = CellType2,
    loc = group1
  )

# Run STARTRAC
tic("Startrac.run")
out2 <- Startrac.run(sce.Tcell2.meta, proj = "GC_AF", n.perm = 1000, cores = 20, verbose = TRUE)

df <- out2@cluster.sig.data %>% filter(aid != "GC_AF", index != "gini")

# Plot
index_colors <- c("expa" = "#87cdcb", "migr" = "#abaad4", "tran" = "#f4b3cb")

ggplot(df, aes(x = majorCluster, y = value)) +
  geom_boxplot(aes(color = index), outlier.shape = NA, fill = NA) +
  geom_jitter(aes(color = index), width = 0.2, size = 2) +
  facet_wrap(~ index, scales = "free_y", ncol = 1, strip.position = "right") +
  scale_color_manual(values = index_colors) +
  labs(x = NULL, y = "STARTRAC index") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "none",
    strip.placement = "outside",
    strip.text.y.right = element_text(angle = 0)
  )


############# Fig2G - pairwise STARTRAC-tran indices

plot(out2,index.type="pairwise.tran",byPatient=T)


############# Fig2H - Pseudotime trajectory analysis for CD8+ T cells using Monocle3

library(Seurat)
library(monocle3)
library(harmony)

# Subset CD8+ T cells
sce.cd8T <- subset(sce.Tcell2, subset = CellType2 %in% c("T01_CD8T", "T02_CD8T", "T03_CD8T"))

# Reinitialize Seurat object
sce.cd8T <- CreateSeuratObject(counts = sce.cd8T@assays$RNA@counts, meta.data = sce.cd8T@meta.data)

# Normalize and correct for confounders
sce.cd8T <- SCTransform(sce.cd8T, vars.to.regress = c("nCount_RNA", "percent_mito"),
                        verbose = TRUE, method = "glmGamPoi", vst.flavor = "v2")

# PCA and batch correction
sce.cd8T <- RunPCA(sce.cd8T, assay = "SCT", npcs = 100)
sce.cd8T <- RunHarmony(sce.cd8T, group.by.vars = "batchid", reduction = "pca", 
                       assay.use = "SCT", reduction.save = "harmony")

# Dimensional reduction and clustering
sce.cd8T <- RunUMAP(sce.cd8T, reduction = "harmony", dims = 1:50)
sce.cd8T <- FindNeighbors(sce.cd8T, reduction = "harmony", dims = 1:50)
sce.cd8T <- FindClusters(sce.cd8T, resolution = 0.5)

saveRDS(sce.cd8T, "scRNA-seq/GC_AF/data/sce.cd8Tcell.rds")

# Convert to monocle3 object
data <- GetAssayData(sce.cd8T, assay = "RNA", slot = "counts")
meta <- sce.cd8T@meta.data
gene_anno <- data.frame(gene_short_name = rownames(data), row.names = rownames(data))

cds <- new_cell_data_set(data, cell_metadata = meta, gene_metadata = gene_anno)

# Monocle3 preprocessing
cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(cds, alignment_group = "batchid")
cds <- reduce_dimension(cds, reduction_method = "UMAP", umap.fast_sgd = TRUE)
cds <- cluster_cells(cds)

# Sync UMAP with Seurat for consistent visualization
cds@int_colData$reducedDims$UMAP <- Embeddings(sce.cd8T, "umap")[colnames(cds), ]

# Learn trajectory
cds <- learn_graph(cds)

# Define root and order cells
root_cells <- colnames(cds)[colData(cds)$CellType == "T01_CD8_Tn_CCR7"]
cds <- order_cells(cds, root_cells = root_cells)

# Plot pseudotime
plot_cells(
    cds,
    color_cells_by = "pseudotime",
    label_cell_groups = FALSE,
    label_leaves = TRUE,
    label_branch_points = TRUE,
    show_trajectory_graph = FALSE,
    graph_label_size = 1.5,
    group_label_size = 4,
    cell_size = 0.8
) +
    scale_color_gradientn(
        colours = c("#3288bd", "#66c2a5", "#ffffbf", "#f46d43", "#9e0142"),
        values = scales::rescale(seq(0, 1, length.out = 5))
    ) +
    labs(x = "UMAP1", y = "UMAP2", color = "Pseudotime") +
    theme(
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 16, hjust = 0.5)
    )


############# Fig2I - CytoTRACE2 Analysis for CD8+ T Cells

library(CytoTRACE2)
library(dplyr)
library(Seurat)

# Prepare input: expression matrix, cell type annotation, UMAP coordinates
expression_data <- as.data.frame(sce.cd8Tcell@assays$RNA@counts)
annotation <- data.frame(annotation = as.character(sce.cd8Tcell$CellType))
rownames(annotation) <- rownames(sce.cd8Tcell@meta.data)
umap <- as.data.frame(Embeddings(sce.cd8Tcell, reduction = "umap"))

# Run CytoTRACE2
cytotrace2_result <- cytotrace2(expression_data, species = "human")

# Merge CytoTRACE2 results with UMAP coordinates
cytotrace2_result$RowNames <- rownames(cytotrace2_result)
umap$RowNames <- rownames(umap)
merged_df <- left_join(cytotrace2_result, umap, by = "RowNames")
rownames(merged_df) <- merged_df$RowNames
merged_df$RowNames <- NULL

# Add CytoTRACE2 result to Seurat object
sce.cd8Tcell <- AddMetaData(sce.cd8Tcell, merged_df)

# Plot UMAP colored by CytoTRACE2_Relative
DimPlot(sce.cd8Tcell, group.by = "CytoTRACE2_Relative", label = TRUE, pt.size = 0.5) +
    scale_color_gradientn(colors = c("#FBE3E8", "#F8CCD9", "#F4B3C5", "#EF8BA6", "#E25778")) +
    labs(x = "UMAP1", y = "UMAP2", color = "CytoTRACE", title = NULL) +
    theme(
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 16, hjust = 0.5)
    )


############# Fig2J - Pseudotime Expression Dynamics (monocle3)

library(monocle3)
library(ggplot2)
library(dplyr)
library(tidyr)

# Define marker genes to track
Track_genes <- c("CCR7", "LEF1", "TCF7", "GNLY", "GZMA", "CX3CR1")

# Subset CDS object
Track_cds <- cds[rowData(cds)$gene_short_name %in% Track_genes, ]

# Extract expression matrix and pseudotime
expr_matrix <- t(exprs(Track_cds)[Track_genes, ])
pseudotime_vals <- pseudotime(Track_cds)

# Convert to long format for ggplot
plot_df <- as.data.frame(expr_matrix)
plot_df$pseudotime <- pseudotime_vals
plot_df_long <- pivot_longer(plot_df, cols = all_of(Track_genes), names_to = "gene", values_to = "expression")

# Plot smoothed gene expression over pseudotime
ggplot(plot_df_long, aes(x = pseudotime, y = expression, color = gene)) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
    scale_color_manual(values = c(
        "CCR7" = my_cols[13],
        "LEF1" = my_cols[12],
        "TCF7" = my_cols[11],
        "GNLY" = my_cols[10],
        "GZMA" = my_cols[9],
        "CX3CR1" = my_cols[8]
    )) +
    labs(x = "Pseudotime", y = "Gene Expression", color = "Gene", title = NULL) +
    theme_classic() +
    theme(
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
        plot.title = element_text(size = 16, hjust = 0.5)
    )

	
############# Fig2K - Differential analysis of T03 vs. T02 in ascites

library(org.Hs.eg.db)
library(clusterProfiler)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(stringr)
library(dplyr)

# Subset for CD8+ T cells from AF group
sce.Tcell2_af <- subset(sce.Tcell2, CellType2 %in% c('T02_CD8T', 'T03_CD8T') & group1 == 'AF')
Idents(sce.Tcell2_af) <- "CellType2"

# Differential expression
deg <- FindMarkers(sce.Tcell2_af, ident.1 = "T03_CD8T", ident.2 = "T02_CD8T", min.pct = 0.25, logfc.threshold = 0.25)
deg <- deg[order(deg$avg_log2FC, decreasing = TRUE), ]
deg_gene <- rownames(deg)

# KEGG enrichment
entrez_id <- mapIds(org.Hs.eg.db, keys = deg_gene, keytype = "SYMBOL", column = "ENTREZID")
kegg_res <- enrichKEGG(gene = entrez_id, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
kegg_res <- setReadable(kegg_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# Extract top enriched pathways
KEGG_top <- kegg_res@result[c(15,24,19,36,91,63,34,4,21,22), ]
KEGG_top$pathway <- factor(KEGG_top$Description, levels = rev(KEGG_top$Description))
KEGG_top$geneID2 <- sapply(strsplit(KEGG_top$geneID, "/"), function(x) paste(x[1:5], collapse = "/"))

# Plot
ggplot(KEGG_top, aes(x = Count, y = pathway, fill = p.adjust)) +
    geom_bar(stat = "identity", width = 0.6) +
    scale_fill_distiller(palette = "Spectral", direction = 1) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = expansion(mult = c(0.1, 0.1)), 
                     labels = function(labels) str_wrap(labels, width = 20)) +
    geom_text(aes(x = 0.1, label = Description), hjust = 0, size = 4.5) +
    geom_text(aes(x = 0.1, label = geneID2), hjust = 0, vjust = 4, size = 4, fontface = "italic") +
    labs(x = "Number of Genes", y = NULL, title = "Up-regulated in T03_CD8T") +
    theme_classic() +
    theme(
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        axis.ticks.y = element_blank(),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        panel.border = element_rect(color = "black", fill = NA, size = 1)
    )


############# Fig2L - CD8 T cell subtype signature survival analysis (TCGA & GSE15459)

library(COSG)
library(GSVA)
library(survival)
library(survminer)
library(clusterProfiler)
library(org.Hs.eg.db)

# Function: Run COSG, GSVA, and survival analysis plot for given data
plot_survival_signature <- function(sce, tpm_path, cl_path, gene_sig_name, plot_title) {
  # Set default assay and cell identities
  DefaultAssay(sce) <- "RNA"
  Idents(sce) <- sce$CellType2
  
  # Run COSG to get top 50 marker genes per subtype
  cosg_res <- cosg(sce, groups = "all", assay = "RNA", slot = "data", mu = 1,
                   remove_lowly_expressed = TRUE, expressed_pct = 0.1, n_genes_user = 50)
  
  # Prepare gene list for GSVA (convert to uppercase)
  gene_list <- lapply(as.list(cosg_res$names), toupper)
  
  # Load expression (TPM) and clinical data
  load(tpm_path) # loads object: tpm
  load(cl_path)  # loads object: cl
  
  # Calculate GSVA enrichment scores using zscore method
  es <- gsva(as.matrix(tpm), gene_list, method = "zscore", kcdf = "Gaussian")
  
  # Extract score for target gene signature
  cl$score <- as.numeric(es[gene_sig_name, ])
  
  # Determine optimal cutpoint to stratify survival groups
  sur.cut <- surv_cutpoint(cl, time = "OS.time", event = "OS", variables = "score")
  sur.cat <- surv_categorize(sur.cut)
  
  # Fit survival model
  sfit <- survfit(Surv(OS, status) ~ score, data = sur.cat)
  
  # Plot Kaplan-Meier curve with p-value and confidence intervals
  p <- ggsurvplot(
    sfit, data = sur.cat, pval = TRUE, conf.int = TRUE,
    pval.coord = c(5, 0.15), pval.size = 6,
    palette = c("#982b2b", "#0074b3"),
    legend = "top", legend.title = "Status",
    legend.labs = c(paste0("High (n = ", sum(sur.cat$score == "high"), ")"),
                    paste0("Low (n = ", sum(sur.cat$score == "low"), ")")),
    font.main = c(12, "bold"),
    font.x = c(14, "bold"),
    font.y = c(14, "bold"),
    font.tickslab = c(12),
    xlab = "Time (Months)",
    title = plot_title,
    risk.table = FALSE
  )
  
  p$plot <- p$plot + theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0.5),
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    legend.position = c(0.75, 0.85),
    legend.background = element_blank()
  )
  
  return(p)
}

# Run and plot for TCGA-STAD dataset
p_tcga <- plot_survival_signature(
  sce = sce.Tcell2,
  tpm_path = "scRNA-seq/GC_AF/data/TCGA_STAD_tpm.Rdata",
  cl_path = "scRNA-seq/GC_AF/data/TCGA_STAD_cl.Rdata",
  gene_sig_name = "T03_CD8T",
  plot_title = "TCGA-STAD"
)

# Run and plot for GSE15459 dataset
p_gse <- plot_survival_signature(
  sce = sce.Tcell2,
  tpm_path = "scRNA-seq/GC_AF/data/GSE15459_tpm.Rdata",
  cl_path = "scRNA-seq/GC_AF/data/GSE15459_cl.Rdata",
  gene_sig_name = "T03_CD8T",
  plot_title = "GSE15459"
)

# Display plots
print(p_tcga)
print(p_gse)

                     
