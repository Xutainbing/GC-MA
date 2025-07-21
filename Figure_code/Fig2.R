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
    theme(
        panel.border = element_rect(fill = NA, color = "black", linewidth = 2),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        text = element_text(size = 20),
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 12, face = "italic")
    ) +
    guides(colour = guide_legend(override.aes = list(size = 4)))


