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


############# Fig4E - PCA plot of average expression of DC subclusters (MA vs PBMC)

library(Seurat)
library(ggplot2)

# Subset DCs
dc_sub <- subset(sce.Myeloid, CellType2 %in% c('M01', 'M02', 'M03', 'M04'))
dc_sub$group1 <- factor(dc_sub$group1, levels = c("MA", "PBMC"))
Idents(dc_sub) <- "group1"
DefaultAssay(dc_sub) <- "RNA"

# Compute average expression per group per celltype
avg_exp <- AverageExpression(dc_sub, add.ident = "CellType2")$RNA
avg_exp <- avg_exp[is.finite(rowSums(avg_exp)), ]

# PCA
pca_res <- prcomp(t(avg_exp))
pca_df <- as.data.frame(pca_res$x[, 1:2])
pca_df$Group <- rep(c("MA", "PBMC"), each = 4)
pca_df$CellType <- rep(c('M01', 'M02', 'M03', 'M04'), times = 2)

# Define plot colors and shapes
group_colors <- c("MA" = "#68A180", "PBMC" = "#3A6963")
cell_shapes <- c(15, 16, 17, 18)  # Square, circle, triangle, diamond

# Plot
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = CellType)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = group_colors) +
  scale_shape_manual(values = cell_shapes) +
  labs(x = "PC1", y = "PC2") +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12)
  ) +
  guides(shape = guide_legend(override.aes = list(size = 4)))


############# Fig4G - GSEA Dot Plot of M03 DC Subset (MA vs PBMC)

library(dplyr)
library(Seurat)
library(clusterProfiler)
library(GseaVis)
library(ggplot2)

# Subset M03 DC
m03_sub = subset(dc_sub ,CellType2 == "M03")

# Set identity class for comparison
Idents(m03_sub) <- "group1"

# Differential expression between MA and PBMC groups
dc_deg <- FindMarkers(m03_sub, ident.1 = "MA", ident.2 = "PBMC", 
                      min.pct = 0.25, logfc.threshold = 0.25)

# Prepare ranked gene list for GSEA
gsea_input <- dc_deg$avg_log2FC
names(gsea_input) <- rownames(dc_deg)
gsea_input <- sort(gsea_input, decreasing = TRUE)

# Load GO biological process gene sets
geneset <- read.gmt("scRNA-seq/GC_AF/data/c5.go.bp.v2025.1.Hs.symbols.gmt")

# Run GSEA analysis
gsea_res <- GSEA(gsea_input, TERM2GENE = geneset, verbose = FALSE,
                 pvalueCutoff = 0.1, pAdjustMethod = "BH")

# Define terms of interest (up- and down-regulated)
go_up <- c(
  "GOBP_MYELOID_LEUKOCYTE_ACTIVATION",
  "GOBP_MACROPHAGE_ACTIVATION",
  "GOBP_IMMUNE_RESPONSE_REGULATING_CELL_SURFACE_RECEPTOR_SIGNALING_PATHWAY",
  "GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY",
  "GOBP_ACTIVATION_OF_IMMUNE_RESPONSE"
)

go_down <- c(
  "GOBP_CYTOPLASMIC_TRANSLATION",
  "GOBP_TRANSLATION",
  "GOBP_RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS",
  "GOBP_OXIDATIVE_PHOSPHORYLATION",
  "GOBP_RNA_SPLICING"
)

# Filter GSEA results to selected terms
go_selected <- c(go_up, go_down)
gsea_res@result <- gsea_res@result[gsea_res@result$Description %in% go_selected, ]

# Create customized color palette for NES values
custom_colors <- colorRampPalette(c("#008bd0", "#eeeeee", "#fc4e00"))(100)

# Generate dot plot of selected GSEA terms
p <- ggplot(gsea_res@result, aes(x = 1, y = Description)) +
  geom_point(aes(size = abs(NES), color = NES)) +  # Size by |NES|, color by NES
  scale_color_gradientn(colors = custom_colors, limits = c(-5, 3), name = "NES") +
  scale_size_continuous(range = c(2, 8)) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    aspect.ratio = 1
  )

p


############# Fig4H - Functional Marker Heatmap of DC Subsets

library(Seurat)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(dplyr)

# Define selected genes related to DC maturation, migration, activation, and chemokines
marker_genes <- c('LAMP3','MARCKSL1','IL15','CD40','CCR7','FSCN1',
                  'SLCO5A1','CCL17','CCL19','CCL22')

# Calculate average expression per DC subtype (M01-M04)
avg_expr <- AverageExpression(dc_sub, features = marker_genes, 
                               group.by = 'CellType2', 
                               slot = 'data', assay = 'RNA')$RNA

# Rename columns for clarity
colnames(avg_expr) <- c('M01','M02','M03','M04')

# Z-score scaling by gene (row-wise)
expr_mat <- t(scale(t(avg_expr)))

# Set up color scale (Z-score from -2 to +2)
col_fun <- colorRamp2(c(-2, 0, 2), c("#0099CC", "white", "#CC0033"))

# Define functional group for annotation and row splitting
gene_groups <- c("Maturation", "Maturation",
                 "Activation", "Activation",
                 "Migration", "Migration", "Migration",
                 "Chemokines", "Chemokines", "Chemokines")

row_anno <- rowAnnotation(
  df = data.frame(Function = gene_groups),
  show_annotation_name = FALSE
)

# Plot heatmap
Heatmap(as.matrix(expr_mat),
        name = "Z-score",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 10, fontface = "italic"),
        column_names_rot = 45,
        column_names_gp = gpar(fontsize = 10),
        border = TRUE,
        rect_gp = gpar(col = "white", lwd = 1),
        col = col_fun,
        right_annotation = row_anno,
        row_split = gene_groups,
        column_gap = unit(2, "mm"))


############# Fig4I - Gene set Score in MA DC Subtypes

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

# Subset MA cells
dc_ma <- subset(dc_sub, subset = group1 == "MA")
Idents(dc_ma) <- "CellType2"

# Define gene sets for each module
gene_sets <- list(
  Antigen_presentation = c("AP1M2", "ACTR1B", "ACTR1A", "BCAP31", "PSME3", "PSMD14", "CLEC4M", 
                           "SEC24B", "IFI30", "SEC23A", "DCTN2", "CENPE", "DCTN6", "SEC24A", 
                           "KIF2C", "KIF3A", "DCTN3", "CHUK", "OSBPL1A", "AP2M1", "AP1S1", 
                           "AP2S1", "CLTA", "CLTC", "PSMB11", "AP1S3", "DYNLL2", "PSMA8", 
                           "CTSD", "CTSE", "CTSL", "CTSV", "CTSS", "CYBA", "CYBB", "AP2A1", 
                           "AP2A2", "AP1B1", "AP2B1", "ACE", "DCTN1", "AP1G1", "DYNC1H1", 
                           "DYNC1I1", "DYNC1I2", "DYNC1LI2", "DNM2", "FCER1G", "FCGR1A", 
                           "FCGR1B", "FCGR2B", "SEC31A", "KIFAP3", "PSME4", "KIF4A", "KIF26A", 
                           "KIF4B", "PYCARD", "RACGAP1", "PDIA3", "HFE", "CD209", "HLA-A", 
                           "HLA-B", "HLA-C", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", 
                           "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", 
                           "HLA-DQB2", "HLA-DRA", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", 
                           "HLA-DRB5", "HLA-E", "HLA-F", "HLA-G", "HLA-H", "MR1", "IKBKB", 
                           "ITGAV", "ITGB5", "ARF1", "KIF2A", "KIF3C", "KIF5A", "KLC1", 
                           "KIF11", "KIF22", "LAG3", "LNPEP", "NCF2", "NCF4", "CD207", 
                           "CLEC4A", "SAR1B", "DYNC1LI1", "DCTN4", "ERAP1", "TREM2", 
                           "MARCHF1", "TAPBPL", "ACTR10", "LGMN", "B2M", "PSMA1", "PSMA2", 
                           "PSMA3", "PSMA4", "PSMA5", "PSMA6", "PSMA7", "PSMB1", "PSMB2", 
                           "PSMB3", "PSMB4", "PSMB5", "PSMB6", "PSMB7", "PSMB8", "PSMB9", 
                           "PSMB10", "KIF15", "PSMC1", "PSMC2", "PSMC3", "PSMC4", "PSMC5", 
                           "PSMC6", "PSMD1", "PSMD2", "PSMD3", "PSMD4", "PSMD5", "PSMD7", 
                           "PSMD8", "PSMD9", "PSMD10", "PSMD11", "PSMD12", "PSMD13", 
                           "PSME1", "PSME2", "SEC13", "ERAP2", "SH3GL2", "KLC2", "NCF1", 
                           "SLC11A1", "SPTBN2", "TAP1", "TAP2", "TAPBP", "TRAF6", "RAB7A", 
                           "CALR", "KIF18A", "CANX", "CAPZA1", "CAPZA2", "CAPZB", "RILP", 
                           "DCTN5", "KIF2B", "IKBKG", "DYNLL1", "VAMP8", "CTSF", "SNAP23", 
                           "AP1S2", "AP1M1", "VAMP3", "CAPZA3", "KIF3B", "CD36", "PSMF1", 
                           "AC027237.1", "SEC22B", "SEC24C", "CD74", "PSMD6", "SEC24D"),
  
  Apoptosis = c("AKT3", "IRAK3", "CHP1", "CHUK", "CSF2RB", "DFFA", "DFFB", "ENDOG", 
                "AKT1", "AKT2", "ENDOD1", "PIK3R5", "APAF1", "BIRC2", "BIRC3", "XIAP", 
                "FAS", "IKBKB", "IL1A", "IL1B", "IL1R1", "IL1RAP", "FASLG", "IL3", 
                "IL3RA", "IRAK1", "IRAK2", "MYD88", "ATM", "NFKB1", "NFKBIA", "NGF", 
                "NTRK1", "IRAK4", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1", 
                "PIK3R2", "CYCS", "PPP3CA", "PPP3CB", "PPP3CC", "PPP3R1", "PPP3R2", 
                "PRKACA", "PRKACB", "PRKACG", "PRKAR1A", "PRKAR1B", "PRKAR2A", 
                "PRKAR2B", "PRKX", "BAD", "BAX", "BCL2", "RELA", "BCL2L1", "BID", 
                "CHP2", "TNF", "TNFRSF1A", "TP53", "TRAF2", "CAPN1", "CAPN2", 
                "CASP3", "CASP6", "CASP7", "CASP8", "CASP9", "CASP10", "PIK3R3", 
                "IKBKG", "TRADD", "RIPK1", "TNFSF10", "FADD", "TNFRSF10D", 
                "TNFRSF10C", "TNFRSF10B", "TNFRSF10A", "CFLAR", "MAP3K14", "AIFM1", 
                "EXOG"),
  
  Differentiation = c("DHRS2", "LILRB2", "CEBPB", "UBD", "BATF", "LILRB1", "BATF2", "CCR7", 
                     "ZBTB46", "CSF2", "AGER", "F2RL1", "FCGR2B", "FLT3", "GAS6", "GATA1", 
                     "TMEM176B", "HLA-G", "HMGB1", "RBPJ", "IL4", "IRF4", "LGALS9", 
                     "AC005840.1", "LYN", "MIR155", "MIR223", "TREM2", "TMEM176A", "BATF3", 
                     "AXL", "PRTN3", "PSEN1", "RELB", "CCL19", "AZI2", "SPI1", "TGFB1", 
                     "TGFBR2", "TRAF6", "TRPM2", "CAMK4", "DCSTAMP"),
  
  Immune_suppressive = c("CD274", "PDCD1LG2", "CD200", "IDO1", "EBI3", "IL4I1", "SOCS1", "SOCS2", "SOCS3")
)

# Filter gene sets to genes present in data
gene_sets <- lapply(gene_sets, function(genes) {
  genes[genes %in% rownames(dc_ma)]
})

# Calculate module scores
for (module_name in names(gene_sets)) {
  dc_ma <- AddModuleScore(dc_ma, features = list(gene_sets[[module_name]]), name = paste0(module_name, "_score"))
}

# Gather scores into a dataframe
scores_df <- data.frame(
  cell_type = dc_ma$CellType2,
  Antigen_presentation = dc_ma$Antigen_presentation_score1,
  Apoptosis = dc_ma$Apoptosis_score1,
  Differentiation = dc_ma$Differentiation_score1,
  Immune_suppressive = dc_ma$Immune_suppressive_score1
)

# Convert to long format for ggplot2
library(tidyr)
scores_long <- scores_df %>%
  pivot_longer(cols = -cell_type, names_to = "Module", values_to = "Score")

# Summarize mean and sd per cell type and module
summary_df <- scores_long %>%
  group_by(cell_type, Module) %>%
  summarise(
    mean = mean(Score),
    sd = sd(Score),
    .groups = "drop"
  )

# Join summary stats back
scores_long <- left_join(scores_long, summary_df, by = c("cell_type", "Module"))

# Define consistent colors for cell types (modify as needed)
cell_colors <- c('#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B')

# Plot function for one module
plot_module <- function(module_name) {
  ggplot(filter(scores_long, Module == module_name), 
         aes(x = cell_type, y = Score, fill = cell_type)) +
    geom_violin(trim = FALSE, color = "white") +
    geom_boxplot(width = 0.2, position = position_dodge(0.9)) +
    geom_point(aes(y = mean), pch = 19, size = 1.5, position = position_dodge(0.9)) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1, size = 0.5, alpha = 0.7, color = "black") +
    scale_fill_manual(values = cell_colors) +
    labs(title = module_name, y = "Score", x = NULL) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = "black", family = "Times"),
      axis.text.y = element_text(size = 12, family = "Times"),
      axis.title.y = element_text(size = 14, family = "Times"),
      axis.line = element_line(color = "black", size = 1),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    )
}

# Generate plots for all modules
plots <- lapply(names(gene_sets), plot_module)

# Combine plots into a 2x2 grid
cowplot::plot_grid(plotlist = plots, ncol = 2, nrow = 2, align = "v", axis = "tlbr")


############# Fig4L - M1/M2 Score in MA Mac Subtypes

mac_sub = subset(sce.Myeloid, subset = CellType%in%c("M08_Mac-IL1B","M09_Mac-C1QC","M10_Mac-PLTP","M11_Mac-SPP1"))
mac_ma = subset(mac_sub, subset = group1==c("MA"))

M1  <- c(
  "IL1B", "INFG", "CSF2", "CD80", "CD86", "IL1R1", "HLA-DR", "TLR2", "TLR4", 
  "IFNR2", "IFNR1", "FCGR1A", "CXCL9", "CXCL10", "CXCL11", "TNF", "IL6", 
  "IL12A", "IL12B", "IL23A", "IL1A", "NOS1", "NOS2", "ROS1", "SLC2A1", 
  "TFEB", "STAT1", "NFKB1", "STAT3", "IFR3", "IFR5", "IFR7", "SOCS3", "JUN"
)

M2 <- c(
  "IL4", "IL13", "CCL17", "CCL22", "CCL24", "IL10", "TGFB1", "CCL13", 
  "VEGFA", "EGF", "PDGFA", "PDGFB", "MMP9", "VEGFB", "VEGFC", "VEGFD", 
  "TGFB2", "TGFB3", "MMP14", "MMP19", "SEPP1", "ARG1", "IDO2", "CD163", 
  "MSR1", "STAB1", "MARCO", "CD36", "FCGR2A", "IL1R2", "IL4R", "CD274", 
  "PDCD1LG2", "PDCD1", "SIRPA", "SIGLEC10", "S100A9", "LILRB1", "LILRB2", 
  "VTCN1", "TEK", "TREM1", "TREM2", "IL1RN", "STAT6", "MAFB", "IFR4", 
  "SOCS1", "CHI3L1", "AXL"
)

### The subsequent code is similar to Fig. 4I.


############# Fig4N - Gene set Score in MA Mac Subtypes

Inflammatory  <- c(
  "ABCA1", "ABI1", "ACVR1B", "ACVR2A", "ADGRE1", "ADM", "ADORA2B", "ADRM1", "AHR", "APLNR", 
  "AQP9", "ATP2A2", "ATP2B1", "ATP2C1", "AXL", "BDKRB1", "BEST1", "BST2", "BTG2", "C3AR1", 
  "C5AR1", "CALCRL", "CCL17", "CCL2", "CCL20", "CCL22", "CCL24", "CCL5", "CCL7", "CCR7", 
  "CCRL2", "CD14", "CD40", "CD48", "CD55", "CD69", "CD70", "CD82", "CDKN1A", "CHST2", 
  "CLEC5A", "CMKLR1", "CSF1", "CSF3", "CSF3R", "CX3CL1", "CXCL10", "CXCL11", "CXCL6", "CXCL8", 
  "CXCL9", "CXCR6", "CYBB", "DCBLD2", "EBI3", "EDN1", "EIF2AK2", "EMP3", "EREG", "F3", "FFAR2", 
  "FPR1", "FZD5", "GABBR1", "GCH1", "GNA15", "GNAI3", "GP1BA", "GPC3", "GPR132", "GPR183", 
  "HAS2", "HBEGF", "HIF1A", "HPN", "HRH1", "ICAM1", "ICAM4", "ICOSLG", "IFITM1", "IFNAR1", 
  "IFNGR2", "IL10", "IL10RA", "IL12B", "IL15", "IL15RA", "IL18", "IL18R1", "IL18RAP", "IL1A", 
  "IL1B", "IL1R1", "IL2RB", "IL4R", "IL6", "IL7R", "INHBA", "IRAK2", "IRF1", "IRF7", "ITGA5", 
  "ITGB3", "ITGB8", "KCNA3", "KCNJ2", "KCNMB2", "KIF1B", "KLF6", "LAMP3", "LCK", "LCP2", 
  "LDLR", "LIF", "LPAR1", "LTA", "LY6E", "LYN", "MARCO", "MEFV", "MEP1A", "MET", "MMP14", 
  "MSR1", "MXD1", "MYC", "NAMPT", "NDP", "NFKB1", "NFKBIA", "NLRP3", "NMI", "NMUR1", "NOD2", 
  "NPFFR2", "OLR1", "OPRK1", "OSM", "OSMR", "P2RX4", "P2RX7", "P2RY2", "PCDH7", "PDE4B", 
  "PDPN", "PIK3R5", "PLAUR", "PROK2", "PSEN1", "PTAFR", "PTGER2", "PTGER4", "PTGIR", "PTPRE", 
  "PVR", "RAF1", "RASGRP1", "RELA", "RGS1", "RGS16", "RHOG", "RIPK2", "RNF144B", "ROS1", 
  "RTP4", "SCARF1", "SCN1B", "SELE", "SELENOS", "SELL", "SEMA4D", "SERPINE1", "SGMS2", 
  "SLAMF1", "SLC11A2", "SLC1A2", "SLC28A2", "SLC31A1", "SLC31A2", "SLC4A4", "SLC7A1", 
  "SLC7A2", "SPHK1", "SRI", "STAB1", "TACR1", "TACR3", "TAPBP", "TIMP1", "TLR1", "TLR2", 
  "TLR3", "TNFAIP6", "TNFRSF1B", "TNFRSF9", "TNFSF10", "TNFSF15", "TNFSF9", "TPBG", "VIP"
)

Angiogenesis <- c(
  "APP", "CCND2", "CXCL6", "ITGAV", "JAG1", "OLR1", 
  "PTK2", "SLCO2A1", "SPP1", "THBD", "TIMP1", 
  "VAV2", "VCAN", "VEGFA"
)


Hypoxia <- c(
  "ADM", "ADORA2B", "AK4", "AKAP12", "ALDOA", "ALDOB", "ALDOC", "AMPD3", "ANGPTL4", 
  "ANKZF1", "ANXA2", "ATF3", "ATP7A", "B3GALT6", "B4GALNT2", "BCAN", "BCL2", "BGN", 
  "BHLHE40", "BNIP3L", "BRS3", "BTG1", "CA12", "CASP6", "CAV1", "CCNG2", "CCRN4L", 
  "CDKN1A", "CDKN1B", "CDKN1C", "CHST2", "CHST3", "CITED2", "COL5A1", "CP", "CSRP2", 
  "CTGF", "CXCR4", "CXCR7", "CYR61", "DCN", "DDIT3", "DDIT4", "DPYSL4", "DTNA", 
  "DUSP1", "EDN2", "EFNA1", "EFNA3", "EGFR", "ENO1", "ENO2", "ENO3", "ERO1L", 
  "ERRFI1", "ETS1", "EXT1", "F3", "FAM162A", "FBP1", "FOS", "FOSL2", "FOXO3", 
  "GAA", "GALK1", "GAPDH", "GAPDHS", "GBE1", "GCK", "GCNT2", "GLRX", "GPC1", 
  "GPC3", "GPC4", "GPI", "GRHPR", "GYS1", "HAS1", "HDLBP", "HEXA", "HK1", "HK2", 
  "HMOX1", "HOXB9", "HS3ST1", "HSPA5", "IDS", "IER3", "IGFBP1", "IGFBP3", "IL6", 
  "ILVBL", "INHA", "IRS2", "ISG20", "JMJD6", "JUN", "KDELR3", "KDM3A", "KIF5A", 
  "KLF6", "KLF7", "KLHL24", "LALBA", "LARGE", "LDHA", "LDHC", "LOX", "LXN", 
  "MAFF", "MAP3K1", "MIF", "MT1E", "MT2A", "MXI1", "MYH9", "NAGK", "NCAN", 
  "NDRG1", "NDST1", "NDST2", "NEDD4L", "NFIL3", "NR3C1", "P4HA1", "P4HA2", "PAM", 
  "PCK1", "PDGFB", "PDK1", "PDK3", "PFKFB3", "PFKL", "PFKP", "PGAM2", "PGF", 
  "PGK1", "PGM1", "PGM2", "PHKG1", "PIM1", "PKLR", "PKP1", "PLAC8", "PLAUR", 
  "PLIN2", "PNRC1", "PPARGC1A", "PPFIA4", "PPP1R15A", "PPP1R3C", "PRDX5", 
  "PRKCA", "PRKCDBP", "PTRF", "PYGM", "RBPJ", "RORA", "RRAGD", "S100A4", "SAP30", 
  "SCARB1", "SDC2", "SDC3", "SDC4"
)

Epithelial_Mesenchymal_Transition <- c(
  "ABI3BP", "ACTA2", "ADAM12", "BGN", "BMP1", "CADM1", "CALD1", "CAP2", "CCN1", "CCN2", 
  "CDH11", "CDH2", "CDH6", "COMP", "CTHRC1", "CXCL1", "CXCL12", "CXCL8", "DAB2", "DCN", 
  "DKK1", "DPYSL3", "DST", "ECM2", "EFEMP2", "ELN", "FAP", "FBLN1", "FBLN2", "FBLN", 
  "FBN1", "FBN2", "FERMT2", "FGF2", "FMOD", "FSTL1", "FSTL3", "FUCA1", "FZD8", "GADD4A", 
  "GADD4B", "GAS1", "GEM", "GJA1", "GLIPR1", "GPC1", "HTRA1", "ID2", "IGFBP3", "IGFBP4", 
  "ITGAV", "JUN", "LAMA2", "LAMA3", "LAMC2", "LOX", "LOXL1", "LOXL2", "LRP1", "LRRC1", 
  "LUM", "MAGEE1", "MATN2", "MATN3", "MCM7", "MEST", "MFAP", "MGP", "MMP14", "MMP2", 
  "MMP3", "MSX1", "MXRA", "MYL9", "MYLK", "NID2", "NNMT", "NTM", "P3H1", "PCOLCE", 
  "PDGFRB", "PFN2", "PMEPA1", "POSTN", "PRRX1", "PRSS2", "PTHLH", "RGS4", "RHOB", "SAT1", 
  "SERPINE2", "SFRP1", "SFRP4", "SGCB", "SGCD", "SGCG", "SLIT2", "SLIT3", "SNAI2", 
  "SPARC", "SPOCK1", "TAGLN", "TFPI2", "TGFB1", "TGFBR3", "THBS2", "THY1", "TIMP3", 
  "TNC", "TNFAIP3", "TPM1", "TPM2", "VCAM1", "VEGFA", "VEGFC", "WIPF1"
)

Extracellular_Matrix_Signature <- c(
  "SDC4", "PLXDC2", "CLEC2D", "LGALS3", "ANXA1", "PLXNA1", "SDC1", "SFTPA2", "SEMA7A", 
  "PLXNB2", "ANXA7", "LGALS9C", "PLXNA3", "SEMA3C", "LGALS8", "LGALS1", "LGALS9", 
  "ANXA11", "ANXA6", "PLXNC1", "CLEC5A", "GREM1", "ELFN1", "CLEC11A", "ANXA4", 
  "ANXA2", "ANXA5", "GPC4", "PLXND1"
)

Phagocytosis <- c(
  "AKT3", "PLA2G4B", "DNM1L", "WASF2", "VAV3", "ARPC1A", "CFL2", "WASF3", "PLA2G4E", 
  "CRK", "CRKL", "DNM1", "DNM2", "DOCK2", "PIKFYVE", "AKT1", "AKT2", "FCGR1A", 
  "FCGR2A", "FCGR2B", "FCGR3A", "PIP5K1C", "PIK3R5", "DNM3", "AMPH", "PLA2G4D", 
  "HCK", "INPP5D", "LIMK2", "LYN", "MARCKS", "MYO10", "PAK1", "ASAP1", "PIK3CA", 
  "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2", "PLA2G4A", "PLCG1", "PLCG2", 
  "PLD1", "PLD2", "PRKCA", "PRKCG", "MAPK1", "MAPK3", "SPHK2", "PTPRC", "RAC1", 
  "RAF1", "RPS6KB1", "NCF1", "SYK", "VAV1", "VAV2", "WAS", "PIP5K1A", "PIP5K1B", 
  "PIP4K2B", "PLA2G6", "PIK3R3", "SCIN", "PLPP1", "ASAP2", "WASL", "FCGR2C", "GAB2"
)

### The subsequent code is similar to Fig. 4I.


############# Fig4R-S: Prepare CellPhoneDB inputs from Epi + Macrophages

library(Seurat)

# Load epithelial and macrophage subsets
epi <- readRDS("scRNA-seq/GC_AF/data/sce.epi.rds")
epi <- subset(epi, CellType == "Epi")
epi <- RenameCells(epi, add.cell.id = "epi")

mac_ma <- RenameCells(mac_ma, add.cell.id = "Mac")
DefaultAssay(mac_ma) <- "RNA"
Idents(mac_ma) <- mac_ma$CellType

# Merge epithelial and macrophage cells
mac_epi <- merge(epi, mac_ma, merge.data = TRUE)
DefaultAssay(mac_epi) <- "RNA"

# Export expression and metadata
write.table(mac_epi@assays$RNA@counts, 
            "scRNA-seq/GC_AF/data/mac_epi_counts.txt", 
            sep = '\t', quote = FALSE)

meta <- data.frame(Cell = colnames(mac_epi), CellType = mac_epi$CellType)
write.table(meta, 
            "scRNA-seq/GC_AF/data/mac_epi_meta.txt", 
            row.names = FALSE, sep = '\t', quote = FALSE)

# Fig4R-S: Run CellPhoneDB statistical analysis

import os
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

os.chdir('scRNA-seq/GC_AF/data')

cpdb_results = cpdb_statistical_analysis_method.call(
    cpdb_file_path='db/v5/cellphonedb.zip',
    meta_file_path='mac_epi_meta.txt',
    counts_file_path='mac_epi_counts.txt',
    counts_data='hgnc_symbol',
    output_path='results/method2_mac_epi',
    threads=20,
    iterations=1000,
    threshold=0.1,
    pvalue=0.05,
    score_interactions=True,
    result_precision=3,
    debug_seed=42,
    subsampling=False,
    separator='|'
)

# Load CellPhoneDB output
library(ktplots)
library(CellChat)
library(tidyr)

pvals <- read.delim("scRNA-seq/GC_AF/data/results/method2_mac_epi/statistical_analysis_pvalues_06_22_2025_152705.txt", check.names = FALSE)
means <- read.delim("scRNA-seq/GC_AF/data/results/method2_mac_epi/statistical_analysis_means_06_22_2025_152705.txt", check.names = FALSE)

# Summarize significant interaction counts
interaction_counts <- do.call(rbind, lapply(14:ncol(pvals), function(i) {
  pair <- strsplit(colnames(pvals)[i], '\\|')[[1]]
  sig_count <- sum(pvals[, i] < 0.05)
  c(Source = pair[1], Target = pair[2], Count = sig_count)
}))

interaction_counts <- as.data.frame(interaction_counts)
interaction_counts$Count <- as.numeric(interaction_counts$Count) / 100  # Normalize

# Pivot to matrix
interaction_matrix <- spread(interaction_counts, Target, Count)
rownames(interaction_matrix) <- interaction_matrix$Source
interaction_matrix <- as.matrix(interaction_matrix[, -1])

# Circle plot for interaction count
netVisual_circle(interaction_matrix, weight.scale = TRUE)

# Plot selected ligand-receptor pairs
plot_cpdb(
  scdata = Myeloid_bclust_epi,
  cell_type1 = "Epi",
  cell_type2 = "M08_Mac-IL1B|M09_Mac-C1QC|M10_Mac-PLTP|M11_Mac-SPP1",
  celltype_key = "CellType",
  genes = c("CD47", "TREM1", "TREM2", "C3", "C5"),
  means = means,
  pvals = pvals
)

