############# Fig6J - Dynamic gene expression along pseudotime trajectory

library(SCP)

# Infer pseudotime trajectory using Slingshot based on UMAP and CellType clusters
epi <- RunSlingshot(
  srt        = epi,
  group.by   = "CellType",
  reduction  = "umap"
)

# Perform differential expression test between CellTypes
epi <- RunDEtest(
  srt          = epi,
  group_by     = "CellType",
  fc.threshold = 0.25,
  only.pos     = FALSE
)

# Filter DEGs with avg_log2FC > 1 and adjusted p-value < 0.05
DEGs <- epi@tools$DEtest_CellType$AllMarkers_wilcox
DEGs <- DEGs[with(DEGs, avg_log2FC > 1 & p_val_adj < 0.05), ]

# Annotate genes with transcription factors (TF) and surface markers (CSPA)
epi <- AnnotateFeatures(
  srt     = epi,
  species = "Homo_sapiens",
  db      = c("TF", "CSPA")
)

# Identify dynamic features along pseudotime trajectory
epi <- RunDynamicFeatures(
  srt         = epi,
  lineages    = c("Lineage"),
  n_candidates = 200
)

# Plot dynamic gene heatmap with annotations
ht <- DynamicHeatmap(
  srt                       = epi2,
  lineages                 = c("Lineage"),
  use_fitted               = TRUE,
  n_split                  = 6,
  reverse_ht               = "Lineage",
  species                  = "Homo_sapiens",
  db                      = "GO_BP",
  anno_terms               = TRUE,
  anno_keys                = TRUE,
  anno_features            = TRUE,
  heatmap_palette          = "viridis",
  cell_annotation          = "CellType",
  separate_annotation      = list("CellType", c("PTPRC", "EPCAM")),
  separate_annotation_palette = c("Paired", "Set1"),
  feature_annotation       = c("TF", "CSPA"),
  feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  pseudotime_label         = 25,
  pseudotime_label_color   = "red",
  height                   = 5,
  width                    = 2
)

# Display the final heatmap
print(ht$plot)

