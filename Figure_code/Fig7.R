############# Fig7B - Feature selection using mRMR on PRJEB25780 dataset

library(limma)
library(mRMRe)

# Load response label and expression matrix
STAD_PRJEB25780_cl  <- readRDS("scRNA-seq/GC_AF/data/STAD_PRJEB25780_cl.Rds")
STAD_PRJEB25780_exp <- readRDS("scRNA-seq/GC_AF/data/STAD_PRJEB25780_exp.Rds")

# Filter genes with low variance (sd <= 0.3)
expr_filt <- STAD_exp[apply(STAD_PRJEB25780_exp, 1, function(x) sd(x) > 0.3), ]

# Create design matrix for differential expression analysis (R = response, N = non-response)
group       <- as.factor(STAD_PRJEB25780_cl$response_NR)
design      <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Fit linear model using limma
fit   <- lmFit(expr_filt, design)
cont  <- makeContrasts(R - N, levels = design)
fit2  <- contrasts.fit(fit, cont)
fit2  <- eBayes(fit2)
deg   <- topTable(fit2, adjust = "fdr", number = Inf)

# Select DEGs with adjusted P < 0.05 and |logFC| > 1
filtered_genes <- rownames(deg[deg$adj.P.Val < 0.05 & abs(deg$logFC) > 1, ])

# Intersect with cluster-specific markers (top 1000 for a specific cluster)
symbols_list   <- epi.markers %>% group_by(cluster) %>% top_n(1000, avg_log2FC)
symbols_list   <- split(symbols_list$gene, symbols_list$cluster)
geneset2       <- symbols_list[[2]]
filtered_genes <- intersect(geneset2, filtered_genes)

# Prepare data for mRMR (transpose to samples × genes, append class label)
X        <- t(STAD_PRJEB25780_exp[filtered_genes, ])
y        <- as.numeric(as.factor(STAD_PRJEB25780_cl$response_NR))
df_mrmr  <- as.data.frame(cbind(X, Class = y))

# Convert to mRMRe data format
data_mrmr    <- mRMR.data(data = df_mrmr)
target_index <- ncol(df_mrmr)

# Run mRMR feature selection, selecting top 10 informative genes
feat_mrmr    <- mRMR.classic(
  data           = data_mrmr,
  target_indices = target_index,
  feature_count  = 10
)

# Extract selected gene names
top_rf_genes <- featureNames(data_mrmr)[solutions(feat_mrmr)[[1]]]


############# Fig7C - TIDE score comparison between risk groups

library(ggpubr)
library(stringr)
library(cowplot)

# Load TIDE result and define risk group
tide <- read.csv("scRNA-seq/GC_AF/data/TCGA_STAD_TIDE_result.csv", row.names = 1, check.names = FALSE)
tide$Risk <- factor(ifelse(str_sub(rownames(tide), 1, 1) == 'L', 'Low', 'High'), levels = c('High', 'Low'))

# Define comparison
my_comparisons <- list(c("Low", "High"))
risk_palette <- c("#E57164", "#A184BC")

# Custom plotting function
plot_violin <- function(df, y_col, y_pos, y_label = NULL) {
  ggviolin(df, x = 'Risk', y = y_col, fill = 'Risk',
           palette = risk_palette, add = 'boxplot',
           add.params = list(fill = "white")) +
    geom_signif(comparisons = my_comparisons, map_signif_level = TRUE,
                y_position = y_pos, textsize = 4, vjust = 0.5, tip_length = 0.02) +
    theme(legend.position = "none", axis.title.x = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 1)) +
    ylab(y_label)
}

# Rename MSI column
colnames(tide)[5] <- "MSI"

# Generate violin plots
pp1 <- plot_violin(tide, "TIDE", c(4, 6))
pp2 <- plot_violin(tide, "Dysfunction", c(3, 6))
pp3 <- plot_violin(tide, "Exclusion", c(4, 6))
pp4 <- plot_violin(tide, "MSI", c(1.25, 2))

# Combine plots
cowplot::plot_grid(pp1, pp2, pp3, pp4, ncol = 2, align = 'v', axis = 'tb')


############# Fig7D - Proportion of responders in each risk group

library(dplyr)
library(ggplot2)

# Calculate response percentage by risk group
dat <- tide %>%
  count(Risk, Responder) %>%
  group_by(Risk) %>%
  mutate(n = n / sum(n)) %>%
  ungroup() %>%
  mutate(Responder = factor(Responder, levels = c("False", "True")))

# Perform Chi-square test
ka <- xtabs(~Responder + Risk, data = tide)
chisq.test(ka)

# Plot response proportion bar plot
ggplot(dat, aes(x = Risk, y = n, fill = Responder)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#A184BC", "#E57164")) +
  geom_label(aes(label = scales::percent(n)),
             color = "white", size = 4, label.size = 0,
             show.legend = FALSE,
             position = position_fill(vjust = 0.5)) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    axis.ticks = element_line(color = "black")
  ) +
  ylab("Percentage") +
  annotate("text", label = "****", x = 1.5, y = 1.05, size = 8) +
  annotate("segment", x = 1, xend = 2, y = 1.04, yend = 1.04)


############# Fig7F - SubMap similarity score between risk groups and responders

library(pheatmap)

# Load SubMap result matrix
df <- read.table("scRNA-seq/GC_AF/data/STAD_PRJEB25780_Risk_Responder_submap.txt", 
                 header = TRUE, row.names = 1, sep = "\t")

# Plot heatmap with annotation and custom color
pheatmap(df,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         number_format = "%.3f",
         number_color = "black",
         border_color = "white",
         cellwidth = 30,
         cellheight = 30,
         fontsize_number = 9,
         color = rev(c("#3b374c", "#44598e", "#64a0c0", "#7ec4b7", "#deebcd")),
         name = "Statistic",
         annotation_row = annotation_row,
         annotation_colors = annotation_colors)


############# Fig7J - Signature score comparison between responders and non-responders

library(ggplot2)
library(ggpubr)

ggplot(STAD_PRJEB25780_cl, aes(x = response_NR, y = signature_score, fill = response_NR)) +
  geom_boxplot(width = 0.6, position = position_dodge(0.9)) +
  stat_compare_means(method = "wilcox.test", label = "p.signif") +
  scale_fill_manual(values = c("#A184BC", "#E57164")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 16),
    axis.line = element_line(color = "black", size = 0.8),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.position = "none"
  ) +
  ylab("Signature score") +
  xlab("") +
  ggtitle("PRJEB25780")


############# Fig7M - ROC curve of signature score for predicting response

library(pROC)

# Generate ROC curve based on signature score and response
roc_obj <- roc(STAD_PRJEB25780_cl$response_NR, STAD_PRJEB25780_cl$signature_score)

# Plot ROC curve
plot.roc(
  roc_obj,
  axes = TRUE,                        # Show axes
  legacy.axes = TRUE,                # Use legacy axes format (0-1)
  col = "steelblue",                 # Curve color
  lty = 1,                            # Line type
  lwd = 3,                            # Line width
  identity = TRUE,                   # Show diagonal line
  identity.col = "grey60",           # Diagonal line color
  identity.lty = 2,                  # Diagonal line type
  identity.lwd = 2,                  # Diagonal line width
  print.auc = TRUE,                  # Display AUC value
  print.auc.pattern = "AUC = 0.929", # AUC label format
  print.auc.x = 0.48,                # AUC label X position
  print.auc.y = 0.13,                # AUC label Y position
  print.auc.cex = 1.2,               # AUC label size
  print.auc.col = "black",           # AUC label color
  auc.polygon = TRUE,                # Fill area under ROC curve
  auc.polygon.col = "skyblue",      
  auc.polygon.border = "darkred",   
  max.auc.polygon = TRUE,           # Highlight max area under curve
  max.auc.polygon.col = "WhiteSmoke",
  max.auc.polygon.lty = 1,
  print.thres = FALSE                # Do not print threshold values
)
                        
                            
