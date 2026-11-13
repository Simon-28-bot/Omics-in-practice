#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("limma")

#Initialize packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyr)
library(dplyr)
library(data.table)
library(limma)
library(Hmisc)
library(ggtext)
library(extrafont)
loadfonts()
library(svglite)
library(ggrepel)
library(DESeq2)
library(pheatmap)
library(org.Hs.eg.db)
library(RColorBrewer)


#Set seed and working directory
set.seed(42)
setwd('C:/School/Medical_epigenomics/Year_1/Q2/Omics_in_practice/Module_2')

#Load the gene-count matrix into R
data.raw <- fread("umi_counts.txt", sep="\t", header = T, data.table = F)
rownames(data.raw) <- data.raw[,1]
data.raw <- data.raw[,-1]

# Compute depth per sample
data.meta <- data.frame(samples =colnames(data.raw), depth = colSums(data.raw))

#Assign column names ranging from 1-24
colnames(data.raw) <- c('D3_R3','D7_R1','D8_R2',
                        'D4_R3', 'D7_R2', '0h_R2',
                        'D7_R3', 'D8_R1', 'D1_R1',
                        'D3_R1', 'D3_R2', 'D2_R3', 
                        '0h_R3', '6h_R3', 'D1_R2',
                        'D2_R1', 'D4_R2', '6h_R1',
                        'D1_R3', 'D2_R2', 'D4_R1',
                        'D8_R3', '0h_R1', '6h_R2')

#PCA############################################################################
#Normalizing using limma-voom
data.norm <- voom(counts = data.raw, normalize.method="quantile")$E

# PCA
pca <- prcomp(t(data.norm),   #Transpose columns and rows
              center = TRUE,  #Centers the data before applying PCA (subtracts the mean). 
              scale. = FALSE) #No scale can be applied due to a constant/zero column to unit variance
data.pca <- data.frame(pca$x)

# Add sample labels
data.pca$Sample <- colnames(data.raw)

#Explain variation
var_expl <- 100 * pca$sdev^2/sum(pca$sdev^2)

#Visualizing the Principal Components
timepoints <- c("0h", "6h", "D1", "D2", "D3", "D4", "D7", "D8")
tint_cols  <- c("#FCC6BB", "#FAA18F", "#F87C63", "#F54927",
                "#F4320B", "#C82909", "#9C2007", "#701705")

colour_map <- setNames(
  rep(tint_cols, times = 3),
  paste0(rep(timepoints, times = 3), "_R", rep(1:3, each = length(timepoints))))

pca.plot <- ggplot(data.pca, aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = Sample),
             shape = 1,
             stroke = 0.9,
             size = 3) +
  xlab(paste0("PC1 (", round(var_expl[1], 1), "%)")) +
  ylab(paste0("PC2 (", round(var_expl[2], 1), "%)")) +
  geom_text_repel(data = subset(data.pca, grepl("_R1$", Sample)),   # only the R1 samples
                  aes(label = sub("_R1$", "", Sample)),             # "0h_R1" -> "0h"
                  family = "Aptos",
                  box.padding = 0.5,
                  point.padding = 0.3,
                  max.overlaps = Inf) +
  scale_colour_manual(
    values = colour_map,
    breaks = paste0(timepoints, "_R3"),
    labels = timepoints) +
  geom_text_repel(data = lab.data,
                aes(label = sub("_R1$", "", Sample)),
                family = "Aptos",
                box.padding = 0.8,       # more space around the label
                point.padding = 0.5,     # more space from points
                force = 3,               # push labels further from points
                seed = 42,               # same layout every time
                max.overlaps = Inf,
                nudge_x = ifelse(lab.data$Sample == "D1_R1", 3, 0),   # shift only D1
                nudge_y = ifelse(lab.data$Sample == "D1_R1", 3, 0))
  theme(
    text                   = element_text(family = "Aptos"),
    axis.title.x           = element_markdown(size = 10), #x-axis label size
    axis.title.y           = element_markdown(size = 10), #y-axis label size
    axis.text.y            = element_markdown(size = 10, colour = "black"), #axis tick label size
    axis.text.x            = element_markdown(size = 10), #axis tick label size
    legend.position        = "right",
    #legend.position.inside = c(.88, .1),
    legend.text            = element_markdown(size = 10), #legend text size
    legend.title           = element_blank(),  #legend title size (if any)
    #plot.title             = element_text(face = "bold", hjust = 0.5, size = 21), # plot title size (if any)
    #plot.title.position    = "panel", # Title position (plot vs panel)
    axis.line              = element_line(colour = "black", linewidth = 0.5), #Axis line (if any)
    axis.ticks             = element_line(colour = "black", linewidth = 0.5), #Axis line tick size (if any)
    panel.background       = element_rect(fill = "transparent", colour = NA), #Panel background (if any)
    plot.background        = element_rect(fill = "transparent", colour = NA), #Plot background (if any)
    #panel.grid.major      = element_line(fill_alpha("black", 0.15)),
    legend.background      = element_rect(fill = "transparent", colour = NA), #Legend background (if any)
    legend.key             = element_rect(fill = "transparent", colour = NA)  #Legend background (if any)
  )
pca.plot

ggsave("pca_plot.svg", # Make sure to not use \ / : * ? " < > |.
       plot = pca.plot, 
       width = 14, 
       height = 11.5, 
       units = "cm", 
       dpi = 300)
ggsave("pca_plot.png", # Make sure to not use \ / : * ? " < > |.
       plot = pca.plot, 
       width = 14, 
       height = 11.5, 
       units = "cm", 
       dpi = 300)


#DESeq2#########################################################################
coldata <- data.frame(condition = factor(c('D3','D7','D8',
                                           'D4','D7','0h',
                                           'D7','D8','D1',
                                           'D3','D3','D2', 
                                           '0h','6h','D1',
                                           'D2','D4','6h',
                                           'D1','D2','D4',
                                           'D8','0h','6h')))
dds <- DESeqDataSetFromMatrix(countData = data.raw, 
                              colData = coldata, 
                              design = ~ condition)
dds <- DESeq(dds)
ntd <- normTransform(dds)
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)

#Volcano plot###################################################################
#To following was used as reference: https://www.biostars.org/p/282295/ 
controls = c('6h','D1','D2','D3','D4','D7','D8')

for (index in controls) {
  
  #Specifying time-points
  res <- results(dds, contrast = c("condition", "0h", index))
  topT <- na.omit(as.data.frame(res)) # Remove genes with NA p-values
  topT$gene <- rownames(topT) # Extract gene names from rownames
  
  #Converting Ensembl IDs in topT to official gene symbols
  clean_ensg <- sub("\\..*", "", topT$gene)
  topT$symbol <- mapIds(org.Hs.eg.db, keys = clean_ensg, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  topT$symbol <- ifelse(is.na(topT$symbol), topT$gene, topT$symbol) # Keep ENSG if no symbol exists
  
  # Classify genes as significantly upregulated, downregulated, or not significant
  topT$status <- "Not significant"
  topT$status[topT$padj < 0.05 & topT$log2FoldChange > 1]  <- "Upregulated"
  topT$status[topT$padj < 0.05 & topT$log2FoldChange < -1] <- "Downregulated"
  
  # Select top 10 upregulated and top 10 downregulated genes (ranked by padj)
  top_up <- head(subset(topT, status == "Upregulated")[order(subset(topT, status == "Upregulated")$padj), ], 10)
  top_down <- head(subset(topT, status == "Downregulated")[order(subset(topT, status == "Downregulated")$padj), ], 10)
  top_genes <- rbind(top_up, top_down)
  
  # Isolate target gene of interest (TFDP1) and prevent duplicate labeling if it's already in top_genes
  tfoi <- subset(topT, symbol == "TFDP1")
  top_genes <- subset(top_genes, symbol != "TFDP1")
  
  # Create the volcano plot (matching the print-friendly style of the PCA)
  volcano.plot <- ggplot(topT, aes(x = log2FoldChange, y = -log10(padj), color = status)) +
    geom_point(alpha = 0.5, size = 1.8) +
    geom_vline(xintercept = c(-1, 1), 
               linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), 
               linetype = "dashed") +
    geom_point(data = tfoi, # Highlight dot for TFDP1
               color = "black", 
               fill = "#FFD700", 
               shape = 21, 
               size = 2, 
               stroke = 1.2) + 
    geom_text_repel(data = top_genes,
                    aes(label = symbol),
                    color = "black",
                    size = 3,
                    box.padding = 0.5,
                    point.padding = 0.3,
                    min.segment.length = 0,
                    segment.color = "grey40",
                    segment.size = 0.3,
                    max.overlaps = Inf) +
    geom_label_repel(data = tfoi,
                     aes(label = symbol),
                     color = "black",
                     fontface = "bold",
                     size = 3,
                     box.padding = 1,
                     point.padding = 0.4,
                     min.segment.length = 0,
                     segment.color = "black",
                     segment.size = 0.5) +
    scale_color_manual(name = NULL, 
                       values = c("Not significant" = "grey70", 
                                  "Upregulated" = "#D95F02", 
                                  "Downregulated" = "#377EB8")) +
    xlab('log2(fold-change)') +
    ylab('-log10(adjusted p-value)') +
    ggtitle(paste0("Volcano plot: 0h vs ", index)) +
    theme(
      text                   = element_text(family = "Aptos"),
      axis.title.x           = element_markdown(size = 10), #x-axis label size
      axis.title.y           = element_markdown(size = 10), #y-axis label size
      axis.text.y            = element_markdown(size = 10, colour = "black"), #axis tick label size
      axis.text.x            = element_markdown(size = 10), #axis tick label size
      plot.title             = element_markdown(size = 12, hjust = 0.5),
      #plot.title.position    = "panel", 
      legend.position        = "right",
      #legend.position.inside = c(.88, .1),
      legend.text            = element_markdown(size = 10), #legend text size
      legend.title           = element_blank(),  #legend title size (if any)
      #plot.title             = element_text(face = "bold", hjust = 0.5, size = 21), # plot title size (if any)
      #plot.title.position    = "panel", # Title position (plot vs panel)
      axis.line              = element_line(colour = "black", linewidth = 0.5), #Axis line (if any)
      axis.ticks             = element_line(colour = "black", linewidth = 0.5), #Axis line tick size (if any)
      panel.background       = element_rect(fill = "transparent", colour = NA), #Panel background (if any)
      plot.background        = element_rect(fill = "transparent", colour = NA), #Plot background (if any)
      panel.grid.major       = element_line(fill_alpha("black", 0.15)),
      legend.background      = element_rect(fill = "transparent", colour = NA), #Legend background (if any)
      legend.key             = element_rect(fill = "transparent", colour = NA)  #Legend background (if any)
    )
  ggsave(paste0("0h_vs_",index, "_volcano_plot.svg"), # Make sure to not use \ / : * ? " < > |.
         plot = volcano.plot, 
         width = 15, 
         height = 11.5, 
         units = "cm", 
         dpi = 300)
  ggsave(paste0("0h_vs_",index, "_volcano_plot.png"), # Make sure to not use \ / : * ? " < > |.
         plot = volcano.plot, 
         width = 15, 
         height = 11.5, 
         units = "cm", 
         dpi = 300)
}



#Heatmap########################################################################
timepoints <- c("6h", "D1", "D2", "D3", "D4", "D7", "D8")
all_res <- do.call(rbind, lapply(timepoints, function(tp) {
  res <- results(dds, contrast = c("condition", "0h", tp))
  df <- na.omit(as.data.frame(res))
  df$gene <- rownames(df)
  df$comparison <- paste0("0h vs ", tp)
  return(df)
}))

all_res$comparison <- factor(all_res$comparison, levels = paste0("0h vs ", timepoints))

##### Creation of heatmaps with pheatmap #####
#Extract variance-stabilized normalized counts from your DESeq2 object
vsd_mat <- assay(vst(dds, blind = T))

#Reuse your existing symbol_map from the volcano plot step to name the rows
symbols <- symbol_map[sub("\\..*", "", rownames(vsd_mat))]
rownames(vsd_mat) <- ifelse(is.na(symbols), rownames(vsd_mat), symbols)


#Collect the top 30 DEGs (lowest padj) from EACH comparison against 0h
timepoints <- c("6h", "D1", "D2", "D3", "D4", "D7", "D8")

top_per_tp <- unlist(lapply(timepoints, function(tp) {
  res_tp <- na.omit(as.data.frame(results(dds, contrast = c("condition", "0h", tp))))
  res_tp$symbol <- symbol_map[sub("\\..*", "", rownames(res_tp))]
  res_tp$symbol <- ifelse(is.na(res_tp$symbol), rownames(res_tp), res_tp$symbol)
  head(res_tp[order(res_tp$padj), "symbol"], 30)
}))

# Combine all unique top DEGs across all time points and include TFDP1
selected_genes <- unique(c('ORC1', 'CDC6', 'TOP2B', 'TOP2A', 'TOP1',"TFDP1", 'RB1', 'RRM2', 'SOD1', 'ATP5MC3', 'CDK1', 'CDK2', 'CDK3', 'CDK4', 'CDK5', 'CDK6', 'CDK7', 'CDK8', 'CDK9', 'CDK10', 'CDK12', 'CDK13', 'CDK14', 'CDK15', 'CDK16', 'CDK17', 'CDK18', 'CDK19', 'CDK20', 'E2F1', 'E2F2', 'E2F3a', 'E2F3b', 'E2F4', 'E2F5', 'E2F6', 'E2F7', 'E2F8'))

# Match against your matrix
present <- intersect(selected_genes, rownames(vsd_mat))

# 4. Compute row z-scores so each gene is shown relative to its own average
heat_mat <- vsd_mat[present, , drop = FALSE]
heat_z <- t(scale(t(heat_mat)))
heat_z <- heat_z[complete.cases(heat_z), , drop = FALSE]

# 5. Sort columns chronologically from 0h to D8 (including R1, R2, R3 order)
time_order <- c("0h", "6h", "D1", "D2", "D3", "D4", "D7", "D8")
sample_time <- factor(sub("_.*", "", colnames(heat_z)), levels = time_order)
col_order <- order(sample_time, colnames(heat_z))

heat_z <- heat_z[, col_order]

# Create custom row labels where only TFDP1 is in bold
custom_row_labels <- ifelse(
  rownames(heat_z) == "TFDP1",
  expression(bold("TFDP1")),
  rownames(heat_z))

# Prepare color-coded time point annotation matching your PCA palette
ann_col <- data.frame(
  Timepoint = sample_time[col_order],
  row.names = colnames(heat_z))

ann_colors <- list(
  Timepoint = c("0h" = "#FCC6BB", "6h" = "#FAA18F", "D1" = "#F87C63", "D2" = "#F54927",
    "D3" = "#F4320B", "D4" = "#C82909", "D7" = "#9C2007", "D8" = "#701705")
)

# 6. Plot the clean, unified heatmap
heatmap.plot <- pheatmap(heat_z,
         color = colorRampPalette(c("#67a9cf", "#f7f7f7", "#ef8a62"))(100), # Blue - White - Red palette
         breaks = seq(-2, 2, length.out = 101),                            # Centers pure white strictly at 0
         border_color = NA, 
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         annotation_col = ann_col,
         annotation_colors = ann_colors,
         show_colnames = TRUE,
         labels_row = custom_row_labels,
         fontsize_row = 6.5, # Slightly smaller font so 60-80 unique gene names fit cleanly
         legend_breaks = c(-2, -1, 0, 1, 2),
         legend_labels = c("-2", "-1", "0", "1", "2"))

ggsave('1heatmap_plot.svg', # Make sure to not use \ / : * ? " < > |.
       plot = heatmap.plot, 
       width = 15, 
       height = 18, 
       units = "cm", 
       dpi = 300)
ggsave('heatmap_plot.png', # Make sure to not use \ / : * ? " < > |.
       plot = heatmap.plot, 
       width = 15, 
       height = 11.5, 
       units = "cm", 
       dpi = 300)























































