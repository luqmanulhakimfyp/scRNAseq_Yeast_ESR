# ==============================================================================
# PROJECT: Non-Parametric Characterization of Single-Cell Gene Expression 
# Heterogeneity and the Environmental Stress Response in S. cerevisiae
# AUTHOR: [Your Library Card Number]
# DESCRIPTION: Master analysis script for non-parametric pipeline and visualizations.
# ==============================================================================

# 1. ENVIRONMENT SETUP & DATA IMPORT
library(Matrix)
library(ggplot2)

raw_matrix <- readMM("E-GEOD-102475.aggregated_filtered_normalised_counts.mtx")
gene_names <- readLines("E-GEOD-102475.aggregated_filtered_normalised_counts.mtx_rows")
gene_names <- sapply(strsplit(gene_names, "\t"), `[`, 2)
rownames(raw_matrix) <- gene_names
full_data <- as.matrix(raw_matrix)

salt_group <- full_data[ , 1:80]
control_group <- full_data[ , 81:163]

message("\n[ INFO ] DATASET LOADED SUCCESSFULLY")
message(paste("[ INFO ] Total Genes Analyzed:", nrow(full_data)))

# 2. NON-PARAMETRIC PIPELINE (Interdecile Range)
analyze_group <- function(data_matrix) { medians <- apply(data_matrix, 1, median); p90 <- apply(data_matrix, 1, function(x) quantile(x, 0.90)); p10 <- apply(data_matrix, 1, function(x) quantile(x, 0.10)); idr <- p90 - p10; disp <- idr / (medians + 1); is_constant <- (medians >= 10) & (disp < 1.0); is_bursting <- (medians < 10) & (p90 > 400); return(list(constant = names(medians[is_constant]), bursting = names(medians[is_bursting]))) }

control_results <- analyze_group(control_group)
salt_results <- analyze_group(salt_group)

message("\n[ INFO ] NON-PARAMETRIC CLASSIFICATION COMPLETE")
message(paste("[ DATA ] Control - Constant Genes:", length(control_results$constant)))
message(paste("[ DATA ] Salt - Constant Genes:", length(salt_results$constant)))

# 3. CROSS-CONDITION INTERSECTION ANALYSIS
core_constant <- intersect(control_results$constant, salt_results$constant)
lost_stability <- setdiff(control_results$constant, salt_results$constant)
gained_stability <- setdiff(salt_results$constant, control_results$constant)

message("\n[ INFO ] CROSS-CONDITION INTERSECTION ANALYSIS")
message(paste("[ DATA ] Core Constant (Stable in both):", length(core_constant)))
message(paste("[ DATA ] Lost Stability (Stable -> Disrupted):", length(lost_stability)))

# 4. FOLD-CHANGE ANALYSIS
ctrl_medians <- apply(control_group, 1, median)
salt_medians <- apply(salt_group, 1, median)
fold_change <- (salt_medians + 1) / (ctrl_medians + 1)
sorted_fc <- sort(fold_change)
top_repressed <- head(names(sorted_fc), 100)
top_induced <- tail(names(sorted_fc), 100)

message("\n[ INFO ] FOLD CHANGE ANALYSIS COMPLETE")

# 5. DATA VISUALIZATION (ggplot2)

# Figure 3A: Baseline Heterogeneity Scatterplot
ctrl_p90 <- apply(control_group, 1, function(x) quantile(x, 0.90))
ctrl_p10 <- apply(control_group, 1, function(x) quantile(x, 0.10))
disp <- (ctrl_p90 - ctrl_p10) / (ctrl_medians + 1)
plot_data <- data.frame(Median = ctrl_medians, Dispersion = disp, Category = "Unclassified Noise")
plot_data$Category[plot_data$Median >= 10 & plot_data$Dispersion < 1.0] <- "Constant (Housekeeping)"
plot_data$Category[plot_data$Median < 10 & ctrl_p90 > 400] <- "Strict Bursting (Cell Cycle)"
plot_data$Category <- factor(plot_data$Category, levels = c("Unclassified Noise", "Constant (Housekeeping)", "Strict Bursting (Cell Cycle)"))
plot_data <- plot_data[order(plot_data$Category), ]

ggplot(plot_data, aes(x = Median + 1, y = Dispersion + 0.01, color = Category)) + geom_point(alpha = 0.8, size = 1.5) + scale_color_manual(values = c("Unclassified Noise" = "gray85", "Constant (Housekeeping)" = "steelblue", "Strict Bursting (Cell Cycle)" = "indianred")) + scale_x_log10(breaks = c(1, 11, 101, 1001), labels = c("0", "10", "100", "1000")) + scale_y_log10() + theme_minimal() + labs(title = "Baseline Gene Expression Heterogeneity", x = "Median Normalized Expression (Log10 Scale)", y = "Relative Dispersion (Log10 Scale)") + theme(text = element_text(size = 14), legend.position = "bottom", legend.title = element_blank())
ggsave("Figure_3A_Scatterplot_Log10.png", width = 8, height = 6, dpi = 300)

# Figure 3B: Optimal Constant Genes
df_baseline <- data.frame(term_name = c("Ribosome & Translation", "Ubiquitin / Isopeptide Bond", "Cytosolic Small Ribosomal Subunit", "rRNA binding", "Ribosome biogenesis"), adjusted_p_value = c(2.13e-106, 1.05e-51, 6.14e-48, 6.84e-11, 4.71e-5))
df_baseline$log_p <- -log10(df_baseline$adjusted_p_value)
ggplot(df_baseline, aes(x = reorder(term_name, log_p), y = log_p)) + geom_bar(stat = "identity", fill = "steelblue") + coord_flip() + theme_minimal() + labs(title = "GO Enrichment: Optimal Constant Genes", x = "Biological Process", y = "-Log10(Benjamini P-value)") + theme(text = element_text(size = 14), axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black"))
ggsave("Figure_3B_Baseline.png", width = 10, height = 5, dpi = 300)

# Figure 4A: Lost Stability
df_lost <- data.frame(term_name = c("Ribosome & Translation", "Ubiquitin / Isopeptide Bond", "Cytosolic Small Ribosomal Subunit", "rRNA binding", "Translational elongation"), adjusted_p_value = c(3.80e-114, 9.36e-42, 6.98e-50, 1.13e-11, 7.44e-10))
df_lost$log_p <- -log10(df_lost$adjusted_p_value)
ggplot(df_lost, aes(x = reorder(term_name, log_p), y = log_p)) + geom_bar(stat = "identity", fill = "indianred") + coord_flip() + theme_minimal() + labs(title = "GO Enrichment: Genes Losing Stability in Salt", x = "Biological Process", y = "-Log10(Benjamini P-value)") + theme(text = element_text(size = 14), axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black"))
ggsave("Figure_4A_LostStability.png", width = 10, height = 5, dpi = 300)

# Figure 4B: Core Constant
df_core <- data.frame(term_name = c("Glycolysis / Gluconeogenesis", "Ubiquitin / Isopeptide Bond"), adjusted_p_value = c(1.99e-10, 9.66e-10))
df_core$log_p <- -log10(df_core$adjusted_p_value)
ggplot(df_core, aes(x = reorder(term_name, log_p), y = log_p)) + geom_bar(stat = "identity", fill = "seagreen") + coord_flip() + theme_minimal() + labs(title = "GO Enrichment: Core Constant Survivors", x = "Biological Process", y = "-Log10(Benjamini P-value)") + theme(text = element_text(size = 14), axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black"))
ggsave("Figure_4B_CoreConstant.png", width = 10, height = 3, dpi = 300)
 # Figure 5A: Fold Change Scatterplot
fc_plot_data <- data.frame(Gene = names(ctrl_medians), Control_Median = ctrl_medians, Salt_Median = salt_medians, Category = "Unaffected / Background"); fc_plot_data$Category[fc_plot_data$Gene %in% top_repressed] <- "Top 100 Repressed (Ribosomal)"; fc_plot_data$Category[fc_plot_data$Gene %in% top_induced] <- "Top 100 Induced (Defense/PPP)"; fc_plot_data$Category <- factor(fc_plot_data$Category, levels = c("Unaffected / Background", "Top 100 Repressed (Ribosomal)", "Top 100 Induced (Defense/PPP)")); fc_plot_data <- fc_plot_data[order(fc_plot_data$Category), ]; ggplot(fc_plot_data, aes(x = Control_Median + 1, y = Salt_Median + 1, color = Category)) + geom_point(alpha = 0.7, size = 1.5) + scale_color_manual(values = c("Unaffected / Background" = "gray85", "Top 100 Repressed (Ribosomal)" = "darkorange", "Top 100 Induced (Defense/PPP)" = "mediumpurple")) + scale_x_log10(breaks = c(1, 11, 101, 1001), labels = c("0", "10", "100", "1000")) + scale_y_log10(breaks = c(1, 11, 101, 1001), labels = c("0", "10", "100", "1000")) + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") + theme_minimal() + labs(title = "Global Transcriptome Shift: Control vs. Salt Stress", x = "Optimal Control Median Expression (Log10)", y = "Salt Stress Median Expression (Log10)") + theme(text = element_text(size = 14), legend.position = "bottom", legend.title = element_blank()); ggsave("Figure_5A_FoldChange_Scatter.png", width = 8, height = 6, dpi = 300)

# Figure 5B: Top 100 Repressed
df_repressed <- data.frame(term_name = c("Ribosomal Large Subunit Biogenesis", "Nucleolus", "Ubiquitin / Isopeptide Bond", "Purine Nucleotide Biosynthesis", "RNA Polymerase I Complex"), adjusted_p_value = c(1.76e-10, 2.16e-16, 6.47e-3, 1.32e-2, 1.58e-3))
df_repressed$log_p <- -log10(df_repressed$adjusted_p_value)
ggplot(df_repressed, aes(x = reorder(term_name, log_p), y = log_p)) + geom_bar(stat = "identity", fill = "darkorange") + coord_flip() + theme_minimal() + labs(title = "GO Enrichment: Top 100 Repressed Genes", x = "Biological Process", y = "-Log10(Benjamini P-value)") + theme(text = element_text(size = 14), axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black"))
ggsave("Figure_5B_Repressed_GO.png", width = 10, height = 5, dpi = 300)

# Figure 5C: Top 100 Induced
df_induced <- data.frame(term_name = c("Carbohydrate Metabolic Process", "Pentose Phosphate Pathway", "Trehalose Biosynthetic Process", "Protein Kinase (S_TKc)"), adjusted_p_value = c(4.53e-2, 9.92e-3, 3.42e-2, 2.23e-2))
df_induced$log_p <- -log10(df_induced$adjusted_p_value)
ggplot(df_induced, aes(x = reorder(term_name, log_p), y = log_p)) + geom_bar(stat = "identity", fill = "mediumpurple") + coord_flip() + theme_minimal() + labs(title = "GO Enrichment: Top 100 Induced Genes", x = "Biological Process", y = "-Log10(Benjamini P-value)") + theme(text = element_text(size = 14), axis.text.x = element_text(color = "black"), axis.text.y = element_text(color = "black"))
ggsave("Figure_5C_Induced_GO.png", width = 10, height = 5, dpi = 300)

message("\n[ INFO ] PIPELINE EXECUTION COMPLETE. ALL FIGURES SAVED.")

