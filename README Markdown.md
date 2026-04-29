# Non-Parametric Characterization of Single-Cell Gene Expression Heterogeneity
This repository contains the computational pipeline used to process, analyze, and visualize single-cell RNA sequencing (scRNA-seq) data from *Saccharomyces cerevisiae* under optimal and hyperosmotic salt stress conditions.
## Project Overview
Traditional parametric statistics often fail to accurately characterize transcriptional noise due to the zero-inflated, right-skewed nature of single-cell data. This pipeline utilizes a non-parametric approach—calculating the Interdecile Range (P90 - P10) relative to the median expression—to isolate true biological bursting from technical dropouts.
### Key Analyses Performed:
1. **Baseline Homeostasis:** Classification of "Constant" (housekeeping) vs. "Strict Bursting" genes in unstressed conditions.
2. **Intersection Analysis:** Tracking the regulatory collapse and stabilization of the transcriptome under 0.6M NaCl salt shock.
3. **Magnitude of the ESR:** Median Fold-Change analysis to quantify the absolute repression of ribosomes and the induction of stress-defense pathways (e.g., Pentose Phosphate Pathway).
## Requirements & Dependencies
This pipeline was executed in **R version 4.5.2 (2025-10-31)**.
* `Matrix` (v1.7-4) - For sparse matrix manipulation.
* `ggplot2` (v4.0.2) - For data visualization and Gene Ontology (GO) plotting.
## Data Source
The pre-processed, SCNorm-normalized count matrix utilized in this pipeline was acquired from the NCBI Gene Expression Omnibus (GEO), Accession **E-GEOD-102475** (Gasch et al., 2017).
## Execution
Run the `main_analysis.R` script. The script autonomously partitions the data, calculates the non-parametric dispersion metrics, performs cross-condition set-theory logic, and outputs all required `ggplot2` `.png` visualizations.
