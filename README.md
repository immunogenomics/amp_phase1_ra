# Accelerating Medicines Partnership Phase 1 Rheumatoid Arthritis (RA)

The preprint version can be viewed and cited:

> Zhang, F., Wei, K., Slowikowski, K., Fonseka, C.Y., Rao, D.A., et al, 2018. Defining Inflammatory Cell States in Rheumatoid Arthritis Joint Synovial Tissues by Integrating Single-cell Transcriptomics and Mass Cytometry. bioRxiv, p.351130. Under review.

## Overview

This repo has the code used to analyze data and make figures in the published manuscript. The files in the repo are organized as follows:

    .
    ├── R
    |── data

`R/` has code for analysis and creating figures.

`data/` has Excel sheets with sample metadata and RData files with processed data ready for analysis.

## Website 

Feel free to use our websites to view the immune and stroma cell populations, and search your favorite genes:
 
1. [Shiny app](https://immunogenomics.io/ampra/): view single-cell RNA-seq, bulk RNA-seq, and mass cytometry data for rheumatoid arthritis
2. [UCSC Cell Browser](https://immunogenomics.io/cellbrowser/): view single-cell RNA-seq datasets for rheumatoid arthritis and lupus.
3. [Broad Institue Single Cell Portal](https://portals.broadinstitute.org/single_cell/study/amp-phase-1): view single-cell RNA-seq datasets for rheumatoid arthritis and lupus. 

## Data Availability

The raw data of this study are available at:

1. ImmPort (study accession code SDY998 and SDY999)
2. dbGAP (study accession: phs001457.v1.p1). 
3. Send me (fanzhang@broadinstitute.org) an email if you have any quesitons or requests for data download.

## Getting Started

Clone this repo:

```bash
cd ~/work/
git clone git@github.com:immunogenomics/amp_phase1_ra.git
cd amp_phase1_ra

# Download the data from Partners
rsync -avh rgs04:/data/srlab/public/srcollab/AMP/amp_phase1_ra/data .
```

## Multiple high-dimensional datasets integration:

#### Optimize leukocyte threshold to classify tissue samples based on leukocyte infiltration by flow cytometry:

    optimal_lymphocyte_threshold.R
    inflamed.html
    
#### CCA analysis of 

1. bulk and single-cell RNA-seq

We use canonical correlation analysis (CCA) to integrate bulk RNA-seq with single-cell RNA-seq. See this file for a detailed walk-through of the analysis:

        R/scRNAseq_bulkRNAseq_integrative_pipeline.R

2. single-cell RNA-seq and mass cytometry

We use the regularized CCA to integrate bulk RNA-seq with mass cytometry. 
See this file for a detailed walk-through of the analysis:

        cytof_bulkRNAseq_integrative_pipeline.R
    
#### Identification of single-cell RNA-seq cluster marker genes:    
    
    cluster_marker_table.R
    
#### Differential analysis for RNA-seq data:

    limma_differential_bulk.R
    
#### Statistical analysis functions for PCA, densitiy analysis, etc:

    pure_functioins.R
    

## Visualization of results for manuscript figures

#### Single cells from RNA-seq in low dimensional space from CCA-based framework (CCA and tSNE):

    cca_bulk_singlecell_visualization

#### Single cells from mass cytometry in low dimensional space (tSNE):

    cytof_results_plot.R
    
#### Plot flow gates for each disease cohort (OA, leukocyte-poor RA, and leukocyte-rich RA):

    Figure2_plots_krenn_flow.R

#### Plot DE markers for each single-cell RNA-seq cluster based on AUC, wilcox p, FC, percent of non-zero expressing cells, etc:

    plot_cluster_markers.R
    
#### Plot post-QC cells based on # genes detected vs percent of Molecules from MT:
    
    plot_genes_detected_mitoch.R
    



