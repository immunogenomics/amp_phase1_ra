# Accelerating Medicines Partnership (AMP) Phase 1 Rheumatoid Arthritis (RA)

## Overview
This repo provides the *data availability*, *code*, *website* for our work on analyzing and integrating single-cell transcriptomics and mass cytometry data to define inflammatory cell states in rheumatoid arthritis joint synovial tissue.

The preprint version can be viewed and cited:

> Zhang, F., Wei, K., Slowikowski, K., Fonseka, C.Y., Rao, D.A., et al, 2018. Defining Inflammatory Cell States in Rheumatoid Arthritis Joint Synovial Tissues by Integrating Single-cell Transcriptomics and Mass Cytometry. bioRxiv, p.351130. In press.


## Data availibility

The raw data of this study are available at:

1. ImmPort (study accession code SDY998 and SDY999): data for single-cell RNA-seq, mass cytometry, bulk RNA-seq, flow cytometry, clinical and histology
2. dbGAP (study accession: phs001457.v1.p1): single-cell RNA-seq and mass cytometry data 

Send us (fanzhang@broadinstitute.org or jmears@broadinstitute.org) an email if you have any quesitons or requests for data download.

## Code 

### Clone this repo:

```bash
cd ~/work/
git clone git@github.com:immunogenomics/amp_phase1_ra.git
cd amp_phase1_ra

# Download the data from Partners
rsync -avh rgs04:/data/srlab/public/srcollab/AMP/amp_phase1_ra/data .


### Structure

The files in the repo are organized as follows:
    .
    ├── R
    |── data

`data/` has Excel sheets with sample metadata and RData files with processed data ready for analysis.

`R/` has code for analysis and creating figures.

+ Classify tissue samples using Mahalanobis distance: `R/optimal_lymphocyte_threshold.R`

+ Integrate bulk with single-cell RNA-seq: `R/scRNAseq_bulkRNAseq_integrative_pipeline.R`

+ Identify cluster marker genes: `R/cluster_marker_table.R`, `R/limma_differential_bulk.R`

+ Functions for PCA, densisty analysis, etc: `R/pure_functioins.R`

+ Visualize results: `R/cytof_results_plot.R`, `plot_cluster_markers.R`, etc


### Clone this repo:

```bash
cd ~/work/
git clone git@github.com:immunogenomics/amp_phase1_ra.git
cd amp_phase1_ra

# Download the data from Partners
rsync -avh rgs04:/data/srlab/public/srcollab/AMP/amp_phase1_ra/data .
``` 

Send us (fanzhang@broadinstitute.org) an email if you have any quesitons for the analysis.


## Website 

Feel free to check out the  websites and search your favorite genes:
 
1. [Shiny app](https://immunogenomics.io/ampra/): view single-cell RNA-seq, bulk RNA-seq, and mass cytometry data for rheumatoid arthritis data.
2. [UCSC Cell Browser](https://immunogenomics.io/cellbrowser/): view single-cell RNA-seq datasets: 1 rheumatoid arthritis dataset and 2 lupus datasets.
3. [Broad Institue Single Cell Portal](https://portals.broadinstitute.org/single_cell/study/amp-phase-1): view single-cell RNA-seq datasets: 1 rheumatoid arthritis datset and 2 lupus datasets. 

Send us (kslowikowski@gmail.com, jmears@broadinstitute.org, and fanzhang@broadinstitute.org) an email if you have any quesitons for the websites. 
