# AMP Phase 1 RA

Accelerating Medicines Partnership Phase 1 Rheumatoid Arthritis

This repo has the code we used to analyze data and make figures in the
published manuscript:

> Zhang et al. TITLE, JOURNAL, DOI

## Overview

The files in the repo are organized as follows:

    .
    ├── R
    ├── data
    └── data-raw

`R/` has code for analysis and creating figures.

`data/` has RData files with processed data ready for analysis.

`data-raw` has other data files including Excel sheets with sample metadata.

## Getting Started

Clone this repo:

```bash
cd ~/work/
git clone git@github.com:immunogenomics/amp_phase1_ra.git
cd amp_phase1_ra

# Download the data from Partners
rsync -avh rgs04:/data/srlab/public/srcollab/AMP/amp_phase1_ra/data .
```

## Analysis

### CCA analysis of bulk and single-cell RNA-seq

We use canonical correlation analysis (CCA) to integrate bulk RNA-seq with
single-cell RNA-seq. See this file for a detailed walk-through of the analysis:

    R/scRNAseq_bulkRNAseq_integrative_pipeline.R

