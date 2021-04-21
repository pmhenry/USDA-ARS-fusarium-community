# Microbial Community Analysis 

## About
This respository contains R workflows used to analyze microbial communities rRNA sequence data. Individual taxa comarisons were performed using the package `DESeq2` while community comparisons were performed via the `microbiome` package. All analyses utlized phylogenetic, taxonomic, and meta data sets to inform comparisons.

## Repository Overview
* `data_raw/`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Original data recieved for analysis.
* `data_prepared/`&nbsp;Formatted data prepared for analysis using the script `format_data_raw.R`.
* `scripts/`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;All R scripts used in this publication analysis.

## Dependencies
* [phyloseq](https://joey711.github.io/phyloseq/) for stacked bar chart and NMDS vizualizations.
* [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) for individual taxa comarisons.
* [microbiome](https://microbiome.github.io/) for community-level comparisons.
* [vegan](https://cran.r-project.org/web/packages/vegan/index.html) is used with `microbiome` to perform PERMANOVA community comparisons.
* [xlsx](https://cran.r-project.org/web/packages/xlsx/index.html) for manipulating excel files.
