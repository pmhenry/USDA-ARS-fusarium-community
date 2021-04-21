# Microbial Community Analysis 

## About
This respository contains R workflows used to analyze microbial communities rRNA count data. Individual taxa comarisons were performed using the package `DESeq2` while community comparisons were performed via the `microbiome` package. All analyses utlized phylogenetic, taxonomic, and meta data sets to inform comparisons.

## Analysis overview
Four analyses were performed to analyze how the fungal genus Fusarium and its associated microbial communities respond to various agricultural crop rotation schemes.
* EF, uses elogation factor oligotype OTU rRNA count data to explore species-level comparisons between treatments.
* ITS, uses interal transcribed spacer rRNA count data to explore genus-level comparisons between treatments.
* comb-wsfs, uses EF rRNA OTU counts to compare the effects of two tillage and two cover cropping treatments.
* RK-T52, uses estimates of *actual* soil microbe EF rRNA OTU abundances to compare species between treatments.
For each analysis, raw data (found in `data_raw/`) was formatted via the script `format_data_raw.R` and exported into the directory 'data_prepared/'. 

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
