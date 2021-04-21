# Setup -------------------------------------------------------------------

# * Globals ---------------------------------------------------------------

## Analysis Prefix, used for naming output directories
vOutDir = "output_comb-wsfs_v3"

vW_ord = 65  ## Figure width,  ordination (mm)
vH_ord = 86  ## Figure height, ordination (mm)

vW_DESeq2 = 60 ## Figure width,  DESeq2 log2FC (mm)
vH_DESeq2 = 89 ## Figure height, DESeq2 log2FC (mm)



# * Files -----------------------------------------------------------------

## Path lists
vPath_current <- as.list(strsplit(dirname(rstudioapi::getSourceEditorContext()$path), "/")[[1]]) ## Windows parsing
vPath_base    <- vPath_current[ 1:length(vPath_current)-1 ]

## Directory strings
vDir_base      <- do.call('file.path', vPath_base)        ## Directory, base
vDir_data_prep <- file.path(vDir_base, "data_prepared")   ## Directory, output (processed) data
vDir_output    <- file.path(vDir_base, "output_analysis") ## Directory, analysis output

## Paths, Data file
vInfile_counts_Bulk   <- file.path(vDir_data_prep, "data_comb-wsfs", "comb-wsfs_count_Bulk.csv")    ## count, Bulk Soil '17 & '18
vInfile_counts_SlRt   <- file.path(vDir_data_prep, "data_comb-wsfs", "comb-wsfs_count_SlRt.csv")    ## count, Soil & Root '17
vInfile_taxa_Bulk     <- file.path(vDir_data_prep, "data_comb-wsfs", "comb-wsfs_taxonomy_Bulk.csv") ## taxonomy, Bulk Soil '17 & '18
vInfile_taxa_SlRt     <- file.path(vDir_data_prep, "data_comb-wsfs", "comb-wsfs_taxonomy_SlRt.csv") ## taxonomy, Soil & Root '17
vInfile_sample_Bulk   <- file.path(vDir_data_prep, "data_comb-wsfs", "comb-wsfs_sample_Bulk.csv")   ## sample, Bulk Soil '17 & '18
vInfile_sample_SlRt   <- file.path(vDir_data_prep, "data_comb-wsfs", "comb-wsfs_sample_SlRt.csv")   ## sample, Soil & Root '17
vInfile_distance      <- file.path(vDir_data_prep, "data_comb-wsfs", "comb-wsfs_distance.phy")      ## distance, used for both analyses


# * Packages --------------------------------------------------------------

library(DESeq2)
library(microbiome)
library(phyloseq)
library(vegan)
library(xlsx)



# * Functions ---------------------------------------------------------------

fExportSVG <- function(vfSVG, vfFolder, vfName, vfW=400, vfH=600) {
  "FUNCTION: Saves an svg object (vfSVG) to a specified location(vfFolder + vfName)."
  library("ggplot2")
  fTryCreatePath(vfFolder)
  vSvgName <- file.path(vfFolder, vfName)
  ggsave(file=vSvgName, plot=vfSVG, width=vfW, height=vfH, units="mm")
}

fExportXLSX <- function(vfXLSX, vfFolder, vfName, vfSheet, vfRow=T,vfCol=T) {
  "FUNCTION: Saves an xlsx object (vfXLSX) to a specified location (vfFolder + vfName)."
  library("xlsx")
  fTryCreatePath(vfFolder)
  vXlsxName <- file.path(vfFolder, vfName)
  write.xlsx(x=vfXLSX, file=vXlsxName, sheetName=vfSheet, append=TRUE, row.names=vfRow, col.names=vfCol)
}

fTryRemoveFile <- function(vPath) {
  "FUNCTION: Tests for the presence of a file (vPath) and deletes it if it exists."
  if (file.exists(vPath)) {file.remove(vPath)}
}

fTryCreatePath <- function(vPath) {
  "FUNCTION: Creates a directory (vPath) if it does not exist."
  dir.create(file.path(vPath), recursive=T, showWarnings=F)
}



# * Initialize Environment --------------------------------------------------

set.seed(64)

## Define output directories
vDirPhyloseq   <- file.path(vDir_output, vOutDir, "Phyloseq")   ## Phyloseq
vDirDESeq2     <- file.path(vDir_output, vOutDir, "DESeq2")     ## DESeq2
vDirMicrobiome <- file.path(vDir_output, vOutDir, "Microbiome") ## Microbiome

## Generate output directories
fTryCreatePath(vDirPhyloseq)
fTryCreatePath(vDirDESeq2)
fTryCreatePath(vDirMicrobiome)

## ggplot
theme_set(theme_bw())

## Import Data Files
dfCoun_Bulk <- read.csv(file =vInfile_counts_Bulk, row.names=1, check.names=FALSE) ## count,    Bulk Soil '17 & '18
dfCoun_SlRt <- read.csv(file =vInfile_counts_SlRt, row.names=1, check.names=FALSE) ## count,    Soil & Root '17
dfTaxa_Bulk <- read.csv(file =vInfile_taxa_Bulk, row.names=1, check.names=FALSE)   ## taxonomy, Bulk Soil '17 & '18
dfTaxa_SlRt <- read.csv(file =vInfile_taxa_SlRt, row.names=1, check.names=FALSE)   ## taxonomy, Soil & Root '17
dfMeta_Bulk <- read.csv(file =vInfile_sample_Bulk, row.names=1, check.names=FALSE) ## sample,   Bulk Soil '17 & '18
dfMeta_SlRt <- read.csv(file =vInfile_sample_SlRt, row.names=1, check.names=FALSE) ## sample,   Soil & Root '17
phyTree     <- read_tree(treefile =vInfile_distance)                               ## distance, for both data sets



# |==================| ----------------------------------------------------
# Phyloseq ----------------------------------------------------------------
# * Initialize  -----------------------------------------------------------

## Phyloseq Object: Bulk, Raw counts
physeq_raw_Bulk <- phyloseq(otu_table(dfCoun_Bulk, taxa_are_rows=F), ## COUNTS
                            tax_table(as.matrix(dfTaxa_Bulk)),       ## TAXONOMY
                            sample_data(dfMeta_Bulk),                ## SAMPLE
                            phyTree)                                 ## DISTANCE

## Phyloseq Object: SlRt, Raw counts
physeq_raw_SlRt <- phyloseq(otu_table(dfCoun_SlRt, taxa_are_rows=F), ## COUNTS
                            tax_table(as.matrix(dfTaxa_SlRt)),       ## TAXONOMY
                            sample_data(dfMeta_SlRt),                ## SAMPLE
                            phyTree)                                 ## DISTANCE


## Format: Remove outlier (low sequencing coverage)
physeq_raw_Bulk <- prune_samples(sample_names(physeq_raw_Bulk) != "CTCC2_2018", physeq_raw_Bulk)

## Phyloseq Object: Bulk, Normalized counts
physeq_norm_Bulk <- transform_sample_counts(physeq_raw_Bulk, function(x) x / sum(x) )

## Phyloseq Object: SlRt, Normalized counts
physeq_norm_SlRt <- transform_sample_counts(physeq_raw_SlRt, function(x) x / sum(x) )



# * Plot: Bulk, Stacked Bar Chart -----------------------------------------

## PLOT: Visualize abundance to determine cutoff
N      <- 30 ## Number of OTUs to display.
cutoff <- 1  ## User-defined % abundance for determining what taxa to include in stacked bar plot
par(mar = c(9, 4, 4, 2) + 0.1) ## make more room on bottom margin
barplot(sort(taxa_sums(physeq_norm_Bulk), TRUE)[1:N]/nsamples(physeq_norm_Bulk)*100, las=2, cex.names=.75, main="Rank-abundance OTU Barplot, 1% cutoff", ylab="relative abundance (%)", mgp = c(2, .5, 0))
lines(c(-50,100),c(cutoff, cutoff), lty = "dotted", ylab="0.01", col="red")

## Create phyloseq object with only top N oligotypes included (14 identified above 1%)
TopNOTUs_Bulk <- names(sort(taxa_sums(physeq_norm_Bulk), TRUE)[1:14]) 
physeq_toptax_Bulk <- prune_taxa(TopNOTUs_Bulk, physeq_norm_Bulk)

## Bar Graph, Individual Replicates
mycolors_14 <- c("#850027", "#c12826", "#aa4f2b", "#ca7c4a", "#feb70e", "#ece372", "#1bd17c", "#628c46", "#664e93", "#ca8096", "#6a1154", "#375670", "#1185c0", "#3e30ff")
vPhy_Bar_Bulk <- plot_bar(physeq_toptax_Bulk, fill="taxa")+
  facet_grid(~full_trtmt+year, scales="free_x")+
  ggtitle("Bulk Soil '17-'18, Taxa Composition")+
  xlab("Repetition")+
  ylab("Relative Abundance")+
  scale_fill_manual(values=mycolors_14)
vPhy_Bar_Bulk
fExportSVG(vPhy_Bar_Bulk, vDirPhyloseq, "barplot_trt-year_Bulk.svg")



# * Plot: SlRt, Stacked Bar Chart -----------------------------------------

## PLOT: Visualize abundance to determine cutoff
N <- 30 ## Number of OTUs to display.
cutoff <- 1
par(mar = c(9, 4, 4, 2) + 0.1) ## make more room on bottom margin
barplot(sort(taxa_sums(physeq_norm_SlRt), TRUE)[1:N]/nsamples(physeq_norm_SlRt)*100, las=2, cex.names=.75, main="Rank-abundance OTU Barplot, 1% cutoff", ylab="relative abundance (%)", mgp = c(2, .5, 0))
lines(c(-50,100),c(cutoff, cutoff), lty = "dotted", ylab="0.01", col="red")

## Create phyloseq object with only top N oligotypes included (21 identified above 1%)
TopNOTUs_SlRt <- names(sort(taxa_sums(physeq_norm_SlRt), TRUE)[1:21]) 
physeq_toptax_SlRt <- prune_taxa(TopNOTUs_SlRt, physeq_norm_SlRt)

## Bar Graph, Individual Replicates
mycolors_21 <- c("#ca8096", "#fa0050", "#9c0750", "#6a1154", "#850027", "#c12826", "#aa4f2b", "#ca7c4a", "#feb70e", "#ece372", "#c2d8a8", "#1bd17c", "#3cea35", "#628c46", "#48534c", "#375670", "#2c1141", "#664e93", "#4da1d5", "#1185c0", "#3e30ff")
vPhy_Bar_SlRt <- plot_bar(physeq_toptax_SlRt, fill="taxa")+
  facet_grid(~full_trtmt+sample_type, scales="free_x")+
  ggtitle("Soil & Root '17, Taxa Composition")+
  xlab("Repetition")+
  ylab("Relative Abundance")+
  scale_fill_manual(values=mycolors_21)
vPhy_Bar_SlRt
fExportSVG(vPhy_Bar_SlRt, vDirPhyloseq, "barplot_trt-year_SlRt.svg")



# * Plot: Bulk, Ordination ------------------------------------------------

## NMDS ordination on the bray-curtis distance.
GP.ord_Bulk <- ordinate(physeq_norm_Bulk, "NMDS", "wunifrac")

## PLOT: Single plot, years pooled (FIGURE 3.1, legend)
mycolors_4 <- c("#C12826", "#CA7C4A", "#628C46", "#1185C0")
vPhy_ord_trt_Bulk <- plot_ordination(physeq_norm_Bulk, GP.ord_Bulk, color="full_trtmt")+
  geom_point(size=1, na.rm=T)+
  stat_ellipse(type="norm", size=.5, na.rm=T)+
  scale_colour_manual(values=mycolors_4)+
  theme(
    text=element_text(size=8),
    legend.background=element_rect(size=0.2,
                                   linetype="solid",
                                   colour ="black"),
    legend.position="right",
    legend.key.size = unit(.1, 'lines'),
    legend.text=element_text(size=8))+
  scale_fill_continuous(guide = guide_legend())
vPhy_ord_trt_Bulk
fExportSVG(vPhy_ord_trt_Bulk, vDirPhyloseq, "fig_3.1_ordination_Bulk_legend.svg", vfW=vW_ord, vfH=vH_ord)

## PLOT: Single plot, years pooled (FIGURE 3.1, no legend)
mycolors_4 <- c("#C12826", "#CA7C4A", "#628C46", "#1185C0")
vPhy_ord_trt_Bulk <- plot_ordination(physeq_norm_Bulk, GP.ord_Bulk, color="full_trtmt")+
  geom_point(size=1, na.rm=T)+
  stat_ellipse(type="norm", size=.5, na.rm=T)+
  scale_colour_manual(values=mycolors_4)+
  theme(
    text=element_text(size=8),
    legend.position="none")+
  scale_fill_continuous(guide = guide_legend())
vPhy_ord_trt_Bulk
fExportSVG(vPhy_ord_trt_Bulk, vDirPhyloseq, "fig_3.1_ordination_Bulk_nolegend.svg", vfW=vW_ord, vfH=vH_ord)



# * Plot: SlRt, Ordination ------------------------------------------------

## NMDS ordination on the bray-curtis distance.
GP.ord_SlRt <- ordinate(physeq_norm_SlRt, "NMDS", "wunifrac")

## PLOT: Single plot, treatment + year (SK FIGURE 3.2, legend)
mycolors_8 <- c("#850027", "#c12826", "#aa4f2b", "#ca7c4a", "#48534b", "#628c46", "#2c1141", "#1185c0")
vPhy_ord_comb_SlRt <- plot_ordination(physeq_norm_SlRt, GP.ord_SlRt, color="comb")+
  geom_point(size=1, na.rm=T)+
  stat_ellipse(type="norm", size=.5, na.rm=T)+
  scale_colour_manual(values=mycolors_8)+
  theme(
    text=element_text(size=8),
    legend.background=element_rect(size=0.2,
                                   linetype="solid",
                                   colour ="black"),
    legend.position="bottom",
    legend.key.size = unit(.1, 'lines'),
    legend.text=element_text(size=8))+
  scale_fill_continuous(guide = guide_legend())
vPhy_ord_comb_SlRt
fExportSVG(vPhy_ord_comb_SlRt, vDirPhyloseq, "fig_3.2_ordination_SlRt_legend.svg", vfW=vW_ord, vfH=vH_ord)

## PLOT: Single plot, treatment + year (SK FIGURE 3.2, no legend)
mycolors_8 <- c("#850027", "#c12826", "#aa4f2b", "#ca7c4a", "#48534b", "#628c46", "#2c1141", "#1185c0")
vPhy_ord_comb_SlRt <- plot_ordination(physeq_norm_SlRt, GP.ord_SlRt, color="comb")+
  geom_point(size=1, na.rm=T)+
  stat_ellipse(type="norm", size=.5, na.rm=T)+
  scale_colour_manual(values=mycolors_8)+
  theme(
    text=element_text(size=8),
    legend.position="none")+
  scale_fill_continuous(guide = guide_legend())
vPhy_ord_comb_SlRt
fExportSVG(vPhy_ord_comb_SlRt, vDirPhyloseq, "fig_3.2_ordination_SlRt_nolegend.svg", vfW=vW_ord, vfH=vH_ord)



# |==================| ----------------------------------------------------
# DESeq2 ------------------------------------------------------------------
# * Initialize ------------------------------------------------------------

## Initialize DESeq2 object (formula is placeholder, updated when "groups" defined)
## Note: Non-normalized counts are used
diagdds_Bulk = phyloseq_to_deseq2(physeq_raw_Bulk, design = ~ full_trtmt)
diagdds_SlRt = phyloseq_to_deseq2(physeq_raw_SlRt, design = ~ sample_type)

## Account for 0 counts 
diagdds_Bulk <- estimateSizeFactors(diagdds_Bulk, type="poscounts")
diagdds_SlRt <- estimateSizeFactors(diagdds_SlRt, type="poscounts")

## Run DESeq2
dds_Bulk = DESeq(diagdds_Bulk, test="Wald", fitType="mean")
dds_SlRt = DESeq(diagdds_SlRt, test="Wald", fitType="mean")



# * Perform Comparisons ---------------------------------------------------

## Function for exracting DESeq2 comparisons
fExtract_comparisons <- function(vComparisons, vLevel, vSub_padj=.05, vSub_log2=1.5, vPhyseq, vDESeq2, vFolder, vLabel) {
  "
  FUNCTION: Extract Comparisons from DESeq2 analysis into sheets of an excel file
  VARIABLES:
    > vComparisons  (matrix) 2-column list of comparisons to extract
    > vLevel        (str) The meta data category that vComparisons corresponds to
    > vSub_padj     (int) Significance cutoff for adjusted p-value (Less than or equal to)
    > vSub_log2     (int) Significance cutoff for log2-fold change (+/- larger/smaller than or equal to, respectively)
    > vPhyseq       Phyloseq object used for DESeq2 analysis containing taxonomy data
    > vDESeq2       DESeq2 object to subset from
    > vFolder       Folder to save ggplot images to
    > vLabel        The name of the comparison being performed
  "
  ## Define Output Folder for individual figures/data sets
  vFolderPath <- file.path(vFolder, vLabel)
  ## Define comparison data file
  vOutfile = paste(vLabel,".xlsx", sep="")
  ## Remove current version of output excel file
  fTryRemoveFile(file.path(vFolderPath,vOutfile))
  ## Initialize data frame for saving fusarium statistics
  vDF_Fus <- data.frame(matrix(nrow=nrow(vComparisons), ncol=6))
  rownames(vDF_Fus) <- paste(vComparisons[,1],vComparisons[,2], sep="_")
  print(rownames(vDF_Fus))
  colnames(vDF_Fus) <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")					
  for (n in 1:nrow(vComparisons)) {
    ## Extract comparison results
    vName          <- paste(vComparisons[n,1], "_", vComparisons[n,2], sep="")
    vFilename_svg  <- paste(vName,".svg", sep="")
    vFilename_xlsx <- paste(vName,".xlsx", sep="")
    res = results(vDESeq2, contrast = c(vLevel, vComparisons[n,1], vComparisons[n,2])) ## Comparison
    ## Subset results
    vSubset <- res[order(res$padj),] ## Sort
    vSubset <- subset(x=vSubset, vSubset$padj           <= vSub_padj) ## Subset Significant p values
    vSubset <- subset(x=vSubset, vSubset$log2FoldChange >= vSub_log2 | vSubset$log2FoldChange <= -vSub_log2 ) ## Subset substantial change
    ## Export to excel
    if (nrow(vSubset)>0) {
      ## Append taxonomy table to results
      vSubset <- cbind(as(vSubset, "data.frame"), as(tax_table(vPhyseq)[rownames(vSubset), ], "matrix"))
      ## Export
      fExportXLSX(vfXLSX=vSubset, vfFolder=vFolderPath, vfName=vOutfile, vfSheet=vName)
      print(paste("Exported:      ",vName))
      ## PLOT
      x = tapply(vSubset$log2FoldChange, vSubset$taxa, function(x) max(x))
      x = sort(x, TRUE)
      vSubset$taxa = factor(as.character(vSubset$taxa), levels=names(x))
      vSVG <- ggplot(vSubset, aes(y=taxa, x=log2FoldChange))+
        geom_vline(xintercept=0, color="blue")+
        geom_point(size=1, shape=23)+ 
        theme(text = element_text(size=8),
              axis.title.x=element_blank())+
        ylab("OTU")
      fExportSVG(vSVG, vFolderPath, vFilename_svg, vfW=vW_DESeq2, vfH=vH_DESeq2)
    } else {
      ## Export
      fExportXLSX(vfXLSX="NONE", vfFolder=vFolderPath, vfName=vOutfile, vfSheet=vName, vfRow=F,vfCol=F)
      print(paste("No Differences:",vName))
    }
  }
}

## Define comparisons
vComp_DESeq2_Bulk <- t(combn(unique(dfMeta_Bulk$full_trtmt), 2))
vComp_DESeq2_SlRt <- data.frame(c1=c("soil"),
                                c2=c("root"))

## Extract Comparisons, Bulk
fExtract_comparisons(vComparisons = vComp_DESeq2_Bulk,
                     vLevel       = "full_trtmt",
                     vSub_padj    = 0.05,
                     vSub_log2    = 1.5,
                     vPhyseq      = physeq_raw_Bulk,
                     vDESeq2      = dds_Bulk,
                     vFolder      = vDirDESeq2,
                     vLabel       = "DESeq2_Bulk")

## Extract Comparisons, SlRt
fExtract_comparisons(vComparisons = vComp_DESeq2_SlRt,
                     vLevel       = "sample_type",
                     vSub_padj    = 0.05,
                     vSub_log2    = 1.5,
                     vPhyseq      = physeq_raw_SlRt,
                     vDESeq2      = dds_SlRt,
                     vFolder      = vDirDESeq2,
                     vLabel       = "DESeq2_SlRt")

## Extract Comparisons, SlRt, NO CUTOFFS (SK FIGURE 3.3)
fExtract_comparisons(vComparisons = vComp_DESeq2_SlRt,
                     vLevel       = "sample_type",
                     vSub_padj    = 1,
                     vSub_log2    = 0,
                     vPhyseq      = physeq_raw_SlRt,
                     vDESeq2      = dds_SlRt,
                     vFolder      = vDirDESeq2,
                     vLabel       = "DESeq2_SlRt_nocut")



# |==================| ----------------------------------------------------
# Microbiome: PERMANOVA ---------------------------------------------------
# * Initialize ------------------------------------------------------------

## Define analysis object as normalized phyloseq data
pseq.rel_Bulk <- microbiome::transform(physeq_raw_Bulk, "compositional")
pseq.rel_SlRt <- microbiome::transform(physeq_raw_SlRt, "compositional")



# * Bulk: Experiment Effects, Permanova ---------------------------------------

fPermanova_Bulk <- function(vfComparisons, vfPhyloseq, vfFolder, vfName) {
  "
  FUNCTION:   Performs permanova comparison between two sample groups.
  VARIABLES:
    > vComparisons, A 2-column matrix of comparisons to analyze
    > vPhyloseq, A normalized phyloseq object from which to subset
    > vfName, String to be appended to output filenames
  OUTPUT:
    > 
  "
  ## Initialize data frame for saving summary statistics
  vDF_Sum <- data.frame(matrix(nrow=nrow(vfComparisons), ncol=5))
  colnames(vDF_Sum) <- c("Treatment_1", "Treatment_2", "Pr(>F)_Beta", "Pr(>F)_Perm", "R2_Perm")
  vDF_Sum$Treatment_1 <- vfComparisons[,1]
  vDF_Sum$Treatment_2 <- vfComparisons[,2]
  
  fTryRemoveFile(file.path(vfFolder,vfName,paste(vfName,".xlsx", sep="")))
  # print(file.path(vfFolder,vfName,paste(vfName,".xlsx", sep="")))
  for (n in 1:nrow(vfComparisons)) {
    ## Names
    vName <- paste(vfComparisons[n,1], vfComparisons[n,2], sep="_")
    vFilename_svg_PCA <- paste(vfName,"_PCA_",     vName,".svg", sep="")
    vFilename_svg_Bar <- paste(vfName,"_TopCoefs_",vName,".svg", sep="")
    vFilename_xlsx    <- paste(vfName,".xlsx", sep="")
    ## Define Current Comparison
    vTitle <- paste(vfName, vfComparisons[n,1], vfComparisons[n,2], sep="_")
    vSheet <- paste(vfComparisons[n,1], vfComparisons[n,2], sep="_")
    vfOut  <- file.path(vfFolder, vfName)
    ## Subset Phyloseq object to selected comparison data
    vTrimKeep <- meta(vfPhyloseq)[,"trt_year"] %in% c(vfComparisons[n,1], vfComparisons[n,2])
    vSubset <- prune_samples(vTrimKeep, vfPhyloseq)
    otu <- abundances(vSubset)
    meta <- meta(vSubset)
    ## Analysis: Permanova
    permanova <- adonis(phyloseq::distance(vSubset, method="wunifrac") ~ trt_year, data = meta, permutations=999)
    ## permanova Results
    vDataPerm <- as.data.frame(permanova$aov.tab)
    # P-value
    ## Check that variance homogeneity assumptions hold (to ensure the reliability of the results). Should be insignificant
    dist <- vegdist(t(otu))
    # print(anova(betadisper(dist, meta$trt_year)))
    vDataPerm <- rbind(PERMANOVA="", statistic=colnames(vDataPerm), vDataPerm)
    vDataDisp <- as.data.frame(anova(betadisper(dist, meta$trt_year, add=T)))
    colnames(vDataDisp) <- c("Df", "SumsOfSqs", "MeanSqs", "F.Model", "Pr(>F)")
    vDataDisp$R2 =c("")
    vDataDisp <- rbind(BETA_DISPERSION="", statistic=colnames(vDataDisp), vDataDisp)
    result <- rbind(vDataPerm,vDataDisp)
    result[is.na(result)] <- ""
    ## Populate summary table
    vDF_Sum[n,"Pr(>F)_Beta"] <- vDataDisp["Groups","Pr(>F)"]
    vDF_Sum[n,"Pr(>F)_Perm"] <- vDataPerm["trt_year","Pr(>F)"]
    vDF_Sum[n,"R2_Perm"]     <- vDataPerm["trt_year","R2"]
    
    ## Export: Permanova comparison and homogeneity check
    fExportXLSX(vfXLSX=result,
                vfFolder=vfOut,
                vfName=vFilename_xlsx,
                vfSheet=vName,
                vfCol=F)
  }
  ## Save results to summary table
  cols.num <- c("Pr(>F)_Beta","Pr(>F)_Perm","R2_Perm")
  vDF_Sum[cols.num] <- sapply(vDF_Sum[cols.num], as.numeric)
  fExportXLSX(vfXLSX=vDF_Sum,
              vfFolder=vfFolder,
              vfName="permanova_summary_Bulk.xlsx",
              vfSheet=vfName,
              vfRow=F)
  print(paste("Exported:",vName))
}

fPermanova_SlRt <- function(vfComparisons, vfPhyloseq, vfFolder, vfName) {
  "
  FUNCTION:   Performs permanova comparison between two sample groups.
  VARIABLES:
    > vComparisons, A 2-column matrix of comparisons to analyze
    > vPhyloseq, A normalized phyloseq object from which to subset
    > vfName, String to be appended to output filenames
  OUTPUT:
    > 
  "
  ## Initialize data frame for saving summary statistics
  vDF_Sum <- data.frame(matrix(nrow=nrow(vfComparisons), ncol=5))
  colnames(vDF_Sum) <- c("Treatment_1", "Treatment_2", "Pr(>F)_Beta", "Pr(>F)_Perm", "R2_Perm")
  vDF_Sum$Treatment_1 <- vfComparisons[,1]
  vDF_Sum$Treatment_2 <- vfComparisons[,2]
  fTryRemoveFile(file.path(vfFolder,vfName,paste(vfName,".xlsx", sep="")))
  # print(file.path(vfFolder,vfName,paste(vfName,".xlsx", sep="")))
  for (n in 1:nrow(vfComparisons)) {
    ## Names
    vName <- paste(vfComparisons[n,1], vfComparisons[n,2], sep="_")
    vFilename_svg_PCA <- paste(vfName,"_PCA_",     vName,".svg", sep="")
    vFilename_svg_Bar <- paste(vfName,"_TopCoefs_",vName,".svg", sep="")
    vFilename_xlsx    <- paste(vfName,".xlsx", sep="")
    ## Define Current Comparison
    vTitle <- paste(vfName, vfComparisons[n,1], vfComparisons[n,2], sep="_")
    vSheet <- paste(vfComparisons[n,1], vfComparisons[n,2], sep="_")
    vfOut  <- file.path(vfFolder, vfName)
    ## Subset Phyloseq object to selected comparison data
    vTrimKeep <- meta(vfPhyloseq)[,"comb"] %in% c(vfComparisons[n,1], vfComparisons[n,2])
    vSubset <- prune_samples(vTrimKeep, vfPhyloseq)
    otu <- abundances(vSubset)
    meta <- meta(vSubset)
    ## Analysis: Permanova
    permanova <- adonis(phyloseq::distance(vSubset, method="wunifrac") ~ comb, data = meta, permutations=999, method = "bray")
    ## permanova Results
    vDataPerm <- as.data.frame(permanova$aov.tab)
    # P-value
    ## Check that variance homogeneity assumptions hold (to ensure the reliability of the results). Should be insignificant
    dist <- vegdist(t(otu))
    vDataPerm <- rbind(PERMANOVA="", statistic=colnames(vDataPerm), vDataPerm)
    vDataDisp <- as.data.frame(anova(betadisper(dist, meta$comb, add=T)))
    colnames(vDataDisp) <- c("Df", "SumsOfSqs", "MeanSqs", "F.Model", "Pr(>F)")
    vDataDisp$R2 =c("")
    vDataDisp <- rbind(BETA_DISPERSION="", statistic=colnames(vDataDisp), vDataDisp)
    result <- rbind(vDataPerm,vDataDisp)
    result[is.na(result)] <- ""
    ## Populate summary table
    vDF_Sum[n,"Pr(>F)_Beta"] <- vDataDisp["Groups","Pr(>F)"]
    vDF_Sum[n,"Pr(>F)_Perm"] <- vDataPerm["comb","Pr(>F)"]
    vDF_Sum[n,"R2_Perm"]     <- vDataPerm["comb","R2"]
    
    ## Export: Permanova comparison and homogeneity check
    fExportXLSX(vfXLSX=result,
                vfFolder=vfOut,
                vfName=vFilename_xlsx,
                vfSheet=vName,
                vfCol=F)
    print(paste("Exported:",vName))
  }
  ## Save results to summary table
  cols.num <- c("Pr(>F)_Beta","Pr(>F)_Perm","R2_Perm")
  vDF_Sum[cols.num] <- sapply(vDF_Sum[cols.num], as.numeric)
  fExportXLSX(vfXLSX=vDF_Sum,
              vfFolder=vfFolder,
              vfName="permanova_summary_SlRt.xlsx",
              vfSheet=vfName,
              vfRow=F)
}

## Define comparisons, Bulk
vComp_Perma_Bulk <- data.frame(c1=c("CTCC2017", "CTNO2017", "STCC2017", "STNO2017"),
                               c2=c("CTCC2018", "CTNO2018", "STCC2018", "STNO2018"))

## Define comparisons, SlRt
vComparisons <- unique(dfMeta_SlRt$comb)
vComp_Perma_SlRt_Soil <- t(combn(vComparisons[grepl("soil", vComparisons, fixed=TRUE)], 2))
vComp_Perma_SlRt_Root <- t(combn(vComparisons[grepl("root", vComparisons, fixed=TRUE)], 2))
vComp_Perma_SlRt_Wthn <- data.frame(c1=c("CTCC_soil", "CTNO_soil", "STCC_soil", "STNO_soil"),
                                    c2=c("CTCC_root", "CTNO_root", "STCC_root", "STNO_root"))

## Remove current version of PERMANOVA summary file
fTryRemoveFile(file.path(vDirMicrobiome,"permanova_summary_Bulk.xlsx"))
fTryRemoveFile(file.path(vDirMicrobiome,"permanova_summary_SlRt.xlsx"))

## Comparisons, Bulk
fPermanova_Bulk(vfComparisons=vComp_Perma_Bulk,
           vfPhyloseq=pseq.rel_Bulk,
           vfFolder=vDirMicrobiome,
           vfName="permanova_Bulk")

## Comparisons, SlRt
fPermanova_SlRt(vfComparisons=vComp_Perma_SlRt_Soil,
           vfPhyloseq=pseq.rel_SlRt,
           vfFolder=vDirMicrobiome,
           vfName="permanova_SlRt_Soil")

fPermanova_SlRt(vfComparisons=vComp_Perma_SlRt_Root,
           vfPhyloseq=pseq.rel_SlRt,
           vfFolder=vDirMicrobiome,
           vfName="permanova_SlRt_Root")

fPermanova_SlRt(vfComparisons=vComp_Perma_SlRt_Wthn,
           vfPhyloseq=pseq.rel_SlRt,
           vfFolder=vDirMicrobiome,
           vfName="permanova_SlRt_Wthn")