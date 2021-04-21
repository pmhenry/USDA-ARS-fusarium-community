# Setup -------------------------------------------------------------------

# * Globals -----------------------------------------------------------------

## Analysis Prefix, used for naming output directories
vOutDir = "output_ITS"

vW_bar = 100 ## Figure width,  bar charts (mm)
vH_bar = 89  ## Figure height, bar charts (mm)

vW_ord = 65  ## Figure width,  ordination (mm)
vH_ord = 86  ## Figure height, ordination (mm) (NOTE: H86 allows space for custom legend)



# * Files -----------------------------------------------------------------

## Path lists
vPath_current <- as.list(strsplit(dirname(rstudioapi::getSourceEditorContext()$path), "/")[[1]]) ## Windows parsing
vPath_base    <- vPath_current[ 1:length(vPath_current)-1 ]

## Directory strings
vDir_base      <- do.call('file.path', vPath_base)        ## Directory, base
vDir_data_prep <- file.path(vDir_base, "data_prepared")   ## Directory, output (processed) data
vDir_output    <- file.path(vDir_base, "output_analysis") ## Directory, analysis output

## Paths, Data file
vInfile_counts   <- file.path(vDir_data_prep, "data_ITS", "ITS_count.csv")    ## count
vInfile_distance <- file.path(vDir_data_prep, "data_ITS", "ITS_distance.phy") ## distance
vInfile_taxa     <- file.path(vDir_data_prep, "data_ITS", "ITS_taxonomy.csv") ## taxa
vInfile_sample   <- file.path(vDir_data_prep, "data_ITS", "ITS_sample.csv")   ## sample



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
dfCounts <- read.csv(file     =vInfile_counts, row.names =1, check.names=FALSE) ## count
phyTree  <- read_tree(treefile =vInfile_distance)                               ## distance
dfTaxono <- read.csv(file     =vInfile_taxa, row.names =1)                      ## taxonomy
dfSample <- read.csv(file     =vInfile_sample, row.names =1)                    ## sample

## Format: Define variable types
dfSample$group <- factor(dfSample$group)



# |==================| ----------------------------------------------------
# Phyloseq ----------------------------------------------------------------
# * Initialize ------------------------------------------------------------

## Phyloseq object: Non-normalized (207 samples x 270 taxa)
physeq_raw <- phyloseq(otu_table(dfCounts, taxa_are_rows=T), ## Counts
                       tax_table(as.matrix(dfTaxono)),       ## Taxonomy
                       sample_data(dfSample),                ## Sample
                       read_tree(vInfile_distance))          ## Distance

## Format: Remove low-abundance samples and taxa
physeq_raw <-prune_samples(sample_sums(physeq_raw)>=1, physeq_raw)
physeq_raw <-prune_taxa(taxa_sums(physeq_raw)>=1, physeq_raw)

## Format: Remove H2O samples
physeq_raw <- subset_samples(physeq_raw, comb != "brswh2o") ## BRSW H2O sample
physeq_raw <- subset_samples(physeq_raw, treatment != "h2") ## BRSW H2O sample

## Format: Remove taxa not classified to the genus level
physeq_raw <- subset_taxa(physeq_raw, genus!="NA")
taxa_names(physeq_raw) <- data.frame(tax_table(physeq_raw))$genus ## Define index as genus

## Phyloseq Object: Normalized counts
physeq_norm <- transform_sample_counts(physeq_raw, function(x) x / sum(x) )



# * Plot: Stacked Bar Chart -----------------------------------------------

## PLOT: Visualize taxa abundances to determine cutoff
N      <- 30 ## Number of taxa to display.
cutoff <- .51 ## User-defined % abundance for determining what taxa to include in stacked bar plot
par(mar = c(9, 4, 4, 2) + 0.1) ## make more room on bottom margin
barplot(sort(taxa_sums(physeq_norm), TRUE)[1:N]/nsamples(physeq_norm)*100, las=2, cex.names=.75, main="Rank-abundance OTU Barplot, 1% cutoff", ylab="relative abundance (%)", mgp = c(2, .5, 0))
lines(c(-50,100),c(cutoff, cutoff), lty = "dotted", ylab="0.01", col="red")

## Create phyloseq object with only top N oligotypes included (21 identified above .51%)
TopNOTUs <- names(sort(taxa_sums(physeq_norm), TRUE)[1:21]) 
physeq_toptax <- prune_taxa(TopNOTUs, physeq_norm)

## Bar Graph, PLANT Individual Replicates (FIGURE 1.3, legend)
mycolors <- c("#ca8096", "#fa0050", "#9c0750", "#6a1154", "#850027", "#c12826", "#aa4f2b", "#ca7c4a", "#feb70e", "#ece372", "#c2d8a8", "#1bd17c", "#3cea35", "#628c46", "#48534c", "#375670", "#2c1141", "#664e93", "#4da1d5", "#1185c0", "#3e30ff")
physeq_toptax_plnt <- prune_samples(sample_data(physeq_toptax)$extracted_from=="plnt", physeq_toptax)
vPhy_Bar <- plot_bar(physeq_toptax_plnt, x="repetition", y="Abundance", fill="genus")+
  facet_grid(experiment~treatment, scales="free_x")+
  ylab("Relative Abundance")+
  scale_fill_manual(values=mycolors)+
  theme(
    text=element_text(size=8),
    axis.text.x = element_text(angle = 0),
    legend.key.size = unit(.75, 'lines'),
    legend.text=element_text(size=8),
    panel.spacing.x=unit(.1, "lines"),
    panel.spacing.y=unit(.2, "lines"),
    axis.title.x=element_blank())+
  guides(fill=guide_legend(ncol=1))+
  scale_y_continuous(labels=function(x)x*100)
vPhy_Bar
fExportSVG(vPhy_Bar, vDirPhyloseq, "fig_1.3_ITS_barplot_plnt_legend.svg", vfW=vW_bar, vfH=vH_bar)

## Bar Graph, PLANT Individual Replicates (FIGURE 1.3, no legend)
mycolors <- c("#ca8096", "#fa0050", "#9c0750", "#6a1154", "#850027", "#c12826", "#aa4f2b", "#ca7c4a", "#feb70e", "#ece372", "#c2d8a8", "#1bd17c", "#3cea35", "#628c46", "#48534c", "#375670", "#2c1141", "#664e93", "#4da1d5", "#1185c0", "#3e30ff")
physeq_toptax_plnt <- prune_samples(sample_data(physeq_toptax)$extracted_from=="plnt", physeq_toptax)
vPhy_Bar <- plot_bar(physeq_toptax_plnt, x="repetition", y="Abundance", fill="genus")+
  facet_grid(experiment~treatment, scales="free_x")+
  ylab("Relative Abundance")+
  scale_fill_manual(values=mycolors)+
  theme(
    text=element_text(size=10),
    axis.text.x = element_text(angle = 0),
    legend.position="none",
    panel.spacing.x=unit(.1, "lines"),
    panel.spacing.y=unit(.2, "lines"),
    axis.title.x=element_blank())+
  guides(fill=guide_legend(ncol=1))+
  scale_y_continuous(labels=function(x)x*100)
vPhy_Bar
fExportSVG(vPhy_Bar, vDirPhyloseq, "fig_1.3_ITS_barplot_plnt_nolegend.svg", vfW=vW_bar, vfH=vH_bar)

## Bar Graph, SOIL Individual Replicates (FIGURE 2.3, legend)
mycolors <- c("#ca8096", "#fa0050", "#9c0750", "#6a1154", "#850027", "#c12826", "#aa4f2b", "#ca7c4a", "#feb70e", "#ece372", "#c2d8a8", "#1bd17c", "#3cea35", "#628c46", "#48534c", "#375670", "#2c1141", "#664e93", "#4da1d5", "#1185c0", "#3e30ff")
physeq_toptax_plnt <- prune_samples(sample_data(physeq_toptax)$extracted_from=="soil", physeq_toptax)
vPhy_Bar <- plot_bar(physeq_toptax_plnt, x="repetition", y="Abundance", fill="genus")+
  facet_grid(experiment~treatment, scales="free_x")+
  xlab("Repetition")+
  ylab("Relative Abundance")+
  scale_fill_manual(values=mycolors)+
  theme(
    text=element_text(size=8),
    axis.text.x = element_text(angle = 0),
    legend.key.size = unit(.75, 'lines'),
    legend.text=element_text(size=8),
    panel.spacing.x=unit(.1, "lines"),
    panel.spacing.y=unit(.2, "lines"),
    axis.title.x=element_blank())+
  guides(fill=guide_legend(ncol=1))+
  scale_y_continuous(labels=function(x)x*100)
vPhy_Bar
fExportSVG(vPhy_Bar, vDirPhyloseq, "fig_2.3_ITS_barplot_soil_legend.svg", vfW=vW_bar, vfH=vH_bar)

## Bar Graph, SOIL Individual Replicates (FIGURE 2.3, no legend)
mycolors <- c("#ca8096", "#fa0050", "#9c0750", "#6a1154", "#850027", "#c12826", "#aa4f2b", "#ca7c4a", "#feb70e", "#ece372", "#c2d8a8", "#1bd17c", "#3cea35", "#628c46", "#48534c", "#375670", "#2c1141", "#664e93", "#4da1d5", "#1185c0", "#3e30ff")
physeq_toptax_plnt <- prune_samples(sample_data(physeq_toptax)$extracted_from=="soil", physeq_toptax)
vPhy_Bar <- plot_bar(physeq_toptax_plnt, x="repetition", y="Abundance", fill="genus")+
  facet_grid(experiment~treatment, scales="free_x")+
  xlab("Repetition")+
  ylab("Relative Abundance")+
  scale_fill_manual(values=mycolors)+
  theme(
    text=element_text(size=10),
    axis.text.x = element_text(angle = 0),
    legend.position="none",
    panel.spacing.x=unit(.1, "lines"),
    panel.spacing.y=unit(.2, "lines"),
    axis.title.x=element_blank())+
  guides(fill=guide_legend(ncol=1))+
  scale_y_continuous(labels=function(x)x*100)
vPhy_Bar
fExportSVG(vPhy_Bar, vDirPhyloseq, "fig_2.3_ITS_barplot_soil_nolegend.svg", vfW=vW_bar, vfH=vH_bar)




# * Plot: Ordination ------------------------------------------------------

## NMDS ordination on the bray-curtis distance.
GP.ord <- ordinate(physeq_norm, "NMDS", "wunifrac")

## PLOT: Plant subset (FIGURE 1.4, legend)
GP.ord_plnt <- subset_samples(physeq_norm, extracted_from=="plnt")
mycolors_7 <- c("#c12826", "#ca7c4a", "#628c46", "#48534c", "#1185c0", "#664e93", "#ca8096")
vPhy_ord_extract <- plot_ordination(GP.ord_plnt, GP.ord, color="treatment")+
  geom_point(size=1, na.rm=T)+
  scale_colour_manual(values=mycolors_7)+
  stat_ellipse(type="norm", na.rm=T)+
  xlim(-.5, .25)+
  ylim(-.2, .4)+
  theme(
    text=element_text(size=8),
    legend.background=element_rect(size=0.2,
                                   linetype="solid",
                                   colour ="black"),
    legend.position="bottom",
    legend.key.size = unit(.1, 'lines'),
    legend.text=element_text(size=8))+
  scale_fill_continuous(guide = guide_legend())
vPhy_ord_extract
fExportSVG(vPhy_ord_extract, vDirPhyloseq, "fig_1.4_ITS_ordination_plnt_legend.svg", vfW=vW_ord, vfH=vH_ord)

## PLOT: Plant subset (FIGURE 1.4, no legend)
GP.ord_plnt <- subset_samples(physeq_norm, extracted_from=="plnt")
mycolors_7 <- c("#c12826", "#ca7c4a", "#628c46", "#48534c", "#1185c0", "#664e93", "#ca8096")
vPhy_ord_extract <- plot_ordination(GP.ord_plnt, GP.ord, color="treatment")+
  geom_point(size=1, na.rm=T)+
  scale_colour_manual(values=mycolors_7)+
  stat_ellipse(type="norm", na.rm=T)+
  xlim(-.5, .25)+
  ylim(-.2, .4)+
  theme(
    text=element_text(size=10),
    legend.background=element_rect(size=0.2,
                                   linetype="solid",
                                   colour ="black"),
    legend.position="none")+
  scale_fill_continuous(guide = guide_legend())
vPhy_ord_extract
fExportSVG(vPhy_ord_extract, vDirPhyloseq, "fig_1.4_ITS_ordination_plnt_nolegend.svg", vfW=vW_ord, vfH=vH_ord)

## PLOT: Single Soil subset (FIGURE 2.4, legend)
GP.ord_soil <- subset_samples(physeq_norm, extracted_from=="soil")
mycolors_7 <- c("#c12826", "#ca7c4a", "#628c46", "#48534c", "#1185c0", "#664e93", "#ca8096")
vPhy_ord_extract <- plot_ordination(GP.ord_soil, GP.ord, color="treatment")+
  geom_point(size=1, na.rm=T)+
  scale_colour_manual(values=mycolors_7)+
  stat_ellipse(type="norm", na.rm=T)+
  xlim(-.15, .5)+
  ylim(-.3, .1)+
  theme(
    text=element_text(size=8),
    legend.background=element_rect(size=0.2,
                                   linetype="solid",
                                   colour ="black"),
    legend.position="bottom",
    legend.key.size = unit(.1, 'lines'),
    legend.text=element_text(size=8))+
  scale_fill_continuous(guide = guide_legend())
vPhy_ord_extract
fExportSVG(vPhy_ord_extract, vDirPhyloseq, "fig_2.4_ITS_ordination_soil_legend.svg", vfW=vW_ord, vfH=vH_ord)

## PLOT: Single Soil subset (FIGURE 2.4, no legend)
GP.ord_soil <- subset_samples(physeq_norm, extracted_from=="soil")
mycolors_7 <- c("#c12826", "#ca7c4a", "#628c46", "#48534c", "#1185c0", "#664e93", "#ca8096")
vPhy_ord_extract <- plot_ordination(GP.ord_soil, GP.ord, color="treatment")+
  geom_point(size=1, na.rm=T)+
  scale_colour_manual(values=mycolors_7)+
  stat_ellipse(type="norm", na.rm=T)+
  xlim(-.15, .5)+
  ylim(-.3, .1)+
  theme(
    text=element_text(size=10),
    legend.background=element_rect(size=0.2,
                                   linetype="solid",
                                   colour ="black"),
    legend.position="none")+
  scale_fill_continuous(guide = guide_legend())
vPhy_ord_extract
fExportSVG(vPhy_ord_extract, vDirPhyloseq, "fig_2.4_ITS_ordination_soil_nolegend.svg", vfW=vW_ord, vfH=vH_ord)



# |==================| ----------------------------------------------------
# DESeq2 ------------------------------------------------------------------

# * Initialize ------------------------------------------------------------

## Initialize DESeq2 object (formula is placeholder, updated when "groups" defined)
## Note: Non-normalized counts are used
diagdds = phyloseq_to_deseq2(physeq_raw, design = ~ group)

## Account for 0 counts
diagdds <- estimateSizeFactors(diagdds, type="poscounts")

## Run DESeq2
dds = DESeq(diagdds, test="Wald", fitType="mean")



# * Perform Comparisons ---------------------------------------------------

fExtract_comparisons <- function(vComparisons, vLevel, vSub_padj=.05, vSub_log2=1.5, vPhyseq, vDESeq2, vFolder, vLabel, vLabel_Fus="fusarium.xlsx") {
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
  ## SUMMARY TABLE: Initialize data frame
  vDF_Sum <- data.frame(matrix(nrow=0, ncol=9))
  colnames(vDF_Sum) <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "OTU", "query", "subject")
  # vDF_Sum$query   <- substr(vComparisons[,1],1,2)
  # vDF_Sum$subject <- substr(vComparisons[,2],1,2)
  
  ## Define Output Folder for individual figures/data sets
  vFolderPath <- file.path(vFolder, vLabel)
  ## Define comparison data file
  vOutfile = paste(vLabel,".xlsx", sep="")
  ## Remove current version of output excel file
  fTryRemoveFile(file.path(vFolderPath,vOutfile))
  ## Initialize data frame for saving fusarium statistics
  vDF_Fus <- data.frame(matrix(nrow=nrow(vComparisons), ncol=6))
  rownames(vDF_Fus) <- paste(vComparisons[,1],vComparisons[,2], sep="_")
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
    ## Collect Fusarium data
    rownames(vDF_Fus[n,]) <- vName
    if (any(rownames(vSubset) %in% "Fusarium")) {
      vDF_Fus[n,] = as.data.frame(subset(vSubset, rownames(vSubset) %in% "Fusarium"))
    } else {
      vDF_Fus[n,] = NA
    }
    ## Export to excel
    if (nrow(vSubset)>0) {
      ## SUMMARY TABLE: Define + Add current comparison
      vSubset2 <- as.data.frame(vSubset)
      vSubset2$query   <- c(rep(substr(vComparisons[n,1],1,2), nrow(vSubset2)))
      vSubset2$subject <- c(rep(substr(vComparisons[n,2],1,2), nrow(vSubset2)))
      vSubset2$OTU <- row.names(vSubset2)
      vDF_Sum <- rbind(vDF_Sum, vSubset2)
      ## Append taxonomy table to results
      vSubset <- cbind(as(vSubset, "data.frame"), as(tax_table(vPhyseq)[rownames(vSubset), ], "matrix"))
      ## Export
      fExportXLSX(vfXLSX=vSubset, vfFolder=vFolderPath, vfName=vOutfile, vfSheet=vName)
      print(paste("Exported:      ",vName))
      ## PLOT
      x = tapply(vSubset$log2FoldChange, vSubset$genus, function(x) max(x))
      x = sort(x, TRUE)
      vSubset$genus = factor(as.character(vSubset$genus), levels=names(x))
      vSVG <- ggplot(vSubset, aes(x=genus, y=log2FoldChange)) + geom_point(size=4) + 
        theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
        ggtitle(vName)+
        geom_hline(yintercept=0, color="red")
      fExportSVG(vSVG, vFolderPath, vFilename_svg)
    } else {
      ## Export
      fExportXLSX(vfXLSX="NONE", vfFolder=vFolderPath, vfName=vOutfile, vfSheet=vName, vfRow=F,vfCol=F)
      print(paste("No Differences:",vName))
    }
  }
  ## SUMMARY TABLE: Export
  fExportXLSX(vfXLSX=vDF_Sum,
              vfFolder=vFolder,
              vfName="ITS_DESeq2_summary.xlsx",
              vfSheet=vLabel,
              vfRow=T)
  ## Export Fusarium Data
  fExportXLSX(vfXLSX=vDF_Fus, vfFolder=vFolder, vfName=vLabel_Fus, vfSheet=vLabel, vfRow=T,vfCol=T)
}

fIdentify_wthn <- function(vfLocations, vfCombos) {
  "
  FUNCTION:
      Identifies treatments with data for 2+ extraction locations and returns a
      data frame of all between-location comparisons.
  VARIABLES:
      > vfLocations   (list) Unique extraction locations
      > vfCombos      (list) All available treatment x location combinations
  "
  vTreatments <- unique(dds$treatment)
  vCombos_wthn <- matrix(nrow=0 ,ncol=2)
  for (vTrt in vTreatments) {
    ## Subset $group based on current $treatment
    vSub <- vGroups[grepl(vTrt, vGroups, fixed=TRUE)]
    ## Initialize list of combinations
    vComb <- t(vSub)
    ## Identify all combinations
    if (length(vSub)>2) {
      vComb <- t(combn(vSub, 2))
    }
    ## Capture DESeq2 analysis for identified $treatment combinations
    if (length(vComb)>1) {
      vCombos_wthn <- rbind(vCombos_wthn,vComb)
      vRange <- length(vComb)/2
    }
  }
  return(vCombos_wthn)
}

## Define comparisons
vGroups <- levels(colData(dds)$group)
vGroups_plnt <- vGroups[grepl("plnt", vGroups, fixed=TRUE)]
vGroups_soil <- vGroups[grepl("soil", vGroups, fixed=TRUE)]
vGroups_brsw <- vGroups[grepl("brsw", vGroups, fixed=TRUE)]
vCombos_plnt <- t(combn(vGroups_plnt, 2))
vCombos_soil <- t(combn(vGroups_soil, 2))
vCombos_brsw <- t(combn(vGroups_brsw, 2))
vCombos_wthn <- fIdentify_wthn(vfLocations=unique(dds$treatment), vfCombos=vGroups)

## Remove current version of summary output files
fTryRemoveFile(file.path(vDirDESeq2,"fusarium.xlsx"))
fTryRemoveFile(file.path(vDirDESeq2,"ITS_DESeq2_summary.xlsx"))

## Extract Comparisons: Compare treatments extracted from plant
fExtract_comparisons(vComparisons = vCombos_plnt,
                     vLevel       = "group",
                     vPhyseq      = physeq_raw,
                     vDESeq2      = dds,
                     vFolder      = vDirDESeq2,
                     vLabel       = "btwn_plnt")

## Extract Comparisons: Compare treatments extracted from soil
fExtract_comparisons(vComparisons = vCombos_soil,
                     vLevel       = "group",
                     vPhyseq      = physeq_raw,
                     vDESeq2      = dds,
                     vFolder      = vDirDESeq2,
                     vLabel       = "btwn_soil")

## Extract Comparisons: Compare treatments extracted from brsw
fExtract_comparisons(vComparisons = vCombos_brsw,
                     vLevel       = "group",
                     vPhyseq      = physeq_raw,
                     vDESeq2      = dds,
                     vFolder      = vDirDESeq2,
                     vLabel       = "btwn_brsw")

## Extract Comparisons: Compare between extraction locations
fExtract_comparisons(vComparisons = vCombos_wthn,
                     vLevel       = "group",
                     vSub_padj    = 0.05,
                     vSub_log2    = 1.5,
                     vPhyseq      = physeq_raw,
                     vDESeq2      = dds,
                     vFolder      = vDirDESeq2,
                     vLabel       = "wthn_extracted")



# |==================| ----------------------------------------------------
# Microbiome: PERMANOVA ---------------------------------------------------

# * Initialize ------------------------------------------------------------

## Define analysis object as normalized phyloseq data
pseq.rel <- microbiome::transform(physeq_raw, transform="compositional")



# * Perform Comparisons ---------------------------------------------------

fPermanova <- function(vfComparisons, vfPhyloseq, vfFolder, vfName) {
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
    vTrimKeep <- meta(vfPhyloseq)[,"group"] %in% c(vfComparisons[n,1], vfComparisons[n,2])
    vSubset <- prune_samples(vTrimKeep, vfPhyloseq)
    otu <- abundances(vSubset)
    meta <- meta(vSubset)
    ## Analysis: Permanova
    permanova <- adonis(phyloseq::distance(vSubset, method="wunifrac") ~ group, data = meta, permutations=999, method = "bray")
    ## permanova Results
    vDataPerm <- as.data.frame(permanova$aov.tab)
    # P-value
    ## Check that variance homogeneity assumptions hold (to ensure the reliability of the results). Should be insignificant
    dist <- vegdist(t(otu))
    # print(anova(betadisper(dist, meta$group)))
    vDataPerm <- rbind(PERMANOVA="", statistic=colnames(vDataPerm), vDataPerm)
    vDataDisp <- as.data.frame(anova(betadisper(dist, meta$group, add=T)))
    colnames(vDataDisp) <- c("Df", "SumsOfSqs", "MeanSqs", "F.Model", "Pr(>F)")
    vDataDisp$R2 =c("")
    vDataDisp <- rbind(BETA_DISPERSION="", statistic=colnames(vDataDisp), vDataDisp)
    result <- rbind(vDataPerm,vDataDisp)
    result[is.na(result)] <- ""
    ## Populate summary table
    vDF_Sum[n,"Pr(>F)_Beta"] <- vDataDisp["Groups","Pr(>F)"]
    vDF_Sum[n,"Pr(>F)_Perm"] <- vDataPerm["group","Pr(>F)"]
    vDF_Sum[n,"R2_Perm"]     <- vDataPerm["group","R2"]
    
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
              vfName="permanova_summary.xlsx",
              vfSheet=vfName,
              vfRow=F)
}

## Remove current version of permanova summary file
fTryRemoveFile(file.path(vDirMicrobiome,"permanova_summary.xlsx"))

fPermanova(vfComparisons=vCombos_plnt, vfPhyloseq=pseq.rel, vfFolder=vDirMicrobiome, vfName="permanova_btwn_plnt")
fPermanova(vfComparisons=vCombos_soil, vfPhyloseq=pseq.rel, vfFolder=vDirMicrobiome, vfName="permanova_btwn_soil")
fPermanova(vfComparisons=vCombos_brsw, vfPhyloseq=pseq.rel, vfFolder=vDirMicrobiome, vfName="permanova_btwn_brsw")
fPermanova(vfComparisons=vCombos_wthn, vfPhyloseq=pseq.rel, vfFolder=vDirMicrobiome, vfName="permanova_wthn_extracted")