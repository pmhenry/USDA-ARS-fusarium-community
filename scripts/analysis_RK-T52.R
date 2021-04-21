# Setup -------------------------------------------------------------------

# * Globals ---------------------------------------------------------------

## Analysis Prefix, used for output folder naming
vOutDir = "output_RK-T52"

vW_bar = 100 ## Figure width,  bar charts (mm)
vH_bar = 89  ## Figure height, bar charts (mm)

vW_ord = 65  ## Figure width,  ordination (mm)
vH_ord = 86  ## Figure height, ordination (mm)



# * Files -----------------------------------------------------------------

## Path lists
vPath_current <- as.list(strsplit(dirname(rstudioapi::getSourceEditorContext()$path), "/")[[1]]) ## Windows parsing
vPath_base    <- vPath_current[ 1:length(vPath_current)-1 ]

## Directory strings
vDir_base      <- do.call('file.path', vPath_base)        ## Directory, base
vDir_data_prep <- file.path(vDir_base, "data_prepared")   ## Directory, output (processed) data
vDir_output    <- file.path(vDir_base, "output_analysis") ## Directory, analysis output

## Paths, Data file
vInfile_counts   <- file.path(vDir_data_prep, "data_RK-T52", "RK-T52_count.csv")    ## count
vInfile_distance <- file.path(vDir_data_prep, "data_RK-T52", "RK-T52_distance.phy") ## distance
vInfile_taxa     <- file.path(vDir_data_prep, "data_RK-T52", "RK-T52_taxonomy.csv") ## taxa
vInfile_sample   <- file.path(vDir_data_prep, "data_RK-T52", "RK-T52_sample.csv")   ## sample



# * Packages --------------------------------------------------------------

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
dfCounts <- read.csv(file      =vInfile_counts, row.names=1, check.names=FALSE) ## count
phyTree  <- read_tree(treefile =vInfile_distance)                               ## distance
dfTaxono <- read.csv(file      =vInfile_taxa, row.names=1)                      ## taxonomy
dfSample <- read.csv(file      =vInfile_sample, row.names=1)                    ## sample



# |==================| ----------------------------------------------------
# Phyloseq ----------------------------------------------------------------
# * Initialize  -----------------------------------------------------------

## Phyloseq object: Non-normalized (125 samples x 34 taxa)
physeq_FOF <- phyloseq(otu_table(dfCounts, taxa_are_rows=T), ## Counts
                       tax_table(as.matrix(dfTaxono)),       ## Taxonomy
                       sample_data(dfSample),                ## Sample
                       phyTree)                              ## Distance

## Format: Remove low-abundance samples and taxa
physeq_FOF <-prune_samples(sample_sums(physeq_FOF)>=1, physeq_FOF)
physeq_FOF <-prune_taxa(taxa_sums(physeq_FOF)>=1, physeq_FOF)



# * Plot: Stacked Bar Chart ----------------------------------------

## Create phyloseq object with only top N oligotypes included (14 identified above 1% in EF analysis)
TopNOTUs <- names(sort(taxa_sums(physeq_FOF), TRUE)[1:14]) 
physeq_toptax <- prune_taxa(TopNOTUs, physeq_FOF)

## Define palette
mycolors <- c("#850027", "#c12826", "#aa4f2b", "#ca7c4a", "#feb70e", "#ece372", "#1bd17c", "#628c46", "#664e93", "#ca8096", "#6a1154", "#375670", "#1185c0", "#3e30ff")

## Bar Graph, Individual Replicates (FIGURE 4.1, legend)
vPhy_Bar <- plot_bar(physeq_toptax, x="REP", y="Abundance", fill="OTU")+
  facet_grid(EXP~CROP, scales="free_y")+
  xlab("Repetition")+
  ylab("CFU/g-estimated abundance (x10^5)")+
  scale_fill_manual(values=mycolors)+
  scale_x_continuous(breaks=seq(0, 6, 1))+
  theme(
    text=element_text(size=8),
    axis.text.x = element_text(angle = 0),
    legend.key.size = unit(.75, 'lines'),
    legend.text=element_text(size=8),
    panel.spacing.x=unit(.1, "lines"),
    panel.spacing.y=unit(.2, "lines"),
    axis.title.x=element_blank())+
  guides(fill=guide_legend(ncol=1))+
  scale_y_continuous(labels=function(x)x/10000)
vPhy_Bar
fExportSVG(vPhy_Bar, vDirPhyloseq, "fig_4.1_RK-T52_barplot_soil_legend.svg", vfW=vW_bar, vfH=vH_bar)

## Bar Graph, Individual Replicates (FIGURE 4.1, no legend)
vPhy_Bar <- plot_bar(physeq_toptax, x="REP", y="Abundance", fill="OTU")+
  facet_grid(EXP~CROP, scales="free_y")+
  xlab("Repetition")+
  ylab("CFU/g-estimated abundance (x10^5)")+
  scale_fill_manual(values=mycolors)+
  scale_x_continuous(breaks=seq(0, 6, 1))+
  theme(
    text=element_text(size=8),
    axis.text.x = element_text(angle = 0),
    legend.position="none",
    panel.spacing.x=unit(.1, "lines"),
    panel.spacing.y=unit(.2, "lines"),
    axis.title.x=element_blank())+
  guides(fill=guide_legend(ncol=1))+
  scale_y_continuous(labels=function(x)x/10000)
vPhy_Bar
fExportSVG(vPhy_Bar, vDirPhyloseq, "fig_4.1_RK-T52_barplot_soil_nolegend.svg", vfW=vW_bar, vfH=vH_bar)

## Bar Graph, Individual Replicates (FIGURE 4.1, no legend, log2)
vPhy_Bar <- plot_bar(physeq_toptax, x="REP", y="Abundance", fill="OTU")+
  facet_grid(EXP~CROP, scales="free_y")+
  xlab("Repetition")+
  ylab("CFU/g-estimated abundance (sqrt)")+
  scale_fill_manual(values=mycolors)+
  scale_x_continuous(breaks=seq(0, 6, 1))+
  theme(
    text=element_text(size=8),
    axis.text.x = element_text(angle = 0),
    legend.position="none",
    panel.spacing.x=unit(.1, "lines"),
    panel.spacing.y=unit(.2, "lines"),
    axis.title.x=element_blank())+
  guides(fill=guide_legend(ncol=1))+
  scale_y_sqrt()
vPhy_Bar
fExportSVG(vPhy_Bar, vDirPhyloseq, "fig_4.1_RK-T52_barplot_soil_nolegend_log2.svg", vfW=vW_bar, vfH=vH_bar)



# * Plot: Ordination ----------------------------------------------

## NMDS ordination on the bray-curtis distance.
GP.ord <- ordinate(physeq_FOF, "NMDS", "wunifrac")

## PLOT: Single Plot (FIGURE 4.2, no legend)
mycolors_7 <- c("#c12826", "#ca7c4a", "#628c46", "#48534c", "#1185c0", "#664e93", "#ca8096")
vPhy_ord_extract <- plot_ordination(physeq_FOF, GP.ord, color="CROP")+
  geom_point(size=1, na.rm=T)+
  scale_colour_manual(values=mycolors_7)+
  stat_ellipse(type="norm", na.rm=T)+
  theme(
    text=element_text(size=10),
    legend.background=element_rect(size=0.2,
                                   linetype="solid",
                                   colour ="black"),
    legend.position="none")+
  scale_fill_continuous(guide = guide_legend())
vPhy_ord_extract
fExportSVG(vPhy_ord_extract, vDirPhyloseq, "fig_4.2_ordination_nolegend.svg", vfW=vW_ord, vfH=vH_ord)



# |==================| ----------------------------------------------------
# Microbiome: PERMANOVA ---------------------------------------------------
# * Initialize ------------------------------------------------------------

vComparisons <- t(combn(unique(dfSample$CROP), 2))



# * Perform Comparisons ---------------------------------------------------

## Define analysis object as normalized phyloseq data
pseq.rel <- microbiome::transform(physeq_FOF, "compositional")

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
    vTrimKeep <- meta(vfPhyloseq)[,"CROP"] %in% c(vfComparisons[n,1], vfComparisons[n,2])
    vSubset <- prune_samples(vTrimKeep, vfPhyloseq)
    otu <- abundances(vSubset)
    meta <- meta(vSubset)
    ## Analysis: Permanova
    permanova <- adonis(phyloseq::distance(vSubset, method="wunifrac") ~ CROP, data = meta, permutations=999)
    ## permanova Results
    vDataPerm <- as.data.frame(permanova$aov.tab)
    # P-value
    ## Check that variance homogeneity assumptions hold (to ensure the reliability of the results). Should be insignificant
    dist <- vegdist(t(otu))
    # print(anova(betadisper(dist, meta$CROP)))
    vDataPerm <- rbind(PERMANOVA="", statistic=colnames(vDataPerm), vDataPerm)
    vDataDisp <- as.data.frame(anova(betadisper(dist, meta$CROP, add=T)))
    colnames(vDataDisp) <- c("Df", "SumsOfSqs", "MeanSqs", "F.Model", "Pr(>F)")
    vDataDisp$R2 =c("")
    vDataDisp <- rbind(BETA_DISPERSION="", statistic=colnames(vDataDisp), vDataDisp)
    result <- rbind(vDataPerm,vDataDisp)
    result[is.na(result)] <- ""
    
    ## Populate summary table
    vDF_Sum[n,"Pr(>F)_Beta"] <- vDataDisp["Groups","Pr(>F)"]
    vDF_Sum[n,"Pr(>F)_Perm"] <- vDataPerm["CROP","Pr(>F)"]
    vDF_Sum[n,"R2_Perm"]     <- vDataPerm["CROP","R2"]
    
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

fPermanova(vfComparisons=vComparisons, vfPhyloseq=pseq.rel, vfFolder=vDirMicrobiome, vfName="permanova_actual_FOF")