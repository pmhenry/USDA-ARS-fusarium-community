# Notes -------------------------------------------------------------------

"
FUNCTION:
    This script is used to process the raw data in order to generate 
consistent sets of data files (counts, sample, taxonomy, distance) for each
analysis.
"



# Functions ---------------------------------------------------------------

fTryCreatePath <- function(vPath) {
  dir.create(file.path(vPath), recursive=T, showWarnings=F)
}



# |==================| ----------------------------------------------------
# Initialize --------------------------------------------------------------


# * Packages --------------------------------------------------------------

library(ape)
library(phyloseq)
library(plyr)
library(dplyr)
library(xlsx)



# * Files -----------------------------------------------------------------

## Path lists
vPath_current <- as.list(strsplit(dirname(rstudioapi::getSourceEditorContext()$path), "/")[[1]]) ## Windows parsing
vPath_base    <- vPath_current[ 1:length(vPath_current)-1 ]

## Directory strings
vDir_base      <- do.call('file.path', vPath_base)      ## Directory, base
vDir_data_raw  <- file.path(vDir_base, "data_raw")      ## Directory, input (raw) data
vDir_data_prep <- file.path(vDir_base, "data_prepared") ## Directory, output (processed) data

## Create Directories
fTryCreatePath(file.path(vDir_data_prep, "data_EF"))
fTryCreatePath(file.path(vDir_data_prep, "data_ITS"))
fTryCreatePath(file.path(vDir_data_prep, "data_RK-T52"))
fTryCreatePath(file.path(vDir_data_prep, "data_comb-wsfs"))




# |==================| ----------------------------------------------------
# Data: EF ----------------------------------------------------------------

# * Import ----------------------------------------------------------------

## Raw Data Files
vInfile_EF_counts   <- file.path(vDir_data_raw, "data_EF", "RK-ef.oligotype-abundance.01-13-20.txt") ## count
vInfile_EF_distance <- file.path(vDir_data_raw, "data_EF", "EF_newick.phy")                          ## distance

## Import: Count file
vRaw_count_EF <- read.table(file=vInfile_EF_counts, header=T, row.names=1)

## Import: Distance object
phyTree_EF <- read_tree(vInfile_EF_distance)



# * Format: Count & Distance ----------------------------------------------

## Remove "levels" column
dfCounts_EF <- vRaw_count_EF[ , !names(vRaw_count_EF) %in% c("Level")]
## Transpose so that rows=samples
dfCounts_EF <- data.frame(t(dfCounts_EF))

## Match taxa names between counts data and tree data
vEFColnames_cnt <- colnames(dfCounts_EF)
vEFColnames_cnt <- gsub("Fusarium_fujikuroi_As", "FfujAs__M00384", vEFColnames_cnt) ## Match "_" (x2) pattern with rest of names
vEFColnames_cnt <- gsub("Fusarium_sambucinum_D", "FsamD__M00384",  vEFColnames_cnt) ## Match "_" (x2) pattern with rest of names
vEFColnames_cnt <- gsub("Fusarium_tricinctum",   "Ftri__M00384",   vEFColnames_cnt) ## Match "_" (x2) pattern with rest of names
vEFColnames_cnt <- gsub("[.]", "-", vEFColnames_cnt)

vEF_Split_Coun <- data.frame(t(data.frame(strsplit(vEFColnames_cnt,      split="_"))))
vEF_Split_Tree <- data.frame(t(data.frame(strsplit(phyTree_EF$tip.label, split="_"))))

vEF_Split_Coun$norm <- paste0(vEF_Split_Coun$X1, "_", vEF_Split_Coun$X3)
vEF_Split_Tree$norm <- paste0(vEF_Split_Tree$X1, "_", vEF_Split_Tree$X2)

## COUNT data ordering
vEF_Comb_Names_Tree <- join(vEF_Split_Coun, vEF_Split_Tree, by="norm") ## Order by original Tree list
vEF_Comb_Names_Tree$final <- paste0(vEF_Comb_Names_Tree$X1, "_", vEF_Comb_Names_Tree$X2)

## TREE data ordering
vEF_Comb_Names_Coun <- join(vEF_Split_Tree, vEF_Split_Coun, by="norm") ## Order by original Counts list
vEF_Comb_Names_Coun$final <- paste0(vEF_Comb_Names_Coun$X1, "_", vEF_Comb_Names_Coun[,5])

## Rename Data with standardized taxa names:
phyTree_EF$tip.label  <- vEF_Comb_Names_Coun$final
colnames(dfCounts_EF) <- vEF_Comb_Names_Tree$final



# * Format: Taxonomy ------------------------------------------------------

## Assign taxonomy at level of "classification": OTU oligotype names
dfTaxono_EF <- data.frame("taxa_OTU"=colnames(dfCounts_EF), row.names=colnames(dfCounts_EF))



# * Format: Sample --------------------------------------------------------

## Create metadata data frame
vEF_Columns_raw <- row.names(dfCounts_EF)
vEF_Columns <- gsub("_ef", "", vEF_Columns_raw)

## Initialize Columns
vEF_Col_ext <- seq(length(vEF_Columns))
vEF_Col_exp <- seq(length(vEF_Columns))
vEF_Col_plt <- seq(length(vEF_Columns))
vEF_Col_rep <- seq(length(vEF_Columns))

## Parse columns names for meta data
for (n in 1:length(vEF_Columns)) {
  vEF_Col_ext[n] <- substr(vEF_Columns[n],1,4)
  vEF_Col_exp[n] <- substr(vEF_Columns[n],5,5)
  vEF_Col_plt[n] <- substr(vEF_Columns[n],6,7)
  vEF_Col_rep[n] <- substr(vEF_Columns[n],8,8)
}

## Create meta data data frame
dfSample_EF <- data.frame(extracted_from=vEF_Col_ext,
                          experiment=vEF_Col_exp,
                          treatment=vEF_Col_plt,
                          repetition=vEF_Col_rep,
                          row.names=vEF_Columns_raw)

## Create column of combined information
dfSample_EF$comb <- paste(dfSample_EF$extracted_from,
                          dfSample_EF$experiment,
                          dfSample_EF$treatment,
                          sep="")

## Create column of combined information (DESeq2)
dfSample_EF$group <- paste(dfSample_EF$treatment,
                           dfSample_EF$extracted_from,
                           sep="")



# * Export ----------------------------------------------------------------

write.csv(x=dfCounts_EF,   file=file.path(vDir_data_prep, "data_EF","EF_count.csv")   ) ## Count
write.tree(phy=phyTree_EF, file=file.path(vDir_data_prep, "data_EF","EF_distance.phy")) ## Distance
write.csv(x=dfTaxono_EF,   file=file.path(vDir_data_prep, "data_EF","EF_taxonomy.csv")) ## Taxonomy
write.csv(x=dfSample_EF,   file=file.path(vDir_data_prep, "data_EF","EF_sample.csv")  ) ## Sample



# |==================| ----------------------------------------------------
# Data: ITS ---------------------------------------------------------------

# * Import ----------------------------------------------------------------

## Raw Data Files
vInfile_ITS_count <- file.path(vDir_data_raw, "data_ITS", "RK-its-unite-genus.abundance.txt")
vInfile_ITS_taxa  <- file.path(vDir_data_raw, "data_ITS","RK-its-unite-genus.taxa_info.txt")
vInfile_ITS_dist  <- file.path(vDir_data_raw, "data_ITS","ITS_newick.phy")

## Import
vRaw_count_ITS <- read.table(vInfile_ITS_count, header=T) ## Count
vRaw_taxa_ITS  <- read.table(vInfile_ITS_taxa,  header=T) ## Taxonomy



# * Format: Count ---------------------------------------------------------

## Identify row order, used (but not yet specified) by the taxonomy file
vITS_Rows <- vRaw_count_ITS$Taxon_Name
## Subset counts file to only include count columns
dfCounts_ITS <- vRaw_count_ITS[3:ncol(vRaw_count_ITS)]
row.names(dfCounts_ITS) <- vITS_Rows

## Remove all instances of genus-level unclassified samples
dfCounts_ITS <- dfCounts_ITS[- grep("unidentified",   rownames(dfCounts_ITS)),]
dfCounts_ITS <- dfCounts_ITS[- grep("Incertae_sedis", rownames(dfCounts_ITS)),]



# * Format: Taxonomy ------------------------------------------------------

## Split taxonomy line by classification level
vITS_Split <- strsplit(vRaw_taxa_ITS$Taxon_Name, split=";")
## Identify longest taxonomy list
vITS_max.length <- max(sapply(vITS_Split, length))
## Add NA values to list elements
vITS_Split <- lapply(vITS_Split, function(v) { c(v, rep(NA, vITS_max.length-length(v)))})
## Convert to df
dfTaxono_ITS <- data.frame(do.call(rbind, vITS_Split))
## Name columns taxonomy levels
colnames(dfTaxono_ITS) <- c("domain","phylum","class", "order", "family", "genus")
## Name rows (matches count file)
rownames(dfTaxono_ITS) <- vITS_Rows

## Remove all instances of genus-level unclassified samples
dfTaxono_ITS <- dfTaxono_ITS[- grep("unidentified",   rownames(dfTaxono_ITS)),]
dfTaxono_ITS <- dfTaxono_ITS[- grep("Incertae_sedis", rownames(dfTaxono_ITS)),]
vITS_TaxRows <- rownames(dfTaxono_ITS)
## Remove "g__" from genus column
dfTaxono_ITS <- data.frame(lapply(dfTaxono_ITS, function(x) {gsub(c("g__"), "", x)}), row.names=vITS_TaxRows)



# * Format: Sample --------------------------------------------------------

## Create metadata data frame
vITS_Columns_raw <- colnames(dfCounts_ITS)
vITS_Columns <- gsub("_its", "", vITS_Columns_raw)

## Initialize Columns
vITS_Col_ext <- seq(length(vITS_Columns))
vITS_Col_exp <- seq(length(vITS_Columns))
vITS_Col_plt <- seq(length(vITS_Columns))
vITS_Col_rep <- seq(length(vITS_Columns))

## Parse columns names for meta data
for (n in 1:length(vITS_Columns)) {
  vITS_Col_ext[n] <- substr(vITS_Columns[n],1,4)
  vITS_Col_exp[n] <- substr(vITS_Columns[n],5,5)
  vITS_Col_plt[n] <- substr(vITS_Columns[n],6,7)
  vITS_Col_rep[n] <- substr(vITS_Columns[n],8,8)
}

## Create meta data data frame
dfSample_ITS <- data.frame(extracted_from=vITS_Col_ext,
                           experiment=vITS_Col_exp,
                           treatment=vITS_Col_plt,
                           repetition=vITS_Col_rep,
                           row.names=vITS_Columns_raw)

## Create column of combined information
dfSample_ITS$comb <- paste(dfSample_ITS$extracted_from,
                           dfSample_ITS$experiment,
                           dfSample_ITS$treatment,
                           sep="")

## Create column of combined information (DESeq2)
dfSample_ITS$group <- paste(dfSample_ITS$treatment,
                            dfSample_ITS$extracted_from,
                            sep="")



# * Format: Distance ------------------------------------------------------

## Nothing to be done, this file will be copied and renamed.



# * Export ----------------------------------------------------------------


write.csv(x=dfCounts_ITS, file=file.path(vDir_data_prep, "data_ITS", "ITS_count.csv")   ) ## Count
write.csv(x=dfTaxono_ITS, file=file.path(vDir_data_prep, "data_ITS", "ITS_taxonomy.csv")) ## Taxonomy
write.csv(x=dfSample_ITS, file=file.path(vDir_data_prep, "data_ITS", "ITS_sample.csv")  ) ## Sample
file.copy(from=vInfile_ITS_dist, to=file.path(vDir_data_prep, "data_ITS","ITS_distance.phy")) ## Distance



# |==================| ----------------------------------------------------
# Data: RK-T52 ---------------------------------------------------------

# * Import ----------------------------------------------------------------

## Data File: Fusarium Abundance
vInfile_RK_FOF      <- file.path(vDir_data_raw, "data_RK-T52", "RKcompiled_T52_soilcfusperg.csv")

## Import
vRaw_count_RK     <- read.table(file=vInfile_EF_counts, header=T, row.names=1) ## Counts (same as EF analysis)
vRaw_count_FOF_RK <- read.csv(file=vInfile_RK_FOF, fileEncoding="UTF-8-BOM")   ## Counts, FOF-adjusted



# * Format: Distance ------------------------------------------------------

## Define Phyloseq Tree object (same as EF analysis)
phyTree_RK <- read_tree(vInfile_EF_distance)



# * Format: Count and Sample ----------------------------------------------

## Remove "levels" column
dfCounts_RK_all <- vRaw_count_RK[ , !names(vRaw_count_RK) %in% c("Level")]
## Transpose so that rows=samples
dfCounts_RK_all <- data.frame(t(dfCounts_RK_all))

## Define FOF index names
vRaw_count_FOF_RK$comb <- paste0("soil",vRaw_count_FOF_RK$EXP,vRaw_count_FOF_RK$CROP,vRaw_count_FOF_RK$REP,"_ef")
row.names(vRaw_count_FOF_RK) <- vRaw_count_FOF_RK$comb

## Filter for soil data by merging with FOF
dfMerge_RK <- merge(x=vRaw_count_FOF_RK,
                    y=dfCounts_RK_all,
                    by.x=0, ## index
                    by.y=0) ## index
row.names(dfMerge_RK) <- dfMerge_RK$comb

## Define Count data
dfCounts_RK <- dfMerge_RK %>% select(8:length(dfMerge_RK))

## Define Meta data
dfSample_RK <- dfMerge_RK %>% select(2:7)

## Match taxa names between counts data and tree data
vRK_Colnames_cnt <- colnames(dfCounts_RK)
vRK_Colnames_cnt <- gsub("Fusarium_fujikuroi_As", "FfujAs__M00384", vRK_Colnames_cnt) ## Match "_" (x2) pattern with rest of names
vRK_Colnames_cnt <- gsub("Fusarium_sambucinum_D", "FsamD__M00384",  vRK_Colnames_cnt) ## Match "_" (x2) pattern with rest of names
vRK_Colnames_cnt <- gsub("Fusarium_tricinctum", "Ftri__M00384",     vRK_Colnames_cnt) ## Match "_" (x2) pattern with rest of names
vRK_Colnames_cnt <- gsub("[.]", "-", vRK_Colnames_cnt)

vRK_Split_Coun <- data.frame(t(data.frame(strsplit(vRK_Colnames_cnt, split="_"))))
vRK_Split_Tree <- data.frame(t(data.frame(strsplit(phyTree_RK$tip.label, split="_"))))

vRK_Split_Coun$norm <- paste0(vRK_Split_Coun$X1, "_", vRK_Split_Coun$X3)
vRK_Split_Tree$norm <- paste0(vRK_Split_Tree$X1, "_", vRK_Split_Tree$X2)

## This sounds reversed, but it works: Use for COUNT data ordering
vRK_Comb_Names_Tree <- join(vRK_Split_Coun, vRK_Split_Tree, by="norm") ## Order by original Tree list
vRK_Comb_Names_Tree$final <- paste0(vRK_Comb_Names_Tree$X1, "_", vRK_Comb_Names_Tree$X2)

## This sounds reversed, but it works: Use for TREE data ordering
vRK_Comb_Names_Coun <- join(vRK_Split_Tree, vRK_Split_Coun, by="norm") ## Order by original Counts list
vRK_Comb_Names_Coun$final <- paste0(vRK_Comb_Names_Coun$X1, "_", vRK_Comb_Names_Coun[,5])

## Rename Data with standardized taxa names:
phyTree_RK$tip.label <-  vRK_Comb_Names_Coun$final
colnames(dfCounts_RK) <- vRK_Comb_Names_Tree$final

## Normalize data
dfCounts_RK_t <- as.data.frame(t(dfCounts_RK))
dfCounts_RK_norm <- sweep(x=dfCounts_RK_t,
                          MARGIN=2,
                          STATS=colSums(dfCounts_RK_t),`/`)
## Transform to FOF abundance
dfCounts_RK_FOF <- dfCounts_RK_norm
for(n in 1:ncol(dfCounts_RK_FOF)){
  dfCounts_RK_FOF[,n] <- dfCounts_RK_FOF[,n] * dfSample_RK$abund_allFus_byavg[n]
}



# * Format: Taxonomy ------------------------------------------------------

## Define taxonomy at level of "classification": OTU oligotype names
dfTaxono_RK <- data.frame("taxa_OTU"=colnames(dfCounts_RK), row.names=colnames(dfCounts_RK))



# * Export ----------------------------------------------------------------

write.csv(x=dfCounts_RK_FOF, file=file.path(vDir_data_prep, "data_RK-T52", "RK-T52_count.csv")   ) ## Count
write.tree(phy=phyTree_RK,   file=file.path(vDir_data_prep, "data_RK-T52", "RK-T52_distance.phy")) ## Distance
write.csv(x=dfTaxono_RK,     file=file.path(vDir_data_prep, "data_RK-T52", "RK-T52_taxonomy.csv")) ## Taxonomy
write.csv(x=dfSample_RK,     file=file.path(vDir_data_prep, "data_RK-T52", "RK-T52_sample.csv")  ) ## Sample



# |==================| ----------------------------------------------------
# Data: comb-wsfs ---------------------------------------------------------

# * Import ----------------------------------------------------------------

## Data File: Counts, Bulk Soil '17 & '18
vInfile_co_counts_Bulk <- file.path(vDir_data_raw, "data_comb-wsfs","comb-wsfs_justsoil_OT.txt")

## Data File: Counts, Soil & Root '17
vInfile_co_counts_SlRt <- file.path(vDir_data_raw, "data_comb-wsfs","comb-wsfs_soilandroot2017_OT.txt")

## Data File: Distance matrix file
vInfile_co_distance    <- file.path(vDir_data_raw, "data_comb-wsfs","wsfs17_repseq_align.phy")

## Import: Counts, Bulk Soil '17 & '18
vRaw_count_coBulk <- read.table(
  file=vInfile_co_counts_Bulk,
  header=T,
  comment.char="",
  fill=T)

## Import: Counts, Soil & Root '17
vRaw_count_coSlRt <- read.table(
  file=vInfile_co_counts_SlRt,
  header=T,
  comment.char="",
  fill=T)

## Import: Distance object
phyTree_co <- read_tree(vInfile_co_distance)



# * Format: Count ---------------------------------------------------------

## Format NAs
vRaw_count_coBulk <- data.frame(lapply(vRaw_count_coBulk, function(x) {gsub("#DIV/0!", 0, x)}))
vRaw_count_coSlRt <- data.frame(lapply(vRaw_count_coSlRt, function(x) {gsub("#DIV/0!", 0, x)}))

## Define Index
row.names(vRaw_count_coBulk) <- vRaw_count_coBulk$ID
row.names(vRaw_count_coSlRt) <- vRaw_count_coSlRt$Barcode_ID

dfCounts_coBulk <- vRaw_count_coBulk %>% select(10:length(vRaw_count_coBulk))           ## Parse count columns
dfCounts_coBulk <- mutate_all(dfCounts_coBulk, function(x) as.numeric(as.character(x))) ## Ensure values are numeric
row.names(dfCounts_coBulk) <- row.names(vRaw_count_coBulk)                              ## Restore row names

dfCounts_coSlRt <- vRaw_count_coSlRt %>% select(11:length(vRaw_count_coSlRt))           ## Parse count columns
dfCounts_coSlRt <- mutate_all(dfCounts_coSlRt, function(x) as.numeric(as.character(x))) ## Ensure values are numeric
row.names(dfCounts_coSlRt) <- row.names(vRaw_count_coSlRt)                              ## Restore row names



# * Format: Taxonomy -------------------------------------------------------

dfTaxono_coBulk <- data.frame(taxa=colnames(dfCounts_coBulk), row.names=colnames(dfCounts_coBulk))
dfTaxono_coSlRt <- data.frame(taxa=colnames(dfCounts_coSlRt), row.names=colnames(dfCounts_coSlRt))



# * Format: Sample --------------------------------------------------------

dfSample_coBulk <- vRaw_count_coBulk %>% select(1:9)
dfSample_coSlRt <- vRaw_count_coSlRt %>% select(1:8,10)

## Define combined column, SlRt
dfSample_coSlRt$comb <- paste0(dfSample_coSlRt$full_trtmt, "_", dfSample_coSlRt$sample_type)

## Define combined column, Bulk
dfSample_coBulk$trt_year <- paste0(dfSample_coBulk$full_trtmt, dfSample_coBulk$year)



# * Format: Distance ------------------------------------------------------

## Manually rename columns
vco_Col_OG <- phyTree_co$tip.label
vco_Col_OG <- gsub("Fusarium_cyanostomum", "cyanostomum_all",       vco_Col_OG)
vco_Col_OG <- gsub("dimerum_all",          "Fusarium_dimerum",      vco_Col_OG)
vco_Col_OG <- gsub("foetens_all",          "Fusarium_foetens",      vco_Col_OG)
vco_Col_OG <- gsub("redolens_all",         "Fusarium_redolens",     vco_Col_OG)
vco_Col_OG <- gsub("sambucinumA_all",      "Fusarium_sambucinum_A", vco_Col_OG)
vco_Col_OG <- gsub("sambucinumD_all",      "Fusarium_sambucinum_D", vco_Col_OG)
vco_Col_OG <- gsub("solani2_all",          "Fusarium_solani_2",     vco_Col_OG)
vco_Col_OG <- gsub("tricinctum_all",       "Fusarium_tricinctum",   vco_Col_OG)
vco_Col_OG <- gsub("unknown1_all",         "Fusarium_unknown_1",    vco_Col_OG)
phyTree_co$tip.label <- vco_Col_OG



# * Export ----------------------------------------------------------------

## Count
write.csv(x=dfCounts_coBulk, file=file.path(vDir_data_prep, "data_comb-wsfs", "comb-wsfs_count_Bulk.csv")) ## Bulk Soil '17 & '18
write.csv(x=dfCounts_coSlRt, file=file.path(vDir_data_prep, "data_comb-wsfs", "comb-wsfs_count_SlRt.csv")) ## Soil & Root '17

## Taxonomy
write.csv(x=dfTaxono_coBulk, file=file.path(vDir_data_prep, "data_comb-wsfs", "comb-wsfs_taxonomy_Bulk.csv")) ## Bulk Soil '17 & '18
write.csv(x=dfTaxono_coSlRt, file=file.path(vDir_data_prep, "data_comb-wsfs", "comb-wsfs_taxonomy_SlRt.csv")) ## Soil & Root '17

## Sample
write.csv(x=dfSample_coBulk, file=file.path(vDir_data_prep, "data_comb-wsfs", "comb-wsfs_sample_Bulk.csv")) ## Bulk Soil '17 & '18
write.csv(x=dfSample_coSlRt, file=file.path(vDir_data_prep, "data_comb-wsfs", "comb-wsfs_sample_SlRt.csv")) ## Soil & Root '17

## Distance
write.tree(phy=phyTree_co, file=file.path(vDir_data_prep, "data_comb-wsfs", "comb-wsfs_distance.phy"))