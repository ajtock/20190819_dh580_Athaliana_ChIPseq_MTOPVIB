#!/applications/R/R-3.4.0/bin/Rscript

# Profile mean coverage around genes and random loci

# Usage via Condor submission system on node7:
# csmit -m 20G -c 1 "./gene_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/REC8_pooled/coverage/REC8_MYC_Rep1_ChIP_norm_allchrs_coverage_coord_tab.bed REC8_MYC_Rep1"

# Source functions to be used in this script
source("/projects/ajt200/Rfunctions/covMatrix_DNAmethMatrix_target_ranLoc_R3.4.0.R")

library(EnrichedHeatmap)
library(genomation)
library(regioneR)

args <- commandArgs(trailingOnly = T)
flankSize <- as.numeric(args[1])
flankName <- as.character(args[2])
winSize <- as.numeric(args[3])
covDatPath <- as.character(args[4])
libName <- as.character(args[5])

matDir <- "./matrices/"
plotDir <- "./plots/"
system(paste0("[ -d ", matDir, " ] || mkdir ", matDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Chromosome definitions
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrStart <- c(rep(1, 5))
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))

# Import genes as GRanges object
representive_genes_uniq <- system("ls /projects/ajt200/TAIR10/representative_genes/representative_genes_uniq_fmt_strand.txt", intern=T)
genesGR <- readGeneric(representive_genes_uniq, header = T, strand = 4, meta.col = list(gene_model = 5))
seqlevels(genesGR) <- sub("", "Chr", seqlevels(genesGR))
print("***********genes***********")
print(genesGR)
# Generate GRanges object containing random loci of same number and
# size distribution as genesGR
set.seed(856291)
ranLocGR <- randomizeRegions(genesGR,
                             genome = genome,
                             per.chromosome = TRUE,
                             allow.overlaps = TRUE)

# Specify locations of normalised per base coverage files
libPath <- system(paste0("ls ", covDatPath), intern = T)

# Import coverage files as GRanges objects and assign to library names
covGR <- readGeneric(libPath, meta.col = list(coverage = 4))
assign(paste0(libName), covGR)

# Define matrix and column mean coverage outfile (mean profiles)
outDF <- list(paste0(matDir, libName,
                     "_norm_cov_feature_smoothed_target_and_", flankName, "_flank_dataframe.txt"),
              paste0(matDir, libName,
                     "_norm_cov_ranLoc_smoothed_target_and_", flankName, "_flank_dataframe.txt"))
outDFcolMeans <- list(paste0(matDir, libName,
                             "_norm_cov_feature_smoothed_target_and_", flankName, "_flank_dataframe_colMeans.txt"),
                      paste0(matDir, libName,
                             "_norm_cov_ranLoc_smoothed_target_and_", flankName, "_flank_dataframe_colMeans.txt"))

# Run covMatrix() function on each coverage GRanges object to obtain matrices
## containing normalised coverage values around target and random loci
covMatrix(signal = covGR,
          feature = genesGR,
          ranLoc = ranLocGR,
          featureSize = mean(width(genesGR)),
          flankSize = flankSize,
          winSize = winSize,
          outDF = outDF,
          outDFcolMeans = outDFcolMeans)
print(paste0(libName, " profile calculation complete"))


