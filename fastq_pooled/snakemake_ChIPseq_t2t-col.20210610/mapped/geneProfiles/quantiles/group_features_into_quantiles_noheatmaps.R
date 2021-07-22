#!/applications/R/R-4.0.0/bin/Rscript

#
# Divide features into quantiles based on mean log2(libName ChIP/control)
# in a given feature region (e.g., promoters).
# Extract and save feature IDs for each quantile for further analyses
# (e.g., GO enrichment and average + 95% CI profile plotting).
# Calculate mean winName-scaled recombination rate (cM/Mb) from
# the promoter to the terminator of each feature.
# Plot feature quantile heatmaps for various genomics datasets
# Plot feature quantile recombination rate densities in a
# heat map or violin plot
#

# Usage:
# /applications/R/R-4.0.0/bin/Rscript group_features_into_quantiles_noheatmaps.R WT_MTOPVIB_HA_Rep1_ChIP '20190819_dh580_Athaliana_ChIPseq_MTOPVIB/fastq_pooled/snakemake_ChIPseq_t2t-col.20210610' both 'Chr1,Chr2,Chr3,Chr4,Chr5' 2200 2000 2kb '2 kb' 10 10bp promoters 4 t2t-col.20210610

#libName <- "WT_MTOPVIB_HA_Rep1_ChIP"
#dirName <- "20190819_dh580_Athaliana_ChIPseq_MTOPVIB/fastq_pooled/snakemake_ChIPseq_t2t-col.20210610"
#align <- "both"
#chrName <- unlist(strsplit("Chr1,Chr2,Chr3,Chr4,Chr5",
#                           split = ","))
#bodyLength <- 2200
#upstream <- 2000
#downstream <- 2000
#flankName <- "2kb"
#flankNamePlot <- "2 kb"
#binSize <- 10
#binName <- "10bp"
#region <- "promoters"
#quantiles <- 4
#refbase <- "t2t-col.20210610"

args <- commandArgs(trailingOnly = T)
libName <- args[1]
dirName <- args[2]
align <- args[3]
chrName <- unlist(strsplit(args[4],
                               split = ","))
bodyLength <- as.numeric(args[5])
upstream <- as.numeric(args[6])
downstream <- as.numeric(args[6])
flankName <- args[7]
flankNamePlot <- args[8]
binSize <- as.numeric(args[9])
binName <- args[10]
region <- args[11]
quantiles <- as.numeric(args[12])
refbase <- args[13]

options(stringsAsFactors = F)
library(EnrichedHeatmap)
library(png)
#library(Cairo)
library(RColorBrewer)
library(circlize)
library(GenomicRanges)
library(dplyr)
library(parallel)
library(doParallel)
registerDoParallel(cores = detectCores())
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

outDir <- paste0("quantiles_by_log2_", libName,
                 "_control_in_", region, "/")
plotDir <- paste0(outDir, "plots/")
system(paste0("[ -d ", outDir, " ] || mkdir ", outDir))
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Genomic definitions
fai <- read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/", refbase, ".fa.fai"), header = F)
chrs <- fai$V1[which(fai$V1 %in% chrName)]
chrLens <- fai$V2[which(fai$V1 %in% chrName)]
regionGR <- GRanges(seqnames = chrs,
                    ranges = IRanges(start = 1,
                                     end = chrLens),
                    strand = "*")

# Define centromere GRanges
CENstart <- c(14841110,3823792,13597188,4203902,11784131)[which(fai$V1 %in% chrName)]
CENend <- c(17559778,6045243,15733925,6977949,14551809)[which(fai$V1 %in% chrName)]
CENGR <- GRanges(seqnames = chrs,
                 ranges = IRanges(start = CENstart,
                                  end = CENend),
                 strand = "*")

# Load ChIP matrix
# feature
ChIP_featureMat <- lapply(seq_along(chrName), function(y) {
  as.matrix(read.table(paste0("/home/ajt200/analysis/",
                              dirName,
                              "/mapped/geneProfiles/matrices/",
                              libName,
                              "_MappedOn_", refbase, "_lowXM_", align, "_sort_norm_genes_in_",
                              chrName[y], "_matrix_bin", binName,
                              "_flank", flankName, ".tab"),
                       header = F, skip = 3))
})
# If features from multiple chromosomes are to be analysed,
# concatenate the corresponding feature coverage matrices
if(length(chrName) > 1) {
 ChIP_featureMat <- do.call(rbind, ChIP_featureMat)
} else {
 ChIP_featureMat <- ChIP_featureMat[[1]]
}

# ranLoc
ChIP_ranLocMat <- lapply(seq_along(chrName), function(y) {
  as.matrix(read.table(paste0("/home/ajt200/analysis/",
                              dirName,
                              "/mapped/geneProfiles/matrices/",
                              libName,
                              "_MappedOn_", refbase, "_lowXM_", align, "_sort_norm_genes_in_",
                              chrName[y], "_ranLoc_matrix_bin", binName,
                              "_flank", flankName, ".tab"),
                       header = F, skip = 3))
})
# If ranLocs from multiple chromosomes are to be analysed,
# concatenate the corresponding ranLoc coverage matrices
if(length(chrName) > 1) {
 ChIP_ranLocMat <- do.call(rbind, ChIP_ranLocMat)
} else {
 ChIP_ranLocMat <- ChIP_ranLocMat[[1]]
}

# Load control matrices
controlNames <- c(
                  "WT_REC8_Myc_Rep1_input"
#                  "WT_CENH3_Rep1_input_SRR4430555",
#                  "WT_H3K9me2_Rep1_input"
#                  "H2AW_input_SRR5298544",
#                  "WT_gDNA_Rep1",
#                  "WT_gDNA_Rep1_R1",
#                  "map_K40_E2",
#                  "map_K45_E2",
#                  "map_K50_E2",
#                  "map_K150_E4",
#                  "map_K200_E4",
#                  "map_K300_E4"
                 )
controlNamesDir <- c(
                     paste0("REC8_pooled/snakemake_ChIPseq_", refbase)
#                    paste0("CENH3_seedlings_Maheshwari_Comai_2017_GenomeRes/snakemake_ChIPseq", refbase),
#                    paste0("170101_Chris_H3K9me2_ChIP/WT/snakemake_ChIPseq", refbase),
#                    paste0("HTA6_HTA7_leaf_Lorkovic_Berger_2017_CurrBiol/snakemake_ChIPseq", refbase),
#                    paste0("150701_Natasha_gDNA/WT/snakemake_ChIPseq", refbase),
#                    paste0("150701_Natasha_gDNA/WT/R1/snakemake_SPO11oligos", refbase),
#                     rep(paste0("nanopore/", refbase, "/genmap_mappability"), 6)
                    )
controlNamesPlot <- c(
                      "Input (REC8)"
#                      "Input (CENH3)",
#                      "Input (H3K9me2)"
#                      "Input (MNase)",
#                      "PE gDNA",
#                      "SE gDNA",
#                      "k=40 e=2",
#                      "k=45 e=2",
#                      "k=50 e=2",
#                      "k=150 e=4",
#                      "k=200 e=4",
#                      "k=300 e=4"
                     )
controlDirs <- sapply(seq_along(controlNamesDir), function(x) {
  paste0("/home/ajt200/analysis/", controlNamesDir[x],
         "/mapped/")
})

## control
# feature
control_featureMats <- mclapply(seq_along(controlNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    if( grepl("map_K", controlNames[x]) ) {
      as.matrix(read.table(paste0(controlDirs[x], "geneProfiles/matrices/",
                                  controlNames[x],
                                  "_MappedOn_", refbase, "_genes_in_",
                                  chrName[y], "_matrix_bin", binName, "_flank", flankName, ".tab"),
                           header = F, skip = 3))
    } else {
      as.matrix(read.table(paste0(controlDirs[x], "geneProfiles/matrices/",
                                  controlNames[x],
                                  "_MappedOn_", refbase, "_lowXM_both_sort_norm_genes_in_",
                                  chrName[y], "_matrix_bin", binName, "_flank", flankName, ".tab"),
                           header = F, skip = 3))
      }
  })
}, mc.cores = length(controlNames))
# If features from multiple chromosomes are to be analysed,
# concatenate the corresponding feature coverage matrices
control_featureMats <- mclapply(seq_along(control_featureMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, control_featureMats[[x]])
  } else {
    control_featureMats[[x]][[1]]
  }
}, mc.cores = length(control_featureMats))

# ranLoc
control_ranLocMats <- mclapply(seq_along(controlNames), function(x) {
  lapply(seq_along(chrName), function(y) {
    if( grepl("map_K", controlNames[x]) ) {
      as.matrix(read.table(paste0(controlDirs[x], "geneProfiles/matrices/",
                                  controlNames[x],
                                  "_MappedOn_", refbase, "_genes_in_",
                                  chrName[y], "_ranLoc_matrix_bin", binName, "_flank", flankName, ".tab"),
                           header = F, skip = 3))
    } else {
      as.matrix(read.table(paste0(controlDirs[x], "geneProfiles/matrices/",
                                  controlNames[x],
                                  "_MappedOn_", refbase, "_lowXM_both_sort_norm_genes_in_",
                                  chrName[y], "_ranLoc_matrix_bin", binName, "_flank", flankName, ".tab"),
                           header = F, skip = 3))
      }
  })
}, mc.cores = length(controlNames))
# If features from multiple chromosomes are to be analysed,
# concatenate the corresponding feature coverage matrices
control_ranLocMats <- mclapply(seq_along(control_ranLocMats), function(x) {
  if(length(chrName) > 1) {
    do.call(rbind, control_ranLocMats[[x]])
  } else {
    control_ranLocMats[[x]][[1]]
  }
}, mc.cores = length(control_ranLocMats))

# Conditionally calculate log2(ChIP/control)
# for each matrix depending on library
# feature
log2ChIP_featureMat <-
if ( grepl("CENH3", libName) ) {
  print(paste0(libName, " library; using ", controlNames[20], " for log2((ChIP+1)/(input+1)) calculation"))
  log2((ChIP_featureMat+1)/(control_featureMats[[20]]+1))
} else if ( grepl("H3K9me2", libName) ) {
  print(paste0(libName, " library; using ", controlNames[30], " for log2((ChIP+1)/(input+1)) calculation"))
  log2((ChIP_featureMat+1)/(control_featureMats[[30]]+1))
} else if ( grepl("MNase", libName) ) {
  print(paste0(libName, " library; using ", controlNames[40], " for log2((MNase+1)/(input+1)) calculation"))
  log2((ChIP_featureMat+1)/(control_featureMats[[40]]+1))
} else if ( grepl("SPO11oligos", libName) ) {
  print(paste0(libName, " library; using ", controlNames[50], " for log2((SPO11-1-oligos+1)/(input+1)) calculation"))
  log2((ChIP_featureMat+1)/(control_featureMats[[50]]+1))
} else {
  print(paste0(libName, " library; using ", controlNames[1], " for log2((ChIP+1)/(input+1)) calculation"))
  log2((ChIP_featureMat+1)/(control_featureMats[[1]]+1))
}

# ranLoc
log2ChIP_ranLocMat <-
if ( grepl("CENH3", libName) ) {
  print(paste0(libName, " library; using ", controlNames[20], " for log2((ChIP+1)/(input+1)) calculation"))
  log2((ChIP_ranLocMat+1)/(control_ranLocMats[[20]]+1))
} else if ( grepl("H3K9me2", libName) ) {
  print(paste0(libName, " library; using ", controlNames[30], " for log2((ChIP+1)/(input+1)) calculation"))
  log2((ChIP_ranLocMat+1)/(control_ranLocMats[[30]]+1))
} else if ( grepl("MNase", libName) ) {
  print(paste0(libName, " library; using ", controlNames[40], " for log2((MNase+1)/(input+1)) calculation"))
  log2((ChIP_ranLocMat+1)/(control_ranLocMats[[40]]+1))
} else if ( grepl("SPO11oligos", libName) ) {
  print(paste0(libName, " library; using ", controlNames[50], " for log2((SPO11-1-oligos+1)/(input+1)) calculation"))
  log2((ChIP_ranLocMat+1)/(control_ranLocMats[[50]]+1))
} else {
  print(paste0(libName, " library; using ", controlNames[1], " for log2((ChIP+1)/(input+1)) calculation"))
  log2((ChIP_ranLocMat+1)/(control_ranLocMats[[1]]+1))
}

# Extract region for ordering of features (adjust promoter/terminator size as necessary)
if( region == "promoters" ) {
  log2ChIP_featureMatRegion <- log2ChIP_featureMat[,(((upstream-1000)/binSize)+1):(upstream/binSize)]
  log2ChIP_ranLocMatRegion <- log2ChIP_ranLocMat[,(((upstream-1000)/binSize)+1):(upstream/binSize)]
} else if ( region == "terminators" ) {
  log2ChIP_featureMatRegion <- log2ChIP_featureMat[,(((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(1000/binSize))]
  log2ChIP_ranLocMatRegion <- log2ChIP_ranLocMat[,(((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(1000/binSize))]
} else if ( region == "bodies" ) {
  log2ChIP_featureMatRegion <- log2ChIP_featureMat[,((upstream/binSize)+1):((upstream+bodyLength)/binSize)]
  log2ChIP_ranLocMatRegion <- log2ChIP_ranLocMat[,((upstream/binSize)+1):((upstream+bodyLength)/binSize)]
} else if ( region == "genes" ) {
  log2ChIP_featureMatRegion <- log2ChIP_featureMat[,(((upstream-1000)/binSize)+1):(((upstream+bodyLength)/binSize)+(1000/binSize))]
  log2ChIP_ranLocMatRegion <- log2ChIP_ranLocMat[,(((upstream-1000)/binSize)+1):(((upstream+bodyLength)/binSize)+(1000/binSize))]
} else {
  print("The region name provided does not match 'promoters', 'terminators', 'bodies', or 'genes'")
}
log2ChIP_featureMatRegionRowMeans <- rowMeans(log2ChIP_featureMatRegion, na.rm = T)
#log2ChIP_featureMatRegionRowMeansSorted <- sort.int(log2ChIP_featureMatRegionRowMeans,
#                                            decreasing = T,
#                                            index.return = T,
#                                            na.last = T)
#log2ChIP_featureMatRegionSorted <- log2ChIP_featureMatRegion[sort.int(log2ChIP_featureMatRegionRowMeans,
#                                                      decreasing = T,
#                                                      index.return = T,
#                                                      na.last = T)$ix,]
#log2ChIP_featureMatSorted <- log2ChIP_featureMat[sort.int(log2ChIP_featureMatRegionRowMeans,
#                                          decreasing = T,
#                                          index.return = T,
#                                          na.last = T)$ix,]
## Replace NAs in log2ChIP_featureMatRegion with 0
#log2ChIP_featureMatRegion[which(is.na(log2ChIP_featureMatRegion))] <- 0

log2ChIP_ranLocMatRegionRowMeans <- rowMeans(log2ChIP_ranLocMatRegion, na.rm = T)
#log2ChIP_ranLocMatRegionRowMeansSorted <- sort.int(log2ChIP_ranLocMatRegionRowMeans,
#                                            decreasing = T,
#                                            index.return = T,
#                                            na.last = T)
#log2ChIP_ranLocMatRegionSorted <- log2ChIP_ranLocMatRegion[sort.int(log2ChIP_ranLocMatRegionRowMeans,
#                                                      decreasing = T,
#                                                      index.return = T,
#                                                      na.last = T)$ix,]
#log2ChIP_ranLocMatSorted <- log2ChIP_ranLocMat[sort.int(log2ChIP_ranLocMatRegionRowMeans,
#                                          decreasing = T,
#                                          index.return = T,
#                                          na.last = T)$ix,]
## Replace NAs in log2ChIP_ranLocMatRegion with 0
#log2ChIP_ranLocMatRegion[which(is.na(log2ChIP_ranLocMatRegion))] <- 0


# Load table of feature coordinates in BED format
features <- lapply(seq_along(chrName), function(x) {
  read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/annotation/genes/", refbase, "_representative_mRNA_",
                    chrName[x], ".bed"),
             header = F)
})
if(length(chrName) > 1) {
  features <- do.call(rbind, features)
} else {
  features <- features[[1]]
}
# Convert 0-based start coordinates (BED)
# into 1-based start coordinates (for output as TSV below)
features[,2] <- features[,2]+1
colnames(features) <- c("chr", "start", "end", "featureID", "score", "strand")
featuresGR <- GRanges(seqnames = features$chr,
                      ranges = IRanges(start = features$start,
                                       end = features$end),
                      strand = features$strand,
                      featureID = features$featureID)

# Load table of ranLoc coordinates in BED format
ranLocs <- lapply(seq_along(chrName), function(x) {
  read.table(paste0("/home/ajt200/analysis/nanopore/", refbase, "/annotation/genes/", refbase, "_representative_mRNA_",
                    chrName[x], "_randomLoci.bed"),
             header = F)
})
if(length(chrName) > 1) {
  ranLocs <- do.call(rbind, ranLocs)
} else {
  ranLocs <- ranLocs[[1]]
}
# Convert 0-based start coordinates (BED)
# into 1-based start coordinates (for output as TSV below)
ranLocs[,2] <- ranLocs[,2]+1
colnames(ranLocs) <- c("chr", "start", "end", "ranLocID", "score", "strand")
ranLocsGR <- GRanges(seqnames = ranLocs$chr,
                     ranges = IRanges(start = ranLocs$start+1,
                                      end = ranLocs$end),
                     strand = ranLocs$strand,
                     ranLocID = ranLocs$ranLocID)

features <- data.frame(features,
                       log2ChIP_featureMatRegionRowMeans = log2ChIP_featureMatRegionRowMeans)
ranLocs <- data.frame(ranLocs,
                      log2ChIP_ranLocMatRegionRowMeans = log2ChIP_ranLocMatRegionRowMeans)

# Group features into quantiles according to decreasing orderingFactor
orderingFactor <- "log2ChIP_featureMatRegionRowMeans"
print(orderingFactor)
# Note that features_DF could be defined on
# line above or below mclapply(), with the same results
features_DF <- data.frame(features,
                          percentile = rank(features[,which(colnames(features) == orderingFactor)]) /
                                       length(features[,which(colnames(features) == orderingFactor)]),
                          quantile = as.character(""))
## Assign 0s to NA values only for coverage data
#if(grepl("_in_", orderingFactor)) {
#  features_DF[,which(colnames(features_DF) == orderingFactor)][
#    which(is.na(features_DF[,which(colnames(features_DF) == orderingFactor)]))] <- 0
#}
quantilesStats <- data.frame()
for(k in 1:quantiles) {
  # First quantile should span 1 to greater than, e.g., 0.75 proportions of features
  if(k < quantiles) {
    features_DF[ !is.na(features_DF[,which(colnames(features_DF) == orderingFactor)]) &
               rank(features_DF[,which(colnames(features_DF) == orderingFactor)]) /
               length(features_DF[,which(colnames(features_DF) == orderingFactor)]) <=
               1-((k-1)/quantiles) &
               rank(features_DF[,which(colnames(features_DF) == orderingFactor)]) /
               length(features_DF[,which(colnames(features_DF) == orderingFactor)]) >
               1-(k/quantiles), ]$quantile <- paste0("Quantile ", k)
  } else {
  # Final quantile should span 0 to, e.g., 0.25 proportions of features
    features_DF[ !is.na(features_DF[,which(colnames(features_DF) == orderingFactor)]) &
               rank(features_DF[,which(colnames(features_DF) == orderingFactor)]) /
               length(features_DF[,which(colnames(features_DF) == orderingFactor)]) <=
               1-((k-1)/quantiles) &
               rank(features_DF[,which(colnames(features_DF) == orderingFactor)]) /
               length(features_DF[,which(colnames(features_DF) == orderingFactor)]) >=
               1-(k/quantiles), ]$quantile <- paste0("Quantile ", k)
  }
  write.table(features_DF[features_DF$quantile == paste0("Quantile ", k),],
              file = paste0(outDir,
                            "quantile", k, "_of_", quantiles,
                            "_by_log2_", libName, "_control_in_", region,
                            "_of_genes_in_t2t-col.20210610_",
                            paste0(chrName, collapse = "_"), ".tsv"),
              quote = FALSE, sep = "\t", row.names = FALSE)
  stats <- data.frame(quantile = as.integer(k),
                      n = as.integer(dim(features_DF[features_DF$quantile == paste0("Quantile ", k),])[1]),
                      mean_width = as.integer(round(mean(
                        (features_DF[features_DF$quantile == paste0("Quantile ", k),]$end -
                         features_DF[features_DF$quantile == paste0("Quantile ", k),]$start) + 1, na.rm = T))),
                      total_width = as.integer(sum(
                        (features_DF[features_DF$quantile == paste0("Quantile ", k),]$end -
                         features_DF[features_DF$quantile == paste0("Quantile ", k),]$start) + 1, na.rm = T)),
                      mean_orderingFactor = as.numeric(mean(features_DF[features_DF$quantile == paste0("Quantile ", k),][,which(colnames(features_DF) == orderingFactor)], na.rm = T)))
  quantilesStats <- rbind(quantilesStats, stats)
}
write.table(quantilesStats,
            file = paste0(outDir,
                          "summary_", quantiles, "quantiles",
                          "_by_log2_", libName, "_control_in_", region,
                          "_of_genes_in_t2t-col.20210610_",
                          paste0(chrName, collapse = "_"), ".tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(features_DF,
            file = paste0(outDir,
                          "features_", quantiles, "quantiles",
                          "_by_log2_", libName, "_control_in_", region,
                          "_of_genes_in_t2t-col.20210610_",
                          paste0(chrName, collapse = "_"), ".tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)

# Divide ranLocs into quantiles based on feature quantile indices
ranLocs_DF <- data.frame(ranLocs,
                         percentile = rank(ranLocs[,which(colnames(ranLocs) == "log2ChIP_ranLocMatRegionRowMeans")]) /
                                      length(ranLocs[,which(colnames(ranLocs) == "log2ChIP_ranLocMatRegionRowMeans")]),
                         random = as.character(""))
# Get row indices for each feature quantile
quantileIndices <- lapply(1:quantiles, function(k) {
  which(features_DF$quantile == paste0("Quantile ", k))
})
for(k in 1:quantiles) {
  ranLocs_DF[quantileIndices[[k]],]$random <- paste0("Random ", k)
}
write.table(ranLocs_DF,
            file = paste0(outDir,
                          "features_", quantiles, "quantiles",
                          "_by_log2_", libName, "_control_in_", region,
                          "_of_genes_in_t2t-col.20210610_",
                          paste0(chrName, collapse = "_"), "_ranLocs.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)

## Order features in each quantile by decreasing log2ChIP_featureMatRegion levels
## to define "row_order" for heatmaps
#combineRowOrders <- function(quantile_bool_list) {
#  do.call("c", lapply(quantile_bool_list, function(x) {
#    quantile_log2ChIP_featureMatRegionRowMeans <- rowMeans(log2ChIP_featureMatRegion[x,], na.rm = T)
#    quantile_log2ChIP_featureMatRegionRowMeans[which(is.na(quantile_log2ChIP_featureMatRegionRowMeans))] <- 0
#    which(x)[order(quantile_log2ChIP_featureMatRegionRowMeans, decreasing = T)]
#  }))
#}
#row_order <- combineRowOrders(quantile_bool_list =
#  lapply(seq_along(1:quantiles), function(k) { 
#    featuresDF$quantile == paste0("Quantile ", k)
#  })
#)
## Confirm row_order is as would be obtained by alternative method
## Note that this alternative 
#stopifnot(identical(row_order,
#                    order(featuresDF$log2ChIP_featureMatRegionRowMeans,
#                          decreasing=T)))
#
## Order feature IDs in each quantile by decreasing log2ChIP_featureMatRegion levels
## for use in GO term enrichment analysis
#listCombineRowOrders <- function(quantile_bool_list) {
#  do.call(list, lapply(quantile_bool_list, function(x) {
#    quantile_log2ChIP_featureMatRegionRowMeans <- rowMeans(log2ChIP_featureMatRegion[x,], na.rm = T)
#    quantile_log2ChIP_featureMatRegionRowMeans[which(is.na(quantile_log2ChIP_featureMatRegionRowMeans))] <- 0
#    which(x)[order(quantile_log2ChIP_featureMatRegionRowMeans, decreasing = T)]
#  }))
#}
#featureIndicesList <- listCombineRowOrders(quantile_bool_list =
#  lapply(seq_along(1:quantiles), function(k) {
#    featuresDF$quantile == paste0("Quantile ", k)
#  })
#)
#stopifnot(identical(row_order,
#                    do.call(c, lapply(featureIndicesList,
#                                      function(x) x))))
## Alternatively, with original ordering:
### Get feature indices for each quantile
##featureIndicesList <- lapply(seq_along(1:quantiles), function(k) {
##  which(featuresDF$quantile == paste0("Quantile ", k))
##})
#
#featureIDsQuantileList <- lapply(seq_along(1:quantiles), function(k) {
#  sub(pattern = "\\.\\d+", replacement = "",
#      x = as.vector(featuresDF[featureIndicesList[[k]],]$featureID))
#})
#sapply(seq_along(featureIDsQuantileList), function(k) {
#  write.table(featureIDsQuantileList[[k]],
#              file = paste0(outDir,
#                            "featureIDs_quantile", k, "_of_", quantiles,
#                            "_by_log2_", libName, "_control_in_",
#                            region, "_of_",
#                            substring(chrName[1][1], first = 1, last = 5), "_in_",
#                            paste0(substring(chrName, first = 10, last = 16),
#                                   collapse = "_"), "_",
#                            substring(chrName[1][1], first = 18), ".txt"),
#              quote = F, row.names = F, col.names = F)
#})

## Load feature matrices for each chromatin dataset, calculate log2(ChIP/control),
## and sort by decreasing log2mat1RegionRowMeans
#ChIPNames <- c(
#               "WT_MTOPVIB_HA_Rep1_ChIP",
#               "DMC1_Rep1_ChIP",
#               "H2AZ_Rep1_ChIP",
#               "H3K4me3_Rep1_ChIP",
#               "H3K4me1_Rep1_ChIP_SRR8126618",
#               "H3K27ac_Rep1_ChIP_SRR8126621",
#               "H3K27me3_ChIP_SRR6350666",
#               "H3K9me2_Rep1_ChIP",
#               "H3K27me1_Rep1_ChIP"
#              )
#ChIPNamesDir <- c(
#                  "20190819_dh580_Athaliana_ChIPseq_MTOPVIB/fastq_pooled/snakemake_ChIPseq_t2t-col.20210610",
#                  "DMC1",
#                  "H2AZ",
#                  "H3K4me3",
#                  "H3K4me1",
#                  "H3K27ac",
#                  "H3K27me3",
#                  "H3K9me2",
#                  "H3K27me1"
#                 )
#ChIPNamesPlot <- c(
#                   "ASY1",
#                   "DMC1",
#                   "H2A.Z",
#                   "H3K4me3",
#                   "H3K4me1",
#                   "H3K27ac",
#                   "H3K27me3",
#                   "H3K9me2",
#                   "H3K27me1"
#                  )
#ChIPColours <- c(
#                 "purple4",
#                 "green2",
#                 "dodgerblue",
#                 "forestgreen",
#                 "goldenrod1",
#                 "orange",
#                 "navy",
#                 "magenta3",
#                 "firebrick1"
#                )
#otherNames <- c(
#                "MNase_Rep1",
#                "DNaseI_Rep1_SRR8447247",
#                "WT_RNAseq_Rep1_ERR2402974",
#                "WT_RNAseq_Rep2_ERR2402973",
#                "WT_RNAseq_Rep3_ERR2402972"
#               )
#otherNamesDir <- c(
#                   "MNase",
#                   "DNaseI",
#                   "RNAseq_meiocyte_Martin_Moore_2018_FrontPlantSci",
#                   "RNAseq_meiocyte_Martin_Moore_2018_FrontPlantSci",
#                   "RNAseq_meiocyte_Martin_Moore_2018_FrontPlantSci"
#                  )
#otherNamesPlot <- c(
#                    "MNase",
#                    "DNaseI",
#                    "RNA-seq Rep1",
#                    "RNA-seq Rep2",
#                    "RNA-seq Rep3"
#                   )
#otherColours <- c(
#                  "darkcyan",
#                  "purple",
#                  "red4",
#                  "red4",
#                  "red4"
#                 )
#sRNANames <- c(
#               "CS+_2_LIB18613_LDI16228"
#              )
#sRNANamesDir <- c(
#                  "sRNAseq_meiocyte_Martin_Moore"
#                 )
#sRNANamesPlot <- c(
#                   "20-nt sRNAs",
#                   "21-nt sRNAs",
#                   "22-nt sRNAs",
#                   "23-nt sRNAs",
#                   "24-nt sRNAs",
#                   "34-nt sRNAs"
#                  )
#sRNAsizes <- c(
#               "20nt",
#               "21nt",
#               "22nt",
#               "23nt",
#               "24nt",
#               "33nt",
#               "34nt"
#              )
#sRNAColours <- c(
#                 "red",
#                 "blue",
#                 "green2",
#                 "darkorange2",
#                 "purple3",
#                 "darkgreen",
#                 "deeppink"
#                )
#DNAmethNames <- c(
#                  "BSseq_Rep8a_SRR6792678"
#                 )
#DNAmethNamesDir <- c(
#                     "BSseq"
#                    )
#DNAmethContexts <- c(
#                     "CpG",
#                     "CHG",
#                     "CHH"
#                    )
#DNAmethNamesPlot <- c(
#                      "mCG",
#                      "mCHG",
#                      "mCHH"
#                     )
#DNAmethColours <- c(
#                    "navy",
#                    "blue",
#                    "deepskyblue1"
#                   )
#
#ChIPDirs <- sapply(seq_along(ChIPNames), function(x) {
#  if(ChIPNames[x] %in% c("H3K4me3_ChIP_SRR6350668",
#                         "H3K27me3_ChIP_SRR6350666",
#                         "H3K36me3_ChIP_SRR6350670",
#                         "H3K9ac_ChIP_SRR6350667",
#                         "CENH3_ChIP_SRR1686799")) {
#    paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
#           ChIPNamesDir[x], "/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/")
#  } else if(ChIPNames[x] %in% c("H3K4me1_Rep1_ChIP_SRR8126618",
#                                "H3K27ac_Rep1_ChIP_SRR8126621")) {
#    paste0("/home/ajt200/analysis/wheat/epigenomics_seedlings_Li_2019_Genome_Biol/",
#           ChIPNamesDir[x], "/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/")
#  } else {
#    paste0("/home/ajt200/analysis/wheat/",
#           ChIPNamesDir[x], "/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/")
#  }
#})
#otherDirs <- sapply(seq_along(otherNames), function(x) {
#  if(otherNames[x] %in% c("MNase_Rep1")) {
#    paste0("/home/ajt200/analysis/wheat/",
#           otherNamesDir[x], "/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/")
#  } else if(otherNames[x] %in% c("DNaseI_Rep1_SRR8447247")) {
#    paste0("/home/ajt200/analysis/wheat/epigenomics_seedlings_Li_2019_Genome_Biol/",
#           otherNamesDir[x], "/snakemake_ChIPseq/mapped/geneProfiles_subgenomes/matrices/")
#  } else if(grepl("RNAseq", otherNames[x])) {
#    paste0("/home/ajt200/analysis/wheat/",
#           otherNamesDir[x], "/snakemake_RNAseq_HISAT2/mapped/geneProfiles_subgenomes/matrices/")
#  } else {
#    stop(paste0("otherNames[", x, "] is not compatible with the specified coverage matrix paths"))
#  }
#})
#sRNADirs <- sapply(seq_along(sRNANames), function(x) {
#  if(sRNANames[x] %in% c("CS+_2_LIB18613_LDI16228")) {
#    paste0("/home/ajt200/analysis/wheat/",
#           sRNANamesDir[x], "/snakemake_sRNAseq/mapped/geneProfiles_subgenomes/matrices/")
#  } else {
#    stop(paste0("sRNANames[", x, "] is not compatible with the specified coverage matrix paths"))
#  }
#})
#DNAmethDirs <- sapply(seq_along(DNAmethNames), function(x) {
#  if(DNAmethNames[x] %in% c("BSseq_Rep8a_SRR6792678")) {
#    paste0("/home/ajt200/analysis/wheat/epigenomics_shoot_leaf_IWGSC_2018_Science/",
#           DNAmethNamesDir[x],
#           "/snakemake_BSseq/coverage/geneProfiles_subgenomes/matrices/")
#  } else {
#    stop(paste0("DNAmethNames[", x, "] is not compatible with the specified coverage matrix paths"))
#  }
#})
#
## ChIP
#ChIP_featureMats <- mclapply(seq_along(ChIPNames), function(x) {
#  lapply(seq_along(chrName), function(y) {
#    as.matrix(read.table(paste0(ChIPDirs[x],
#                                ChIPNames[x],
#                                "_MappedOn_wheat_v1.0_lowXM_", align, "_sort_norm_",
#                                chrName[y], "_matrix_bin", binName,
#                                "_flank", flankName, ".tab"),
#                         header = F, skip = 3))
#  })
#}, mc.cores = length(ChIPNames))
## If features from all 3 subgenomes are to be analysed,
## concatenate the 3 corresponding feature coverage matrices
#ChIP_featureMats <- mclapply(seq_along(ChIP_featureMats), function(x) {
#  if(length(chrName) == 3) {
#    do.call(rbind, ChIP_featureMats[[x]])
#  } else {
#    ChIP_featureMats[[x]][[1]]
#  }
#}, mc.cores = length(ChIP_featureMats))
#
## Conditionally calculate log2(ChIP/input) or log2(ChIP/MNase)
## for each matrix depending on library
#log2ChIP_featureMats <- mclapply(seq_along(ChIP_featureMats), function(x) {
#  if(ChIPNames[x] %in% c(
#                         "WT_MTOPVIB_HA_Rep1_ChIP",
#                         "DMC1_Rep1_ChIP",
#                         "H3K4me3_ChIP_SRR6350668",
#                         "H3K27me3_ChIP_SRR6350666",
#                         "H3K36me3_ChIP_SRR6350670",
#                         "H3K9ac_ChIP_SRR6350667",
#                         "H3K4me1_Rep1_ChIP_SRR8126618",
#                         "H3K27ac_Rep1_ChIP_SRR8126621"
#                        )) {
#    print(paste0(ChIPNames[x], " was sonication-based; using ", controlNames[1], " for log2((ChIP+1)/(control+1)) calculation"))
#    log2((ChIP_featureMats[[x]]+1)/(control_ranLocMats[[1]]+1))
#  } else {
#    print(paste0(ChIPNames[x], " was MNase-based; using ", controlNames[2], " for log2((ChIP+1)/(control+1)) calculation"))
#    log2((ChIP_featureMats[[x]]+1)/(control_ranLocMats[[2]]+1))
#  }
#}, mc.cores = length(ChIP_featureMats))
#
#for(x in seq_along(log2ChIP_featureMats)) {
#  attr(log2ChIP_featureMats[[x]], "upstream_index") = 1:(upstream/binSize)
#  attr(log2ChIP_featureMats[[x]], "target_index") = ((upstream/binSize)+1):((upstream+bodyLength)/binSize)
#  attr(log2ChIP_featureMats[[x]], "downstream_index") = (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))
#  attr(log2ChIP_featureMats[[x]], "extend") = c(upstream, downstream)
#  attr(log2ChIP_featureMats[[x]], "smooth") = FALSE
#  attr(log2ChIP_featureMats[[x]], "signal_name") = ChIPNamesPlot[x]
#  attr(log2ChIP_featureMats[[x]], "target_name") = chrName
#  attr(log2ChIP_featureMats[[x]], "target_is_single_point") = FALSE
#  attr(log2ChIP_featureMats[[x]], "background") = 0
#  attr(log2ChIP_featureMats[[x]], "signal_is_categorical") = FALSE
#  class(log2ChIP_featureMats[[x]]) = c("normalizedMatrix", "matrix")
#}
#
#for(x in seq_along(control_ranLocMats)) {
#  attr(control_ranLocMats[[x]], "upstream_index") = 1:(upstream/binSize)
#  attr(control_ranLocMats[[x]], "target_index") = ((upstream/binSize)+1):((upstream+bodyLength)/binSize)
#  attr(control_ranLocMats[[x]], "downstream_index") = (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))
#  attr(control_ranLocMats[[x]], "extend") = c(upstream, downstream)
#  attr(control_ranLocMats[[x]], "smooth") = FALSE
#  attr(control_ranLocMats[[x]], "signal_name") = controlNamesPlot[x]
#  attr(control_ranLocMats[[x]], "target_name") = chrName
#  attr(control_ranLocMats[[x]], "target_is_single_point") = FALSE
#  attr(control_ranLocMats[[x]], "background") = 0
#  attr(control_ranLocMats[[x]], "signal_is_categorical") = FALSE
#  class(control_ranLocMats[[x]]) = c("normalizedMatrix", "matrix")
#}
#
## other
#othermats <- mclapply(seq_along(otherNames), function(x) {
#  lapply(seq_along(chrName), function(y) {
#    otherFile <- system(paste0("ls ", otherDirs[x],
#                               otherNames[x],
#                               "_MappedOn_wheat_v1.0*", align, "_sort_norm_",
#                               chrName[y], "_matrix_bin", binName,
#                               "_flank", flankName, ".tab"),
#                        intern = T)
#    as.matrix(read.table(otherFile,
#                         header = F, skip = 3))
#  })
#}, mc.cores = length(otherNames))
## If features from all 3 subgenomes are to be analysed,
## concatenate the 3 corresponding feature coverage matrices
#othermats <- mclapply(seq_along(othermats), function(x) {
#  if(length(chrName) == 3) {
#    do.call(rbind, othermats[[x]])
#  } else {
#    othermats[[x]][[1]]
#  }
#}, mc.cores = length(othermats))
#
#for(x in seq_along(othermats)) {
#  attr(othermats[[x]], "upstream_index") = 1:(upstream/binSize)
#  attr(othermats[[x]], "target_index") = ((upstream/binSize)+1):((upstream+bodyLength)/binSize)
#  attr(othermats[[x]], "downstream_index") = (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))
#  attr(othermats[[x]], "extend") = c(upstream, downstream)
#  attr(othermats[[x]], "smooth") = FALSE
#  attr(othermats[[x]], "signal_name") = otherNamesPlot[x]
#  attr(othermats[[x]], "target_name") = chrName
#  attr(othermats[[x]], "target_is_single_point") = FALSE
#  attr(othermats[[x]], "background") = 0
#  attr(othermats[[x]], "signal_is_categorical") = FALSE
#  class(othermats[[x]]) = c("normalizedMatrix", "matrix")
#}
#
## sRNA
#sRNAmats <- mclapply(seq_along(sRNAsizes), function(x) {
#  lapply(seq_along(chrName), function(y) {
#    as.matrix(read.table(paste0(sRNADirs,
#                                sRNANames,
#                                "_MappedOn_wheat_v1.0_", align, "_", sRNAsizes[x], "_sort_norm_",
#                                chrName[y], "_matrix_bin", binName,
#                                "_flank", flankName, ".tab"),
#                         header = F, skip = 3))
#  })
#}, mc.cores = length(sRNAsizes))
## If features from all 3 subgenomes are to be analysed,
## concatenate the 3 corresponding feature coverage matrices
#sRNAmats <- mclapply(seq_along(sRNAmats), function(x) {
#  if(length(chrName) == 3) {
#    do.call(rbind, sRNAmats[[x]])
#  } else {
#    sRNAmats[[x]][[1]]
#  }
#}, mc.cores = length(sRNAmats))
#
#for(x in seq_along(sRNAmats)) {
#  attr(sRNAmats[[x]], "upstream_index") = 1:(upstream/binSize)
#  attr(sRNAmats[[x]], "target_index") = ((upstream/binSize)+1):((upstream+bodyLength)/binSize)
#  attr(sRNAmats[[x]], "downstream_index") = (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))
#  attr(sRNAmats[[x]], "extend") = c(upstream, downstream)
#  attr(sRNAmats[[x]], "smooth") = FALSE
#  attr(sRNAmats[[x]], "signal_name") = sRNANamesPlot[x]
#  attr(sRNAmats[[x]], "target_name") = chrName
#  attr(sRNAmats[[x]], "target_is_single_point") = FALSE
#  attr(sRNAmats[[x]], "background") = 0
#  attr(sRNAmats[[x]], "signal_is_categorical") = FALSE
#  class(sRNAmats[[x]]) = c("normalizedMatrix", "matrix")
#}
#
## DNAmeth
#DNAmethmats <- mclapply(seq_along(DNAmethContexts), function(x) {
#  lapply(seq_along(chrName), function(y) {
#    as.matrix(read.table(paste0(DNAmethDirs,
#                                DNAmethNames,
#                                "_MappedOn_wheat_v1.0_incl_organelles_controls_dedup_", DNAmethContexts[x], "_",
#                                chrName[y], "_matrix_bin", binName,
#                                "_flank", flankName, ".tab"),
#                         header = F, skip = 3))
#  })
#}, mc.cores = length(DNAmethContexts))
## If features from all 3 subgenomes are to be analysed,
## concatenate the 3 corresponding feature coverage matrices
#DNAmethmats <- mclapply(seq_along(DNAmethmats), function(x) {
#  if(length(chrName) == 3) {
#    do.call(rbind, DNAmethmats[[x]])
#  } else {
#    DNAmethmats[[x]][[1]]
#  }
#}, mc.cores = length(DNAmethmats))
#
#for(x in seq_along(DNAmethmats)) {
#  attr(DNAmethmats[[x]], "upstream_index") = 1:(upstream/binSize)
#  attr(DNAmethmats[[x]], "target_index") = ((upstream/binSize)+1):((upstream+bodyLength)/binSize)
#  attr(DNAmethmats[[x]], "downstream_index") = (((upstream+bodyLength)/binSize)+1):(((upstream+bodyLength)/binSize)+(downstream/binSize))
#  attr(DNAmethmats[[x]], "extend") = c(upstream, downstream)
#  attr(DNAmethmats[[x]], "smooth") = FALSE
#  attr(DNAmethmats[[x]], "signal_name") = DNAmethNamesPlot[x]
#  attr(DNAmethmats[[x]], "target_name") = chrName
#  attr(DNAmethmats[[x]], "target_is_single_point") = FALSE
#  attr(DNAmethmats[[x]], "background") = "NA"
#  attr(DNAmethmats[[x]], "signal_is_categorical") = FALSE
#  class(DNAmethmats[[x]]) = c("normalizedMatrix", "matrix")
#}
#
#
#if(grepl("genes", chrName)) {
#  featureStartLab <- "TSS"
#  featureEndLab <- "TTS"
#} else {
#  featureStartLab <- "Start"
#  featureEndLab <- "End"
#}
#
## Heatmap plotting function
## Note that for plotting heatmaps for individual datasets in separate PDFs,
## must edit this function - print(EnrichedHeatmap(...))
#featureHeatmap <- function(mat,
#                           col_fun,
#                           colour,
#                           datName) {
#  EnrichedHeatmap(mat = mat,
#                  col = col_fun,
#                  column_title = datName,
#                  top_annotation = HeatmapAnnotation(enriched = anno_enriched(gp = gpar(col = colour,
#                                                                                        lwd = 2),
#                                                                              yaxis_side = "left",
#                                                                              yaxis_facing = "left",
#                                                                              yaxis_gp = gpar(fontsize = 10),
#                                                                              pos_line_gp = gpar(col = "black",
#                                                                                                 lty = 2,
#                                                                                                 lwd = 2))),
#                  top_annotation_height = unit(2, "cm"),
#                  width = unit(6, "cm"),
#                  name = datName,
#                  heatmap_legend_param = list(title = datName,
#                                              title_position = "topcenter",
#                                              title_gp = gpar(font = 2, fontsize = 12),
#                                              legend_direction = "horizontal",
#                                              labels_gp = gpar(fontsize = 10)),
#                  axis_name = c(paste0("-", flankNamePlot),
#                                featureStartLab, featureEndLab,
#                                paste0("+", flankNamePlot)),
#                  axis_name_gp = gpar(fontsize = 12),
#                  border = FALSE,
#                  pos_line_gp = gpar(col = "white", lty = 2, lwd = 2),
#                  # If converting into png with pdfTotiffTopng.sh,
#                  # set use_raster to FALSE
#                  #use_raster = FALSE)
#                  use_raster = TRUE, raster_device = "png", raster_quality = 10)
#}
#
## Define heatmap colours
#rich8to6equal <- c("#0000CB", "#0081FF", "#87CEFA", "#FDEE02", "#FFAB00", "#FF3300")
##quantileColours <- c("darkorange1", "green2", "purple3", "deepskyblue")
##quantileColours <- colorRampPalette(c("red", "blue"))(4)
#quantileColours <- c("red", "purple", "blue", "navy")
#
## Create quantile colour block "heatmap"
#quantileBlockhtmp <- Heatmap(featuresDF$quantile,
#                             col = structure(quantileColours,
#                                             names = paste0("Quantile ", 1:quantiles)),
#                             show_row_names = FALSE, show_heatmap_legend = FALSE,
#                             width = unit(3, "mm"), name = "") 
#quantilecMMbheatmap <- Heatmap(featuresDF$cMMb,
#                               cluster_rows = FALSE,
#                               col = colorRamp2(quantile(featuresDF$cMMb,
#                                                         c(0.60, 0.50, 0.40, 0.30),
#                                                         na.rm = T),
#                                                quantileColours), 
#                               na_col = "grey40",
#                               show_row_names = FALSE, show_heatmap_legend = TRUE,
#                               heatmap_legend_param = list(title = "cM/Mb",
#                                                           title_position = "topcenter",
#                                                           title_gp = gpar(font = 2, fontsize = 12),
#                                                           legend_direction = "horizontal",
#                                                           labels_gp = gpar(fontsize = 10)),
#                               width = unit(3, "cm"), name = "")
#quantilecMMbrowAnno <- rowAnnotation(cMMb = anno_points(featuresDF$cMMb,
#                                                        which = "row",
#                                                        size = unit(1, "mm"),
#                                                        gp = gpar(col = "black"),
#                                                        axis_param = list(at = c(0, 1.5), labels = c("0", "1.5")),
#                                                        width = unit(3, "cm")))
## Plot together
#log2ChIPhtmpList <- mclapply(seq_along(ChIPNames), function(x) {
#  ChIP_col_fun <- colorRamp2(quantile(log2ChIP_featureMats[[x]],
#                                      c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
#                                      na.rm = T),
#                             rich8to6equal)
#  featureHeatmap(mat = log2ChIP_featureMats[[x]],
#                 col_fun = ChIP_col_fun,
#                 colour = quantileColours,
#                 datName = ChIPNamesPlot[x])
#}, mc.cores = length(log2ChIP_featureMats))
#otherhtmpList <- mclapply(seq_along(otherNames), function(x) {
#  ChIP_col_fun <- colorRamp2(quantile(othermats[[x]],
#                                      c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
#                                      na.rm = T),
#                             rich8to6equal)
#  featureHeatmap(mat = othermats[[x]],
#                 col_fun = ChIP_col_fun,
#                 colour = quantileColours,
#                 datName = otherNamesPlot[x])
#}, mc.cores = length(othermats))
#controlhtmpList <- mclapply(seq_along(controlNames), function(x) {
#  ChIP_col_fun <- colorRamp2(quantile(control_ranLocMats[[x]],
#                                      c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95),
#                                      na.rm = T),
#                             rich8to6equal)
#  featureHeatmap(mat = control_ranLocMats[[x]],
#                 col_fun = ChIP_col_fun,
#                 colour = quantileColours,
#                 datName = controlNamesPlot[x])
#}, mc.cores = length(control_ranLocMats))
##sRNAhtmpList <- mclapply(seq_along(sRNANamesPlot), function(x) {
##  ChIP_col_fun <- colorRamp2(quantile(sRNAmats[[x]],
##                                      c(0.998, 0.9982, 0.9984, 0.9986, 0.9988, 0.999),
##                                      na.rm = T),
##                             rich8to6equal)
##  featureHeatmap(mat = sRNAmats[[x]],
##                 col_fun = ChIP_col_fun,
##                 colour = quantileColours,
##                 datName = sRNANamesPlot[x])
##}, mc.cores = length(sRNAmats))
#DNAmethhtmpList <- mclapply(seq_along(DNAmethNamesPlot), function(x) {
#  ChIP_col_fun <- colorRamp2(quantile(DNAmethmats[[x]],
#                                      c(0.50, 0.60, 0.70, 0.80, 0.90, 0.95),
#                                      na.rm = T),
#                             rich8to6equal)
#  featureHeatmap(mat = DNAmethmats[[x]],
#                 col_fun = ChIP_col_fun,
#                 colour = quantileColours,
#                 datName = DNAmethNamesPlot[x])
#}, mc.cores = length(DNAmethmats))
#
#htmpList <- c(quantileBlockhtmp,
#              quantilecMMbheatmap,
#              log2ChIPhtmpList,
#              controlhtmpList[[1]],
#              otherhtmpList,
##              sRNAhtmpList,
#              DNAmethhtmpList)
#htmps <- NULL
#for(x in 1:length(htmpList)) {
#  htmps <- htmps + htmpList[[x]]
#}
#pdf(paste0(plotDir, "log2ChIPcontrol_around_",
#           substring(chrName[1][1], first = 1, last = 5), "_in_",
#           paste0(substring(chrName, first = 10, last = 16),
#                  collapse = "_"), "_",
#           substring(chrName[1][1], first = 18),
#           "_heatmaps_quantiled_by_log2_", libName, "_control_in_", region, ".pdf"),
#    width = 3*length(htmpList),
#    height = 10)
#draw(htmps,
#     split = featuresDF$quantile,
#     row_order = row_order,
#     heatmap_legend_side = "bottom",
#     gap = unit(c(1, 1, rep(14, length(htmpList)-2)), "mm")
#    )
#dev.off()
