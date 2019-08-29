#!/applications/R/R-3.5.0/bin/Rscript

###########################################################################
# Calculate log2 ratio of ChIP and input coverage values in 10-kb windows #
# for generating chromosome-scale plots                                   #
# and apply genome-wide Z-score standardisation                           #
###########################################################################

# Usage on hydrogen node7
# csmit -m 20G -c 1 "Rscript ./log2_ChIPinput_perwin_commandArgs_v220819.R '/home/ajt200/analysis/20190819_dh580_Athaliana_ChIPseq_MTOPVIB/fastq_pooled/coverage/' '/home/ajt200/analysis/REC8_pooled/coverage/' MTOPVIB_HA_Rep1_ChIP REC8_MYC_Rep1_input 10000 10kb 101"

library(segmentSeq)

args <- commandArgs(trailingOnly = TRUE)
ChIPDir <- args[1]
inputDir <- args[2]
ChIPname <- args[3]
inputname <- args[4]
winSize <- as.numeric(args[5])
winName <- as.character(args[6])
N <- as.numeric(args[7])

ChIPfile <- paste0(ChIPDir, ChIPname, "_norm_allchrs_coverage_coord_tab.bed")
inputfile <- paste0(inputDir, inputname, "_norm_allchrs_coverage_coord_tab.bed")
inputnameChIP <- paste0(inputname, "_", ChIPname)
outDir <- "./" 

ChIPcovDat <- read.table(ChIPfile,
                         colClasses = c(NA, rep("NULL", 2), NA))
inputcovDat <- read.table(inputfile, 
                          colClasses = c(NA, rep("NULL", 2), NA))

# Genomic definitions
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
# Pericentromeric regions are as defined in Supplemental Table S26
# of Ziolkowski et al. (2017) Genes Dev. 31
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)

# Make cumulative genomes
sumchr <- cumsum(c(0, chrLens))
print(sumchr)
sumchr_tot <- sumchr[length(sumchr)]
print(sumchr_tot)

centromeres <- sapply(seq_along(centromeres), function(x) {
  centromeres[x] + sumchr[x]
})
print(centromeres)
pericenStart <- sapply(seq_along(pericenStart), function(x) {
  pericenStart[x] + sumchr[x]
})
print(pericenStart)
pericenEnd <- sapply(seq_along(pericenEnd), function(x) {
  pericenEnd[x] + sumchr[x]
})
print(pericenEnd)

# Calculate mean coverage values within windows and log2(ChIP/input) ratios
print(winName)
cumChIPwinDat <- NULL
cuminputwinDat <- NULL
for(i in 1:5) {
  print(i)
  # Define coverage coordinates
  chrChIP <- ChIPcovDat[ChIPcovDat[,1] == chrs[i],]
  chrinput <- inputcovDat[inputcovDat[,1] == chrs[i],]
  ChIPcov <- chrChIP[,2]
  inputcov <- chrinput[,2]
  covCoords <- seq(1, length(ChIPcov), by = 1)
  covIRcoords <- IRanges(start = covCoords, width = 1)
  covGRcoords <- GRanges(seqnames = chrs[i], strand = "*",
                         ranges = covIRcoords)
  # Define adjacent windows
  seqWindows <- seq(1, chrLens[i], by = winSize)
  cumWindows <- seqWindows + sumchr[i]
  windowsIRanges <- IRanges(start = seqWindows,
                            width = winSize)
  windowsIRanges <- windowsIRanges[-length(windowsIRanges)]
  windowsIRanges <- append(windowsIRanges,
                           IRanges(start = seqWindows[length(seqWindows)],
                                   end = chrLens[i]))
  windowsGRanges <- GRanges(seqnames = chrs[i], strand = "*",
                            ranges = windowsIRanges)
  print(windowsGRanges)
  # Identify overlapping window and coverage coordinates
  overlaps <- getOverlaps(coordinates = windowsGRanges,
                          segments = covGRcoords,
                          overlapType = "overlapping",
                          whichOverlaps = TRUE)
  #+1 is offset to avoid infinite values
  ChIPcovWinVals <- sapply(overlaps, function(x) mean(ChIPcov[x]))
  inputcovWinVals <- sapply(overlaps, function(x) mean(inputcov[x]))
  ChIPwinDat <- data.frame(chr = rep(chrs[i], times = length(ChIPcovWinVals)),
                           cumWindows = as.integer(cumWindows), 
                           coverage = ChIPcovWinVals)
  inputwinDat <- data.frame(chr = rep(chrs[i], times = length(inputcovWinVals)),
                            cumWindows = as.integer(cumWindows), 
                            coverage = inputcovWinVals)
  cumChIPwinDat <- rbind(cumChIPwinDat, ChIPwinDat)
  cuminputwinDat <- rbind(cuminputwinDat, inputwinDat)
}
write.table(cumChIPwinDat,
            file = paste0(outDir, ChIPname, "_genome_norm_coverage_",
                          winName, ".txt"))
write.table(cuminputwinDat,
            file = paste0(outDir, inputnameChIP, "_genome_norm_coverage_",
                          winName, ".txt"))
cumlog2trans <- log2((cumChIPwinDat$coverage+1)/(cuminputwinDat$coverage+1))
cumlog2trans <- data.frame(chr = cumChIPwinDat$chr,
                           cumWindows = as.integer(cumChIPwinDat$cumWindows),
                           log2cov = as.numeric(cumlog2trans))
write.table(cumlog2trans,
            file = paste0(outDir, "log2_", ChIPname, "_", inputname,
                          "_genome_norm_coverage_noZscore_",
                          winName, ".txt"))

chrlog2trans <- lapply(seq_along(chrs), function(i) {
  cumlog2trans[cumlog2trans$chr == chrs[i],]
})


# Calculate moving average of current window,
#### (N/2) previous windows (where N is even) OR
# (N/2)-0.5 previous windows (where N is odd),
# and
#### (N/2) subsequent windows (where N is even) OR
# (N/2)-0.5 subsequent windows (where N is odd)
# (the higher N is, the greater the smoothing)
stopifnot(N %% 2 != 0)
flank <- (N/2)-0.5
# Define MA filter coefficients
f <- rep(1/N, N)

filt_log2trans <- NULL
for(i in seq_along(chrs)) {
  filt_chrlog2trans <- stats::filter(chrlog2trans[[i]]$log2cov,
                                     filter = f,
                                     sides = 2)
  filt_chrlog2trans[1:flank] <- filt_chrlog2trans[flank+1]
  filt_chrlog2trans[(length(filt_chrlog2trans)-flank+1):length(filt_chrlog2trans)] <- filt_chrlog2trans[(length(filt_chrlog2trans)-flank)]
  filt_chrlog2transDF <- data.frame(chr = chrlog2trans[[i]]$chr,
                                    cumWindows = as.integer(chrlog2trans[[i]]$cumWindows),
                                    filt_log2cov = as.numeric(filt_chrlog2trans))
  filt_log2trans <- rbind(filt_log2trans, filt_chrlog2transDF)
}
write.table(filt_log2trans,
            file = paste0(outDir, "filt_log2_", ChIPname, "_", inputname,
                          "_genome_norm_coverage_noZscore_",
                          winName, ".txt"))
