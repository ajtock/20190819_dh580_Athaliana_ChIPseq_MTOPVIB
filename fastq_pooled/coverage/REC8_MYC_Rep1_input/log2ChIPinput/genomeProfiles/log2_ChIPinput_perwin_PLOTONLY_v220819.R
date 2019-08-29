#!/applications/R/R-3.5.0/bin/Rscript

# Generate genome-wide plots of chromosome profiles

# MTOPVIB (magenta3)
# SPO11-1-oligos (dodgerblue2)

# Usage:
# Rscript ./log2_ChIPinput_perwin_PLOTONLY_v220819.R 10kb

args <- commandArgs(trailingOnly = T)
winName <- as.character(args[1])

library(parallel)

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

plotDir <- "./plots/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

MTOPVIBDir <- "/home/ajt200/analysis/20190819_dh580_Athaliana_ChIPseq_MTOPVIB/fastq_pooled/coverage/REC8_MYC_Rep1_input/log2ChIPinput/genomeProfiles/"
MTOPVIBprefix <- "MTOPVIB_HA_Rep1_ChIP_REC8_MYC_Rep1_input" 
MTOPVIBname <- "MTOPVIB"
MTOPVIBcolour <- "firebrick1"

inDir1 <- "/home/ajt200/analysis/20190819_dh580_Athaliana_ChIPseq_MTOPVIB/fastq_pooled/coverage/REC8_MYC_Rep1_input/log2ChIPinput/genomeProfiles/"
prefix1 <- "WT_SPO11oligo_RPI1_WT_nakedDNA_R1"
name1 <- "SPO11-1"
colour1 <- "dodgerblue2"

inDir2 <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/genomeProfiles/genomewideZscore/"
prefix2 <- "REC8_HA_Rep2_ChIP_REC8_MYC_Rep1_input"
name2 <- "REC8-HA"
colour2 <- "forestgreen"

# MTOPVIB
inDirMTOPVIB <- c(MTOPVIBDir)
prefixMTOPVIB <- c(MTOPVIBprefix)
nameMTOPVIB <- c(MTOPVIBname)
colourMTOPVIB <- c(MTOPVIBcolour)

filt_MTOPVIBlog2trans_list <- mclapply(seq_along(nameMTOPVIB), function(x) {
  read.table(paste0(inDirMTOPVIB[x], "filt_log2_", prefixMTOPVIB[x],
                    "_genome_norm_coverage_noZscore_", winName, ".txt"))
}, mc.cores = length(nameMTOPVIB))

# All
inDirAll <- c(inDir1)
prefixAll <- c(prefix1)
nameAll <- c(name1)
colourAll <- c(colour1)

filt_log2trans_list <- mclapply(seq_along(nameAll), function(x) {
  read.table(paste0(inDirAll[x], "filt_log2_", prefixAll[x],
                    "_genome_norm_coverage_noZscore_", winName, ".txt"))
}, mc.cores = length(nameAll))

filt_log2trans_list_Z <- read.table(paste0(inDir2, "filt_log2_", prefix2,
                                           "_genome_norm_coverage_Zscore_", winName, ".txt"))

# Feature frequency
inDirFeatures <- "/home/ajt200/analysis/REC8_pooled/coverage/common_input_MYC_Rep1/log2ChIPinput/genomeProfiles/genomewideZscore/"

# TEs
filt_TEs <- read.table(paste0(inDirFeatures, "filt_TE_frequency_genome_",
                              winName, ".txt"))
# genes
filt_genes <- read.table(paste0(inDirFeatures, "filt_gene_frequency_genome_",
                                winName, ".txt"))
# COs
filt_COs <- read.table(paste0(inDirFeatures, "filt_CO_frequency_genome_",
                              winName, ".txt"))
# DNA methylation in 200-kb windows
filt_meth <- read.table(paste0(inDirFeatures, "filt_DNAmeth_GSM980986_WT_rep2_genome_200kb.txt"))


# Function to plot MTOPVIB genome-scale coverage overlaid with other datasets
dat1to2diffY <- function(xplot,
                         dat1, dat2,
                         dat1Lab, dat2Lab,
                         dat1Colour, dat2Colour) {
  plot(xplot, dat1, type = "l", lwd = 2, col = dat1Colour,
       ylim = c(min(dat1),
                max(dat1)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  axis(side = 2, lwd.tick = 2, cex.axis = 2.00)
  mtext(side = 2, line = 3.50, cex = 2,
        text = dat1Lab, col = dat1Colour)
  par(new = T)
  plot(xplot, dat2, type = "l", lwd = 2, col = dat2Colour,
       ylim = c(min(dat2),
                max(dat2)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 3.00, adj = c(0.5, -2.0), labels = dat2Lab, xpd = NA, srt = -90,
       col = dat2Colour)
  axis(side = 4, at = pretty(dat2), lwd.tick = 2, cex.axis = 2.00)
  axis(side = 1, lwd.tick = 2,
       at = c(0, 2e7, 4e7, 6e7, 8e7, 1e8, 1.2e8),
       labels = c("0", "20", "40", "60", "80", "100", "120"), cex.axis = 2.00)
  mtext(side = 1, line = 3.50, cex = 2,
        text = "Coordinates (Mb)")
  abline(v = sumchr, lty = 1, lwd = 2)
  abline(v = centromeres, lty = 5, lwd = 2)
  rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 1, lwd = 1, col = "slategray1")
  rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 3, lwd = 1, col = "slategray1")
  box(lwd = 2)
}

# Function to plot MTOPVIB genome-scale coverage overlaid with other datasets 
dat1to2GenomePlot <- function(xplot,
                              dat1, dat2,
                              legendNames, plotColours) {
  plot(xplot, dat2, type = "l", lwd = 2, col = plotColours[2],
       ylim = c(min(dat1, dat2),
                max(dat1, dat2)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  lines(xplot, dat1, type = "l", lwd = 2, col = plotColours[1])
  axis(side = 1, lwd.tick = 2,
       at = c(0, 2e7, 4e7, 6e7, 8e7, 1e8, 1.2e8),
       labels = c("0", "20", "40", "60", "80", "100", "120"), cex.axis = 2.00)
  mtext(side = 1, line = 3.50, cex = 2,
        text = "Coordinates (Mb)")
  axis(side = 2, lwd.tick = 2, cex.axis = 2.00)
  mtext(side = 2, line = 3.00, cex = 2,
        text = expression("Log"[2]*"(ChIP/input)"))
  abline(v = sumchr, lty = 1, lwd = 2)
  abline(v = centromeres, lty = 5, lwd = 2)
  rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 1, lwd = 1, col = "slategray1")
  rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 3, lwd = 1, col = "slategray1")
  box(lwd = 2)
  legend("topleft",
         legend = legendNames,
         col = "white",
         text.col = plotColours,
         ncol = 1, cex = 1.2, lwd = 1.5, bty = "n")
}

# Function to plot MTOPVIB genome-scale coverage overlaid with other datasets
dat1to2vFeatureGenomePlot <- function(xplot,
                                      dat1, dat2, feature1,
                                      legendNames, legendPos, feature1Lab,
                                      plotColours, feature1Colour) {
  plot(xplot, feature1, type = "l", lwd = 2, col = feature1Colour,
       ylim = c(0, max(feature1)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 3.00, adj = c(0.5, -2.0), labels = feature1Lab, xpd = NA, srt = -90,
       col = feature1Colour)
  axis(side = 4, at = pretty(feature1), lwd.tick = 2, cex.axis = 2.00)
  par(new = T)
  plot(xplot, dat2, type = "l", lwd = 2, col = plotColours[2],
       ylim = c(min(dat1, dat2),
                max(dat1, dat2)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  lines(xplot, dat1, type = "l", lwd = 2, col = plotColours[1])
  axis(side = 1, lwd.tick = 2,
       at = c(0, 2e7, 4e7, 6e7, 8e7, 1e8, 1.2e8),
       labels = c("0", "20", "40", "60", "80", "100", "120"), cex.axis = 2.00)
  mtext(side = 1, line = 3.50, cex = 2,
        text = "Coordinates (Mb)")
  axis(side = 2, lwd.tick = 2, cex.axis = 2.00)
  mtext(side = 2, line = 3.00, cex = 2,
        text = expression("Log"[2]*"(ChIP/input)"))
  abline(v = sumchr, lty = 1, lwd = 2)
  abline(v = centromeres, lty = 5, lwd = 2)
  rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 1, lwd = 1, col = "slategray1")
  rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 3, lwd = 1, col = "slategray1")
  box(lwd = 2)
  legend(legendPos,
         legend = legendNames,
         col = "white",
         text.col = plotColours,
         ncol = 1, cex = 1.2, lwd = 1.5, bty = "n")
}

# Function to plot MTOPVIB genome-scale coverage overlaid with other datasets
dat1vFeatureGenomePlot <- function(xplot,
                                   dat1, feature1,
                                   legendNames, feature1Lab,
                                   plotColours, feature1Colour) {
  plot(xplot, feature1, type = "l", lwd = 2, col = feature1Colour,
       ylim = c(0, max(feature1)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 3.00, adj = c(0.5, -2.0), labels = feature1Lab, xpd = NA, srt = -90,
       col = feature1Colour)
  axis(side = 4, at = pretty(feature1), lwd.tick = 2, cex.axis = 2.00)
  par(new = T)
  plot(xplot, dat1, type = "l", lwd = 2, col = plotColours[1],
       ylim = c(min(dat1),
                max(dat1)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  axis(side = 1, lwd.tick = 2,
       at = c(0, 2e7, 4e7, 6e7, 8e7, 1e8, 1.2e8),
       labels = c("0", "20", "40", "60", "80", "100", "120"), cex.axis = 2.00)
  mtext(side = 1, line = 3.50, cex = 2,
        text = "Coordinates (Mb)")
  axis(side = 2, lwd.tick = 2, cex.axis = 2.00)
  mtext(side = 2, line = 3.00, cex = 2,
        text = expression("Log"[2]*"(ChIP/input)"))
  abline(v = sumchr, lty = 1, lwd = 2)
  abline(v = centromeres, lty = 5, lwd = 2)
  rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 1, lwd = 1, col = "slategray1")
  rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 3, lwd = 1, col = "slategray1")
  box(lwd = 2)
  legend("topleft",
         legend = legendNames,
         col = "white",
         text.col = plotColours,
         ncol = 1, cex = 1.2, lwd = 1.5, bty = "n")
}

# Function to plot MTOPVIB genome-scale coverage overlaid with other datasets 
dat1to4GenomePlot <- function(xplot,
                              dat1, dat2, dat3, dat4,
                              legendNames, plotColours) {
  plot(xplot, dat4, type = "l", lwd = 2, col = plotColours[4],
       ylim = c(min(dat1, dat2, dat3, dat4),
                max(dat1, dat2, dat3, dat4)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  lines(xplot, dat3, type = "l", lwd = 2, col = plotColours[3])
  lines(xplot, dat2, type = "l", lwd = 2, col = plotColours[2])
  lines(xplot, dat1, type = "l", lwd = 2, col = plotColours[1])
  axis(side = 1, lwd.tick = 2,
       at = c(0, 2e7, 4e7, 6e7, 8e7, 1e8, 1.2e8),
       labels = c("0", "20", "40", "60", "80", "100", "120"), cex.axis = 2.00)
  mtext(side = 1, line = 3.50, cex = 2,
        text = "Coordinates (Mb)")
  axis(side = 2, lwd.tick = 2, cex.axis = 2.00)
  mtext(side = 2, line = 3.00, cex = 2,
        text = expression("Log"[2]*"(ChIP/input)"))
  abline(v = sumchr, lty = 1, lwd = 2)
  abline(v = centromeres, lty = 5, lwd = 2)
  rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 1, lwd = 1, col = "slategray1")
  rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 3, lwd = 1, col = "slategray1")
  box(lwd = 2)
  legend("topleft",
         legend = legendNames,
         col = "white",
         text.col = plotColours,
         ncol = 1, cex = 1.2, lwd = 1.5, bty = "n")
}

# Function to plot MTOPVIB genome-scale coverage overlaid with DNA methylation (each context)
MTOPVIBvsDNAmethSepGenomePlot <- function(xplot1,
                                       xplot2,
                                       dat1A,
                                       dat2,
                                       dat1ALab,
                                       legendNames, plotColours) {
  plot(xplot1, dat1A, type = "l", lwd = 2, col = plotColours[1],
       ylim = c(min(dat1A), max(dat1A)),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  axis(side = 2, lwd.tick = 2, cex.axis = 2.00)
  mtext(side = 2, line = 3.50, cex = 2,
        #text = expression("Log"[2]*"(ChIP/input)"))
        text = dat1ALab, col = plotColours[1])
  par(new = T)
  plot(xplot2, dat2[,3], type = "l", lwd = 2, col = plotColours[2],
       ylim = c(min(dat2[,3], dat2[,4], dat2[,5]), max(dat2[,3], dat2[,4], dat2[,5])),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  p <- par('usr')
  text(p[2], mean(p[3:4]), cex = 3.00, adj = c(0.5, -2.0), labels = "DNA methylation", xpd = NA, srt = -90,
       col = "blue")
  lines(xplot2, dat2[,4], type = "l", lwd = 2, col = plotColours[3])
  lines(xplot2, dat2[,5], type = "l", lwd = 2, col = plotColours[4])
  axis(side = 4, at = pretty(c(dat2[,3], dat2[,4], dat2[,5])), lwd.tick = 2, cex.axis = 2.00)
  axis(side = 1, lwd.tick = 2,
       at = c(0, 2e7, 4e7, 6e7, 8e7, 1e8, 1.2e8),
       labels = c("0", "20", "40", "60", "80", "100", "120"), cex.axis = 2.00)
  mtext(side = 1, line = 3.50, cex = 2,
        text = "Coordinates (Mb)")
  abline(v = sumchr, lty = 1, lwd = 2)
  abline(v = centromeres, lty = 5, lwd = 2)
  rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 1, lwd = 1, col = "slategray1")
  rug(x = c(pericenStart, pericenEnd), ticksize = 0.03, side = 3, lwd = 1, col = "slategray1")
  box(lwd = 2)
  legend("topleft",
         legend = legendNames,
         col = "white",
         text.col = plotColours[2:4],
         ncol = 1, cex = 1.2, lwd = 1.5, bty = "n")
}


#
pdf(paste0(plotDir, "log2_", MTOPVIBprefix, "_", winName, "_genomeProfiles_v220819.pdf"),
    height = 28, width = 12)
par(mfcol = c(7, 1))
par(mar = c(5.1, 6.1, 2.1, 6.1))
par(mgp = c(3, 1.5, 0))
dat1to2diffY(xplot = filt_MTOPVIBlog2trans_list[[1]]$cumWindows,
             dat1 = filt_MTOPVIBlog2trans_list[[1]]$filt_log2cov,
             dat2 = filt_log2trans_list[[1]]$filt_log2cov,
             dat1Lab = nameMTOPVIB,
             dat2Lab = nameAll[1],
             dat1Colour = colourMTOPVIB,
             dat2Colour = colourAll[1])
dat1to2diffY(xplot = filt_MTOPVIBlog2trans_list[[1]]$cumWindows,
             dat1 = filt_MTOPVIBlog2trans_list[[1]]$filt_log2cov,
             dat2 = filt_log2trans_list_Z$filt_ZscoreLog2cov,
             dat1Lab = nameMTOPVIB,
             dat2Lab = name2,
             dat1Colour = colourMTOPVIB,
             dat2Colour = colour2)
dat1to2diffY(xplot = filt_MTOPVIBlog2trans_list[[1]]$cumWindows,
             dat1 = filt_MTOPVIBlog2trans_list[[1]]$filt_log2cov,
             dat2 = filt_TEs$filt_TEs,
             dat1Lab = nameMTOPVIB,
             dat2Lab = "TEs",
             dat1Colour = colourMTOPVIB,
             dat2Colour = "darkgreen")
MTOPVIBvsDNAmethSepGenomePlot(xplot1 = filt_MTOPVIBlog2trans_list[[1]]$cumWindows,
                           xplot2 = filt_meth$cumWindows,
                           dat1A = filt_MTOPVIBlog2trans_list[[1]]$filt_log2cov,
                           dat2 = filt_meth,
                           dat1ALab = nameMTOPVIB,
                           legendNames = c("mCG", "mCHG", "mCHH"),
                           plotColours = c(colourMTOPVIB[1], "navy", "blue", "deepskyblue1"))
dat1to2diffY(xplot = filt_MTOPVIBlog2trans_list[[1]]$cumWindows,
             dat1 = filt_MTOPVIBlog2trans_list[[1]]$filt_log2cov,
             dat2 = filt_genes$filt_genes,
             dat1Lab = nameMTOPVIB,
             dat2Lab = "Genes",
             dat1Colour = colourMTOPVIB,
             dat2Colour = "goldenrod4")
dat1to2diffY(xplot = filt_MTOPVIBlog2trans_list[[1]]$cumWindows,
             dat1 = filt_MTOPVIBlog2trans_list[[1]]$filt_log2cov,
             dat2 = filt_COs$filt_COs,
             dat1Lab = nameMTOPVIB,
             dat2Lab = "Crossovers",
             dat1Colour = colourMTOPVIB,
             dat2Colour = "darkblue")
dat1to2diffY(xplot = filt_MTOPVIBlog2trans_list[[1]]$cumWindows,
             dat1 = filt_genes$filt_genes,
             dat2 = filt_COs$filt_COs,
             dat1Lab = "Genes",
             dat2Lab = "Crossovers",
             dat1Colour = "goldenrod4",
             dat2Colour = "darkblue")
dev.off()
