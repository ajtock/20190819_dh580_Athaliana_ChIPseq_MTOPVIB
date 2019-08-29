#!/applications/R/R-3.5.0/bin/Rscript

# Plot feature average coverage profiles with 95% CIs

# Usage:
# /applications/R/R-3.5.0/bin/Rscript features_avgProfileRibbon_noCorr_4profiles.R Genes 'Genes' 2000 2kb '2 kb' 20 20bp 'MTOPVIB_HA_Rep1,SPO11_1_oligos_RPI1,REC8_HA_Rep2,MNase' 'MTOPVIB,SPO11-1,REC8-HA,MNase' 'firebrick1,dodgerblue2,forestgreen,purple4' '0.69,0.30'

# wt,kss colours
# REC8: 'red,red4'
# SPO11-1-oligos: 'dodgerblue2,navy'
# H3K9me2: 'green2,darkgreen'
# mC*: 'orange,orange4'

#featureName <- "Genes"
#featureNamePlot <- "Genes"
#upstream <- 2000
#downstream <- 2000
#flankName <- "2kb"
#flankNamePlot <- "2 kb"
#binSize <- 20
#binName <- "20bp"
#libNames <- unlist(strsplit("REC8_HA_Rep2,SPO11_1_oligos_RPI1,MNase,H3K9me2",
#                            split = ","))
#libNamesPlot <- unlist(strsplit("REC8-HA,SPO11-1,Nucleosomes,H3K9me2",
#                                split = ","))
#colours <- unlist(strsplit("red,dodgerblue2,aquamarine,green2",
#                           split = ","))
### Custom annotation legends using textGrob()
## top left
#legendPos <- as.numeric(unlist(strsplit("0.02,0.94",
#                                        split = ",")))
## middle left upper
#legendPos <- as.numeric(unlist(strsplit("0.02,0.70",
#                                        split = ",")))
## middle left 
#legendPos <- as.numeric(unlist(strsplit("0.02,0.52",
#                                        split = ",")))
## middle left lower
#legendPos <- as.numeric(unlist(strsplit("0.02,0.40",
#                                        split = ",")))
## middle right
#legendPos <- as.numeric(unlist(strsplit("0.70,0.54",
#                                        split = ",")))
## bottom left
#legendPos <- as.numeric(unlist(strsplit("0.02,0.30",
#                                        split = ",")))

### ggplot2 legends
## top left
#legendPos <- as.numeric(unlist(strsplit("0.165,0.92",
#                                        split = ",")))
## OR bottom right
#legendPos <- as.numeric(unlist(strsplit("0.84,0.08",
#                                        split = ",")))

args <- commandArgs(trailingOnly = T)
featureName <- args[1]
featureNamePlot <- args[2]
upstream <- as.numeric(args[3])
downstream <- as.numeric(args[3])
flankName <- args[4]
flankNamePlot <- args[5]
binSize <- as.numeric(args[6])
binName <- args[7]
libNames <- unlist(strsplit(args[8],
                            split = ","))
libNamesPlot <- unlist(strsplit(args[9],
                                split = ","))
colours <- unlist(strsplit(args[10],
                           split = ","))
legendPos <- as.numeric(unlist(strsplit(args[11],
                                        split = ",")))

library(parallel)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)
library(extrafont)

matDir <- "./matrices/"
plotDir <- "./plots/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Load feature coverage matrix for each ChIP-seq dataset
featureMats <- mclapply(seq_along(libNames), function(x) {
  as.matrix(read.table(paste0(matDir,
                              libNames[x],                            
                              "_norm_cov_feature_smoothed_target_and_",
                              flankName, "_flank_dataframe.txt"),
                       header = T))
}, mc.cores = length(libNames))

# Load ranLoc coverage matrix for each ChIP-seq dataset
ranLocMats <- mclapply(seq_along(libNames), function(x) {
  as.matrix(read.table(paste0(matDir,
                              libNames[x],                            
                              "_norm_cov_ranLoc_smoothed_target_and_",
                              flankName, "_flank_dataframe.txt"),
                       header = T))
}, mc.cores = length(libNames))

## feature
# Transpose matrix and convert to dataframe
# in which first column is window name
wideDFfeature_list <- mclapply(seq_along(featureMats), function(x) {
  data.frame(window = colnames(featureMats[[x]]),
             t(featureMats[[x]]))
}, mc.cores = length(featureMats))

# Convert into tidy data.frame (long format)
tidyDFfeature_list  <- mclapply(seq_along(wideDFfeature_list), function(x) {
  gather(data  = wideDFfeature_list[[x]],
         key   = feature,
         value = coverage,
         -window)
}, mc.cores = length(wideDFfeature_list))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFfeature_list)) {
  tidyDFfeature_list[[x]]$window <- factor(tidyDFfeature_list[[x]]$window,
                                           levels = as.character(wideDFfeature_list[[x]]$window))
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (features) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFfeature_list  <- mclapply(seq_along(tidyDFfeature_list), function(x) {
  data.frame(window = as.character(wideDFfeature_list[[x]]$window),
             n      = tapply(X     = tidyDFfeature_list[[x]]$coverage,
                             INDEX = tidyDFfeature_list[[x]]$window,
                             FUN   = length),
             mean   = tapply(X     = tidyDFfeature_list[[x]]$coverage,
                             INDEX = tidyDFfeature_list[[x]]$window,
                             FUN   = mean),
             sd     = tapply(X     = tidyDFfeature_list[[x]]$coverage,
                             INDEX = tidyDFfeature_list[[x]]$window,
                             FUN   = sd))
}, mc.cores = length(tidyDFfeature_list))

for(x in seq_along(summaryDFfeature_list)) {
  summaryDFfeature_list[[x]]$window <- factor(summaryDFfeature_list[[x]]$window,
                                              levels = as.character(wideDFfeature_list[[x]]$window))
  summaryDFfeature_list[[x]]$winNo <- factor(1:dim(summaryDFfeature_list[[x]])[1])
  summaryDFfeature_list[[x]]$sem <- summaryDFfeature_list[[x]]$sd/sqrt(summaryDFfeature_list[[x]]$n-1)
  summaryDFfeature_list[[x]]$CI_lower <- summaryDFfeature_list[[x]]$mean -
    qt(0.975, df = summaryDFfeature_list[[x]]$n-1)*summaryDFfeature_list[[x]]$sem
  summaryDFfeature_list[[x]]$CI_upper <- summaryDFfeature_list[[x]]$mean +
    qt(0.975, df = summaryDFfeature_list[[x]]$n-1)*summaryDFfeature_list[[x]]$sem
}

names(summaryDFfeature_list) <- libNamesPlot

# Convert list summaryDFfeature_list into a single data.frame for plotting
summaryDFfeature <- bind_rows(summaryDFfeature_list, .id = "libName")
summaryDFfeature$libName <- factor(summaryDFfeature$libName,
                                   levels = names(summaryDFfeature_list))


## ranLoc
# Transpose matrix and convert to dataframe
# in which first column is window name
wideDFranLoc_list <- mclapply(seq_along(ranLocMats), function(x) {
  data.frame(window = colnames(ranLocMats[[x]]),
             t(ranLocMats[[x]]))
}, mc.cores = length(ranLocMats))

# Convert into tidy data.frame (long format)
tidyDFranLoc_list  <- mclapply(seq_along(wideDFranLoc_list), function(x) {
  gather(data  = wideDFranLoc_list[[x]],
         key   = ranLoc,
         value = coverage,
         -window)
}, mc.cores = length(wideDFranLoc_list))

# Order levels of factor "window" so that sequential levels
# correspond to sequential windows
for(x in seq_along(tidyDFranLoc_list)) {
  tidyDFranLoc_list[[x]]$window <- factor(tidyDFranLoc_list[[x]]$window,
                                          levels = as.character(wideDFranLoc_list[[x]]$window))
}

# Create summary data.frame in which each row corresponds to a window (Column 1),
# Column2 is the number of coverage values (ranLocs) per window,
# Column3 is the mean of coverage values per window,
# Column4 is the standard deviation of coverage values per window,
# Column5 is the standard error of the mean of coverage values per window,
# Column6 is the lower bound of the 95% confidence interval, and
# Column7 is the upper bound of the 95% confidence interval
summaryDFranLoc_list  <- mclapply(seq_along(tidyDFranLoc_list), function(x) {
  data.frame(window = as.character(wideDFranLoc_list[[x]]$window),
             n      = tapply(X     = tidyDFranLoc_list[[x]]$coverage,
                             INDEX = tidyDFranLoc_list[[x]]$window,
                             FUN   = length),
             mean   = tapply(X     = tidyDFranLoc_list[[x]]$coverage,
                             INDEX = tidyDFranLoc_list[[x]]$window,
                             FUN   = mean),
             sd     = tapply(X     = tidyDFranLoc_list[[x]]$coverage,
                             INDEX = tidyDFranLoc_list[[x]]$window,
                             FUN   = sd))
}, mc.cores = length(tidyDFranLoc_list))

for(x in seq_along(summaryDFranLoc_list)) {
  summaryDFranLoc_list[[x]]$window <- factor(summaryDFranLoc_list[[x]]$window,
                                             levels = as.character(wideDFranLoc_list[[x]]$window))
  summaryDFranLoc_list[[x]]$winNo <- factor(1:dim(summaryDFranLoc_list[[x]])[1])
  summaryDFranLoc_list[[x]]$sem <- summaryDFranLoc_list[[x]]$sd/sqrt(summaryDFranLoc_list[[x]]$n-1)
  summaryDFranLoc_list[[x]]$CI_lower <- summaryDFranLoc_list[[x]]$mean -
    qt(0.975, df = summaryDFranLoc_list[[x]]$n-1)*summaryDFranLoc_list[[x]]$sem
  summaryDFranLoc_list[[x]]$CI_upper <- summaryDFranLoc_list[[x]]$mean +
    qt(0.975, df = summaryDFranLoc_list[[x]]$n-1)*summaryDFranLoc_list[[x]]$sem
}

names(summaryDFranLoc_list) <- libNamesPlot

# Convert list summaryDFranLoc_list into a single data.frame for plotting
summaryDFranLoc <- bind_rows(summaryDFranLoc_list, .id = "libName")
summaryDFranLoc$libName <- factor(summaryDFranLoc$libName,
                                  levels = names(summaryDFranLoc_list))

if(featureName == "Genes") {
  featureStartLab <- "TSS"
  featureEndLab <- "TTS"
} else {
  featureStartLab <- "Start"
  featureEndLab <- "End"
}

# Define y-axis limits
#ymin <- min(c(summaryDFfeature$mean-summaryDFfeature$sem,
#              summaryDFranLoc$mean-summaryDFranLoc$sem))
#ymax <- max(c(summaryDFfeature$mean+summaryDFfeature$sem,
#              summaryDFranLoc$mean+summaryDFranLoc$sem))
ymin <- min(c(summaryDFfeature$CI_lower,
              summaryDFranLoc$CI_lower))
ymax <- max(c(summaryDFfeature$CI_upper,
              summaryDFranLoc$CI_upper))

# Function for formatting y-axis labels
# with a given number of decimals
fmt_decimals <- function(decimals) {
  function(x) format(x, nsmall = decimals, scientific = FALSE)
}

# Define legend labels
legendLabs <- lapply(seq_along(libNamesPlot), function(x) {
  grobTree(textGrob(bquote(.(libNamesPlot[x])),
                    x = legendPos[1], y = legendPos[2]-((x-1)*0.08), just = "left",
                    gp = gpar(col = colours[x], fontsize = 18))) 
})

# Plot average coverage profiles with 95% CI ribbon
## feature
ggObjGA <- NULL
ggObj1 <- NULL
ggObj1 <- ggplot(data = summaryDFfeature,
                 mapping = aes(x = winNo,
                               y = mean,
                               group = libName),
                ) +
  geom_line(data = summaryDFfeature,
            mapping = aes(colour = libName),
            size = 1) +
  scale_colour_manual(values = colours) +
  geom_ribbon(data = summaryDFfeature,
              #mapping = aes(ymin = mean-sem,
              #              ymax = mean+sem,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = libName),
              alpha = 0.4) +
  scale_fill_manual(values = colours) +
  scale_y_continuous(limits = c(ymin, ymax),
                     labels = function(x) sprintf("%5.2f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFfeature_list[[1]])[1])-(downstream/binSize),
                              dim(summaryDFfeature_list[[1]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              featureStartLab,
                              featureEndLab,
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFfeature_list[[1]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = expression("Log"[2]*"(ChIP/input)")) +
  annotation_custom(legendLabs[[1]]) +
  annotation_custom(legendLabs[[2]]) +
  annotation_custom(legendLabs[[3]]) +
  annotation_custom(legendLabs[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = "black"),
        legend.position = "none",
        #legend.text = element_text(size = 10),
        #legend.background = element_rect(fill = "transparent"),
        #legend.key = element_rect(colour = "transparent",
        #                          fill = "transparent"),
        #legend.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,0.9,0.0,0.3), "cm"), 
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote(.(featureNamePlot) ~ "(" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFfeature$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
## ranLoc
ggObj2 <- NULL
ggObj2 <- ggplot(data = summaryDFranLoc,
                 mapping = aes(x = winNo,
                               y = mean,
                               group = libName),
                ) +
  geom_line(data = summaryDFranLoc,
            mapping = aes(colour = libName),
            size = 1) +
  scale_colour_manual(values = colours) +
  geom_ribbon(data = summaryDFranLoc,
              #mapping = aes(ymin = mean-sem,
              #              ymax = mean+sem,
              mapping = aes(ymin = CI_lower,
                            ymax = CI_upper,
                            fill = libName),
              alpha = 0.4) +
  scale_fill_manual(values = colours) +
  scale_y_continuous(limits = c(ymin, ymax),
                     labels = function(x) sprintf("%5.2f", x)) +
  scale_x_discrete(breaks = c(1,
                              (upstream/binSize)+1,
                              (dim(summaryDFranLoc_list[[1]])[1])-(downstream/binSize),
                              dim(summaryDFranLoc_list[[1]])[1]),
                   labels = c(paste0("-", flankNamePlot),
                              "Start",
                              "End",
                              paste0("+", flankNamePlot))) +
  geom_vline(xintercept = c((upstream/binSize)+1,
                            (dim(summaryDFranLoc_list[[1]])[1])-(downstream/binSize)),
             linetype = "dashed",
             size = 1) +
  labs(x = "",
       y = "") +
  annotation_custom(legendLabs[[1]]) +
  annotation_custom(legendLabs[[2]]) +
  annotation_custom(legendLabs[[3]]) +
  annotation_custom(legendLabs[[4]]) +
  theme_bw() +
  theme(
        axis.ticks = element_line(size = 1.0, colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.text.x = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black", family = "Luxi Mono"),
        axis.title = element_text(size = 30, colour = "black"),
        legend.position = "none",
        #legend.text = element_text(size = 10),
        #legend.background = element_rect(fill = "transparent"),
        #legend.key = element_rect(colour = "transparent",
        #                          fill = "transparent"),
        #legend.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 3.5, colour = "black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,0.9,0.0,0.3), "cm"), 
        plot.title = element_text(hjust = 0.5, size = 30)) +
  ggtitle(bquote("Random loci (" * italic("n") ~ "=" ~
                 .(prettyNum(summaryDFranLoc$n[1],
                             big.mark = ",", trim = T)) *
                 ")"))
ggObjGA <- grid.arrange(ggObj1, ggObj2, nrow = 1, ncol = 2)
ggsave(paste0(plotDir,
              paste(libNames, collapse = "_"),
              "_around_", featureName, ".pdf"),
       plot = ggObjGA,
       height = 6.5, width = 14)

