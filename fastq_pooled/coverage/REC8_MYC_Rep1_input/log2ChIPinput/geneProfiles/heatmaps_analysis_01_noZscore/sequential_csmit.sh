#!/bin/bash

csmit -m 30G -c 1 "/applications/R/R-3.4.0/bin/Rscript ./gene_Profiles_commandArgs.R 2000 2kb 20 /home/ajt200/analysis/20190819_dh580_Athaliana_ChIPseq_MTOPVIB/fastq_pooled/coverage/REC8_MYC_Rep1_input/log2ChIPinput/noZscore/log2_MTOPVIB_HA_Rep1_ChIP_REC8_MYC_Rep1_input_norm_allchrs_coverage_coord_tab_noZscore.bed MTOPVIB_HA_Rep1" & sleep 10;
wait

