# Angela Oill
# XY Variant Calling Project
# Plotting results comparing default alignment to SCC alignment

#----------------#
# Load Libraries #
#----------------#

library("dplyr")
library("ggplot2")
library(cowplot)
library(ggpubr)
library(viridis)

library(tidyverse)
library(ggpubr)
library(rstatix)

library("RColorBrewer")



# Read in number of simulated variants for the same regions

##                          ##
## Read in data - simulated ##
##                          ##

# FEMALES #
chrX_PARs_sim_EUR <- read.table("../../simulated_variants/20_samples/EUR_females_chrX_PARs_diploid_simulated_counts.txt",
                       header = T)

chrX_nonPARs_sim_EUR <- read.table("../../simulated_variants/20_samples/EUR_females_chrX_nonPARs_diploid_simulated_counts.txt",
                       header = T)

chrX_PARs_sim_ASN <- read.table("../../simulated_variants/20_samples/ASN_females_chrX_PARs_diploid_simulated_counts.txt",
                               header = T)

chrX_nonPARs_sim_ASN <- read.table("../../simulated_variants/20_samples/ASN_females_chrX_nonPARs_diploid_simulated_counts.txt",
                                  header = T)

chrX_PARs_sim_AFR <- read.table("../../simulated_variants/20_samples/AFR_females_chrX_PARs_diploid_simulated_counts.txt",
                                header = T)

chrX_nonPARs_sim_AFR <- read.table("../../simulated_variants/20_samples/AFR_females_chrX_nonPARs_diploid_simulated_counts.txt",
                                   header = T)

# MALES #
chrX_PARs_sim_males_EUR <- read.table("../../simulated_variants/20_samples/EUR_males_chrX_PARs_diploid_simulated_counts.txt",
                                  header = T)

chrX_nonPARs_sim_males_EUR <- read.table("../../simulated_variants/20_samples/EUR_males_chrX_nonPARs_haploid_simulated_counts.txt",
                                                         header = T)

chrY_sim_males_EUR <- read.table("../../simulated_variants/20_samples/EUR_males_chrY_haploid_simulated_counts.txt",
                                                         header = T)


chrX_PARs_sim_males_ASN <- read.table("../../simulated_variants/20_samples/ASN_males_chrX_PARs_diploid_simulated_counts.txt",
                                     header = T)

chrX_nonPARs_sim_males_ASN <- read.table("../../simulated_variants/20_samples/ASN_males_chrX_nonPARs_haploid_simulated_counts.txt",
                                        header = T)

chrY_sim_males_ASN <- read.table("../../simulated_variants/20_samples/ASN_males_chrY_haploid_simulated_counts.txt",
                                header = T)

chrX_PARs_sim_males_AFR <- read.table("../../simulated_variants/20_samples/AFR_males_chrX_PARs_diploid_simulated_counts.txt",
                                      header = T)

chrX_nonPARs_sim_males_AFR <- read.table("../../simulated_variants/20_samples/AFR_males_chrX_nonPARs_haploid_simulated_counts.txt",
                                         header = T)

chrY_sim_males_AFR <- read.table("../../simulated_variants/20_samples/AFR_males_chrY_haploid_simulated_counts.txt",
                                 header = T)

##              ##
## Read in data ##
##              ##
# FEMALES #
## Chromosome X PARs - EUR ##
chrX_PARs_dat_EUR <- read.table(
  "../../performance_metrics/EUR/females_20_samples/data/EUR_females_chrX_PARs_diploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrX_PARs_dat_EUR <- chrX_PARs_dat_EUR[,!(names(chrX_PARs_dat_EUR) %in% drop)]


# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_PARs_dat_EUR_merge <- merge(chrX_PARs_dat_EUR, chrX_PARs_sim_EUR, by = "Sample")

chrX_PARs_dat_EUR_merge$Chromosome <- "chrX_PARs"
chrX_PARs_dat_EUR_merge$SampleSet <- "EUR"
chrX_PARs_dat_EUR_merge$FPtoTP <- chrX_PARs_dat_EUR_merge$FP / chrX_PARs_dat_EUR_merge$TP
chrX_PARs_dat_EUR_merge$FNtoTP <- chrX_PARs_dat_EUR_merge$FN / chrX_PARs_dat_EUR_merge$TP
chrX_PARs_dat_EUR_merge$FPtoSimVars <- chrX_PARs_dat_EUR_merge$FP / chrX_PARs_dat_EUR_merge$Num_sim_variants
chrX_PARs_dat_EUR_merge$FNtoSimVars <- chrX_PARs_dat_EUR_merge$FN / chrX_PARs_dat_EUR_merge$Num_sim_variants
chrX_PARs_dat_EUR_merge$TPtoSimVars <- chrX_PARs_dat_EUR_merge$TP / chrX_PARs_dat_EUR_merge$Num_sim_variants

## Chromosome X PARs - ASN ##
chrX_PARs_dat_ASN <- read.table(
  "../../performance_metrics/ASN/females/data/ASN_females_chrX_PARs_diploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrX_PARs_dat_ASN <- chrX_PARs_dat_ASN[,!(names(chrX_PARs_dat_ASN) %in% drop)]


# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_PARs_dat_ASN_merge <- merge(chrX_PARs_dat_ASN, chrX_PARs_sim_ASN, by = "Sample")

chrX_PARs_dat_ASN_merge$Chromosome <- "chrX_PARs"
chrX_PARs_dat_ASN_merge$SampleSet <- "ASN"
chrX_PARs_dat_ASN_merge$FPtoTP <- chrX_PARs_dat_ASN_merge$FP / chrX_PARs_dat_ASN_merge$TP
chrX_PARs_dat_ASN_merge$FNtoTP <- chrX_PARs_dat_ASN_merge$FN / chrX_PARs_dat_ASN_merge$TP
chrX_PARs_dat_ASN_merge$FPtoSimVars <- chrX_PARs_dat_ASN_merge$FP / chrX_PARs_dat_ASN_merge$Num_sim_variants
chrX_PARs_dat_ASN_merge$FNtoSimVars <- chrX_PARs_dat_ASN_merge$FN / chrX_PARs_dat_ASN_merge$Num_sim_variants
chrX_PARs_dat_ASN_merge$TPtoSimVars <- chrX_PARs_dat_ASN_merge$TP / chrX_PARs_dat_ASN_merge$Num_sim_variants

## Chromosome X PARs - AFR ##
chrX_PARs_dat_AFR <- read.table(
  "../../performance_metrics/AFR/females/data/AFR_females_chrX_PARs_diploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrX_PARs_dat_AFR <- chrX_PARs_dat_AFR[,!(names(chrX_PARs_dat_AFR) %in% drop)]


# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_PARs_dat_AFR_merge <- merge(chrX_PARs_dat_AFR, chrX_PARs_sim_AFR, by = "Sample")

chrX_PARs_dat_AFR_merge$Chromosome <- "chrX_PARs"
chrX_PARs_dat_AFR_merge$SampleSet <- "AFR"
chrX_PARs_dat_AFR_merge$FPtoTP <- chrX_PARs_dat_AFR_merge$FP / chrX_PARs_dat_AFR_merge$TP
chrX_PARs_dat_AFR_merge$FNtoTP <- chrX_PARs_dat_AFR_merge$FN / chrX_PARs_dat_AFR_merge$TP
chrX_PARs_dat_AFR_merge$FPtoSimVars <- chrX_PARs_dat_AFR_merge$FP / chrX_PARs_dat_AFR_merge$Num_sim_variants
chrX_PARs_dat_AFR_merge$FNtoSimVars <- chrX_PARs_dat_AFR_merge$FN / chrX_PARs_dat_AFR_merge$Num_sim_variants
chrX_PARs_dat_AFR_merge$TPtoSimVars <- chrX_PARs_dat_AFR_merge$TP / chrX_PARs_dat_AFR_merge$Num_sim_variants


## Chromosome X non-PARs - EUR ##
chrX_nonPARs_dat_EUR <- read.table(
  "../../performance_metrics/EUR/females_20_samples/data/EUR_females_chrX_nonPARs_diploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrX_nonPARs_dat_EUR <- chrX_nonPARs_dat_EUR[,!(names(chrX_nonPARs_dat_EUR) %in% drop)]


# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_nonPARs_dat_EUR_merge <- merge(chrX_nonPARs_dat_EUR, chrX_nonPARs_sim_EUR, by = "Sample")

chrX_nonPARs_dat_EUR_merge$Chromosome <- "chrX_nonPARs"
chrX_nonPARs_dat_EUR_merge$SampleSet <- "EUR"
chrX_nonPARs_dat_EUR_merge$FPtoTP <- chrX_nonPARs_dat_EUR_merge$FP / chrX_nonPARs_dat_EUR_merge$TP
chrX_nonPARs_dat_EUR_merge$FNtoTP <- chrX_nonPARs_dat_EUR_merge$FN / chrX_nonPARs_dat_EUR_merge$TP
chrX_nonPARs_dat_EUR_merge$FPtoSimVars <- chrX_nonPARs_dat_EUR_merge$FP / chrX_nonPARs_dat_EUR_merge$Num_sim_variants
chrX_nonPARs_dat_EUR_merge$FNtoSimVars <- chrX_nonPARs_dat_EUR_merge$FN / chrX_nonPARs_dat_EUR_merge$Num_sim_variants
chrX_nonPARs_dat_EUR_merge$TPtoSimVars <- chrX_nonPARs_dat_EUR_merge$TP / chrX_nonPARs_dat_EUR_merge$Num_sim_variants

## Chromosome X non-PARs - ASN ##
chrX_nonPARs_dat_ASN <- read.table(
  "../../performance_metrics/ASN/females/data/ASN_females_chrX_nonPARs_diploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrX_nonPARs_dat_ASN <- chrX_nonPARs_dat_ASN[,!(names(chrX_nonPARs_dat_ASN) %in% drop)]


# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_nonPARs_dat_ASN_merge <- merge(chrX_nonPARs_dat_ASN, chrX_nonPARs_sim_ASN, by = "Sample")

chrX_nonPARs_dat_ASN_merge$Chromosome <- "chrX_nonPARs"
chrX_nonPARs_dat_ASN_merge$SampleSet <- "ASN"
chrX_nonPARs_dat_ASN_merge$FPtoTP <- chrX_nonPARs_dat_ASN_merge$FP / chrX_nonPARs_dat_ASN_merge$TP
chrX_nonPARs_dat_ASN_merge$FNtoTP <- chrX_nonPARs_dat_ASN_merge$FN / chrX_nonPARs_dat_ASN_merge$TP
chrX_nonPARs_dat_ASN_merge$FPtoSimVars <- chrX_nonPARs_dat_ASN_merge$FP / chrX_nonPARs_dat_ASN_merge$Num_sim_variants
chrX_nonPARs_dat_ASN_merge$FNtoSimVars <- chrX_nonPARs_dat_ASN_merge$FN / chrX_nonPARs_dat_ASN_merge$Num_sim_variants
chrX_nonPARs_dat_ASN_merge$TPtoSimVars <- chrX_nonPARs_dat_ASN_merge$TP / chrX_nonPARs_dat_ASN_merge$Num_sim_variants

## Chromosome X non-PARs - AFR ##
chrX_nonPARs_dat_AFR <- read.table(
  "../../performance_metrics/AFR/females/data/AFR_females_chrX_nonPARs_diploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrX_nonPARs_dat_AFR <- chrX_nonPARs_dat_AFR[,!(names(chrX_nonPARs_dat_AFR) %in% drop)]


# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_nonPARs_dat_AFR_merge <- merge(chrX_nonPARs_dat_AFR, chrX_nonPARs_sim_AFR, by = "Sample")

chrX_nonPARs_dat_AFR_merge$Chromosome <- "chrX_nonPARs"
chrX_nonPARs_dat_AFR_merge$SampleSet <- "AFR"
chrX_nonPARs_dat_AFR_merge$FPtoTP <- chrX_nonPARs_dat_AFR_merge$FP / chrX_nonPARs_dat_AFR_merge$TP
chrX_nonPARs_dat_AFR_merge$FNtoTP <- chrX_nonPARs_dat_AFR_merge$FN / chrX_nonPARs_dat_AFR_merge$TP
chrX_nonPARs_dat_AFR_merge$FPtoSimVars <- chrX_nonPARs_dat_AFR_merge$FP / chrX_nonPARs_dat_AFR_merge$Num_sim_variants
chrX_nonPARs_dat_AFR_merge$FNtoSimVars <- chrX_nonPARs_dat_AFR_merge$FN / chrX_nonPARs_dat_AFR_merge$Num_sim_variants
chrX_nonPARs_dat_AFR_merge$TPtoSimVars <- chrX_nonPARs_dat_AFR_merge$TP / chrX_nonPARs_dat_AFR_merge$Num_sim_variants


# MALES #
## Chromosome X PARs - EUR ##
chrX_PARs_dat_males_EUR <- read.table(
  "../../performance_metrics/EUR/males_20_samples/data/EUR_males_chrX_PARs_diploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrX_PARs_dat_males_EUR <- chrX_PARs_dat_males_EUR[,!(names(chrX_PARs_dat_males_EUR) %in% drop)]


# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_PARs_dat_males_EUR_merge <- merge(chrX_PARs_dat_males_EUR, chrX_PARs_sim_males_EUR, by = "Sample")

chrX_PARs_dat_males_EUR_merge$Chromosome <- "chrX_PARs"
chrX_PARs_dat_males_EUR_merge$SampleSet <- "EUR"
chrX_PARs_dat_males_EUR_merge$FPtoTP <- chrX_PARs_dat_males_EUR_merge$FP / chrX_PARs_dat_males_EUR_merge$TP
chrX_PARs_dat_males_EUR_merge$FNtoTP <- chrX_PARs_dat_males_EUR_merge$FN / chrX_PARs_dat_males_EUR_merge$TP
chrX_PARs_dat_males_EUR_merge$FPtoSimVars <- chrX_PARs_dat_males_EUR_merge$FP / chrX_PARs_dat_males_EUR_merge$Num_sim_variants
chrX_PARs_dat_males_EUR_merge$FNtoSimVars <- chrX_PARs_dat_males_EUR_merge$FN / chrX_PARs_dat_males_EUR_merge$Num_sim_variants
chrX_PARs_dat_males_EUR_merge$TPtoSimVars <- chrX_PARs_dat_males_EUR_merge$TP / chrX_PARs_dat_males_EUR_merge$Num_sim_variants

## Chromosome X PARs - ASN ##
chrX_PARs_dat_males_ASN <- read.table(
  "../../performance_metrics/ASN/males/data/ASN_males_chrX_PARs_diploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrX_PARs_dat_males_ASN <- chrX_PARs_dat_males_ASN[,!(names(chrX_PARs_dat_males_ASN) %in% drop)]


# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_PARs_dat_males_ASN_merge <- merge(chrX_PARs_dat_males_ASN, chrX_PARs_sim_males_ASN, by = "Sample")

chrX_PARs_dat_males_ASN_merge$Chromosome <- "chrX_PARs"
chrX_PARs_dat_males_ASN_merge$SampleSet <- "ASN"
chrX_PARs_dat_males_ASN_merge$FPtoTP <- chrX_PARs_dat_males_ASN_merge$FP / chrX_PARs_dat_males_ASN_merge$TP
chrX_PARs_dat_males_ASN_merge$FNtoTP <- chrX_PARs_dat_males_ASN_merge$FN / chrX_PARs_dat_males_ASN_merge$TP
chrX_PARs_dat_males_ASN_merge$FPtoSimVars <- chrX_PARs_dat_males_ASN_merge$FP / chrX_PARs_dat_males_ASN_merge$Num_sim_variants
chrX_PARs_dat_males_ASN_merge$FNtoSimVars <- chrX_PARs_dat_males_ASN_merge$FN / chrX_PARs_dat_males_ASN_merge$Num_sim_variants
chrX_PARs_dat_males_ASN_merge$TPtoSimVars <- chrX_PARs_dat_males_ASN_merge$TP / chrX_PARs_dat_males_ASN_merge$Num_sim_variants

## Chromosome X PARs - AFR ##
chrX_PARs_dat_males_AFR <- read.table(
  "../../performance_metrics/AFR/males/data/AFR_males_chrX_PARs_diploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrX_PARs_dat_males_AFR <- chrX_PARs_dat_males_AFR[,!(names(chrX_PARs_dat_males_AFR) %in% drop)]


# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_PARs_dat_males_AFR_merge <- merge(chrX_PARs_dat_males_AFR, chrX_PARs_sim_males_AFR, by = "Sample")

chrX_PARs_dat_males_AFR_merge$Chromosome <- "chrX_PARs"
chrX_PARs_dat_males_AFR_merge$SampleSet <- "AFR"
chrX_PARs_dat_males_AFR_merge$FPtoTP <- chrX_PARs_dat_males_AFR_merge$FP / chrX_PARs_dat_males_AFR_merge$TP
chrX_PARs_dat_males_AFR_merge$FNtoTP <- chrX_PARs_dat_males_AFR_merge$FN / chrX_PARs_dat_males_AFR_merge$TP
chrX_PARs_dat_males_AFR_merge$FPtoSimVars <- chrX_PARs_dat_males_AFR_merge$FP / chrX_PARs_dat_males_AFR_merge$Num_sim_variants
chrX_PARs_dat_males_AFR_merge$FNtoSimVars <- chrX_PARs_dat_males_AFR_merge$FN / chrX_PARs_dat_males_AFR_merge$Num_sim_variants
chrX_PARs_dat_males_AFR_merge$TPtoSimVars <- chrX_PARs_dat_males_AFR_merge$TP / chrX_PARs_dat_males_AFR_merge$Num_sim_variants


## Chromosome X non-PARs - EUR ##
chrX_nonPARs_dat_males_EUR <- read.table(
  "../../performance_metrics/EUR/males_20_samples/data/EUR_males_chrX_nonPARs_haploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrX_nonPARs_dat_males_EUR <- chrX_nonPARs_dat_males_EUR[,!(names(chrX_nonPARs_dat_males_EUR) %in% drop)]


# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_nonPARs_dat_males_EUR_merge <- merge(chrX_nonPARs_dat_males_EUR, chrX_nonPARs_sim_males_EUR, by = "Sample")

chrX_nonPARs_dat_males_EUR_merge$Chromosome <- "chrX_nonPARs"
chrX_nonPARs_dat_males_EUR_merge$SampleSet <- "EUR"
chrX_nonPARs_dat_males_EUR_merge$FPtoTP <- chrX_nonPARs_dat_males_EUR_merge$FP / chrX_nonPARs_dat_males_EUR_merge$TP
chrX_nonPARs_dat_males_EUR_merge$FNtoTP <- chrX_nonPARs_dat_males_EUR_merge$FN / chrX_nonPARs_dat_males_EUR_merge$TP
chrX_nonPARs_dat_males_EUR_merge$FPtoSimVars <- chrX_nonPARs_dat_males_EUR_merge$FP / chrX_nonPARs_dat_males_EUR_merge$Num_sim_variants
chrX_nonPARs_dat_males_EUR_merge$FNtoSimVars <- chrX_nonPARs_dat_males_EUR_merge$FN / chrX_nonPARs_dat_males_EUR_merge$Num_sim_variants
chrX_nonPARs_dat_males_EUR_merge$TPtoSimVars <- chrX_nonPARs_dat_males_EUR_merge$TP / chrX_nonPARs_dat_males_EUR_merge$Num_sim_variants

## Chromosome X non-PARs - ASN ##
chrX_nonPARs_dat_males_ASN <- read.table(
  "../../performance_metrics/ASN/males/data/ASN_males_chrX_nonPARs_haploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrX_nonPARs_dat_males_ASN <- chrX_nonPARs_dat_males_ASN[,!(names(chrX_nonPARs_dat_males_ASN) %in% drop)]


# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_nonPARs_dat_males_ASN_merge <- merge(chrX_nonPARs_dat_males_ASN, chrX_nonPARs_sim_males_ASN, by = "Sample")

chrX_nonPARs_dat_males_ASN_merge$Chromosome <- "chrX_nonPARs"
chrX_nonPARs_dat_males_ASN_merge$SampleSet <- "ASN"
chrX_nonPARs_dat_males_ASN_merge$FPtoTP <- chrX_nonPARs_dat_males_ASN_merge$FP / chrX_nonPARs_dat_males_ASN_merge$TP
chrX_nonPARs_dat_males_ASN_merge$FNtoTP <- chrX_nonPARs_dat_males_ASN_merge$FN / chrX_nonPARs_dat_males_ASN_merge$TP
chrX_nonPARs_dat_males_ASN_merge$FPtoSimVars <- chrX_nonPARs_dat_males_ASN_merge$FP / chrX_nonPARs_dat_males_ASN_merge$Num_sim_variants
chrX_nonPARs_dat_males_ASN_merge$FNtoSimVars <- chrX_nonPARs_dat_males_ASN_merge$FN / chrX_nonPARs_dat_males_ASN_merge$Num_sim_variants
chrX_nonPARs_dat_males_ASN_merge$TPtoSimVars <- chrX_nonPARs_dat_males_ASN_merge$TP / chrX_nonPARs_dat_males_ASN_merge$Num_sim_variants

## Chromosome X non-PARs - AFR ##
chrX_nonPARs_dat_males_AFR <- read.table(
  "../../performance_metrics/AFR/males/data/AFR_males_chrX_nonPARs_haploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrX_nonPARs_dat_males_AFR <- chrX_nonPARs_dat_males_AFR[,!(names(chrX_nonPARs_dat_males_AFR) %in% drop)]


# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_nonPARs_dat_males_AFR_merge <- merge(chrX_nonPARs_dat_males_AFR, chrX_nonPARs_sim_males_AFR, by = "Sample")

chrX_nonPARs_dat_males_AFR_merge$Chromosome <- "chrX_nonPARs"
chrX_nonPARs_dat_males_AFR_merge$SampleSet <- "AFR"
chrX_nonPARs_dat_males_AFR_merge$FPtoTP <- chrX_nonPARs_dat_males_AFR_merge$FP / chrX_nonPARs_dat_males_AFR_merge$TP
chrX_nonPARs_dat_males_AFR_merge$FNtoTP <- chrX_nonPARs_dat_males_AFR_merge$FN / chrX_nonPARs_dat_males_AFR_merge$TP
chrX_nonPARs_dat_males_AFR_merge$FPtoSimVars <- chrX_nonPARs_dat_males_AFR_merge$FP / chrX_nonPARs_dat_males_AFR_merge$Num_sim_variants
chrX_nonPARs_dat_males_AFR_merge$FNtoSimVars <- chrX_nonPARs_dat_males_AFR_merge$FN / chrX_nonPARs_dat_males_AFR_merge$Num_sim_variants
chrX_nonPARs_dat_males_AFR_merge$TPtoSimVars <- chrX_nonPARs_dat_males_AFR_merge$TP / chrX_nonPARs_dat_males_AFR_merge$Num_sim_variants

## Chromosome Y - EUR ##
chrY_dat_males_EUR <- read.table(
  "../../performance_metrics/EUR/males_20_samples/data/EUR_males_chrY_haploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrY_dat_males_EUR <- chrY_dat_males_EUR[,!(names(chrY_dat_males_EUR) %in% drop)]


# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrY_dat_males_EUR_merge <- merge(chrY_dat_males_EUR, chrY_sim_males_EUR, by = "Sample")

chrY_dat_males_EUR_merge$Chromosome <- "chrY"
chrY_dat_males_EUR_merge$SampleSet <- "EUR"
chrY_dat_males_EUR_merge$FPtoTP <- chrY_dat_males_EUR_merge$FP / chrY_dat_males_EUR_merge$TP
chrY_dat_males_EUR_merge$FNtoTP <- chrY_dat_males_EUR_merge$FN / chrY_dat_males_EUR_merge$TP
chrY_dat_males_EUR_merge$FPtoSimVars <- chrY_dat_males_EUR_merge$FP / chrY_dat_males_EUR_merge$Num_sim_variants
chrY_dat_males_EUR_merge$FNtoSimVars <- chrY_dat_males_EUR_merge$FN / chrY_dat_males_EUR_merge$Num_sim_variants
chrY_dat_males_EUR_merge$TPtoSimVars <- chrY_dat_males_EUR_merge$TP / chrY_dat_males_EUR_merge$Num_sim_variants

## Chromosome Y - ASN ##
chrY_dat_males_ASN <- read.table(
  "../../performance_metrics/ASN/males/data/ASN_males_chrY_haploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrY_dat_males_ASN <- chrY_dat_males_ASN[,!(names(chrY_dat_males_ASN) %in% drop)]


# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrY_dat_males_ASN_merge <- merge(chrY_dat_males_ASN, chrY_sim_males_ASN, by = "Sample")

chrY_dat_males_ASN_merge$Chromosome <- "chrY"
chrY_dat_males_ASN_merge$SampleSet <- "ASN"
chrY_dat_males_ASN_merge$FPtoTP <- chrY_dat_males_ASN_merge$FP / chrY_dat_males_ASN_merge$TP
chrY_dat_males_ASN_merge$FNtoTP <- chrY_dat_males_ASN_merge$FN / chrY_dat_males_ASN_merge$TP
chrY_dat_males_ASN_merge$FPtoSimVars <- chrY_dat_males_ASN_merge$FP / chrY_dat_males_ASN_merge$Num_sim_variants
chrY_dat_males_ASN_merge$FNtoSimVars <- chrY_dat_males_ASN_merge$FN / chrY_dat_males_ASN_merge$Num_sim_variants
chrY_dat_males_ASN_merge$TPtoSimVars <- chrY_dat_males_ASN_merge$TP / chrY_dat_males_ASN_merge$Num_sim_variants

## Chromosome Y - AFR ##
chrY_dat_males_AFR <- read.table(
  "../../performance_metrics/AFR/males/data/AFR_males_chrY_haploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrY_dat_males_AFR <- chrY_dat_males_AFR[,!(names(chrY_dat_males_AFR) %in% drop)]


# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrY_dat_males_AFR_merge <- merge(chrY_dat_males_AFR, chrY_sim_males_AFR, by = "Sample")

chrY_dat_males_AFR_merge$Chromosome <- "chrY"
chrY_dat_males_AFR_merge$SampleSet <- "AFR"
chrY_dat_males_AFR_merge$FPtoTP <- chrY_dat_males_AFR_merge$FP / chrY_dat_males_AFR_merge$TP
chrY_dat_males_AFR_merge$FNtoTP <- chrY_dat_males_AFR_merge$FN / chrY_dat_males_AFR_merge$TP
chrY_dat_males_AFR_merge$FPtoSimVars <- chrY_dat_males_AFR_merge$FP / chrY_dat_males_AFR_merge$Num_sim_variants
chrY_dat_males_AFR_merge$FNtoSimVars <- chrY_dat_males_AFR_merge$FN / chrY_dat_males_AFR_merge$Num_sim_variants
chrY_dat_males_AFR_merge$TPtoSimVars <- chrY_dat_males_AFR_merge$TP / chrY_dat_males_AFR_merge$Num_sim_variants





##                        ##
## Prep data for plotting ##
##                        ##
dat_X_females <- rbind(
  chrX_PARs_dat_EUR_merge, chrX_PARs_dat_ASN_merge, chrX_PARs_dat_AFR_merge,
  chrX_nonPARs_dat_EUR_merge, chrX_nonPARs_dat_ASN_merge, chrX_nonPARs_dat_AFR_merge
  )


dat_XY_males <- rbind(
  chrX_PARs_dat_males_EUR_merge, chrX_PARs_dat_males_ASN_merge, chrX_PARs_dat_males_AFR_merge,
  chrX_nonPARs_dat_males_EUR_merge, chrX_nonPARs_dat_males_ASN_merge, chrX_nonPARs_dat_males_AFR_merge,
  chrY_dat_males_EUR_merge, chrY_dat_males_ASN_merge, chrY_dat_males_AFR_merge
)




##      ##
## Plot ##
##      ##

FPtoSimVars_females <- dat_X_females %>%
  arrange(FPtoSimVars) %>%
  mutate(name = factor(Chromosome, 
                       levels=c(
                                "chrX_PARs",
                                "chrX_nonPARs"
                       ))) %>%
  ggplot( aes(x=name, y=FPtoSimVars, color=SampleSet)) +
  ylim(0,max(c(dat_XY_males$FPtoSimVars, dat_X_females$FPtoSimVars))) +
  geom_boxplot(outlier.shape=NA) +
  geom_point(position=position_jitterdodge(.2)) +
  #geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  scale_color_discrete(name = "Ancestry", labels = c("EUR"="European", 
                                                   "ASN"="Asian",
                                                   "AFR"="African")) +
  ylab("False Positives/Total # Simulated") +
  #theme(legend.position = "none") + 
  xlab("Chromosome") +
  scale_x_discrete(labels=c( 
                            "chrX_PARs" = "PARs",
                            "chrX_nonPARs" = expression("X"["non-PARs"])
                            )) +
  ggtitle("False Positives/Total # Simulated") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1))



# MALES #
FPtoSimVars_males <- dat_XY_males %>%
  arrange(FPtoSimVars) %>%
  mutate(name = factor(Chromosome, 
                       levels=c(
                         "chrX_PARs",
                         "chrX_nonPARs",
                         "chrY"
                       ))) %>%
  ggplot( aes(x=name, y=FPtoSimVars, color=SampleSet)) +
  ylim(0,max(c(dat_XY_males$FPtoSimVars, dat_X_females$FPtoSimVars))) +
  geom_boxplot(outlier.shape=NA) +
  geom_point(position=position_jitterdodge(.2)) +
  #geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  scale_color_discrete(name = "Ancestry", labels = c("EUR"="European", 
                                                     "ASN"="Asian",
                                                     "AFR"="African")) +
  ylab("False Positives/Total # Simulated") +
  #theme(legend.position = "none") + 
  xlab("Chromosome") +
  scale_x_discrete(labels=c( 
    "chrX_PARs" = "PARs",
    "chrX_nonPARs" = expression("X"["non-PARs"]),
    "chrY" = expression("Y"["non-PARs"])
  )) +
  ggtitle("False Positives/Total # Simulated") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1))




FNtoSimVars_females <- dat_X_females %>%
  arrange(FNtoSimVars) %>%
  mutate(name = factor(Chromosome, 
                       levels=c(
                         "chrX_PARs",
                         "chrX_nonPARs"
                       ))) %>%
  ggplot( aes(x=name, y=FNtoSimVars, color=SampleSet)) +
  ylim(0,max(c(dat_XY_males$FNtoSimVars, dat_X_females$FNtoSimVars))) +
  geom_boxplot(outlier.shape=NA) +
  geom_point(position=position_jitterdodge(.2)) +
  #geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  scale_color_discrete(name = "Ancestry", labels = c("EUR"="European", 
                                                     "ASN"="Asian",
                                                     "AFR"="African")) +
  ylab("False Negatives/Total # Simulated") +
  #theme(legend.position = "none") + 
  xlab("Chromosome") +
  scale_x_discrete(labels=c( 
    "chrX_PARs" = "PARs",
    "chrX_nonPARs" = expression("X"["non-PARs"])
  )) +
  ggtitle("False Negatives/Total # Simulated") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1))



# MALES #
FNtoSimVars_males <- dat_XY_males %>%
  arrange(FNtoSimVars) %>%
  mutate(name = factor(Chromosome, 
                       levels=c(
                         "chrX_PARs",
                         "chrX_nonPARs",
                         "chrY"
                       ))) %>%
  ggplot( aes(x=name, y=FNtoSimVars, color=SampleSet)) +
  ylim(0,max(c(dat_XY_males$FNtoSimVars, dat_X_females$FNtoSimVars))) +
  geom_boxplot(outlier.shape=NA) +
  geom_point(position=position_jitterdodge(.2)) +
  #geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  scale_color_discrete(name = "Ancestry", labels = c("EUR"="European", 
                                                     "ASN"="Asian",
                                                     "AFR"="African")) +
  ylab("False Negatives/Total # Simulated") +
  #theme(legend.position = "none") + 
  xlab("Chromosome") +
  scale_x_discrete(labels=c( 
    "chrX_PARs" = "PARs",
    "chrX_nonPARs" = expression("X"["non-PARs"]),
    "chrY" = expression("Y"["non-PARs"])
  )) +
  ggtitle("False Negatives/Total # Simulated") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1))





##              ##
## Output plots ##
##              ##
pdf("../plots/SCC_jitter_females_males_by_ancestry.pdf",
    width = 10, height = 8)

ggarrange(
  FPtoSimVars_females, FPtoSimVars_males,
  FNtoSimVars_females, FNtoSimVars_males,
  ncol = 2,
  nrow = 2,
  labels = c("A", "B", "C", "D"),
  common.legend = T,
  legend = "top"
  )

dev.off()

