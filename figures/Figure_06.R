# Angela Oill
# XY Variant Calling Project

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
chrX_PARs_sim_10 <- read.table("../../simulated_variants/EUR_females_chrX_PARs_diploid_simulated_counts.txt",
                       header = T)

chrX_nonPARs_sim_10 <- read.table("../../simulated_variants/EUR_females_chrX_nonPARs_diploid_simulated_counts.txt",
                       header = T)

chrX_PARs_sim_20 <- read.table("../../simulated_variants/20_samples/EUR_females_chrX_PARs_diploid_simulated_counts.txt",
                               header = T)

chrX_nonPARs_sim_20 <- read.table("../../simulated_variants/20_samples/EUR_females_chrX_nonPARs_diploid_simulated_counts.txt",
                                  header = T)

# MALES #
chrX_PARs_sim_males_10 <- read.table("../../simulated_variants/EUR_males_chrX_PARs_diploid_simulated_counts.txt",
                                  header = T)

chrX_nonPARs_sim_males_10 <- read.table("../../simulated_variants/EUR_males_chrX_nonPARs_haploid_simulated_counts.txt",
                                                         header = T)

chrY_sim_males_10 <- read.table("../../simulated_variants/EUR_males_chrY_haploid_simulated_counts.txt",
                                                         header = T)


chrX_PARs_sim_males_20 <- read.table("../../simulated_variants/20_samples/EUR_males_chrX_PARs_diploid_simulated_counts.txt",
                                     header = T)

chrX_nonPARs_sim_males_20 <- read.table("../../simulated_variants/20_samples/EUR_males_chrX_nonPARs_haploid_simulated_counts.txt",
                                        header = T)

chrY_sim_males_20 <- read.table("../../simulated_variants/20_samples/EUR_males_chrY_haploid_simulated_counts.txt",
                                header = T)

##              ##
## Read in data ##
##              ##
# FEMALES #
## Chromosome X PARs - 10 samples ##
chrX_PARs_dat_10 <- read.table(
  "../../performance_metrics/EUR/females/data/EUR_females_chrX_PARs_diploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrX_PARs_dat_10 <- chrX_PARs_dat_10[,!(names(chrX_PARs_dat_10) %in% drop)]


# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_PARs_dat_10_merge <- merge(chrX_PARs_dat_10, chrX_PARs_sim_10, by = "Sample")

chrX_PARs_dat_10_merge$Chromosome <- "chrX_PARs"
chrX_PARs_dat_10_merge$SampleSet <- "10_samples"
chrX_PARs_dat_10_merge$FPtoTP <- chrX_PARs_dat_10_merge$FP / chrX_PARs_dat_10_merge$TP
chrX_PARs_dat_10_merge$FNtoTP <- chrX_PARs_dat_10_merge$FN / chrX_PARs_dat_10_merge$TP
chrX_PARs_dat_10_merge$FPtoSimVars <- chrX_PARs_dat_10_merge$FP / chrX_PARs_dat_10_merge$Num_sim_variants
chrX_PARs_dat_10_merge$FNtoSimVars <- chrX_PARs_dat_10_merge$FN / chrX_PARs_dat_10_merge$Num_sim_variants
chrX_PARs_dat_10_merge$TPtoSimVars <- chrX_PARs_dat_10_merge$TP / chrX_PARs_dat_10_merge$Num_sim_variants

## Chromosome X PARs - 20 samples ##
chrX_PARs_dat_20 <- read.table(
  "../../performance_metrics/EUR/females_20_samples/data/EUR_females_chrX_PARs_diploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrX_PARs_dat_20 <- chrX_PARs_dat_20[,!(names(chrX_PARs_dat_20) %in% drop)]


# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_PARs_dat_20_merge <- merge(chrX_PARs_dat_20, chrX_PARs_sim_20, by = "Sample")

chrX_PARs_dat_20_merge$Chromosome <- "chrX_PARs"
chrX_PARs_dat_20_merge$SampleSet <- "20_samples"
chrX_PARs_dat_20_merge$FPtoTP <- chrX_PARs_dat_20_merge$FP / chrX_PARs_dat_20_merge$TP
chrX_PARs_dat_20_merge$FNtoTP <- chrX_PARs_dat_20_merge$FN / chrX_PARs_dat_20_merge$TP
chrX_PARs_dat_20_merge$FPtoSimVars <- chrX_PARs_dat_20_merge$FP / chrX_PARs_dat_20_merge$Num_sim_variants
chrX_PARs_dat_20_merge$FNtoSimVars <- chrX_PARs_dat_20_merge$FN / chrX_PARs_dat_20_merge$Num_sim_variants
chrX_PARs_dat_20_merge$TPtoSimVars <- chrX_PARs_dat_20_merge$TP / chrX_PARs_dat_20_merge$Num_sim_variants

# merge same samples and do paired t-test
head(chrX_PARs_dat_20_merge)
head(chrX_PARs_dat_10_merge)

chrX_PARs_dat_10_20_merge <- merge(chrX_PARs_dat_10_merge, chrX_PARs_dat_20_merge, by = "Sample")

t.test(chrX_PARs_dat_10_20_merge$FPtoSimVars.x, chrX_PARs_dat_10_20_merge$FPtoSimVars.y, paired = TRUE, alternative = "two.sided")
# t = 1.7392, df = 9, p-value = 0.116
t.test(chrX_PARs_dat_10_20_merge$FNtoSimVars.x, chrX_PARs_dat_10_20_merge$FNtoSimVars.y, paired = TRUE, alternative = "two.sided")
# t = -0.2591, df = 9, p-value = 0.8014


## Chromosome X non-PARs - 10 samples ##
chrX_nonPARs_dat_10 <- read.table(
  "../../performance_metrics/EUR/females/data/EUR_females_chrX_nonPARs_diploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrX_nonPARs_dat_10 <- chrX_nonPARs_dat_10[,!(names(chrX_nonPARs_dat_10) %in% drop)]


# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_nonPARs_dat_10_merge <- merge(chrX_nonPARs_dat_10, chrX_nonPARs_sim_10, by = "Sample")

chrX_nonPARs_dat_10_merge$Chromosome <- "chrX_nonPARs"
chrX_nonPARs_dat_10_merge$SampleSet <- "10_samples"
chrX_nonPARs_dat_10_merge$FPtoTP <- chrX_nonPARs_dat_10_merge$FP / chrX_nonPARs_dat_10_merge$TP
chrX_nonPARs_dat_10_merge$FNtoTP <- chrX_nonPARs_dat_10_merge$FN / chrX_nonPARs_dat_10_merge$TP
chrX_nonPARs_dat_10_merge$FPtoSimVars <- chrX_nonPARs_dat_10_merge$FP / chrX_nonPARs_dat_10_merge$Num_sim_variants
chrX_nonPARs_dat_10_merge$FNtoSimVars <- chrX_nonPARs_dat_10_merge$FN / chrX_nonPARs_dat_10_merge$Num_sim_variants
chrX_nonPARs_dat_10_merge$TPtoSimVars <- chrX_nonPARs_dat_10_merge$TP / chrX_nonPARs_dat_10_merge$Num_sim_variants

## Chromosome X non-PARs - 20 samples ##
chrX_nonPARs_dat_20 <- read.table(
  "../../performance_metrics/EUR/females_20_samples/data/EUR_females_chrX_nonPARs_diploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrX_nonPARs_dat_20 <- chrX_nonPARs_dat_20[,!(names(chrX_nonPARs_dat_20) %in% drop)]


# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_nonPARs_dat_20_merge <- merge(chrX_nonPARs_dat_20, chrX_nonPARs_sim_20, by = "Sample")

chrX_nonPARs_dat_20_merge$Chromosome <- "chrX_nonPARs"
chrX_nonPARs_dat_20_merge$SampleSet <- "20_samples"
chrX_nonPARs_dat_20_merge$FPtoTP <- chrX_nonPARs_dat_20_merge$FP / chrX_nonPARs_dat_20_merge$TP
chrX_nonPARs_dat_20_merge$FNtoTP <- chrX_nonPARs_dat_20_merge$FN / chrX_nonPARs_dat_20_merge$TP
chrX_nonPARs_dat_20_merge$FPtoSimVars <- chrX_nonPARs_dat_20_merge$FP / chrX_nonPARs_dat_20_merge$Num_sim_variants
chrX_nonPARs_dat_20_merge$FNtoSimVars <- chrX_nonPARs_dat_20_merge$FN / chrX_nonPARs_dat_20_merge$Num_sim_variants
chrX_nonPARs_dat_20_merge$TPtoSimVars <- chrX_nonPARs_dat_20_merge$TP / chrX_nonPARs_dat_20_merge$Num_sim_variants



# merge same samples and do paired t-test
head(chrX_nonPARs_dat_20_merge)
head(chrX_nonPARs_dat_10_merge)

chrX_nonPARs_dat_10_20_merge <- merge(chrX_nonPARs_dat_10_merge, chrX_nonPARs_dat_20_merge, by = "Sample")

t.test(chrX_nonPARs_dat_10_20_merge$FPtoSimVars.x, chrX_nonPARs_dat_10_20_merge$FPtoSimVars.y, paired = TRUE, alternative = "two.sided")
# t = 4.5659, df = 9, p-value = 0.001355
t.test(chrX_nonPARs_dat_10_20_merge$FNtoSimVars.x, chrX_nonPARs_dat_10_20_merge$FNtoSimVars.y, paired = TRUE, alternative = "two.sided")
# t = -0.70244, df = 9, p-value = 0.5002

mean(chrX_nonPARs_dat_10_20_merge$FPtoSimVars.x)
mean(chrX_nonPARs_dat_10_20_merge$FPtoSimVars.y)

mean(chrX_nonPARs_dat_10_20_merge$FP.x)
mean(chrX_nonPARs_dat_10_20_merge$FP.y)

# MALES #
## Chromosome X PARs - 10 samples ##
chrX_PARs_dat_males_10 <- read.table(
  "../../performance_metrics/EUR/males/data/EUR_males_chrX_PARs_diploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrX_PARs_dat_males_10 <- chrX_PARs_dat_males_10[,!(names(chrX_PARs_dat_males_10) %in% drop)]


# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_PARs_dat_males_10_merge <- merge(chrX_PARs_dat_males_10, chrX_PARs_sim_males_10, by = "Sample")

chrX_PARs_dat_males_10_merge$Chromosome <- "chrX_PARs"
chrX_PARs_dat_males_10_merge$SampleSet <- "10_samples"
chrX_PARs_dat_males_10_merge$FPtoTP <- chrX_PARs_dat_males_10_merge$FP / chrX_PARs_dat_males_10_merge$TP
chrX_PARs_dat_males_10_merge$FNtoTP <- chrX_PARs_dat_males_10_merge$FN / chrX_PARs_dat_males_10_merge$TP
chrX_PARs_dat_males_10_merge$FPtoSimVars <- chrX_PARs_dat_males_10_merge$FP / chrX_PARs_dat_males_10_merge$Num_sim_variants
chrX_PARs_dat_males_10_merge$FNtoSimVars <- chrX_PARs_dat_males_10_merge$FN / chrX_PARs_dat_males_10_merge$Num_sim_variants
chrX_PARs_dat_males_10_merge$TPtoSimVars <- chrX_PARs_dat_males_10_merge$TP / chrX_PARs_dat_males_10_merge$Num_sim_variants

## Chromosome X PARs - 20 samples ##
chrX_PARs_dat_males_20 <- read.table(
  "../../performance_metrics/EUR/males_20_samples/data/EUR_males_chrX_PARs_diploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrX_PARs_dat_males_20 <- chrX_PARs_dat_males_20[,!(names(chrX_PARs_dat_males_20) %in% drop)]


# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_PARs_dat_males_20_merge <- merge(chrX_PARs_dat_males_20, chrX_PARs_sim_males_20, by = "Sample")

chrX_PARs_dat_males_20_merge$Chromosome <- "chrX_PARs"
chrX_PARs_dat_males_20_merge$SampleSet <- "20_samples"
chrX_PARs_dat_males_20_merge$FPtoTP <- chrX_PARs_dat_males_20_merge$FP / chrX_PARs_dat_males_20_merge$TP
chrX_PARs_dat_males_20_merge$FNtoTP <- chrX_PARs_dat_males_20_merge$FN / chrX_PARs_dat_males_20_merge$TP
chrX_PARs_dat_males_20_merge$FPtoSimVars <- chrX_PARs_dat_males_20_merge$FP / chrX_PARs_dat_males_20_merge$Num_sim_variants
chrX_PARs_dat_males_20_merge$FNtoSimVars <- chrX_PARs_dat_males_20_merge$FN / chrX_PARs_dat_males_20_merge$Num_sim_variants
chrX_PARs_dat_males_20_merge$TPtoSimVars <- chrX_PARs_dat_males_20_merge$TP / chrX_PARs_dat_males_20_merge$Num_sim_variants


# merge same samples and do paired t-test
head(chrX_PARs_dat_males_20_merge)
head(chrX_PARs_dat_males_10_merge)

chrX_PARs_dat_males_10_20_merge <- merge(chrX_PARs_dat_males_10_merge, chrX_PARs_dat_males_20_merge, by = "Sample")

t.test(chrX_PARs_dat_males_10_20_merge$FPtoSimVars.x, chrX_PARs_dat_males_10_20_merge$FPtoSimVars.y, paired = TRUE, alternative = "two.sided")
# t = 2.2213, df = 9, p-value = 0.05345
t.test(chrX_PARs_dat_males_10_20_merge$FNtoSimVars.x, chrX_PARs_dat_males_10_20_merge$FNtoSimVars.y, paired = TRUE, alternative = "two.sided")
# t = -18.899, df = 9, p-value = 1.495e-08

mean(chrX_PARs_dat_males_10_20_merge$FNtoSimVars.x)
mean(chrX_PARs_dat_males_10_20_merge$FNtoSimVars.y)

mean(chrX_PARs_dat_males_10_20_merge$FN.x)
mean(chrX_PARs_dat_males_10_20_merge$FN.y)

## Chromosome X non-PARs - 10 samples ##
chrX_nonPARs_dat_males_10 <- read.table(
  "../../performance_metrics/EUR/males/data/EUR_males_chrX_nonPARs_haploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrX_nonPARs_dat_males_10 <- chrX_nonPARs_dat_males_10[,!(names(chrX_nonPARs_dat_males_10) %in% drop)]


# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_nonPARs_dat_males_10_merge <- merge(chrX_nonPARs_dat_males_10, chrX_nonPARs_sim_males_10, by = "Sample")

chrX_nonPARs_dat_males_10_merge$Chromosome <- "chrX_nonPARs"
chrX_nonPARs_dat_males_10_merge$SampleSet <- "10_samples"
chrX_nonPARs_dat_males_10_merge$FPtoTP <- chrX_nonPARs_dat_males_10_merge$FP / chrX_nonPARs_dat_males_10_merge$TP
chrX_nonPARs_dat_males_10_merge$FNtoTP <- chrX_nonPARs_dat_males_10_merge$FN / chrX_nonPARs_dat_males_10_merge$TP
chrX_nonPARs_dat_males_10_merge$FPtoSimVars <- chrX_nonPARs_dat_males_10_merge$FP / chrX_nonPARs_dat_males_10_merge$Num_sim_variants
chrX_nonPARs_dat_males_10_merge$FNtoSimVars <- chrX_nonPARs_dat_males_10_merge$FN / chrX_nonPARs_dat_males_10_merge$Num_sim_variants
chrX_nonPARs_dat_males_10_merge$TPtoSimVars <- chrX_nonPARs_dat_males_10_merge$TP / chrX_nonPARs_dat_males_10_merge$Num_sim_variants

## Chromosome X non-PARs - 20 samples ##
chrX_nonPARs_dat_males_20 <- read.table(
  "../../performance_metrics/EUR/males_20_samples/data/EUR_males_chrX_nonPARs_haploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrX_nonPARs_dat_males_20 <- chrX_nonPARs_dat_males_20[,!(names(chrX_nonPARs_dat_males_20) %in% drop)]


# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_nonPARs_dat_males_20_merge <- merge(chrX_nonPARs_dat_males_20, chrX_nonPARs_sim_males_20, by = "Sample")

chrX_nonPARs_dat_males_20_merge$Chromosome <- "chrX_nonPARs"
chrX_nonPARs_dat_males_20_merge$SampleSet <- "20_samples"
chrX_nonPARs_dat_males_20_merge$FPtoTP <- chrX_nonPARs_dat_males_20_merge$FP / chrX_nonPARs_dat_males_20_merge$TP
chrX_nonPARs_dat_males_20_merge$FNtoTP <- chrX_nonPARs_dat_males_20_merge$FN / chrX_nonPARs_dat_males_20_merge$TP
chrX_nonPARs_dat_males_20_merge$FPtoSimVars <- chrX_nonPARs_dat_males_20_merge$FP / chrX_nonPARs_dat_males_20_merge$Num_sim_variants
chrX_nonPARs_dat_males_20_merge$FNtoSimVars <- chrX_nonPARs_dat_males_20_merge$FN / chrX_nonPARs_dat_males_20_merge$Num_sim_variants
chrX_nonPARs_dat_males_20_merge$TPtoSimVars <- chrX_nonPARs_dat_males_20_merge$TP / chrX_nonPARs_dat_males_20_merge$Num_sim_variants

# merge same samples and do paired t-test
head(chrX_nonPARs_dat_males_20_merge)
head(chrX_nonPARs_dat_males_10_merge)

chrX_nonPARs_dat_males_10_20_merge <- merge(chrX_nonPARs_dat_males_10_merge, chrX_nonPARs_dat_males_20_merge, by = "Sample")

t.test(chrX_nonPARs_dat_males_10_20_merge$FPtoSimVars.x, chrX_nonPARs_dat_males_10_20_merge$FPtoSimVars.y, paired = TRUE, alternative = "two.sided")
# t = 2.2213, df = 9, p-value = 0.7689
t.test(chrX_nonPARs_dat_males_10_20_merge$FNtoSimVars.x, chrX_nonPARs_dat_males_10_20_merge$FNtoSimVars.y, paired = TRUE, alternative = "two.sided")
# t = -18.899, df = 9, p-value = 2.672e-11

mean(chrX_nonPARs_dat_males_10_20_merge$FNtoSimVars.x)
mean(chrX_nonPARs_dat_males_10_20_merge$FNtoSimVars.y)

mean(chrX_nonPARs_dat_males_10_20_merge$FN.x)
mean(chrX_nonPARs_dat_males_10_20_merge$FN.y)


## Chromosome Y - 10 samples ##
chrY_dat_males_10 <- read.table(
  "../../performance_metrics/EUR/males/data/EUR_males_chrY_haploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrY_dat_males_10 <- chrY_dat_males_10[,!(names(chrY_dat_males_10) %in% drop)]


# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrY_dat_males_10_merge <- merge(chrY_dat_males_10, chrY_sim_males_10, by = "Sample")

chrY_dat_males_10_merge$Chromosome <- "chrY"
chrY_dat_males_10_merge$SampleSet <- "10_samples"
chrY_dat_males_10_merge$FPtoTP <- chrY_dat_males_10_merge$FP / chrY_dat_males_10_merge$TP
chrY_dat_males_10_merge$FNtoTP <- chrY_dat_males_10_merge$FN / chrY_dat_males_10_merge$TP
chrY_dat_males_10_merge$FPtoSimVars <- chrY_dat_males_10_merge$FP / chrY_dat_males_10_merge$Num_sim_variants
chrY_dat_males_10_merge$FNtoSimVars <- chrY_dat_males_10_merge$FN / chrY_dat_males_10_merge$Num_sim_variants
chrY_dat_males_10_merge$TPtoSimVars <- chrY_dat_males_10_merge$TP / chrY_dat_males_10_merge$Num_sim_variants

## Chromosome Y - 20 samples ##
chrY_dat_males_20 <- read.table(
  "../../performance_metrics/EUR/males_20_samples/data/EUR_males_chrY_haploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrY_dat_males_20 <- chrY_dat_males_20[,!(names(chrY_dat_males_20) %in% drop)]


# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrY_dat_males_20_merge <- merge(chrY_dat_males_20, chrY_sim_males_20, by = "Sample")

chrY_dat_males_20_merge$Chromosome <- "chrY"
chrY_dat_males_20_merge$SampleSet <- "20_samples"
chrY_dat_males_20_merge$FPtoTP <- chrY_dat_males_20_merge$FP / chrY_dat_males_20_merge$TP
chrY_dat_males_20_merge$FNtoTP <- chrY_dat_males_20_merge$FN / chrY_dat_males_20_merge$TP
chrY_dat_males_20_merge$FPtoSimVars <- chrY_dat_males_20_merge$FP / chrY_dat_males_20_merge$Num_sim_variants
chrY_dat_males_20_merge$FNtoSimVars <- chrY_dat_males_20_merge$FN / chrY_dat_males_20_merge$Num_sim_variants
chrY_dat_males_20_merge$TPtoSimVars <- chrY_dat_males_20_merge$TP / chrY_dat_males_20_merge$Num_sim_variants



# merge same samples and do paired t-test
head(chrY_dat_males_20_merge)
head(chrY_dat_males_10_merge)

chrY_dat_males_10_20_merge <- merge(chrY_dat_males_10_merge, chrY_dat_males_20_merge, by = "Sample")

t.test(chrY_dat_males_10_20_merge$FPtoSimVars.x, chrY_dat_males_10_20_merge$FPtoSimVars.y, paired = TRUE, alternative = "two.sided")
# t = 0.59399, df = 9, p-value = 0.5671
t.test(chrY_dat_males_10_20_merge$FNtoSimVars.x, chrY_dat_males_10_20_merge$FNtoSimVars.y, paired = TRUE, alternative = "two.sided")
# t = 4.186, df = 9, p-value = 0.002355

mean(chrY_dat_males_10_20_merge$FNtoSimVars.x)
mean(chrY_dat_males_10_20_merge$FNtoSimVars.y)

mean(chrY_dat_males_10_20_merge$FN.x)
mean(chrY_dat_males_10_20_merge$FN.y)


##                        ##
## Prep data for plotting ##
##                        ##
dat_X_females <- rbind(
  chrX_PARs_dat_10_merge, chrX_PARs_dat_20_merge,
  chrX_nonPARs_dat_10_merge, chrX_nonPARs_dat_20_merge
  )


dat_XY_males <- rbind(
  chrX_PARs_dat_males_10_merge, chrX_PARs_dat_males_20_merge,
  chrX_nonPARs_dat_males_10_merge, chrX_nonPARs_dat_males_20_merge,
  chrY_dat_males_10_merge, chrY_dat_males_20_merge
)


# Subset weight data before treatment
#before <- subset(chrX_PARs_dat_males_10_20_merge,  SampleSet.x == "10_samples", FNtoSimVars.x,
#                 drop = TRUE)
# subset weight data after treatment
#after <- subset(chrX_PARs_dat_males_10_20_merge,  SampleSet.y == "20_samples", FNtoSimVars.y,
#                drop = TRUE)
# Plot paired data
#library(PairedData)
#pd <- paired(before, after)
#plot(pd, type = "profile") + theme_bw()

dat_X_females_10_20 <- rbind(
  chrX_PARs_dat_10_20_merge,
  chrX_nonPARs_dat_10_20_merge
)

myvars <- c("Sample", "SampleSet.x", "Chromosome.x", "FPtoSimVars.x", "FNtoSimVars.x")
dat_X_females_10_tmp <- dat_X_females_10_20[myvars]
colnames(dat_X_females_10_tmp) <- c("Sample", "SampleSet", "Chromosome", "FPtoSimVars", "FNtoSimVars")

myvars <- c("Sample", "SampleSet.y", "Chromosome.y", "FPtoSimVars.y", "FNtoSimVars.y")
dat_X_females_20_tmp <- dat_X_females_10_20[myvars]
colnames(dat_X_females_20_tmp) <- c("Sample", "SampleSet", "Chromosome", "FPtoSimVars", "FNtoSimVars")


dat_X_females_10_20_rows <- rbind(dat_X_females_10_tmp, dat_X_females_20_tmp)

dat_XY_males <- rbind(
  chrX_PARs_dat_males_10_20_merge,
  chrX_nonPARs_dat_males_10_20_merge,
  chrY_dat_males_10_20_merge
)

myvars <- c("Sample", "SampleSet.x", "Chromosome.x", "FPtoSimVars.x", "FNtoSimVars.x")
dat_X_males_10_tmp <- dat_XY_males[myvars]
colnames(dat_X_males_10_tmp) <- c("Sample", "SampleSet", "Chromosome", "FPtoSimVars", "FNtoSimVars")

myvars <- c("Sample", "SampleSet.y", "Chromosome.y", "FPtoSimVars.y", "FNtoSimVars.y")
dat_X_males_20_tmp <- dat_XY_males[myvars]
colnames(dat_X_males_20_tmp) <- c("Sample", "SampleSet", "Chromosome", "FPtoSimVars", "FNtoSimVars")


dat_X_males_10_20_rows <- rbind(dat_X_males_10_tmp, dat_X_males_20_tmp)


# Plot paired data #
# FEMALES #
femalesFPtoSimVars1020matched <- dat_X_females_10_20_rows %>%
  arrange(FPtoSimVars) %>%
  mutate(name = factor(Chromosome, 
                       levels=c(
                         "chrX_PARs",
                         "chrX_nonPARs"
                       ))) %>%
  ggplot( aes(x=name, y=FPtoSimVars, color=SampleSet)) +
  ylim(0,max(c(dat_X_males_10_20_rows$FPtoSimVars, dat_X_females_10_20_rows$FPtoSimVars))) +
  geom_boxplot(outlier.shape=NA) +
  geom_point(position=position_jitterdodge(.2)) +
  #geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  scale_color_manual(values = c("blue", "lightblue"),
                     name = "# Joint Genotyped", labels = c("10_samples"="10 Samples", 
                                                            "20_samples"="20 Samples")) +
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

femalesFNtoSimVars1020matched <- dat_X_females_10_20_rows %>%
  arrange(FNtoSimVars) %>%
  mutate(name = factor(Chromosome, 
                       levels=c(
                         "chrX_PARs",
                         "chrX_nonPARs"
                       ))) %>%
  ggplot( aes(x=name, y=FNtoSimVars, color=SampleSet)) +
  ylim(0,max(c(dat_X_males_10_20_rows$FNtoSimVars, dat_X_females_10_20_rows$FNtoSimVars))) +
  geom_boxplot(outlier.shape=NA) +
  geom_point(position=position_jitterdodge(.2)) +
  #geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  scale_color_manual(values = c("blue", "lightblue"),
                     name = "# Joint Genotyped", labels = c("10_samples"="10 Samples", 
                                                            "20_samples"="20 Samples")) +
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
malesFPtoSimVars1020matched <- dat_X_males_10_20_rows %>%
  arrange(FPtoSimVars) %>%
  mutate(name = factor(Chromosome, 
                       levels=c(
                         "chrX_PARs",
                         "chrX_nonPARs",
                         "chrY"
                       ))) %>%
  ggplot( aes(x=name, y=FPtoSimVars, color=SampleSet)) +
  ylim(0,max(c(dat_X_males_10_20_rows$FPtoSimVars, dat_X_females_10_20_rows$FPtoSimVars))) +
  geom_boxplot(outlier.shape=NA) +
  geom_point(position=position_jitterdodge(.2)) +
  #geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  scale_color_manual(values = c("blue", "lightblue"),
                     name = "# Joint Genotyped", labels = c("10_samples"="10 Samples", 
                                                            "20_samples"="20 Samples")) +
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

malesFNtoSimVars1020matched <- dat_X_males_10_20_rows %>%
  arrange(FNtoSimVars) %>%
  mutate(name = factor(Chromosome, 
                       levels=c(
                         "chrX_PARs",
                         "chrX_nonPARs",
                         "chrY"
                       ))) %>%
  ggplot( aes(x=name, y=FNtoSimVars, color=SampleSet)) +
  ylim(0,max(c(dat_X_males_10_20_rows$FNtoSimVars, dat_X_females_10_20_rows$FNtoSimVars))) +
  geom_boxplot(outlier.shape=NA) +
  geom_point(position=position_jitterdodge(.2)) +
  #geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  scale_color_manual(values = c("blue", "lightblue"),
                     name = "# Joint Genotyped", labels = c("10_samples"="10 Samples", 
                                                            "20_samples"="20 Samples")) +
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
  scale_color_manual(values = c("blue", "lightblue"),
                     name = "# Joint Genotyped", labels = c("10_samples"="10 Samples", 
                                                            "20_samples"="20 Samples")) +
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
  scale_color_manual(values = c("blue", "lightblue"),
                     name = "# Joint Genotyped", labels = c("10_samples"="10 Samples", 
                                                            "20_samples"="20 Samples")) +
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
  scale_color_manual(values = c("blue", "lightblue"),
                     name = "# Joint Genotyped", labels = c("10_samples"="10 Samples", 
                                                            "20_samples"="20 Samples")) +
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
FNtoSimVars_males <-dat_XY_males %>%
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
  scale_color_manual(values = c("blue", "lightblue"),
                     name = "# Joint Genotyped", labels = c("10_samples"="10 Samples", 
                                                            "20_samples"="20 Samples")) +
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
pdf("../plots/SCC_jitter_females_males_10vs20.pdf",
    width = 7.5, height = 7.5)

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


pdf("../../plots/SCC_jitter_females_males_10vs20_10matched.pdf",
    width = 7.5, height = 7.5)
ggarrange(
  femalesFPtoSimVars1020matched, malesFPtoSimVars1020matched,
  femalesFNtoSimVars1020matched, malesFNtoSimVars1020matched,
  ncol = 2,
  nrow = 2,
  labels = c("A", "B", "C", "D"),
  common.legend = T,
  legend = "top"
)
dev.off()
