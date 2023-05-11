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
chrX_PARs_sim <- read.table("../simulated_variants/EUR_females_chrX_PARs_diploid_simulated_counts.txt",
                       header = T)

chrX_nonPARs_no_XTR_no_amplicons_sim <- read.table("../simulated_variants/EUR_females_chrX_nonPARs_diploid_simulated_counts_nonPARsminusXTRminusAmplicons.txt",
                       header = T)

chrX_XTR_sim <- read.table("../simulated_variants/EUR_females_chrX_nonPARs_diploid_simulated_counts_XTR.txt",
                                      header = T)

chrX_Amplicons_sim <- read.table("../simulated_variants/EUR_females_chrX_nonPARs_diploid_simulated_counts_Amplicons.txt",
                                              header = T)


# MALES #
chrX_PARs_sim_males <- read.table("../simulated_variants/EUR_males_chrX_PARs_diploid_simulated_counts.txt",
                                  header = T)

chrX_nonPARs_no_XTR_no_amplicons_sim_males <- read.table("../simulated_variants/EUR_males_chrX_nonPARs_haploid_simulated_counts_nonPARsminusXTRminusAmplicons.txt",
                                                         header = T)

chrX_XTR_sim_males <- read.table("../simulated_variants/EUR_males_chrX_nonPARs_haploid_simulated_counts_XTR.txt",
                                 header = T)

chrX_Amplicons_sim_males <- read.table("../simulated_variants/EUR_males_chrX_nonPARs_haploid_simulated_counts_Amplicons.txt",
                                       header = T)

chrY_nonPARs_no_XTR_no_amplicons_sim_males <- read.table("../simulated_variants/EUR_males_chrY_haploid_simulated_counts_nonPARsminusXTRminusAmplicons.txt",
                                                         header = T)

chrY_XTR_sim_males <- read.table("../simulated_variants/EUR_males_chrY_haploid_simulated_counts_XTR.txt",
                                 header = T)

chrY_Amplicons_sim_males <- read.table("../simulated_variants/EUR_males_chrY_haploid_simulated_counts_Amplicons.txt",
                                       header = T)


##              ##
## Read in data ##
##              ##
# FEMALES #
## Chromosome X PARs ##
chrX_PARs_dat <- read.table(
  "../performance_metrics/EUR/females/data/EUR_females_chrX_PARs_diploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrX_PARs_dat <- chrX_PARs_dat[,!(names(chrX_PARs_dat) %in% drop)]


# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_PARs_dat_merge <- merge(chrX_PARs_dat, chrX_PARs_sim, by = "Sample")

chrX_PARs_dat_merge$Chromosome <- "chrX_PARs"
chrX_PARs_dat_merge$FPtoTP <- chrX_PARs_dat_merge$FP / chrX_PARs_dat_merge$TP
chrX_PARs_dat_merge$FNtoTP <- chrX_PARs_dat_merge$FN / chrX_PARs_dat_merge$TP
chrX_PARs_dat_merge$FPtoSimVars <- chrX_PARs_dat_merge$FP / chrX_PARs_dat_merge$Num_sim_variants
chrX_PARs_dat_merge$FNtoSimVars <- chrX_PARs_dat_merge$FN / chrX_PARs_dat_merge$Num_sim_variants
chrX_PARs_dat_merge$TPtoSimVars <- chrX_PARs_dat_merge$TP / chrX_PARs_dat_merge$Num_sim_variants


## Chromosome X nonPARs minus XTR minus amplicons ##
chrX_nonPARs_no_XTR_no_amplicons <- read.table(
  "../performance_metrics/EUR/females/by_region/EUR_females_chrX_nonPARs_diploid_golden_vs_called_performance_metrics_nonPARsminusXTRminusAmplicons.txt",
  header = T
)

# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_nonPARs_no_XTR_no_amplicons_merge <- merge(chrX_nonPARs_no_XTR_no_amplicons, 
                                                           chrX_nonPARs_no_XTR_no_amplicons_sim, by = "Sample")

chrX_nonPARs_no_XTR_no_amplicons_merge$Chromosome <- "chrX_nonPARs_no_XTR_no_amplicons"
chrX_nonPARs_no_XTR_no_amplicons_merge$FPtoTP <- chrX_nonPARs_no_XTR_no_amplicons_merge$FP / chrX_nonPARs_no_XTR_no_amplicons_merge$TP
chrX_nonPARs_no_XTR_no_amplicons_merge$FNtoTP <- chrX_nonPARs_no_XTR_no_amplicons_merge$FN / chrX_nonPARs_no_XTR_no_amplicons_merge$TP
chrX_nonPARs_no_XTR_no_amplicons_merge$FPtoSimVars <- chrX_nonPARs_no_XTR_no_amplicons_merge$FP / chrX_nonPARs_no_XTR_no_amplicons_merge$Num_sim_variants
chrX_nonPARs_no_XTR_no_amplicons_merge$FNtoSimVars <- chrX_nonPARs_no_XTR_no_amplicons_merge$FN / chrX_nonPARs_no_XTR_no_amplicons_merge$Num_sim_variants
chrX_nonPARs_no_XTR_no_amplicons_merge$TPtoSimVars <- chrX_nonPARs_no_XTR_no_amplicons_merge$TP / chrX_nonPARs_no_XTR_no_amplicons_merge$Num_sim_variants


## Chromosome X XTR  ##
chrX_XTR <- read.table(
  "../performance_metrics/EUR/females/by_region/EUR_females_chrX_nonPARs_diploid_golden_vs_called_performance_metrics_XTR.txt",
  header = T
)

# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_XTR_merge <- merge(chrX_XTR, chrX_XTR_sim, by = "Sample")

chrX_XTR_merge$Chromosome <- "chrX_XTR"
chrX_XTR_merge$FPtoTP <- chrX_XTR_merge$FP / chrX_XTR_merge$TP
chrX_XTR_merge$FNtoTP <- chrX_XTR_merge$FN / chrX_XTR_merge$TP
chrX_XTR_merge$FPtoSimVars <- chrX_XTR_merge$FP / chrX_XTR_merge$Num_sim_variants
chrX_XTR_merge$FNtoSimVars <- chrX_XTR_merge$FN / chrX_XTR_merge$Num_sim_variants
chrX_XTR_merge$TPtoSimVars <- chrX_XTR_merge$TP / chrX_XTR_merge$Num_sim_variants

## Chromosome X Amplicons ##
chrX_Amplicons <- read.table(
  "../performance_metrics/EUR/females/by_region/EUR_females_chrX_nonPARs_diploid_golden_vs_called_performance_metrics_Amplicons.txt",
  header = T
)

# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_Amplicons_merge <- merge(chrX_Amplicons, chrX_Amplicons_sim, by = "Sample")

chrX_Amplicons_merge$Chromosome <- "chrX_Amplicons"
chrX_Amplicons_merge$FPtoTP <- chrX_Amplicons_merge$FP / chrX_Amplicons_merge$TP
chrX_Amplicons_merge$FNtoTP <- chrX_Amplicons_merge$FN / chrX_Amplicons_merge$TP
chrX_Amplicons_merge$FPtoSimVars <- chrX_Amplicons_merge$FP / chrX_Amplicons_merge$Num_sim_variants
chrX_Amplicons_merge$FNtoSimVars <- chrX_Amplicons_merge$FN / chrX_Amplicons_merge$Num_sim_variants
chrX_Amplicons_merge$TPtoSimVars <- chrX_Amplicons_merge$TP / chrX_Amplicons_merge$Num_sim_variants


# MALES #
## Chromosome X PARs ##
chrX_PARs_dat_males <- read.table(
  "../performance_metrics/EUR/males/data/EUR_males_chrX_PARs_diploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrX_PARs_dat_males <- chrX_PARs_dat_males[,!(names(chrX_PARs_dat_males) %in% drop)]


# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_PARs_dat_males_merge <- merge(chrX_PARs_dat_males, chrX_PARs_sim_males, by = "Sample")

chrX_PARs_dat_males_merge$Chromosome <- "chrX_PARs"
chrX_PARs_dat_males_merge$FPtoTP <- chrX_PARs_dat_males_merge$FP / chrX_PARs_dat_males_merge$TP
chrX_PARs_dat_males_merge$FNtoTP <- chrX_PARs_dat_males_merge$FN / chrX_PARs_dat_males_merge$TP
chrX_PARs_dat_males_merge$FPtoSimVars <- chrX_PARs_dat_males_merge$FP / chrX_PARs_dat_males_merge$Num_sim_variants
chrX_PARs_dat_males_merge$FNtoSimVars <- chrX_PARs_dat_males_merge$FN / chrX_PARs_dat_males_merge$Num_sim_variants
chrX_PARs_dat_males_merge$TPtoSimVars <- chrX_PARs_dat_males_merge$TP / chrX_PARs_dat_males_merge$Num_sim_variants


## Chromosome X nonPARs minus XTR minus amplicons ##
chrX_nonPARs_no_XTR_no_amplicons_males <- read.table(
  "../performance_metrics/EUR/males/by_region/EUR_males_chrX_nonPARs_haploid_golden_vs_called_performance_metrics_nonPARsminusXTRminusAmplicons.txt",
  header = T
)

# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_nonPARs_no_XTR_no_amplicons_males_merge <- merge(chrX_nonPARs_no_XTR_no_amplicons_males, 
                                                      chrX_nonPARs_no_XTR_no_amplicons_sim_males, by = "Sample")

chrX_nonPARs_no_XTR_no_amplicons_males_merge$Chromosome <- "chrX_nonPARs_no_XTR_no_amplicons"
chrX_nonPARs_no_XTR_no_amplicons_males_merge$FPtoTP <- chrX_nonPARs_no_XTR_no_amplicons_males_merge$FP / chrX_nonPARs_no_XTR_no_amplicons_males_merge$TP
chrX_nonPARs_no_XTR_no_amplicons_males_merge$FNtoTP <- chrX_nonPARs_no_XTR_no_amplicons_males_merge$FN / chrX_nonPARs_no_XTR_no_amplicons_males_merge$TP
chrX_nonPARs_no_XTR_no_amplicons_males_merge$FPtoSimVars <- chrX_nonPARs_no_XTR_no_amplicons_males_merge$FP / chrX_nonPARs_no_XTR_no_amplicons_males_merge$Num_sim_variants
chrX_nonPARs_no_XTR_no_amplicons_males_merge$FNtoSimVars <- chrX_nonPARs_no_XTR_no_amplicons_males_merge$FN / chrX_nonPARs_no_XTR_no_amplicons_males_merge$Num_sim_variants
chrX_nonPARs_no_XTR_no_amplicons_males_merge$TPtoSimVars <- chrX_nonPARs_no_XTR_no_amplicons_males_merge$TP / chrX_nonPARs_no_XTR_no_amplicons_males_merge$Num_sim_variants


## Chromosome X XTR  ##
chrX_XTR_males <- read.table(
  "../performance_metrics/EUR/males/by_region/EUR_males_chrX_nonPARs_haploid_golden_vs_called_performance_metrics_XTR.txt",
  header = T
)

# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_XTR_males_merge <- merge(chrX_XTR_males, chrX_XTR_sim_males, by = "Sample")

chrX_XTR_males_merge$Chromosome <- "chrX_XTR"
chrX_XTR_males_merge$FPtoTP <- chrX_XTR_males_merge$FP / chrX_XTR_males_merge$TP
chrX_XTR_males_merge$FNtoTP <- chrX_XTR_males_merge$FN / chrX_XTR_males_merge$TP
chrX_XTR_males_merge$FPtoSimVars <- chrX_XTR_males_merge$FP / chrX_XTR_males_merge$Num_sim_variants
chrX_XTR_males_merge$FNtoSimVars <- chrX_XTR_males_merge$FN / chrX_XTR_males_merge$Num_sim_variants
chrX_XTR_males_merge$TPtoSimVars <- chrX_XTR_males_merge$TP / chrX_XTR_males_merge$Num_sim_variants

## Chromosome X Amplicons ##
chrX_Amplicons_males <- read.table(
  "../performance_metrics/EUR/males/by_region/EUR_males_chrX_nonPARs_haploid_golden_vs_called_performance_metrics_Amplicons.txt",
  header = T
)

# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_Amplicons_males_merge <- merge(chrX_Amplicons_males, chrX_Amplicons_sim_males, by = "Sample")

chrX_Amplicons_males_merge$Chromosome <- "chrX_Amplicons"
chrX_Amplicons_males_merge$FPtoTP <- chrX_Amplicons_males_merge$FP / chrX_Amplicons_males_merge$TP
chrX_Amplicons_males_merge$FNtoTP <- chrX_Amplicons_males_merge$FN / chrX_Amplicons_males_merge$TP
chrX_Amplicons_males_merge$FPtoSimVars <- chrX_Amplicons_males_merge$FP / chrX_Amplicons_males_merge$Num_sim_variants
chrX_Amplicons_males_merge$FNtoSimVars <- chrX_Amplicons_males_merge$FN / chrX_Amplicons_males_merge$Num_sim_variants
chrX_Amplicons_males_merge$TPtoSimVars <- chrX_Amplicons_males_merge$TP / chrX_Amplicons_males_merge$Num_sim_variants


## Chromosome X nonPARs minus XTR minus amplicons ##
chrY_nonPARs_no_XTR_no_amplicons_males <- read.table(
  "../performance_metrics/EUR/males/by_region/EUR_males_chrY_haploid_golden_vs_called_performance_metrics_nonPARsminusXTRminusAmplicons.txt",
  header = T
)

# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrY_nonPARs_no_XTR_no_amplicons_males_merge <- merge(chrY_nonPARs_no_XTR_no_amplicons_males, 
                                                      chrY_nonPARs_no_XTR_no_amplicons_sim_males, by = "Sample")

chrY_nonPARs_no_XTR_no_amplicons_males_merge$Chromosome <- "chrY_nonPARs_no_XTR_no_amplicons"
chrY_nonPARs_no_XTR_no_amplicons_males_merge$FPtoTP <- chrY_nonPARs_no_XTR_no_amplicons_males_merge$FP / chrY_nonPARs_no_XTR_no_amplicons_males_merge$TP
chrY_nonPARs_no_XTR_no_amplicons_males_merge$FNtoTP <- chrY_nonPARs_no_XTR_no_amplicons_males_merge$FN / chrY_nonPARs_no_XTR_no_amplicons_males_merge$TP
chrY_nonPARs_no_XTR_no_amplicons_males_merge$FPtoSimVars <- chrY_nonPARs_no_XTR_no_amplicons_males_merge$FP / chrY_nonPARs_no_XTR_no_amplicons_males_merge$Num_sim_variants
chrY_nonPARs_no_XTR_no_amplicons_males_merge$FNtoSimVars <- chrY_nonPARs_no_XTR_no_amplicons_males_merge$FN / chrY_nonPARs_no_XTR_no_amplicons_males_merge$Num_sim_variants
chrY_nonPARs_no_XTR_no_amplicons_males_merge$TPtoSimVars <- chrY_nonPARs_no_XTR_no_amplicons_males_merge$TP / chrY_nonPARs_no_XTR_no_amplicons_males_merge$Num_sim_variants


## Chromosome X XTR  ##
chrY_XTR_males <- read.table(
  "../performance_metrics/EUR/males/by_region/EUR_males_chrY_haploid_golden_vs_called_performance_metrics_XTR.txt",
  header = T
)

# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrY_XTR_males_merge <- merge(chrY_XTR_males, chrY_XTR_sim_males, by = "Sample")

chrY_XTR_males_merge$Chromosome <- "chrY_XTR"
chrY_XTR_males_merge$FPtoTP <- chrY_XTR_males_merge$FP / chrY_XTR_males_merge$TP
chrY_XTR_males_merge$FNtoTP <- chrY_XTR_males_merge$FN / chrY_XTR_males_merge$TP
chrY_XTR_males_merge$FPtoSimVars <- chrY_XTR_males_merge$FP / chrY_XTR_males_merge$Num_sim_variants
chrY_XTR_males_merge$FNtoSimVars <- chrY_XTR_males_merge$FN / chrY_XTR_males_merge$Num_sim_variants
chrY_XTR_males_merge$TPtoSimVars <- chrY_XTR_males_merge$TP / chrY_XTR_males_merge$Num_sim_variants

## Chromosome X Amplicons ##
chrY_Amplicons_males <- read.table(
  "../performance_metrics/EUR/males/by_region/EUR_males_chrY_haploid_golden_vs_called_performance_metrics_Amplicons.txt",
  header = T
)

# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrY_Amplicons_males_merge <- merge(chrY_Amplicons_males, chrY_Amplicons_sim_males, by = "Sample")

chrY_Amplicons_males_merge$Chromosome <- "chrY_Amplicons"
chrY_Amplicons_males_merge$FPtoTP <- chrY_Amplicons_males_merge$FP / chrY_Amplicons_males_merge$TP
chrY_Amplicons_males_merge$FNtoTP <- chrY_Amplicons_males_merge$FN / chrY_Amplicons_males_merge$TP
chrY_Amplicons_males_merge$FPtoSimVars <- chrY_Amplicons_males_merge$FP / chrY_Amplicons_males_merge$Num_sim_variants
chrY_Amplicons_males_merge$FNtoSimVars <- chrY_Amplicons_males_merge$FN / chrY_Amplicons_males_merge$Num_sim_variants
chrY_Amplicons_males_merge$TPtoSimVars <- chrY_Amplicons_males_merge$TP / chrY_Amplicons_males_merge$Num_sim_variants

##                        ##
## Prep data for plotting ##
##                        ##
dat_X_females <- rbind(
  chrX_PARs_dat_merge, 
  chrX_nonPARs_no_XTR_no_amplicons_merge,
  chrX_XTR_merge,
  chrX_Amplicons_merge
  )


dat_X_males <- rbind(
  chrX_PARs_dat_males_merge, 
  chrX_nonPARs_no_XTR_no_amplicons_males_merge,
  chrX_XTR_males_merge,
  chrX_Amplicons_males_merge
)

dat_Y_males <- rbind(
  chrY_nonPARs_no_XTR_no_amplicons_males_merge,
  chrY_XTR_males_merge,
  chrY_Amplicons_males_merge
)

# Get means #
# females
mean(chrX_PARs_dat_merge$FPtoSimVars)
mean(chrX_PARs_dat_merge$FNtoSimVars)

mean(chrX_nonPARs_no_XTR_no_amplicons_merge$FPtoSimVars)
mean(chrX_nonPARs_no_XTR_no_amplicons_merge$FNtoSimVars)

mean(chrX_XTR_merge$FPtoSimVars)
mean(chrX_XTR_merge$FNtoSimVars)

mean(chrX_Amplicons_merge$FPtoSimVars)
mean(chrX_Amplicons_merge$FNtoSimVars)

# males
mean(chrX_PARs_dat_males_merge$FPtoSimVars)
mean(chrX_PARs_dat_males_merge$FNtoSimVars)

mean(chrX_nonPARs_no_XTR_no_amplicons_males_merge$FPtoSimVars)
mean(chrX_nonPARs_no_XTR_no_amplicons_males_merge$FNtoSimVars)

mean(chrX_XTR_males_merge$FPtoSimVars)
mean(chrX_XTR_males_merge$FNtoSimVars)

mean(chrX_Amplicons_males_merge$FPtoSimVars)
mean(chrX_Amplicons_males_merge$FNtoSimVars)

mean(chrY_nonPARs_no_XTR_no_amplicons_males_merge$FPtoSimVars)
mean(chrY_nonPARs_no_XTR_no_amplicons_males_merge$FNtoSimVars)

mean(chrY_XTR_males_merge$FPtoSimVars)
mean(chrY_XTR_males_merge$FNtoSimVars)

mean(chrY_Amplicons_males_merge$FPtoSimVars)
mean(chrY_Amplicons_males_merge$FNtoSimVars)

# means continued
mean(chrX_PARs_dat_merge$FP)
mean(chrX_PARs_dat_merge$FN)
mean(chrX_PARs_dat_merge$Num_sim_variants)

mean(chrX_nonPARs_no_XTR_no_amplicons_merge$FP)
mean(chrX_nonPARs_no_XTR_no_amplicons_merge$FN)
mean(chrX_nonPARs_no_XTR_no_amplicons_merge$Num_sim_variants)

mean(chrX_XTR_merge$FP)
mean(chrX_XTR_merge$FN)
mean(chrX_XTR_merge$Num_sim_variants)

mean(chrX_Amplicons_merge$FP)
mean(chrX_Amplicons_merge$FN)
mean(chrX_Amplicons_merge$Num_sim_variants)

# males
mean(chrX_PARs_dat_males_merge$FP)
mean(chrX_PARs_dat_males_merge$FN)
mean(chrX_PARs_dat_males_merge$Num_sim_variants)

mean(chrX_nonPARs_no_XTR_no_amplicons_males_merge$FP)
mean(chrX_nonPARs_no_XTR_no_amplicons_males_merge$FN)
mean(chrX_nonPARs_no_XTR_no_amplicons_males_merge$Num_sim_variants)

mean(chrX_XTR_males_merge$FP)
mean(chrX_XTR_males_merge$FN)
mean(chrX_XTR_males_merge$Num_sim_variants)

mean(chrX_Amplicons_males_merge$FP)
mean(chrX_Amplicons_males_merge$FN)
mean(chrX_Amplicons_males_merge$Num_sim_variants)

mean(chrY_nonPARs_no_XTR_no_amplicons_males_merge$FP)
mean(chrY_nonPARs_no_XTR_no_amplicons_males_merge$FN)
mean(chrY_nonPARs_no_XTR_no_amplicons_males_merge$Num_sim_variants)

mean(chrY_XTR_males_merge$FP)
mean(chrY_XTR_males_merge$FN)
mean(chrY_XTR_males_merge$Num_sim_variants)

mean(chrY_Amplicons_males_merge$FP)
mean(chrY_Amplicons_males_merge$FN)
mean(chrY_Amplicons_males_merge$Num_sim_variants)


##      ##
## Plot ##
##      ##

FPtoSimVars_females <- dat_X_females %>%
  arrange(FPtoSimVars) %>%
  mutate(name = factor(Chromosome, 
                       levels=c(
                                "chrX_PARs",
                                "chrX_nonPARs_no_XTR_no_amplicons", 
                                "chrX_XTR", 
                                "chrX_Amplicons"
                       ))) %>%
  ggplot( aes(x=name, y=FPtoSimVars, color=name)) +
  ylim(0,max(c(dat_X_males$FPtoSimVars, dat_X_females$FPtoSimVars, dat_Y_males$FPtoSimVars))) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  scale_color_manual(values = c("darkolivegreen4", "#B2DF8A", "steelblue", "darkmagenta")) +
  ylab("False Positives/Total # Simulated") +
  #scale_color_manual(values = c("#A6CEE3", "#FB9A99", "steelblue", "#B2DF8A", "darkmagenta",
  #                              "darkolivegreen4", "darkolivegreen4")) +
  #scale_color_manual(values = c("#A6CEE3", "#FB9A99", "darkolivegreen4", "#B2DF8A", "#2ECC71", "darkmagenta")) +
  theme(legend.position = "none") + 
  xlab("Chromosome X") +
  scale_x_discrete(labels=c( 
                            "chrX_PARs" = "PARs",
                            "chrX_nonPARs_no_XTR_no_amplicons" = expression("X"["non-PARs"]),
                            "chrX_XTR" = expression("X"["XTR"]),
                            "chrX_Amplicons" = expression("X"["Amplicons"])
                            )) +
  #ggtitle("False Positives/True Positives") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1))


FNtoSimVars_females <- dat_X_females %>%
  arrange(FNtoSimVars) %>%
  mutate(name = factor(Chromosome, 
                       levels=c(
                         "chrX_PARs",
                         "chrX_nonPARs_no_XTR_no_amplicons", 
                         "chrX_XTR", 
                         "chrX_Amplicons"
                       ))) %>%
  ggplot( aes(x=name, y=FNtoSimVars, color=name)) +
  ylim(0,1) +
  #ylim(0,max(c(dat_X_males$FNtoSimVars, dat_X_females$FNtoSimVars, dat_Y_males$FNtoSimVars))) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  scale_color_manual(values = c("darkolivegreen4", "#B2DF8A", "steelblue", "darkmagenta")) +
  ylab("False Negatives/Total # Simulated") +
  theme(legend.position = "none") + 
  xlab("Chromosome X") +
  scale_x_discrete(labels=c( 
    "chrX_PARs" = "PARs",
    "chrX_nonPARs_no_XTR_no_amplicons" = expression("X"["non-PARs"]),
    "chrX_XTR" = expression("X"["XTR"]),
    "chrX_Amplicons" = expression("X"["Amplicons"])
  )) +
  #ggtitle("False Negatives/True Positives") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1))



FPtoSimVars_males_X <- dat_X_males %>%
  arrange(FPtoSimVars) %>%
  mutate(name = factor(Chromosome, 
                       levels=c(
                         "chrX_PARs",
                         "chrX_nonPARs_no_XTR_no_amplicons", 
                         "chrX_XTR", 
                         "chrX_Amplicons"
                       ))) %>%
  ggplot( aes(x=name, y=FPtoSimVars, color=name)) +
  ylim(0,max(c(dat_X_males$FPtoSimVars, dat_X_females$FPtoSimVars, dat_Y_males$FPtoSimVars))) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  scale_color_manual(values = c("darkolivegreen4", "#B2DF8A", "steelblue", "darkmagenta")) +
  ylab("False Positives/Total # Simulated") +
  #scale_color_manual(values = c("#A6CEE3", "#FB9A99", "steelblue", "#B2DF8A", "darkmagenta",
  #                              "darkolivegreen4", "darkolivegreen4")) +
  #scale_color_manual(values = c("#A6CEE3", "#FB9A99", "darkolivegreen4", "#B2DF8A", "#2ECC71", "darkmagenta")) +
  theme(legend.position = "none") + 
  xlab("Chromosome X") +
  scale_x_discrete(labels=c( 
    "chrX_PARs" = "PARs",
    "chrX_nonPARs_no_XTR_no_amplicons" = expression("X"["non-PARs"]),
    "chrX_XTR" = expression("X"["XTR"]),
    "chrX_Amplicons" = expression("X"["Amplicons"])
  )) +
  #ggtitle("False Positives/True Positives") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1))


FNtoSimVars_males_X <- dat_X_males %>%
  arrange(FNtoSimVars) %>%
  mutate(name = factor(Chromosome, 
                       levels=c(
                         "chrX_PARs",
                         "chrX_nonPARs_no_XTR_no_amplicons", 
                         "chrX_XTR", 
                         "chrX_Amplicons"
                       ))) %>%
  ggplot( aes(x=name, y=FNtoSimVars, color=name)) +
  ylim(0,1) +
  #ylim(0,max(c(dat_X_males$FNtoSimVars, dat_X_females$FNtoSimVars, dat_Y_males$FNtoSimVars))) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  scale_color_manual(values = c("darkolivegreen4", "#B2DF8A", "steelblue", "darkmagenta")) +
  ylab("False Negatives/Total # Simulated") +
  theme(legend.position = "none") + 
  xlab("Chromosome X") +
  scale_x_discrete(labels=c( 
    "chrX_PARs" = "PARs",
    "chrX_nonPARs_no_XTR_no_amplicons" = expression("X"["non-PARs"]),
    "chrX_XTR" = expression("X"["XTR"]),
    "chrX_Amplicons" = expression("X"["Amplicons"])
  )) +
  #ggtitle("False Negatives/True Positives") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1))


FPtoSimVars_males_Y <- dat_Y_males %>%
  arrange(FPtoSimVars) %>%
  mutate(name = factor(Chromosome, 
                       levels=c(
                         "chrY_nonPARs_no_XTR_no_amplicons", 
                         "chrY_XTR", 
                         "chrY_Amplicons"
                       ))) %>%
  ggplot( aes(x=name, y=FPtoSimVars, color=name)) +
  #ylim(0,0.075) +
  ylim(0,max(c(dat_X_males$FPtoSimVars, dat_X_females$FPtoSimVars, dat_Y_males$FPtoSimVars))) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  scale_color_manual(values = c("#FF7F00", "steelblue", "darkmagenta")) +
  ylab("False Positives/Total # Simulated") +
  #scale_color_manual(values = c("#A6CEE3", "#FB9A99", "steelblue", "#B2DF8A", "darkmagenta",
  #                              "darkolivegreen4", "darkolivegreen4")) +
  #scale_color_manual(values = c("#A6CEE3", "#FB9A99", "darkolivegreen4", "#B2DF8A", "#2ECC71", "darkmagenta")) +
  theme(legend.position = "none") + 
  xlab("Chromosome Y") +
  scale_x_discrete(labels=c( 
    "chrY_nonPARs_no_XTR_no_amplicons" = expression("Y"["non-PARs"]),
    "chrY_XTR" = expression("Y"["XTR"]),
    "chrY_Amplicons" = expression("Y"["Amplicons"])
  )) +
  #ggtitle("False Positives/True Positives") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1))


FNtoSimVars_males_Y <- dat_Y_males %>%
  arrange(FNtoSimVars) %>%
  mutate(name = factor(Chromosome, 
                       levels=c(
                         "chrY_nonPARs_no_XTR_no_amplicons", 
                         "chrY_XTR", 
                         "chrY_Amplicons"
                       ))) %>%
  ggplot( aes(x=name, y=FNtoSimVars, color=name)) +
  ylim(0,1) +
  #ylim(0,max(c(dat_X_males$FNtoSimVars, dat_X_females$FNtoSimVars, dat_Y_males$FNtoSimVars))) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  scale_color_manual(values = c("#FF7F00", "steelblue", "darkmagenta")) +
  ylab("False Negatives/Total # Simulated") +
  theme(legend.position = "none") + 
  xlab("Chromosome Y") +
  scale_x_discrete(labels=c( 
    "chrY_PARs" = "PARs",
    "chrY_nonPARs_no_XTR_no_amplicons" = expression("Y"["non-PARs"]),
    "chrY_XTR" = expression("Y"["XTR"]),
    "chrY_Amplicons" = expression("Y"["Amplicons"])
  )) +
  #ggtitle("False Negatives/True Positives") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1))



##              ##
## Output plots ##
##              ##

plot_X_females <- ggarrange(
  FPtoSimVars_females, FNtoSimVars_females,
  ncol = 1,
  nrow = 2,
  labels = c("A", "D")
)
plot_X_females_annoated <- annotate_figure(plot_X_females, top = text_grob("Chromosome X", face = "bold", size = 14))
#print(plot_X_females_annoated)


plot_X_males <- ggarrange(
  FPtoSimVars_males_X, FNtoSimVars_males_X,
  ncol = 1,
  nrow = 2,
  labels = c("B", "E")
)
plot_X_males_annoated <- annotate_figure(plot_X_males, top = text_grob("Chromosome X", face = "bold", size = 14))
#print(plot_X_males_annoated)


plot_Y_males <- ggarrange(
  FPtoSimVars_males_Y, FNtoSimVars_males_Y,
  ncol = 1,
  nrow = 2,
  labels = c("C", "F")
)
plot_Y_males_annoated <- annotate_figure(plot_Y_males, top = text_grob("Chromosome Y", face = "bold", size = 14))
#print(plot_Y_males_annoated)


pdf("../plots/SCC_jitter_females_males_by_region_02.pdf",
    width = 10, height = 8)

ggarrange(
  plot_X_females_annoated, plot_X_males_annoated, plot_Y_males_annoated,
  ncol = 3,
  nrow = 1
)

dev.off()

pdf("../plots/SCC_jitter_females_males_by_region_same_axis_as_default.pdf",
    width = 10, height = 8)

ggarrange(
  plot_X_females_annoated, plot_X_males_annoated, plot_Y_males_annoated,
  ncol = 3,
  nrow = 1
)

dev.off()




###############################################################################
###############################################################################
# Complex vs not-complex
###############################################################################
###############################################################################
# Read in number of simulated variants for the same regions

##                          ##
## Read in data - simulated ##
##                          ##

# FEMALES #
chrX_PARs_complex_sim <- read.table("../simulated_variants/EUR_females_chrX_PARs_diploid_simulated_counts_complex.txt",
                                    header = T)

chrX_PARs_not_complex_sim <- read.table("../simulated_variants/EUR_females_chrX_PARs_diploid_simulated_counts_not_complex.txt",
                                        header = T)

chrX_complex_sim <- read.table("../simulated_variants/EUR_females_chrX_nonPARs_diploid_simulated_counts_complex.txt",
                               header = T)

chrX_not_complex_sim <- read.table("../simulated_variants/EUR_females_chrX_nonPARs_diploid_simulated_counts_not_complex.txt",
                                   header = T)


# MALES #
chrX_PARs_complex_sim_males <- read.table("../simulated_variants/EUR_males_chrX_PARs_diploid_simulated_counts_complex.txt",
                                          header = T)

chrX_PARs_not_complex_sim_males <- read.table("../simulated_variants/EUR_males_chrX_PARs_diploid_simulated_counts_not_complex.txt",
                                              header = T)

chrX_complex_sim_males <- read.table("../simulated_variants/EUR_males_chrX_nonPARs_haploid_simulated_counts_complex.txt",
                                     header = T)

chrX_not_complex_sim_males <- read.table("../simulated_variants/EUR_males_chrX_nonPARs_haploid_simulated_counts_not_complex.txt",
                                         header = T)

chrY_complex_sim_males <- read.table("../simulated_variants/EUR_males_chrY_haploid_simulated_counts_complex.txt",
                                     header = T)

chrY_not_complex_sim_males <- read.table("../simulated_variants/EUR_males_chrY_haploid_simulated_counts_not_complex.txt",
                                         header = T)


##              ##
## Read in data ##
##              ##
# FEMALES #
## Chromosome X PARs - complex ##
chrX_PARs_complex_dat <- read.table(
  "../performance_metrics/EUR/females/by_region/EUR_females_chrX_PARs_diploid_golden_vs_called_performance_metrics_complex.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrX_PARs_complex_dat <- chrX_PARs_complex_dat[,!(names(chrX_PARs_complex_dat) %in% drop)]


# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_PARs_complex_dat_merge <- merge(chrX_PARs_complex_dat, chrX_PARs_complex_sim, by = "Sample")

chrX_PARs_complex_dat_merge$Chromosome <- "chrX_PARs_complex"
chrX_PARs_complex_dat_merge$FPtoTP <- chrX_PARs_complex_dat_merge$FP / chrX_PARs_complex_dat_merge$TP
chrX_PARs_complex_dat_merge$FNtoTP <- chrX_PARs_complex_dat_merge$FN / chrX_PARs_complex_dat_merge$TP
chrX_PARs_complex_dat_merge$FPtoSimVars <- chrX_PARs_complex_dat_merge$FP / chrX_PARs_complex_dat_merge$Num_sim_variants
chrX_PARs_complex_dat_merge$FNtoSimVars <- chrX_PARs_complex_dat_merge$FN / chrX_PARs_complex_dat_merge$Num_sim_variants
chrX_PARs_complex_dat_merge$TPtoSimVars <- chrX_PARs_complex_dat_merge$TP / chrX_PARs_complex_dat_merge$Num_sim_variants

## Chromosome X PARs - complex ##
chrX_PARs_not_complex_dat <- read.table(
  "../performance_metrics/EUR/females/by_region/EUR_females_chrX_PARs_diploid_golden_vs_called_performance_metrics_not_complex.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrX_PARs_not_complex_dat <- chrX_PARs_not_complex_dat[,!(names(chrX_PARs_not_complex_dat) %in% drop)]


# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_PARs_not_complex_dat_merge <- merge(chrX_PARs_not_complex_dat, chrX_PARs_not_complex_sim, by = "Sample")

chrX_PARs_not_complex_dat_merge$Chromosome <- "chrX_PARs_not_complex"
chrX_PARs_not_complex_dat_merge$FPtoTP <- chrX_PARs_not_complex_dat_merge$FP / chrX_PARs_not_complex_dat_merge$TP
chrX_PARs_not_complex_dat_merge$FNtoTP <- chrX_PARs_not_complex_dat_merge$FN / chrX_PARs_not_complex_dat_merge$TP
chrX_PARs_not_complex_dat_merge$FPtoSimVars <- chrX_PARs_not_complex_dat_merge$FP / chrX_PARs_not_complex_dat_merge$Num_sim_variants
chrX_PARs_not_complex_dat_merge$FNtoSimVars <- chrX_PARs_not_complex_dat_merge$FN / chrX_PARs_not_complex_dat_merge$Num_sim_variants
chrX_PARs_not_complex_dat_merge$TPtoSimVars <- chrX_PARs_not_complex_dat_merge$TP / chrX_PARs_not_complex_dat_merge$Num_sim_variants


## Chromosome X complex  ##
chrX_complex <- read.table(
  "../performance_metrics/EUR/females/by_region/EUR_females_chrX_nonPARs_diploid_golden_vs_called_performance_metrics_complex.txt",
  header = T
)

# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_complex_merge <- merge(chrX_complex, chrX_complex_sim, by = "Sample")

chrX_complex_merge$Chromosome <- "chrX_complex"
chrX_complex_merge$FPtoTP <- chrX_complex_merge$FP / chrX_complex_merge$TP
chrX_complex_merge$FNtoTP <- chrX_complex_merge$FN / chrX_complex_merge$TP
chrX_complex_merge$FPtoSimVars <- chrX_complex_merge$FP / chrX_complex_merge$Num_sim_variants
chrX_complex_merge$FNtoSimVars <- chrX_complex_merge$FN / chrX_complex_merge$Num_sim_variants
chrX_complex_merge$TPtoSimVars <- chrX_complex_merge$TP / chrX_complex_merge$Num_sim_variants

## Chromosome X not_complex ##
chrX_not_complex <- read.table(
  "../performance_metrics/EUR/females/by_region/EUR_females_chrX_nonPARs_diploid_golden_vs_called_performance_metrics_not_complex.txt",
  header = T
)

# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_not_complex_merge <- merge(chrX_not_complex, chrX_not_complex_sim, by = "Sample")

chrX_not_complex_merge$Chromosome <- "chrX_not_complex"
chrX_not_complex_merge$FPtoTP <- chrX_not_complex_merge$FP / chrX_not_complex_merge$TP
chrX_not_complex_merge$FNtoTP <- chrX_not_complex_merge$FN / chrX_not_complex_merge$TP
chrX_not_complex_merge$FPtoSimVars <- chrX_not_complex_merge$FP / chrX_not_complex_merge$Num_sim_variants
chrX_not_complex_merge$FNtoSimVars <- chrX_not_complex_merge$FN / chrX_not_complex_merge$Num_sim_variants
chrX_not_complex_merge$TPtoSimVars <- chrX_not_complex_merge$TP / chrX_not_complex_merge$Num_sim_variants


# MALES #
## Chromosome X PARs - complex ##
chrX_PARs_complex_dat_males <- read.table(
  "../performance_metrics/EUR/males/by_region/EUR_males_chrX_PARs_diploid_golden_vs_called_performance_metrics_complex.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrX_PARs_complex_dat_males <- chrX_PARs_complex_dat_males[,!(names(chrX_PARs_complex_dat_males) %in% drop)]


# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_PARs_complex_dat_males_merge <- merge(chrX_PARs_complex_dat_males, chrX_PARs_complex_sim_males, by = "Sample")

chrX_PARs_complex_dat_males_merge$Chromosome <- "chrX_PARs_complex"
chrX_PARs_complex_dat_males_merge$FPtoTP <- chrX_PARs_complex_dat_males_merge$FP / chrX_PARs_complex_dat_males_merge$TP
chrX_PARs_complex_dat_males_merge$FNtoTP <- chrX_PARs_complex_dat_males_merge$FN / chrX_PARs_complex_dat_males_merge$TP
chrX_PARs_complex_dat_males_merge$FPtoSimVars <- chrX_PARs_complex_dat_males_merge$FP / chrX_PARs_complex_dat_males_merge$Num_sim_variants
chrX_PARs_complex_dat_males_merge$FNtoSimVars <- chrX_PARs_complex_dat_males_merge$FN / chrX_PARs_complex_dat_males_merge$Num_sim_variants
chrX_PARs_complex_dat_males_merge$TPtoSimVars <- chrX_PARs_complex_dat_males_merge$TP / chrX_PARs_complex_dat_males_merge$Num_sim_variants

## Chromosome X PARs - complex ##
chrX_PARs_not_complex_dat_males <- read.table(
  "../performance_metrics/EUR/males/by_region/EUR_males_chrX_PARs_diploid_golden_vs_called_performance_metrics_not_complex.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrX_PARs_not_complex_dat_males <- chrX_PARs_not_complex_dat_males[,!(names(chrX_PARs_not_complex_dat_males) %in% drop)]


# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_PARs_not_complex_dat_males_merge <- merge(chrX_PARs_not_complex_dat_males, chrX_PARs_not_complex_sim_males, by = "Sample")

chrX_PARs_not_complex_dat_males_merge$Chromosome <- "chrX_PARs_not_complex"
chrX_PARs_not_complex_dat_males_merge$FPtoTP <- chrX_PARs_not_complex_dat_males_merge$FP / chrX_PARs_not_complex_dat_males_merge$TP
chrX_PARs_not_complex_dat_males_merge$FNtoTP <- chrX_PARs_not_complex_dat_males_merge$FN / chrX_PARs_not_complex_dat_males_merge$TP
chrX_PARs_not_complex_dat_males_merge$FPtoSimVars <- chrX_PARs_not_complex_dat_males_merge$FP / chrX_PARs_not_complex_dat_males_merge$Num_sim_variants
chrX_PARs_not_complex_dat_males_merge$FNtoSimVars <- chrX_PARs_not_complex_dat_males_merge$FN / chrX_PARs_not_complex_dat_males_merge$Num_sim_variants
chrX_PARs_not_complex_dat_males_merge$TPtoSimVars <- chrX_PARs_not_complex_dat_males_merge$TP / chrX_PARs_not_complex_dat_males_merge$Num_sim_variants


## Chromosome X complex  ##
chrX_complex_males <- read.table(
  "../performance_metrics/EUR/males/by_region/EUR_males_chrX_nonPARs_haploid_golden_vs_called_performance_metrics_complex.txt",
  header = T
)

# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_complex_males_merge <- merge(chrX_complex_males, chrX_complex_sim_males, by = "Sample")

chrX_complex_males_merge$Chromosome <- "chrX_complex"
chrX_complex_males_merge$FPtoTP <- chrX_complex_males_merge$FP / chrX_complex_males_merge$TP
chrX_complex_males_merge$FNtoTP <- chrX_complex_males_merge$FN / chrX_complex_males_merge$TP
chrX_complex_males_merge$FPtoSimVars <- chrX_complex_males_merge$FP / chrX_complex_males_merge$Num_sim_variants
chrX_complex_males_merge$FNtoSimVars <- chrX_complex_males_merge$FN / chrX_complex_males_merge$Num_sim_variants
chrX_complex_males_merge$TPtoSimVars <- chrX_complex_males_merge$TP / chrX_complex_males_merge$Num_sim_variants

## Chromosome X not_complex ##
chrX_not_complex_males <- read.table(
  "../performance_metrics/EUR/males/by_region/EUR_males_chrX_nonPARs_haploid_golden_vs_called_performance_metrics_not_complex.txt",
  header = T
)

# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_not_complex_males_merge <- merge(chrX_not_complex_males, chrX_not_complex_sim_males, by = "Sample")

chrX_not_complex_males_merge$Chromosome <- "chrX_not_complex"
chrX_not_complex_males_merge$FPtoTP <- chrX_not_complex_males_merge$FP / chrX_not_complex_males_merge$TP
chrX_not_complex_males_merge$FNtoTP <- chrX_not_complex_males_merge$FN / chrX_not_complex_males_merge$TP
chrX_not_complex_males_merge$FPtoSimVars <- chrX_not_complex_males_merge$FP / chrX_not_complex_males_merge$Num_sim_variants
chrX_not_complex_males_merge$FNtoSimVars <- chrX_not_complex_males_merge$FN / chrX_not_complex_males_merge$Num_sim_variants
chrX_not_complex_males_merge$TPtoSimVars <- chrX_not_complex_males_merge$TP / chrX_not_complex_males_merge$Num_sim_variants


## Chromosome X complex  ##
chrY_complex_males <- read.table(
  "../performance_metrics/EUR/males/by_region/EUR_males_chrY_haploid_golden_vs_called_performance_metrics_complex.txt",
  header = T
)

# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrY_complex_males_merge <- merge(chrY_complex_males, chrY_complex_sim_males, by = "Sample")

chrY_complex_males_merge$Chromosome <- "chrY_complex"
chrY_complex_males_merge$FPtoTP <- chrY_complex_males_merge$FP / chrY_complex_males_merge$TP
chrY_complex_males_merge$FNtoTP <- chrY_complex_males_merge$FN / chrY_complex_males_merge$TP
chrY_complex_males_merge$FPtoSimVars <- chrY_complex_males_merge$FP / chrY_complex_males_merge$Num_sim_variants
chrY_complex_males_merge$FNtoSimVars <- chrY_complex_males_merge$FN / chrY_complex_males_merge$Num_sim_variants
chrY_complex_males_merge$TPtoSimVars <- chrY_complex_males_merge$TP / chrY_complex_males_merge$Num_sim_variants

## Chromosome X not_complex ##
chrY_not_complex_males <- read.table(
  "../performance_metrics/EUR/males/by_region/EUR_males_chrY_haploid_golden_vs_called_performance_metrics_not_complex.txt",
  header = T
)

# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrY_not_complex_males_merge <- merge(chrY_not_complex_males, chrY_not_complex_sim_males, by = "Sample")

chrY_not_complex_males_merge$Chromosome <- "chrY_not_complex"
chrY_not_complex_males_merge$FPtoTP <- chrY_not_complex_males_merge$FP / chrY_not_complex_males_merge$TP
chrY_not_complex_males_merge$FNtoTP <- chrY_not_complex_males_merge$FN / chrY_not_complex_males_merge$TP
chrY_not_complex_males_merge$FPtoSimVars <- chrY_not_complex_males_merge$FP / chrY_not_complex_males_merge$Num_sim_variants
chrY_not_complex_males_merge$FNtoSimVars <- chrY_not_complex_males_merge$FN / chrY_not_complex_males_merge$Num_sim_variants
chrY_not_complex_males_merge$TPtoSimVars <- chrY_not_complex_males_merge$TP / chrY_not_complex_males_merge$Num_sim_variants

##                        ##
## Prep data for plotting ##
##                        ##
dat_X_females <- rbind(
  #chrX_PARs_complex_dat_merge, 
  #chrX_PARs_not_complex_dat_merge, 
  chrX_complex_merge,
  chrX_not_complex_merge
)


dat_X_males <- rbind(
  #chrX_PARs_complex_dat_males_merge, 
  #chrX_PARs_not_complex_dat_males_merge, 
  chrX_complex_males_merge,
  chrX_not_complex_males_merge
)

dat_Y_males <- rbind(
  chrY_complex_males_merge,
  chrY_not_complex_males_merge
)


##      ##
## Plot ##
##      ##

FPtoSimVars_females <- dat_X_females %>%
  arrange(FPtoSimVars) %>%
  mutate(name = factor(Chromosome, 
                       levels=c(
                         #"chrX_PARs_complex",
                         #"chrX_PARs_not_complex",
                         "chrX_complex", 
                         "chrX_not_complex"
                       ))) %>%
  ggplot( aes(x=name, y=FPtoSimVars, color=name)) +
  #ylim(0,max(c(dat_X_males$FPtoSimVars, dat_X_females$FPtoSimVars))) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  #scale_color_manual(values = c("darkolivegreen4", "darkolivegreen4", "#B2DF8A", "#B2DF8A")) +
  scale_color_manual(values = c("black","black")) +
  ylab("False Positives/Total # Simulated") +
  #scale_color_manual(values = c("#A6CEE3", "#FB9A99", "steelblue", "#B2DF8A", "darkmagenta",
  #                              "darkolivegreen4", "darkolivegreen4")) +
  #scale_color_manual(values = c("#A6CEE3", "#FB9A99", "darkolivegreen4", "#B2DF8A", "#2ECC71", "darkmagenta")) +
  theme(legend.position = "none") + 
  xlab("Chromosome X non-PARs") +
  scale_x_discrete(labels=c( 
    #"chrX_PARs_complex" = expression("PARs"["complex"]),
    #"chrX_PARs_not_complex" = expression("PARs"["not complex"]),
    "chrX_complex" = expression("X"["complex"]),
    "chrX_not_complex" = expression("X"["not_complex"])
  )) +
  #ggtitle("False Positives/True Positives") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1))


FNtoSimVars_females <- dat_X_females %>%
  arrange(FNtoSimVars) %>%
  mutate(name = factor(Chromosome, 
                       levels=c(
                         #"chrX_PARs_complex",
                         #"chrX_PARs_not_complex",
                         "chrX_complex", 
                         "chrX_not_complex"
                       ))) %>%
  ggplot( aes(x=name, y=FNtoSimVars, color=name)) +
  #ylim(0,max(c(dat_X_males$FNtoSimVars, dat_X_females$FNtoSimVars))) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  scale_color_manual(values = c("black","black")) +
  ylab("False Negatives/Total # Simulated") +
  theme(legend.position = "none") + 
  xlab("Chromosome X non-PARs") +
  scale_x_discrete(labels=c( 
    "chrX_PARs_complex" = expression("PARs"["complex"]),
    "chrX_PARs_not_complex" = expression("PARs"["not complex"]),
    "chrX_complex" = expression("X"["complex"]),
    "chrX_not_complex" = expression("X"["not_complex"])
  )) +
  #ggtitle("False Negatives/True Positives") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1))



FPtoSimVars_males_X <- dat_X_males %>%
  arrange(FPtoSimVars) %>%
  mutate(name = factor(Chromosome, 
                       levels=c(
                         #"chrX_PARs_complex",
                         #"chrX_PARs_not_complex",
                         "chrX_complex", 
                         "chrX_not_complex"
                       ))) %>%
  ggplot( aes(x=name, y=FPtoSimVars, color=name)) +
  #ylim(0,max(c(dat_X_males$FPtoSimVars, dat_X_females$FPtoSimVars))) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  scale_color_manual(values = c("black","black")) +
  ylab("False Positives/Total # Simulated") +
  theme(legend.position = "none") + 
  xlab("Chromosome X non-PARs") +
  scale_x_discrete(labels=c( 
    #"chrX_PARs_complex" = expression("PARs"["complex"]),
    #"chrX_PARs_not_complex" = expression("PARs"["not complex"]),
    "chrX_complex" = expression("X"["complex"]),
    "chrX_not_complex" = expression("X"["not_complex"])
  )) +
  #ggtitle("False Positives/True Positives") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1))


FNtoSimVars_males_X <- dat_X_males %>%
  arrange(FNtoSimVars) %>%
  mutate(name = factor(Chromosome, 
                       levels=c(
                         #"chrX_PARs_complex",
                         #"chrX_PARs_not_complex",
                         "chrX_complex", 
                         "chrX_not_complex"
                       ))) %>%
  ggplot( aes(x=name, y=FNtoSimVars, color=name)) +
  #ylim(0,max(c(dat_X_males$FNtoSimVars, dat_X_females$FNtoSimVars))) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  #scale_color_manual(values = c("#B2DF8A","#B2DF8A")) +
  scale_color_manual(values = c("black","black")) +
  ylab("False Negatives/Total # Simulated") +
  theme(legend.position = "none") + 
  xlab("Chromosome X non-PARs") +
  scale_x_discrete(labels=c( 
    #"chrX_PARs_complex" = expression("PARs"["complex"]),
    #"chrX_PARs_not_complex" = expression("PARs"["not complex"]),
    "chrX_complex" = expression("X"["complex"]),
    "chrX_not_complex" = expression("X"["not_complex"])
  )) +
  #ggtitle("False Negatives/True Positives") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1))


#FPtoSimVars_males_Y <- 
dat_Y_males %>%
  arrange(FPtoSimVars) %>%
  mutate(name = factor(Chromosome, 
                       levels=c(
                         "chrY_complex", 
                         "chrY_not_complex"
                       ))) %>%
  ggplot( aes(x=name, y=FPtoSimVars, color=name)) +
  #ylim(0,0.075) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  scale_color_manual(values = c("#FF7F00", "#FF7F00")) +
  ylab("False Positives/Total # Simulated") +
  theme(legend.position = "none") + 
  xlab("Chromosome Y") +
  scale_x_discrete(labels=c( 
    "chrY_complex" = expression("Y"["complex"]),
    "chrY_not_complex" = expression("Y"["not_complex"])
  )) +
  #ggtitle("False Positives/True Positives") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1))


#FNtoSimVars_males_Y <- 
dat_Y_males %>%
  arrange(FNtoSimVars) %>%
  mutate(name = factor(Chromosome, 
                       levels=c(
                         "chrY_complex", 
                         "chrY_not_complex"
                       ))) %>%
  ggplot( aes(x=name, y=FNtoSimVars, color=name)) +
  #ylim(0,0.075) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2)) +
  theme_classic() +
  scale_color_manual(values = c("#FF7F00", "#FF7F00")) +
  ylab("False Negatives/Total # Simulated") +
  theme(legend.position = "none") + 
  xlab("Chromosome Y") +
  scale_x_discrete(labels=c( 
    "chrY_PARs" = "PARs",
    "chrY_complex" = expression("Y"["complex"]),
    "chrY_not_complex" = expression("Y"["not_complex"])
  )) +
  #ggtitle("False Negatives/True Positives") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1))



##              ##
## Output plots ##
##              ##

plot_X_females <- ggarrange(
  FPtoSimVars_females, FNtoSimVars_females,
  ncol = 1,
  nrow = 2,
  labels = c("A", "C")
)
plot_X_females_annoated <- annotate_figure(plot_X_females, top = text_grob("Simulated females (XX)", face = "bold", size = 14))
#print(plot_X_females_annoated)


plot_X_males <- ggarrange(
  FPtoSimVars_males_X, FNtoSimVars_males_X,
  ncol = 1,
  nrow = 2,
  labels = c("B", "D")
)
plot_X_males_annoated <- annotate_figure(plot_X_males, top = text_grob("Simulated males (XY)", face = "bold", size = 14))
#print(plot_X_males_annoated)


#plot_Y_males <- ggarrange(
#  FPtoSimVars_males_Y, FNtoSimVars_males_Y,
#  ncol = 1,
#  nrow = 2,
#  labels = c("C", "F")
#)
#plot_Y_males_annoated <- annotate_figure(plot_Y_males, top = text_grob("Chromosome Y", face = "bold", size = 14))
#print(plot_Y_males_annoated)


pdf("../plots/SCC_jitter_females_males_by_region_complex_notcomplex.pdf",
    width = 10, height = 8)

ggarrange(
  plot_X_females_annoated, plot_X_males_annoated,
  ncol = 2,
  nrow = 1
)

dev.off()


