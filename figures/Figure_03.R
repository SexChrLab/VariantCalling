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
#brewer.pal(n = 8, name = "Paired")

# color names for chr 8, PARs, nonPARs diploid, nonPARs haploid, mtDNA diploid,
# mtDNA haploid, Y diploid, Y haploid
# "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00"
# "chr8_dip" = "#A6CEE3"
# "chrX_PARs_dip" = "#1F78B4"
# "chrX_nonPARs_dip" = "#B2DF8A"
# "chrX_nonPARs_hap" = "#33A02C"
# "chrY_dip" = "#FB9A99"
# "chrY_hap" = "#E31A1C"
# "chrM_dip" = "#FDBF6F"
# "chrM_hap" = "#FF7F00"


#-----------------------------------------------------------------------------#
## PREP DATA FOR PLOTTING
#-----------------------------------------------------------------------------#

#########
# Males #
#########
##                          ##
## Read in data - simulated ##
##                          ##
chrM_sim_males <- read.table("../simulated_variants/EUR_males_chrM_haploid_simulated_counts.txt",
                             header = T)

chrX_PARs_sim_males <- read.table("../simulated_variants/EUR_males_chrX_PARs_diploid_simulated_counts.txt",
                                  header = T)

chrX_nonPARs_sim_males <- read.table("../simulated_variants/EUR_males_chrX_nonPARs_haploid_simulated_counts.txt",
                                            header = T)


chrY_sim_males <- read.table("../simulated_variants/EUR_males_chrY_haploid_simulated_counts.txt",
                                    header = T)

##                    ##
## Read in data - SCC ##
##                    ##
dat_all_SCC_males <- c()

## Autosomes
for(i in 1:22) {
  print(paste("Processing chromosome", i))
  chri_sim <- read.table(paste("../simulated_variants/EUR_males_chr", i, 
                               "_diploid_simulated_counts.txt", sep = ""),
                         header = T)
  chri_dat_SCC_males <- read.table(
    paste("../performance_metrics/EUR/males/data/EUR_males_chr", i,
          "_diploid_golden_vs_called_performance_metrics.txt", 
          sep = ""),
    header = T
  )
  drop <- c("Validated_Genotypes")
  chri_dat_SCC_males <- chri_dat_SCC_males[,!(names(chri_dat_SCC_males) %in% drop)]
  
  # Merge with simulated counts and then do FPs and FNs over total simulated variants
  chri_dat_SCC_males_merge <- merge(chri_dat_SCC_males, chri_sim, by = "Sample")
  
  chri_dat_SCC_males_merge$Chromosome <- paste("chr", i, sep = "")
  chri_dat_SCC_males_merge$alignment <- "SCC"
  chri_dat_SCC_males_merge$FPtoTP <- chri_dat_SCC_males_merge$FP / chri_dat_SCC_males_merge$TP
  chri_dat_SCC_males_merge$FNtoTP <- chri_dat_SCC_males_merge$FN / chri_dat_SCC_males_merge$TP
  chri_dat_SCC_males_merge$FPtoSimVars <- chri_dat_SCC_males_merge$FP / chri_dat_SCC_males_merge$Num_sim_variants
  chri_dat_SCC_males_merge$FNtoSimVars <- chri_dat_SCC_males_merge$FN / chri_dat_SCC_males_merge$Num_sim_variants
  chri_dat_SCC_males_merge$TPtoSimVars <- chri_dat_SCC_males_merge$TP / chri_dat_SCC_males_merge$Num_sim_variants
  
  dat_all_SCC_males <- rbind(dat_all_SCC_males, chri_dat_SCC_males_merge)
}

males_default_max_FPtoSimVars <- max(dat_all_SCC_males$FPtoSimVars)
males_default_min_FPtoSimVars <- min(dat_all_SCC_males$FPtoSimVars)
males_default_max_FNtoSimVars <- max(dat_all_SCC_males$FNtoSimVars)
males_default_min_FNtoSimVars <- min(dat_all_SCC_males$FNtoSimVars)

## Chromosome M - haploid 
chrM_dat_SCC_males <- read.table(
  "../performance_metrics/EUR/males/data/EUR_males_chrM_haploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrM_dat_SCC_males <- chrM_dat_SCC_males[,!(names(chrM_dat_SCC_males) %in% drop)]

# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrM_dat_SCC_males_merge <- merge(chrM_dat_SCC_males, chrM_sim_males, by = "Sample")

chrM_dat_SCC_males_merge$Chromosome <- "chrM_hap"
chrM_dat_SCC_males_merge$alignment <- "SCC"
chrM_dat_SCC_males_merge$FPtoTP <- chrM_dat_SCC_males_merge$FP / chrM_dat_SCC_males_merge$TP
chrM_dat_SCC_males_merge$FNtoTP <- chrM_dat_SCC_males_merge$FN / chrM_dat_SCC_males_merge$TP
chrM_dat_SCC_males_merge$FPtoSimVars <- chrM_dat_SCC_males_merge$FP / chrM_dat_SCC_males_merge$Num_sim_variants
chrM_dat_SCC_males_merge$FNtoSimVars <- chrM_dat_SCC_males_merge$FN / chrM_dat_SCC_males_merge$Num_sim_variants
chrM_dat_SCC_males_merge$TPtoSimVars <- chrM_dat_SCC_males_merge$TP / chrM_dat_SCC_males_merge$Num_sim_variants

dat_all_SCC_males <- rbind(dat_all_SCC_males, chrM_dat_SCC_males_merge)

## Chromosome M - diploid 
chrM_dat_SCC_males <- read.table(
  "../performance_metrics/EUR/males/data/EUR_males_chrM_diploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrM_dat_SCC_males <- chrM_dat_SCC_males[,!(names(chrM_dat_SCC_males) %in% drop)]

# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrM_dat_SCC_males_merge <- merge(chrM_dat_SCC_males, chrM_sim_males, by = "Sample")

chrM_dat_SCC_males_merge$Chromosome <- "chrM_dip"
chrM_dat_SCC_males_merge$alignment <- "SCC"
chrM_dat_SCC_males_merge$FPtoTP <- chrM_dat_SCC_males_merge$FP / chrM_dat_SCC_males_merge$TP
chrM_dat_SCC_males_merge$FNtoTP <- chrM_dat_SCC_males_merge$FN / chrM_dat_SCC_males_merge$TP
chrM_dat_SCC_males_merge$FPtoSimVars <- chrM_dat_SCC_males_merge$FP / chrM_dat_SCC_males_merge$Num_sim_variants
chrM_dat_SCC_males_merge$FNtoSimVars <- chrM_dat_SCC_males_merge$FN / chrM_dat_SCC_males_merge$Num_sim_variants
chrM_dat_SCC_males_merge$TPtoSimVars <- chrM_dat_SCC_males_merge$TP / chrM_dat_SCC_males_merge$Num_sim_variants

dat_all_SCC_males <- rbind(dat_all_SCC_males, chrM_dat_SCC_males_merge)


## Chromosome X nonPARs - haploid
chrX_nonPARs_dat_SCC_males <- read.table(
  "../performance_metrics/EUR/males/data/EUR_males_chrX_nonPARs_haploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrX_nonPARs_dat_SCC_males <- chrX_nonPARs_dat_SCC_males[,!(names(chrX_nonPARs_dat_SCC_males) %in% drop)]

# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_nonPARs_dat_SCC_males_merge <- merge(chrX_nonPARs_dat_SCC_males, chrX_nonPARs_sim_males, by = "Sample")

chrX_nonPARs_dat_SCC_males_merge$Chromosome <- "chrX_nonPARs_hap"
chrX_nonPARs_dat_SCC_males_merge$alignment <- "SCC"
chrX_nonPARs_dat_SCC_males_merge$FPtoTP <- chrX_nonPARs_dat_SCC_males_merge$FP / chrX_nonPARs_dat_SCC_males_merge$TP
chrX_nonPARs_dat_SCC_males_merge$FNtoTP <- chrX_nonPARs_dat_SCC_males_merge$FN / chrX_nonPARs_dat_SCC_males_merge$TP
chrX_nonPARs_dat_SCC_males_merge$FPtoSimVars <- chrX_nonPARs_dat_SCC_males_merge$FP / chrX_nonPARs_dat_SCC_males_merge$Num_sim_variants
chrX_nonPARs_dat_SCC_males_merge$FNtoSimVars <- chrX_nonPARs_dat_SCC_males_merge$FN / chrX_nonPARs_dat_SCC_males_merge$Num_sim_variants
chrX_nonPARs_dat_SCC_males_merge$TPtoSimVars <- chrX_nonPARs_dat_SCC_males_merge$TP / chrX_nonPARs_dat_SCC_males_merge$Num_sim_variants

dat_all_SCC_males <- rbind(dat_all_SCC_males, chrX_nonPARs_dat_SCC_males_merge)


## Chromosome X nonPARs - diploid
chrX_nonPARs_dat_SCC_males <- read.table(
  "../performance_metrics/EUR/males/data/EUR_males_chrX_nonPARs_diploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrX_nonPARs_dat_SCC_males <- chrX_nonPARs_dat_SCC_males[,!(names(chrX_nonPARs_dat_SCC_males) %in% drop)]

# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_nonPARs_dat_SCC_males_merge <- merge(chrX_nonPARs_dat_SCC_males, chrX_nonPARs_sim_males, by = "Sample")

chrX_nonPARs_dat_SCC_males_merge$Chromosome <- "chrX_nonPARs_dip"
chrX_nonPARs_dat_SCC_males_merge$alignment <- "SCC"
chrX_nonPARs_dat_SCC_males_merge$FPtoTP <- chrX_nonPARs_dat_SCC_males_merge$FP / chrX_nonPARs_dat_SCC_males_merge$TP
chrX_nonPARs_dat_SCC_males_merge$FNtoTP <- chrX_nonPARs_dat_SCC_males_merge$FN / chrX_nonPARs_dat_SCC_males_merge$TP
chrX_nonPARs_dat_SCC_males_merge$FPtoSimVars <- chrX_nonPARs_dat_SCC_males_merge$FP / chrX_nonPARs_dat_SCC_males_merge$Num_sim_variants
chrX_nonPARs_dat_SCC_males_merge$FNtoSimVars <- chrX_nonPARs_dat_SCC_males_merge$FN / chrX_nonPARs_dat_SCC_males_merge$Num_sim_variants
chrX_nonPARs_dat_SCC_males_merge$TPtoSimVars <- chrX_nonPARs_dat_SCC_males_merge$TP / chrX_nonPARs_dat_SCC_males_merge$Num_sim_variants

dat_all_SCC_males <- rbind(dat_all_SCC_males, chrX_nonPARs_dat_SCC_males_merge)


## Chromosome Y - haploid
chrY_dat_SCC_males <- read.table(
  "../performance_metrics/EUR/males/data/EUR_males_chrY_haploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrY_dat_SCC_males <- chrY_dat_SCC_males[,!(names(chrY_dat_SCC_males) %in% drop)]

# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrY_dat_SCC_males_merge <- merge(chrY_dat_SCC_males, chrY_sim_males, by = "Sample")

chrY_dat_SCC_males_merge$Chromosome <- "chrY_hap"
chrY_dat_SCC_males_merge$alignment <- "SCC"
chrY_dat_SCC_males_merge$FPtoTP <- chrY_dat_SCC_males_merge$FP / chrY_dat_SCC_males_merge$TP
chrY_dat_SCC_males_merge$FNtoTP <- chrY_dat_SCC_males_merge$FN / chrY_dat_SCC_males_merge$TP
chrY_dat_SCC_males_merge$FPtoSimVars <- chrY_dat_SCC_males_merge$FP / chrY_dat_SCC_males_merge$Num_sim_variants
chrY_dat_SCC_males_merge$FNtoSimVars <- chrY_dat_SCC_males_merge$FN / chrY_dat_SCC_males_merge$Num_sim_variants
chrY_dat_SCC_males_merge$TPtoSimVars <- chrY_dat_SCC_males_merge$TP / chrY_dat_SCC_males_merge$Num_sim_variants

dat_all_SCC_males <- rbind(dat_all_SCC_males, chrY_dat_SCC_males_merge)

## Chromosome Y - diploid
chrY_dat_SCC_males <- read.table(
  "../performance_metrics/EUR/males/data/EUR_males_chrY_diploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrY_dat_SCC_males <- chrY_dat_SCC_males[,!(names(chrY_dat_SCC_males) %in% drop)]

# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrY_dat_SCC_males_merge <- merge(chrY_dat_SCC_males, chrY_sim_males, by = "Sample")

chrY_dat_SCC_males_merge$Chromosome <- "chrY_dip"
chrY_dat_SCC_males_merge$alignment <- "SCC"
chrY_dat_SCC_males_merge$FPtoTP <- chrY_dat_SCC_males_merge$FP / chrY_dat_SCC_males_merge$TP
chrY_dat_SCC_males_merge$FNtoTP <- chrY_dat_SCC_males_merge$FN / chrY_dat_SCC_males_merge$TP
chrY_dat_SCC_males_merge$FPtoSimVars <- chrY_dat_SCC_males_merge$FP / chrY_dat_SCC_males_merge$Num_sim_variants
chrY_dat_SCC_males_merge$FNtoSimVars <- chrY_dat_SCC_males_merge$FN / chrY_dat_SCC_males_merge$Num_sim_variants
chrY_dat_SCC_males_merge$TPtoSimVars <- chrY_dat_SCC_males_merge$TP / chrY_dat_SCC_males_merge$Num_sim_variants

dat_all_SCC_males <- rbind(dat_all_SCC_males, chrY_dat_SCC_males_merge)

## Chromosome X PARs
chrX_PARs_dat_SCC_males <- read.table(
  "../performance_metrics/EUR/males/data/EUR_males_chrX_PARs_diploid_golden_vs_called_performance_metrics.txt",
  header = T
)

drop <- c("Validated_Genotypes")
chrX_PARs_dat_SCC_males <- chrX_PARs_dat_SCC_males[,!(names(chrX_PARs_dat_SCC_males) %in% drop)]

# Merge with simulated counts and then do FPs and FNs over total simulated variants
chrX_PARs_dat_SCC_males_merge <- merge(chrX_PARs_dat_SCC_males, chrX_nonPARs_sim_males, by = "Sample")

chrX_PARs_dat_SCC_males_merge$Chromosome <- "chrX_PARs"
chrX_PARs_dat_SCC_males_merge$alignment <- "SCC"
chrX_PARs_dat_SCC_males_merge$FPtoTP <- chrX_PARs_dat_SCC_males_merge$FP / chrX_PARs_dat_SCC_males_merge$TP
chrX_PARs_dat_SCC_males_merge$FNtoTP <- chrX_PARs_dat_SCC_males_merge$FN / chrX_PARs_dat_SCC_males_merge$TP
chrX_PARs_dat_SCC_males_merge$FPtoSimVars <- chrX_PARs_dat_SCC_males_merge$FP / chrX_PARs_dat_SCC_males_merge$Num_sim_variants
chrX_PARs_dat_SCC_males_merge$FNtoSimVars <- chrX_PARs_dat_SCC_males_merge$FN / chrX_PARs_dat_SCC_males_merge$Num_sim_variants
chrX_PARs_dat_SCC_males_merge$TPtoSimVars <- chrX_PARs_dat_SCC_males_merge$TP / chrX_PARs_dat_SCC_males_merge$Num_sim_variants

dat_all_SCC_males <- rbind(dat_all_SCC_males, chrX_PARs_dat_SCC_males_merge)






#-----------------------------------------------------------------------------#
## PLOT ##
#-----------------------------------------------------------------------------#
# Remove chrM for now and PARs
dat_all_SCC_males_noM <- dat_all_SCC_males[dat_all_SCC_males$Chromosome != "chrM_hap", ]  
dat_all_SCC_males_noM <- dat_all_SCC_males_noM[dat_all_SCC_males_noM$Chromosome != "chrM_dip", ]  
dat_all_SCC_males_noM_PARs <- dat_all_SCC_males_noM[dat_all_SCC_males_noM$Chromosome != "chrX_PARs", ]  
#dat_all_SCC_males_noM_PARs_Y <- dat_all_SCC_males_noM_PARs[dat_all_SCC_males_noM_PARs$Chromosome != "chrY", ]  
#dat_all_SCC_males_noM_Y <- dat_all_SCC_males_noM[dat_all_SCC_males_noM$Chromosome != "chrY", ]  


males_FPtoSimVars_X_Y_d_h <- dat_all_SCC_males_noM_PARs %>%
  arrange(FPtoSimVars) %>%
  mutate(name = factor(Chromosome, 
                       levels=c("chr1", "chr2", "chr3", "chr4", 
                                "chr5", "chr6", "chr7", "chr8", 
                                "chr9", "chr10", "chr11", "chr12", 
                                "chr13", "chr14", "chr15", "chr16", 
                                "chr17", "chr18", "chr19", "chr20", 
                                "chr21", "chr22",
                                #"chrM",
                                #"chrX_PARs",
                                "chrX_nonPARs_dip", 
                                "chrX_nonPARs_hap", 
                                "chrY_dip",
                                "chrY_hap"
                       ))) %>%
  ggplot( aes(x=name, y=FPtoSimVars, color=Chromosome)) +
  #ggplot( aes(x=name, y=FPtoSimVars, color=name, fill=alignment)) +
  #ylim(0,0.004) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_dodge(.75), size =.5) +
  theme_classic() +
  # color names for chr 8, PARs, nonPARs diploid, nonPARs haploid, mtDNA diploid,
  # mtDNA haploid, Y diploid, Y haploid
  # "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00"
  ylab("False Positives/Total # Simulated") +
  #theme(legend.position = "none") + 
  xlab("Chromosome") +
  scale_x_discrete(labels=c("chr1" = "1", "chr2" = "2", "chr3" = "3", "chr4" = "4", 
                            "chr5" = "5", "chr6" = "6", "chr7" = "7", "chr8" = "8", 
                            "chr9" = "9", "chr10" = "10", "chr11" = "11", "chr12" = "12", 
                            "chr13" = "13", "chr14" = "14", "chr15" = "15", "chr16" = "16", 
                            "chr17" = "17", "chr18" = "18", "chr19" = "19", "chr20" = "20", 
                            "chr21" = "21", "chr22" = "22",
                            "chrM" = "M", 
                            "chrX_PARs" = "PARs",
                            "chrX_nonPARs_dip" = expression("X"["non-PARs (Dip)"]),
                            "chrX_nonPARs_hap" = expression("X"["non-PARs (Hap)"]),
                            #"chrY" = "Y non-PARs"
                            "chrY_dip" = expression("Y"["non-PARs (Dip)"]),
                            "chrY_hap" = expression("Y"["non-PARs (Hap)"])
  )) +
  #scale_color_manual(values = c(rep("grey", 22), "#B2DF8A", "#FF7F00")) +
  scale_color_manual(values = c(rep("grey", 22), "#B2DF8A", "#33A02C", "#FDBF6F", "#FF7F00")) +
  ggtitle("False Positives/Total # Simulated") +
  theme(
    plot.title = element_text(hjust = 0.5),
    #plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    legend.position="none",
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_hline(yintercept=males_default_max_FPtoSimVars, linetype="dashed", color = "red") +
  geom_hline(yintercept=males_default_min_FPtoSimVars, linetype="dashed", color = "red")



males_FNtoSimVars_X_Y_d_h <- dat_all_SCC_males_noM_PARs %>%
  arrange(FNtoSimVars) %>%
  mutate(name = factor(Chromosome, 
                       levels=c("chr1", "chr2", "chr3", "chr4", 
                                "chr5", "chr6", "chr7", "chr8", 
                                "chr9", "chr10", "chr11", "chr12", 
                                "chr13", "chr14", "chr15", "chr16", 
                                "chr17", "chr18", "chr19", "chr20", 
                                "chr21", "chr22",
                                #"chrM",
                                #"chrX_PARs",
                                "chrX_nonPARs_dip", 
                                "chrX_nonPARs_hap", 
                                "chrY_dip",
                                "chrY_hap"
                       ))) %>%
  ggplot( aes(x=name, y=FNtoSimVars, color=Chromosome)) +
  #ggplot( aes(x=name, y=FPtoSimVars, color=name, fill=alignment)) +
  #ylim(0,0.004) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_dodge(.75), size =.5) +
  theme_classic() +
  # color names for chr 8, PARs, nonPARs diploid, nonPARs haploid, mtDNA diploid,
  # mtDNA haploid, Y diploid, Y haploid
  # "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00"
  ylab("False Negatives/Total # Simulated") +
  #theme(legend.position = "none") + 
  xlab("Chromosome") +
  scale_x_discrete(labels=c("chr1" = "1", "chr2" = "2", "chr3" = "3", "chr4" = "4", 
                            "chr5" = "5", "chr6" = "6", "chr7" = "7", "chr8" = "8", 
                            "chr9" = "9", "chr10" = "10", "chr11" = "11", "chr12" = "12", 
                            "chr13" = "13", "chr14" = "14", "chr15" = "15", "chr16" = "16", 
                            "chr17" = "17", "chr18" = "18", "chr19" = "19", "chr20" = "20", 
                            "chr21" = "21", "chr22" = "22",
                            "chrM" = "M", 
                            "chrX_PARs" = "PARs",
                            "chrX_nonPARs_dip" = expression("X"["non-PARs (Dip)"]),
                            "chrX_nonPARs_hap" = expression("X"["non-PARs (Hap)"]),
                            #"chrY" = "Y non-PARs"
                            "chrY_dip" = expression("Y"["non-PARs (Dip)"]),
                            "chrY_hap" = expression("Y"["non-PARs (Hap)"])
  )) +
  #scale_color_manual(values = c(rep("grey", 22), "#B2DF8A", "#FF7F00")) +
  scale_color_manual(values = c(rep("grey", 22), "#B2DF8A", "#33A02C", "#FDBF6F", "#FF7F00")) +
  ggtitle("False Negatives/Total # Simulated") +
  theme(
    plot.title = element_text(hjust = 0.5),
    #plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    legend.position="none",
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_hline(yintercept=males_default_max_FNtoSimVars, linetype="dashed", color = "red") +
  geom_hline(yintercept=males_default_min_FNtoSimVars, linetype="dashed", color = "red")





# Plot hap vs dip on X and Y no autosomes 
dat_all_SCC_males_XYonly <- dat_all_SCC_males[ which(dat_all_SCC_males$Chromosome=="chrX_nonPARs_hap" |
                                                       dat_all_SCC_males$Chromosome=="chrX_nonPARs_dip" |
                                                       dat_all_SCC_males$Chromosome=="chrY_dip" |
                                                       dat_all_SCC_males$Chromosome=="chrY_hap"),]

                
males_FPtoSimVars_X_Y_d_h_02 <- dat_all_SCC_males_XYonly %>%
  arrange(FPtoSimVars) %>%
  mutate(name = factor(Chromosome, 
                       levels=c(
                                #"chrM",
                                #"chrX_PARs",
                                "chrX_nonPARs_dip", 
                                "chrX_nonPARs_hap", 
                                "chrY_dip",
                                "chrY_hap"
                       ))) %>%
  ggplot( aes(x=name, y=FPtoSimVars, color=Chromosome)) +
  #ggplot( aes(x=name, y=FPtoSimVars, color=name, fill=alignment)) +
  #ylim(0,0.004) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_dodge(.75)) +
  theme_classic() +
  # color names for chr 8, PARs, nonPARs diploid, nonPARs haploid, mtDNA diploid,
  # mtDNA haploid, Y diploid, Y haploid
  # "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00"
  ylab("False Positives/Total # Simulated") +
  #theme(legend.position = "none") + 
  xlab("Chromosome") +
  scale_x_discrete(labels=c("chr1" = "1", "chr2" = "2", "chr3" = "3", "chr4" = "4", 
                            "chr5" = "5", "chr6" = "6", "chr7" = "7", "chr8" = "8", 
                            "chr9" = "9", "chr10" = "10", "chr11" = "11", "chr12" = "12", 
                            "chr13" = "13", "chr14" = "14", "chr15" = "15", "chr16" = "16", 
                            "chr17" = "17", "chr18" = "18", "chr19" = "19", "chr20" = "20", 
                            "chr21" = "21", "chr22" = "22",
                            "chrM" = "M", 
                            "chrX_PARs" = "PARs",
                            "chrX_nonPARs_dip" = expression("X"["non-PARs (Dip)"]),
                            "chrX_nonPARs_hap" = expression("X"["non-PARs (Hap)"]),
                            #"chrY" = "Y non-PARs"
                            "chrY_dip" = expression("Y"["non-PARs (Dip)"]),
                            "chrY_hap" = expression("Y"["non-PARs (Hap)"])
  )) +
  #scale_color_manual(values = c(rep("grey", 22), "#B2DF8A", "#FF7F00")) +
  scale_color_manual(values = c("#B2DF8A", "#33A02C", "#FDBF6F", "#FF7F00")) +
  ggtitle("False Positives/Total # Simulated") +
  theme(
    plot.title = element_text(hjust = 0.5),
    #plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    legend.position="none",
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_hline(yintercept=males_default_max_FPtoSimVars, linetype="dashed", color = "red") +
  geom_hline(yintercept=males_default_min_FPtoSimVars, linetype="dashed", color = "red")



males_FNtoSimVars_X_Y_d_h_02 <- dat_all_SCC_males_XYonly %>%
  arrange(FNtoSimVars) %>%
  mutate(name = factor(Chromosome, 
                       levels=c(
                                #"chrM",
                                #"chrX_PARs",
                                "chrX_nonPARs_dip", 
                                "chrX_nonPARs_hap", 
                                "chrY_dip",
                                "chrY_hap"
                       ))) %>%
  ggplot( aes(x=name, y=FNtoSimVars, color=Chromosome)) +
  #ggplot( aes(x=name, y=FPtoSimVars, color=name, fill=alignment)) +
  #ylim(0,0.004) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_dodge(.75)) +
  theme_classic() +
  # color names for chr 8, PARs, nonPARs diploid, nonPARs haploid, mtDNA diploid,
  # mtDNA haploid, Y diploid, Y haploid
  # "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00"
  ylab("False Negatives/Total # Simulated") +
  #theme(legend.position = "none") + 
  xlab("Chromosome") +
  scale_x_discrete(labels=c("chr1" = "1", "chr2" = "2", "chr3" = "3", "chr4" = "4", 
                            "chr5" = "5", "chr6" = "6", "chr7" = "7", "chr8" = "8", 
                            "chr9" = "9", "chr10" = "10", "chr11" = "11", "chr12" = "12", 
                            "chr13" = "13", "chr14" = "14", "chr15" = "15", "chr16" = "16", 
                            "chr17" = "17", "chr18" = "18", "chr19" = "19", "chr20" = "20", 
                            "chr21" = "21", "chr22" = "22",
                            "chrM" = "M", 
                            "chrX_PARs" = "PARs",
                            "chrX_nonPARs_dip" = expression("X"["non-PARs (Dip)"]),
                            "chrX_nonPARs_hap" = expression("X"["non-PARs (Hap)"]),
                            #"chrY" = "Y non-PARs"
                            "chrY_dip" = expression("Y"["non-PARs (Dip)"]),
                            "chrY_hap" = expression("Y"["non-PARs (Hap)"])
  )) +
  #scale_color_manual(values = c(rep("grey", 22), "#B2DF8A", "#FF7F00")) +
  scale_color_manual(values = c("#B2DF8A", "#33A02C", "#FDBF6F", "#FF7F00")) +
  ggtitle("False Negatives/Total # Simulated") +
  theme(
    plot.title = element_text(hjust = 0.5),
    #plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    legend.position="none",
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_hline(yintercept=males_default_max_FNtoSimVars, linetype="dashed", color = "red") +
  geom_hline(yintercept=males_default_min_FNtoSimVars, linetype="dashed", color = "red")





#-----------------------------------------------------------------------------#
## OUTPUT PLOTS ##
#-----------------------------------------------------------------------------#
#males_FPtoSimVars_X_Y_d_h
#males_FNtoSimVars_X_Y_d_h

male_XY_plot <- ggarrange(
  males_FPtoSimVars_X_Y_d_h, 
  males_FNtoSimVars_X_Y_d_h, 
  ncol = 1,
  nrow = 2,
  align = c("h"),
  labels = c("A", "B")
)

male_XY_plot_annoated <- annotate_figure(male_XY_plot, top = text_grob("Simulated males (XY)", face = "bold", size = 14))


pdf("../plots/potential_figure_SCC_jitter_all_autosomes_XY_males_hap_dip.pdf",
    width = 6,
    height = 7.5)
print(male_XY_plot_annoated)
dev.off()

# NEW
male_XY_plot <- ggarrange(
  males_FPtoSimVars_X_Y_d_h, 
  males_FNtoSimVars_X_Y_d_h, 
  nrow = 1,
  ncol = 2,
  align = c("h"),
  labels = c("A", "B")
)

male_XY_plot_annoated <- annotate_figure(male_XY_plot, top = text_grob("Simulated males (XY)", face = "bold", size = 14))


pdf("../plots/potential_figure_SCC_jitter_all_autosomes_XY_males_hap_dip_NEW.pdf",
    width = 10,
    height = 4.5)
print(male_XY_plot_annoated)
dev.off()



male_XY_plot <- ggarrange(
  males_FPtoSimVars_X_Y_d_h_02, 
  males_FNtoSimVars_X_Y_d_h_02, 
  ncol = 2,
  nrow = 1,
  align = c("h"),
  labels = c("A", "B")
)

male_XY_plot_annoated <- annotate_figure(male_XY_plot, top = text_grob("Simulated males (XY)", face = "bold", size = 14))


pdf("../plots/potential_figure_SCC_jitter_XY_males_hap_dip_auto_ranges.pdf",
    width = 7.5,
    height = 4)
print(male_XY_plot_annoated)
dev.off()
