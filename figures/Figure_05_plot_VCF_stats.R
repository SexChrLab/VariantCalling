# Plot VCF stats

################
# Load libraries
################
library(ggplot2)
library(cowplot)
library(ggpubr)
library(dplyr)

###############################################################################
###############################################################################
# MALES #
###############################################################################
###############################################################################

##############
# Read in data
##############

males_dat_autosomes <- c()


# Autosomes
for(i in 1:22) {
  print(paste("Processing chromosome", i))
  chri_dat <- read.table(paste("../performance_metrics/EUR/males/vcf_stats/chr", i, 
                               "_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_gatkHardFilter_stats.txt", sep = ""),
                         header = T)
  chri_dat$Chromosome <- "Autosomes"

  males_dat_autosomes <- rbind(males_dat_autosomes, chri_dat)
}

# PARs
males_dat_PARs <- read.table(
  "../performance_metrics/EUR/males/vcf_stats/chrX_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_gatkHardFilter_PARs_stats.txt",
  header = T)
males_dat_PARs$Chromosome <- "PARs"

# X non PARs
males_dat_chrX_nonPARs <- read.table(
  "../performance_metrics/EUR/males/vcf_stats/chrX_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs_gatkHardFilter_nonPARs_stats.txt",
  header = T)
males_dat_chrX_nonPARs$Chromosome <- "X (non-PARs)"

# Y non PARs
males_dat_chrY <- read.table(
  "../performance_metrics/EUR/males/vcf_stats/chrY_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs_gatkHardFilter_stats.txt",
  header = T)
males_dat_chrY$Chromosome <- "Y (non-PARs)"


#####################################
# Merge all data into one data frame
#####################################
males_all_data_merged <- rbind(males_dat_autosomes, males_dat_PARs,
                         males_dat_chrX_nonPARs, males_dat_chrY)

males_all_data_merged <- rbind(males_dat_autosomes,
                               males_dat_chrX_nonPARs, males_dat_chrY)

###############################################
# Plot autosomes vs sex chromosomes on one plot 
###############################################
#----#
# QD #
#----#
male_QD <- ggplot(males_all_data_merged,aes(x=QD, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated males (XY)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=2, linetype="dashed", color = "red")


#-----#
# SOR #
#-----#
male_SOR <- ggplot(males_all_data_merged,aes(x=SOR, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated males (XY)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=3, linetype="dashed", color = "red")


#----#
# FS #
#----#
male_FS <- ggplot(males_all_data_merged,aes(x=FS, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated males (XY)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=60, linetype="dashed", color = "red")


#----#
# MQ #
#----#
male_MQ <- ggplot(males_all_data_merged,aes(x=MQ, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated males (XY)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=40, linetype="dashed", color = "red")


#-----------#
# MQRankSum #
#-----------#
male_MQRankSum <- ggplot(males_all_data_merged,aes(x=MQRankSum, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated males (XY)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=-12.5, linetype="dashed", color = "red")


#----------------#
# ReadPosRankSum #
#----------------#
male_ReadPosRankSum <- ggplot(males_all_data_merged,aes(x=ReadPosRankSum, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated males (XY)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=-8, linetype="dashed", color = "red")



#----#
# AN #
#----#
male_AN <- ggplot(males_all_data_merged,aes(x=AN, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated males (XY)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=20, linetype="dashed", color = "red") # annotating autosomal expectation


#----#
# DP #
#----#
mean(males_dat_autosomes$DP)
mean(males_dat_chrX_nonPARs$DP)
mean(males_dat_chrY$DP)
# mean auto: 133.9931
# mean X: 46.17204
# mean Y: 46.44892

male_DP <- ggplot(males_all_data_merged,aes(x=DP, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated males (XY)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=67, linetype="dashed", color = "red") + # annotating autosomal expectation
  geom_vline(xintercept=201, linetype="dashed", color = "red") # annotating autosomal expectation






###############################################################################
###############################################################################
# FEMALES #
###############################################################################
###############################################################################

##############
# Read in data
##############

females_dat_autosomes <- c()


# Autosomes
for(i in 1:22) {
  print(paste("Processing chromosome", i))
  chri_dat <- read.table(paste("../performance_metrics/EUR/females/vcf_stats/chr", i, 
                               "_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs_gatkHardFilter_stats.txt", sep = ""),
                         header = T)
  chri_dat$Chromosome <- "Autosomes"
  
  females_dat_autosomes <- rbind(females_dat_autosomes, chri_dat)
}

# PARs
females_dat_PARs <- read.table(
  "../performance_metrics/EUR/females/vcf_stats/chrX_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs_gatkHardFilter_PARs_stats.txt",
  header = T)
females_dat_PARs$Chromosome <- "PARs"

# X non PARs
females_dat_chrX_nonPARs <- read.table(
  "../performance_metrics/EUR/females/vcf_stats/chrX_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs_gatkHardFilter_nonPARs_stats.txt",
  header = T)
females_dat_chrX_nonPARs$Chromosome <- "X (non-PARs)"


#####################################
# Merge all data into one data frame
#####################################
females_all_data_merged <- rbind(females_dat_autosomes, females_dat_PARs,
                                 females_dat_chrX_nonPARs)

females_all_data_merged <- rbind(females_dat_autosomes,
                                 females_dat_chrX_nonPARs)

###############################################
# Plot autosomes vs sex chromosomes on one plot 
###############################################
#----#
# QD #
#----#
female_QD <- ggplot(females_all_data_merged,aes(x=QD, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated females (XX)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=2, linetype="dashed", color = "red")


#-----#
# SOR #
#-----#
female_SOR <- ggplot(females_all_data_merged,aes(x=SOR, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated females (XX)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=3, linetype="dashed", color = "red")


#----#
# FS #
#----#
female_FS <- ggplot(females_all_data_merged,aes(x=FS, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated females (XX)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=60, linetype="dashed", color = "red")



#----#
# MQ #
#----#
female_MQ <- ggplot(females_all_data_merged,aes(x=MQ, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated females (XX)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=40, linetype="dashed", color = "red")



#-----------#
# MQRankSum #
#-----------#
female_MQRankSum <- ggplot(females_all_data_merged,aes(x=MQRankSum, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated females (XX)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=-12.5, linetype="dashed", color = "red")


#----------------#
# ReadPosRankSum #
#----------------#
female_ReadPosRankSum <- ggplot(females_all_data_merged,aes(x=ReadPosRankSum, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated females (XX)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=-8, linetype="dashed", color = "red")


#----#
# AN #
#----#
female_AN <- ggplot(females_all_data_merged,aes(x=AN, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated females (XX)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=20, linetype="dashed", color = "red") # annotating autosomal expectation


#----#
# DP #
#----#
mean(females_dat_autosomes$DP)
mean(females_dat_chrX_nonPARs$DP)
# mean auto: 134.0021
# mean X: 135.8596

female_DP <- ggplot(females_all_data_merged,aes(x=DP, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated females (XX)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=67, linetype="dashed", color = "red") + # annotating autosomal expectation
  geom_vline(xintercept=201, linetype="dashed", color = "red") # annotating autosomal expectation


###############################################################################
###############################################################################
# PLOT ALL DATA - MALES FEMALES SIDE-BY-SIDE #
###############################################################################
###############################################################################
# QD
QD_plot <- ggarrange(
  female_QD, 
  male_QD, 
  ncol = 2,
  nrow = 1,
  align = c("h"),
  labels = c("A", "B"),
  legend = "top"
)

# SOR
SOR_plot <- ggarrange(
  female_SOR, 
  male_SOR, 
  ncol = 2,
  nrow = 1,
  align = c("h"),
  labels = c("A", "B"),
  legend = "top"
)

# FS
FS_plot <- ggarrange(
  female_FS, 
  male_FS, 
  ncol = 2,
  nrow = 1,
  align = c("h"),
  labels = c("A", "B"),
  legend = "top"
)

# MQ
MQ_plot <- ggarrange(
  female_MQ, 
  male_MQ, 
  ncol = 2,
  nrow = 1,
  align = c("h"),
  labels = c("A", "B"),
  legend = "top"
)

# MQRankSum
MQRankSum_plot <- ggarrange(
  female_MQRankSum, 
  male_MQRankSum, 
  ncol = 2,
  nrow = 1,
  align = c("h"),
  labels = c("A", "B"),
  legend = "top"
)

# ReadPosRankSum
ReadPosRankSum_plot <- ggarrange(
  female_ReadPosRankSum, 
  male_ReadPosRankSum, 
  ncol = 2,
  nrow = 1,
  align = c("h"),
  labels = c("A", "B"),
  legend = "top"
)

# AN
AN_plot <- ggarrange(
  female_AN, 
  male_AN, 
  ncol = 2,
  nrow = 1,
  align = c("h"),
  labels = c("A", "B"),
  legend = "top"
)

# DP
DP_plot <- ggarrange(
  female_DP, 
  male_DP, 
  ncol = 2,
  nrow = 1,
  align = c("h"),
  labels = c("A", "B"),
  legend = "top"
)


pdf("../plots/SCC_female_male_VCF_stats.pdf",
    height = 5, width = 9)

print(QD_plot)
print(SOR_plot)
print(FS_plot)
print(MQ_plot)
print(MQRankSum_plot)
print(ReadPosRankSum_plot)
print(AN_plot)
print(DP_plot)

dev.off()

pdf("../plots/SCC_female_male_VCF_stats_on_one.pdf",
    height = 11, width = 13)
ggarrange(female_QD, male_QD, female_SOR, male_SOR, 
          female_FS, male_FS, female_MQ, male_MQ,
          female_MQRankSum, male_MQRankSum, female_ReadPosRankSum, male_ReadPosRankSum,
          female_AN, male_AN, female_DP, male_DP,
          ncol = 4,
          nrow = 4,
          align = c("h"),
          labels = c("A", "", "B", "", 
                     "C", "", "D", "",
                     "E", "", "F", "",
                     "G", "", "H", ""),
          legend = "top"
          )
dev.off()



###############################################################################
#######################
#######################
# DEFAULT DATA 
#######################
#######################
###############################################################################


###############################################################################
###############################################################################
# MALES #
###############################################################################
###############################################################################

##############
# Read in data
##############

males_dat_autosomes <- c()


# Autosomes
for(i in 1:22) {
  print(paste("Processing chromosome", i))
  chri_dat <- read.table(paste("../performance_metrics/EUR/males_default/vcf_stats/chr", i, 
                               "_GRCh38_default_gatk_diploid_called_raw_SNPs_gatkHardFilter_stats.txt", sep = ""),
                         header = T)
  chri_dat$Chromosome <- "Autosomes"
  
  males_dat_autosomes <- rbind(males_dat_autosomes, chri_dat)
}

# PARs
males_dat_PARs <- read.table(
  "../performance_metrics/EUR/males_default/vcf_stats/chrX_GRCh38_default_gatk_diploid_called_raw_SNPs_gatkHardFilter_PARs_stats.txt",
  header = T)
males_dat_PARs$Chromosome <- "PARs"

# X non PARs
males_dat_chrX_nonPARs <- read.table(
  "../performance_metrics/EUR/males_default/vcf_stats/chrX_GRCh38_default_gatk_haploid_called_raw_SNPs_gatkHardFilter_nonPARs_stats.txt",
  header = T)
males_dat_chrX_nonPARs$Chromosome <- "X (non-PARs)"

# Y non PARs
males_dat_chrY <- read.table(
  "../performance_metrics/EUR/males_default/vcf_stats/chrY_GRCh38_default_gatk_haploid_called_raw_SNPs_gatkHardFilter_stats.txt",
  header = T)
males_dat_chrY$Chromosome <- "Y (non-PARs)"


#####################################
# Merge all data into one data frame
#####################################
males_all_data_merged <- rbind(males_dat_autosomes, males_dat_PARs,
                               males_dat_chrX_nonPARs, males_dat_chrY)

males_all_data_merged <- rbind(males_dat_autosomes,
                               males_dat_chrX_nonPARs, males_dat_chrY)

###############################################
# Plot autosomes vs sex chromosomes on one plot 
###############################################
#----#
# QD #
#----#
male_QD <- ggplot(males_all_data_merged,aes(x=QD, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated males (XY)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=2, linetype="dashed", color = "red")


#-----#
# SOR #
#-----#
male_SOR <- ggplot(males_all_data_merged,aes(x=SOR, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated males (XY)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=3, linetype="dashed", color = "red")


#----#
# FS #
#----#
male_FS <- ggplot(males_all_data_merged,aes(x=FS, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated males (XY)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=60, linetype="dashed", color = "red")


#----#
# MQ #
#----#
male_MQ <- ggplot(males_all_data_merged,aes(x=MQ, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated males (XY)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=40, linetype="dashed", color = "red")


#-----------#
# MQRankSum #
#-----------#
male_MQRankSum <- ggplot(males_all_data_merged,aes(x=MQRankSum, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated males (XY)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=-12.5, linetype="dashed", color = "red")


#----------------#
# ReadPosRankSum #
#----------------#
male_ReadPosRankSum <- ggplot(males_all_data_merged,aes(x=ReadPosRankSum, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated males (XY)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=-8, linetype="dashed", color = "red")



#----#
# AN #
#----#
male_AN <- ggplot(males_all_data_merged,aes(x=AN, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated males (XY)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=20, linetype="dashed", color = "red") # annotating autosomal expectation


#----#
# DP #
#----#
mean(males_dat_autosomes$DP)
mean(males_dat_chrX_nonPARs$DP)
mean(males_dat_chrY$DP)
# mean auto: 133.9932
# mean X: 46.17181
# mean Y: 46.44863

male_DP <- ggplot(males_all_data_merged,aes(x=DP, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated males (XY)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=67, linetype="dashed", color = "red") + # annotating autosomal expectation
  geom_vline(xintercept=201, linetype="dashed", color = "red") # annotating autosomal expectation






###############################################################################
###############################################################################
# FEMALES #
###############################################################################
###############################################################################

##############
# Read in data
##############

females_dat_autosomes <- c()


# Autosomes
for(i in 1:22) {
  print(paste("Processing chromosome", i))
  chri_dat <- read.table(paste("../performance_metrics/EUR/females_default/vcf_stats/chr", i, 
                               "_GRCh38_default_gatk_diploid_called_raw_SNPs_gatkHardFilter_stats.txt", sep = ""),
                         header = T)
  chri_dat$Chromosome <- "Autosomes"
  
  females_dat_autosomes <- rbind(females_dat_autosomes, chri_dat)
}

# PARs
females_dat_PARs <- read.table(
  "../performance_metrics/EUR/females_default/vcf_stats/chrX_GRCh38_default_gatk_diploid_called_raw_SNPs_gatkHardFilter_PARs_stats.txt",
  header = T)
females_dat_PARs$Chromosome <- "PARs"

# X non PARs
females_dat_chrX_nonPARs <- read.table(
  "../performance_metrics/EUR/females_default/vcf_stats/chrX_GRCh38_default_gatk_diploid_called_raw_SNPs_gatkHardFilter_nonPARs_stats.txt",
  header = T)
females_dat_chrX_nonPARs$Chromosome <- "X (non-PARs)"


#####################################
# Merge all data into one data frame
#####################################
females_all_data_merged <- rbind(females_dat_autosomes, females_dat_PARs,
                                 females_dat_chrX_nonPARs)

females_all_data_merged <- rbind(females_dat_autosomes,
                                 females_dat_chrX_nonPARs)

###############################################
# Plot autosomes vs sex chromosomes on one plot 
###############################################
#----#
# QD #
#----#
female_QD <- ggplot(females_all_data_merged,aes(x=QD, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated females (XX)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=2, linetype="dashed", color = "red")


#-----#
# SOR #
#-----#
female_SOR <- ggplot(females_all_data_merged,aes(x=SOR, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated females (XX)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=3, linetype="dashed", color = "red")


#----#
# FS #
#----#
female_FS <- ggplot(females_all_data_merged,aes(x=FS, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated females (XX)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=60, linetype="dashed", color = "red")



#----#
# MQ #
#----#
female_MQ <- ggplot(females_all_data_merged,aes(x=MQ, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated females (XX)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=40, linetype="dashed", color = "red")



#-----------#
# MQRankSum #
#-----------#
female_MQRankSum <- ggplot(females_all_data_merged,aes(x=MQRankSum, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated females (XX)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=-12.5, linetype="dashed", color = "red")


#----------------#
# ReadPosRankSum #
#----------------#
female_ReadPosRankSum <- ggplot(females_all_data_merged,aes(x=ReadPosRankSum, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated females (XX)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=-8, linetype="dashed", color = "red")


#----#
# AN #
#----#
female_AN <- ggplot(females_all_data_merged,aes(x=AN, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated females (XX)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=20, linetype="dashed", color = "red") # annotating autosomal expectation


#----#
# DP #
#----#
mean(females_dat_autosomes$DP)
mean(females_dat_chrX_nonPARs$DP)
# mean auto: 133.9986
# mean X: 135.4546

female_DP <- ggplot(females_all_data_merged,aes(x=DP, fill=Chromosome)) + geom_density(alpha=0.25) +
  theme_classic() +
  guides(fill=guide_legend(title="Key")) +
  ggtitle("Simulated females (XX)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9)) + 
  geom_vline(xintercept=67, linetype="dashed", color = "red") + # annotating autosomal expectation
  geom_vline(xintercept=201, linetype="dashed", color = "red") # annotating autosomal expectation


###############################################################################
###############################################################################
# PLOT ALL DATA - MALES FEMALES SIDE-BY-SIDE #
###############################################################################
###############################################################################
# QD
QD_plot <- ggarrange(
  female_QD, 
  male_QD, 
  ncol = 2,
  nrow = 1,
  align = c("h"),
  labels = c("A", "B"),
  legend = "top"
)

# SOR
SOR_plot <- ggarrange(
  female_SOR, 
  male_SOR, 
  ncol = 2,
  nrow = 1,
  align = c("h"),
  labels = c("A", "B"),
  legend = "top"
)

# FS
FS_plot <- ggarrange(
  female_FS, 
  male_FS, 
  ncol = 2,
  nrow = 1,
  align = c("h"),
  labels = c("A", "B"),
  legend = "top"
)

# MQ
MQ_plot <- ggarrange(
  female_MQ, 
  male_MQ, 
  ncol = 2,
  nrow = 1,
  align = c("h"),
  labels = c("A", "B"),
  legend = "top"
)

# MQRankSum
MQRankSum_plot <- ggarrange(
  female_MQRankSum, 
  male_MQRankSum, 
  ncol = 2,
  nrow = 1,
  align = c("h"),
  labels = c("A", "B"),
  legend = "top"
)

# ReadPosRankSum
ReadPosRankSum_plot <- ggarrange(
  female_ReadPosRankSum, 
  male_ReadPosRankSum, 
  ncol = 2,
  nrow = 1,
  align = c("h"),
  labels = c("A", "B"),
  legend = "top"
)

# AN
AN_plot <- ggarrange(
  female_AN, 
  male_AN, 
  ncol = 2,
  nrow = 1,
  align = c("h"),
  labels = c("A", "B"),
  legend = "top"
)

# DP
DP_plot <- ggarrange(
  female_DP, 
  male_DP, 
  ncol = 2,
  nrow = 1,
  align = c("h"),
  labels = c("A", "B"),
  legend = "top"
)


pdf("../plots/default_female_male_VCF_stats.pdf",
    height = 5, width = 9)

print(QD_plot)
print(SOR_plot)
print(FS_plot)
print(MQ_plot)
print(MQRankSum_plot)
print(ReadPosRankSum_plot)
print(AN_plot)
print(DP_plot)

dev.off()
