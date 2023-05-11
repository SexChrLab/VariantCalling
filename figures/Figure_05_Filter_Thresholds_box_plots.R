# Angela Oill
# Box plots of TP and FP for each filter

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# PLOT DIFFERENT THRESHOLDS
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

#----------------#
# Load Libraries #
#----------------#

library("dplyr")
library("ggplot2")
library(cowplot)
library(ggpubr)
library(viridis)
library("RColorBrewer")
#brewer.pal(n = 8, name = "Paired")

library(tidyverse)
library(ggpubr)
library(rstatix)

library("Rmisc")



#-----------------------------------------------------------------------------#
# QD #
#-----------------------------------------------------------------------------#
# Plot the different thresholds
QD_thresholds <- c("1.0", "1.5","2.0", "12.0", "16.0", "20.0", "28.0")

########
# Chr8 #
########
QD_data_all_thresholds_8 <- c()
for (val in QD_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/males/per_filter/QD/ALL_males_chr8_autos_diploid_QD_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chr8"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  QD_data_all_thresholds_8 <- rbind(QD_data_all_thresholds_8, dati)
}


QD_data_all_thresholds_8_reordered <- QD_data_all_thresholds_8 %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=QD_thresholds)) 

QD_8_TP <- ggplot(QD_data_all_thresholds_8_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("QD Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")

QD_8_FP <- ggplot(QD_data_all_thresholds_8_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("QD Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")


QD_8_plot <- ggarrange(
  QD_8_TP, QD_8_FP,
  ncol = 2,
  nrow = 1#,
  #labels = c("A", "B")
)
QD_8_plot_annoated <- annotate_figure(QD_8_plot, top = text_grob("Chromosome 8", face = "bold", size = 14))
#print(QD_8_plot_annoated)


########
# PARs # 
########
QD_data_all_thresholds_PARs <- c()
for (val in QD_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/males/per_filter/QD/ALL_males_chrX_PARs_diploid_QD_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chrXPARs"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  QD_data_all_thresholds_PARs <- rbind(QD_data_all_thresholds_PARs, dati)
}

QD_data_all_thresholds_PARs_reordered <- QD_data_all_thresholds_PARs %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=QD_thresholds)) 

QD_PARs_TP <- ggplot(QD_data_all_thresholds_PARs_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("QD Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")

QD_PARs_FP <- ggplot(QD_data_all_thresholds_PARs_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("QD Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")


QD_PARs_plot <- ggarrange(
  QD_PARs_TP, QD_PARs_FP,
  ncol = 2,
  nrow = 1,
  labels = c("C", "D")
)
QD_PARs_plot_annoated <- annotate_figure(QD_PARs_plot, top = text_grob("PARs", face = "bold", size = 14))
#print(QD_PARs_plot_annoated)


########################
# Chromosome X nonPARs # 
########################
# TODO test out what this would look like with diploid called data
QD_data_all_thresholds_X_nonPARs <- c()
for (val in QD_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/males/per_filter/QD/ALL_males_chrX_nonPARs_haploid_QD_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chrXnonPARs"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  QD_data_all_thresholds_X_nonPARs <- rbind(QD_data_all_thresholds_X_nonPARs, dati)
}


QD_data_all_thresholds_X_nonPARs_reordered <- QD_data_all_thresholds_X_nonPARs %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=QD_thresholds)) 

QD_X_nonPARs_TP <- ggplot(QD_data_all_thresholds_X_nonPARs_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("QD Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")

QD_X_nonPARs_FP <- ggplot(QD_data_all_thresholds_X_nonPARs_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("QD Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")


QD_X_nonPARs_plot <- ggarrange(
  QD_X_nonPARs_TP, QD_X_nonPARs_FP,
  ncol = 2,
  nrow = 1#,
  #labels = c("E", "F")
)
QD_X_nonPARs_plot_annoated <- annotate_figure(QD_X_nonPARs_plot, top = text_grob("Chromosome X non-PARs", face = "bold", size = 14))
#print(QD_X_nonPARs_plot_annoated)

################
# Y chromosome # 
################
QD_data_all_thresholds_Y <- c()
for (val in QD_thresholds) {
  #print(paste("../performance_metrics/EUR/males/per_filter/QD/ALL_males_chrY_haploid_QD_",
  #            val, "_golden_vs_called_performance_metrics.txt", sep = ""))
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/males/per_filter/QD/ALL_males_chrY_haploid_QD_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chrY"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  QD_data_all_thresholds_Y <- rbind(QD_data_all_thresholds_Y, dati)
}

QD_data_all_thresholds_Y_reordered <- QD_data_all_thresholds_Y %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=QD_thresholds)) 

QD_Y_TP <- ggplot(QD_data_all_thresholds_Y_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("QD Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")

QD_Y_FP <- ggplot(QD_data_all_thresholds_Y_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("QD Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")


QD_Y_plot <- ggarrange(
  QD_Y_TP, QD_Y_FP,
  ncol = 2,
  nrow = 1#,
  #labels = c("G", "H")
)
QD_Y_plot_annoated <- annotate_figure(QD_Y_plot, top = text_grob("Chromosome Y", face = "bold", size = 14))
#print(QD_Y_plot_annoated)


#########
# mtDNA # 
#########
QD_data_all_thresholds_M <- c()
for (val in QD_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/males/per_filter/QD/ALL_males_chrM_haploid_QD_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chrM"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  QD_data_all_thresholds_M <- rbind(QD_data_all_thresholds_M, dati)
}


QD_data_all_thresholds_M_reordered <- QD_data_all_thresholds_M %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=QD_thresholds)) 

QD_M_TP <- ggplot(QD_data_all_thresholds_M_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("QD Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")

QD_M_FP <- ggplot(QD_data_all_thresholds_M_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(
    height = .1, width = .2, 
    #position=position_jitter(0.2), 
    color = "#440154FF"
  ) +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("QD Threshold") +
  scale_y_continuous(limits = c(0,15)) +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")


QD_M_plot <- ggarrange(
  QD_M_TP, QD_M_FP,
  ncol = 2,
  nrow = 1,
  labels = c("I", "J")
)
QD_M_plot_annoated <- annotate_figure(QD_M_plot, top = text_grob("mtDNA", face = "bold", size = 14))
#print(QD_M_plot_annoated)


########
# PLOT #
########
pdf("../performance_metrics/EUR/males/per_filter/QD/QD_FP_TP_all_thresholds_haploid_boxplot.pdf",
    width = 6, height = 12)

ggarrange(
  QD_8_plot_annoated, QD_PARs_plot_annoated, 
  QD_X_nonPARs_plot_annoated, QD_Y_plot_annoated,
  QD_M_plot_annoated,
  ncol = 1,
  nrow = 5
)

dev.off()



pdf("../performance_metrics/EUR/males/per_filter/QD/QD_FP_TP_all_thresholds_haploid_boxplot_2.pdf",
    width = 10, height = 8)

ggarrange(
  QD_8_plot_annoated, QD_PARs_plot_annoated, 
  QD_X_nonPARs_plot_annoated, QD_Y_plot_annoated,
  QD_M_plot_annoated,
  ncol = 2,
  nrow = 3
)

dev.off()

pdf("../performance_metrics/EUR/males/per_filter/QD/QD_FP_TP_all_thresholds_haploid_boxplot_3.pdf",
    width = 6, height = 11)

ggarrange(
  QD_8_plot_annoated, 
  QD_X_nonPARs_plot_annoated, QD_Y_plot_annoated,
  ncol = 1,
  nrow = 3
)

dev.off()

#----------------------------------------------------------#
# NEW - prep for one large plot (potential main text figure)
#----------------------------------------------------------#
QD_X_nonPARs_plot <- ggarrange(
  QD_X_nonPARs_TP, QD_X_nonPARs_FP,
  ncol = 2,
  nrow = 1
)
QD_X_nonPARs_plot_annoated <- annotate_figure(QD_X_nonPARs_plot, top = text_grob("Chromosome X non-PARs", face = "bold", size = 14))

QD_Y_plot <- ggarrange(
  QD_Y_TP, QD_Y_FP,
  ncol = 2,
  nrow = 1
)
QD_Y_plot_annoated <- annotate_figure(QD_Y_plot, top = text_grob("Chromosome Y", face = "bold", size = 14))


QD_plot_males <- ggarrange(
  QD_X_nonPARs_plot_annoated, QD_Y_plot_annoated,
  ncol = 2,
  nrow = 1
)

#-----------------------------------------------------------------------------#
# MQ #
#-----------------------------------------------------------------------------#

# Plot the different thresholds
MQ_thresholds <- c("20.0", "30.0", "40.0", "50.0", "60.0")

########
# Chr8 #
########
MQ_data_all_thresholds_8 <- c()
for (val in MQ_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/males/per_filter/MQ/ALL_males_chr8_autos_diploid_MQ_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chr8"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  MQ_data_all_thresholds_8 <- rbind(MQ_data_all_thresholds_8, dati)
}


MQ_data_all_thresholds_8_reordered <- MQ_data_all_thresholds_8 %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=MQ_thresholds)) 

MQ_8_TP <- ggplot(MQ_data_all_thresholds_8_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("MQ Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")

MQ_8_FP <- ggplot(MQ_data_all_thresholds_8_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("MQ Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")


MQ_8_plot <- ggarrange(
  MQ_8_TP, MQ_8_FP,
  ncol = 2,
  nrow = 1,
  labels = c("A", "B")
)
MQ_8_plot_annoated <- annotate_figure(MQ_8_plot, top = text_grob("Chromosome 8", face = "bold", size = 14))
#print(MQ_8_plot_annoated)


########
# PARs # 
########
MQ_data_all_thresholds_PARs <- c()
for (val in MQ_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/males/per_filter/MQ/ALL_males_chrX_PARs_diploid_MQ_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chrXPARs"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  MQ_data_all_thresholds_PARs <- rbind(MQ_data_all_thresholds_PARs, dati)
}

MQ_data_all_thresholds_PARs_reordered <- MQ_data_all_thresholds_PARs %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=MQ_thresholds)) 

MQ_PARs_TP <- ggplot(MQ_data_all_thresholds_PARs_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("MQ Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")

MQ_PARs_FP <- ggplot(MQ_data_all_thresholds_PARs_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("MQ Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")


MQ_PARs_plot <- ggarrange(
  MQ_PARs_TP, MQ_PARs_FP,
  ncol = 2,
  nrow = 1,
  labels = c("C", "D")
)
MQ_PARs_plot_annoated <- annotate_figure(MQ_PARs_plot, top = text_grob("PARs", face = "bold", size = 14))
#print(MQ_PARs_plot_annoated)


########################
# Chromosome X nonPARs # 
########################
# TODO test out what this would look like with diploid called data
MQ_data_all_thresholds_X_nonPARs <- c()
for (val in MQ_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/males/per_filter/MQ/ALL_males_chrX_nonPARs_haploid_MQ_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chrXnonPARs"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  MQ_data_all_thresholds_X_nonPARs <- rbind(MQ_data_all_thresholds_X_nonPARs, dati)
}


MQ_data_all_thresholds_X_nonPARs_reordered <- MQ_data_all_thresholds_X_nonPARs %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=MQ_thresholds)) 

MQ_X_nonPARs_TP <- ggplot(MQ_data_all_thresholds_X_nonPARs_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("MQ Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")

MQ_X_nonPARs_FP <- ggplot(MQ_data_all_thresholds_X_nonPARs_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("MQ Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")


MQ_X_nonPARs_plot <- ggarrange(
  MQ_X_nonPARs_TP, MQ_X_nonPARs_FP,
  ncol = 2,
  nrow = 1,
  labels = c("E", "F")
)
MQ_X_nonPARs_plot_annoated <- annotate_figure(MQ_X_nonPARs_plot, top = text_grob("Chromosome X non-PARs", face = "bold", size = 14))
#print(MQ_X_nonPARs_plot_annoated)

################
# Y chromosome # 
################
MQ_data_all_thresholds_Y <- c()
for (val in MQ_thresholds) {
  #print(paste("../performance_metrics/EUR/males/per_filter/MQ/ALL_males_chrY_haploid_MQ_",
  #            val, "_golden_vs_called_performance_metrics.txt", sep = ""))
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/males/per_filter/MQ/ALL_males_chrY_haploid_MQ_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chrY"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  MQ_data_all_thresholds_Y <- rbind(MQ_data_all_thresholds_Y, dati)
}

MQ_data_all_thresholds_Y_reordered <- MQ_data_all_thresholds_Y %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=MQ_thresholds)) 

MQ_Y_TP <- ggplot(MQ_data_all_thresholds_Y_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("MQ Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")

MQ_Y_FP <- ggplot(MQ_data_all_thresholds_Y_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("MQ Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")


MQ_Y_plot <- ggarrange(
  MQ_Y_TP, MQ_Y_FP,
  ncol = 2,
  nrow = 1,
  labels = c("G", "H")
)
MQ_Y_plot_annoated <- annotate_figure(MQ_Y_plot, top = text_grob("Chromosome Y", face = "bold", size = 14))
#print(MQ_Y_plot_annoated)


#########
# mtDNA # 
#########
MQ_data_all_thresholds_M <- c()
for (val in MQ_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/males/per_filter/MQ/ALL_males_chrM_haploid_MQ_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chrM"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  MQ_data_all_thresholds_M <- rbind(MQ_data_all_thresholds_M, dati)
}


MQ_data_all_thresholds_M_reordered <- MQ_data_all_thresholds_M %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=MQ_thresholds)) 

MQ_M_TP <- ggplot(MQ_data_all_thresholds_M_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("MQ Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")

MQ_M_FP <- ggplot(MQ_data_all_thresholds_M_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(
    height = .1, width = .2, 
    #position=position_jitter(0.2), 
    color = "#440154FF"
  ) +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("MQ Threshold") +
  scale_y_continuous(limits = c(0,15)) +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")


MQ_M_plot <- ggarrange(
  MQ_M_TP, MQ_M_FP,
  ncol = 2,
  nrow = 1,
  labels = c("I", "J")
)
MQ_M_plot_annoated <- annotate_figure(MQ_M_plot, top = text_grob("mtDNA", face = "bold", size = 14))
#print(MQ_M_plot_annoated)


########
# PLOT #
########
pdf("../performance_metrics/EUR/males/per_filter/MQ/MQ_FP_TP_all_thresholds_haploid_boxplot.pdf",
    width = 6, height = 12)

ggarrange(
  MQ_8_plot_annoated, MQ_PARs_plot_annoated, 
  MQ_X_nonPARs_plot_annoated, MQ_Y_plot_annoated,
  MQ_M_plot_annoated,
  ncol = 1,
  nrow = 5
)

dev.off()



pdf("../performance_metrics/EUR/males/per_filter/MQ/MQ_FP_TP_all_thresholds_haploid_boxplot_2.pdf",
    width = 10, height = 8)

ggarrange(
  MQ_8_plot_annoated, MQ_PARs_plot_annoated, 
  MQ_X_nonPARs_plot_annoated, MQ_Y_plot_annoated,
  MQ_M_plot_annoated,
  ncol = 2,
  nrow = 3
)

dev.off()


#----------------------------------------------------------#
# NEW - prep for one large plot (potential main text figure)
#----------------------------------------------------------#
MQ_X_nonPARs_plot <- ggarrange(
  MQ_X_nonPARs_TP, MQ_X_nonPARs_FP,
  ncol = 2,
  nrow = 1
)
MQ_X_nonPARs_plot_annoated <- annotate_figure(MQ_X_nonPARs_plot, top = text_grob("Chromosome X non-PARs", face = "bold", size = 14))

MQ_Y_plot <- ggarrange(
  MQ_Y_TP, MQ_Y_FP,
  ncol = 2,
  nrow = 1
)
MQ_Y_plot_annoated <- annotate_figure(MQ_Y_plot, top = text_grob("Chromosome Y", face = "bold", size = 14))


MQ_plot_males <- ggarrange(
  MQ_X_nonPARs_plot_annoated, MQ_Y_plot_annoated,
  ncol = 2,
  nrow = 1
)



#-----------------------------------------------------------------------------#
# DP 
#-----------------------------------------------------------------------------#
DP_thresholds <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "15", "20")

########
# Chr8 #
########
DP_data_all_thresholds_8 <- c()
for (val in DP_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/males/per_filter/DP/ALL_males_chr8_autos_diploid_DP_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chr8"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  DP_data_all_thresholds_8 <- rbind(DP_data_all_thresholds_8, dati)
}


DP_data_all_thresholds_8_reordered <- DP_data_all_thresholds_8 %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=DP_thresholds)) 

DP_8_TP <- ggplot(DP_data_all_thresholds_8_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("DP Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 6, linetype="dotted", 
             color = "red")

DP_8_FP <- ggplot(DP_data_all_thresholds_8_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("DP Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 6, linetype="dotted", 
             color = "red")
#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")


DP_8_plot <- ggarrange(
  DP_8_TP, DP_8_FP,
  ncol = 2,
  nrow = 1,
  labels = c("A", "B")
)
DP_8_plot_annoated <- annotate_figure(DP_8_plot, top = text_grob("Chromosome 8", face = "bold", size = 14))
#print(DP_8_plot_annoated)


########
# PARs # 
########
DP_data_all_thresholds_PARs <- c()
for (val in DP_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/males/per_filter/DP/ALL_males_chrX_PARs_diploid_DP_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chrXPARs"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  DP_data_all_thresholds_PARs <- rbind(DP_data_all_thresholds_PARs, dati)
}

DP_data_all_thresholds_PARs_reordered <- DP_data_all_thresholds_PARs %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=DP_thresholds)) 

DP_PARs_TP <- ggplot(DP_data_all_thresholds_PARs_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("DP Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 6, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")

DP_PARs_FP <- ggplot(DP_data_all_thresholds_PARs_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("DP Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 6, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")


DP_PARs_plot <- ggarrange(
  DP_PARs_TP, DP_PARs_FP,
  ncol = 2,
  nrow = 1,
  labels = c("C", "D")
)
DP_PARs_plot_annoated <- annotate_figure(DP_PARs_plot, top = text_grob("PARs", face = "bold", size = 14))
#print(DP_PARs_plot_annoated)


########################
# Chromosome X nonPARs # 
########################
# TODO test out what this would look like with diploid called data
DP_data_all_thresholds_X_nonPARs <- c()
for (val in DP_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/males/per_filter/DP/ALL_males_chrX_nonPARs_haploid_DP_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chrXnonPARs"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  DP_data_all_thresholds_X_nonPARs <- rbind(DP_data_all_thresholds_X_nonPARs, dati)
}


DP_data_all_thresholds_X_nonPARs_reordered <- DP_data_all_thresholds_X_nonPARs %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=DP_thresholds)) 

DP_X_nonPARs_TP <- ggplot(DP_data_all_thresholds_X_nonPARs_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("DP Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 6, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")

DP_X_nonPARs_FP <- ggplot(DP_data_all_thresholds_X_nonPARs_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("DP Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 6, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")


DP_X_nonPARs_plot <- ggarrange(
  DP_X_nonPARs_TP, DP_X_nonPARs_FP,
  ncol = 2,
  nrow = 1,
  labels = c("E", "F")
)
DP_X_nonPARs_plot_annoated <- annotate_figure(DP_X_nonPARs_plot, top = text_grob("Chromosome X non-PARs", face = "bold", size = 14))
#print(DP_X_nonPARs_plot_annoated)

################
# Y chromosome # 
################
DP_data_all_thresholds_Y <- c()
for (val in DP_thresholds) {
  #print(paste("../performance_metrics/EUR/males/per_filter/DP/ALL_males_chrY_haploid_DP_",
  #            val, "_golden_vs_called_performance_metrics.txt", sep = ""))
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/males/per_filter/DP/ALL_males_chrY_haploid_DP_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chrY"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  DP_data_all_thresholds_Y <- rbind(DP_data_all_thresholds_Y, dati)
}

DP_data_all_thresholds_Y_reordered <- DP_data_all_thresholds_Y %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=DP_thresholds)) 

DP_Y_TP <- ggplot(DP_data_all_thresholds_Y_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("DP Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 6, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")

DP_Y_FP <- ggplot(DP_data_all_thresholds_Y_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("DP Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 6, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")


DP_Y_plot <- ggarrange(
  DP_Y_TP, DP_Y_FP,
  ncol = 2,
  nrow = 1,
  labels = c("G", "H")
)
DP_Y_plot_annoated <- annotate_figure(DP_Y_plot, top = text_grob("Chromosome Y", face = "bold", size = 14))
#print(DP_Y_plot_annoated)


#########
# mtDNA # 
#########
DP_data_all_thresholds_M <- c()
for (val in DP_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/males/per_filter/DP/ALL_males_chrM_haploid_DP_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chrM"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  DP_data_all_thresholds_M <- rbind(DP_data_all_thresholds_M, dati)
}


DP_data_all_thresholds_M_reordered <- DP_data_all_thresholds_M %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=DP_thresholds)) 

DP_M_TP <- ggplot(DP_data_all_thresholds_M_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("DP Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 6, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")

DP_M_FP <- ggplot(DP_data_all_thresholds_M_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(
    height = .1, width = .2, 
    #position=position_jitter(0.2), 
    color = "#440154FF"
  ) +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("DP Threshold") +
  scale_y_continuous(limits = c(0,15)) +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 6, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")


DP_M_plot <- ggarrange(
  DP_M_TP, DP_M_FP,
  ncol = 2,
  nrow = 1,
  labels = c("I", "J")
)
DP_M_plot_annoated <- annotate_figure(DP_M_plot, top = text_grob("mtDNA", face = "bold", size = 14))
#print(DP_M_plot_annoated)


########
# PLOT #
########
pdf("../performance_metrics/EUR/males/per_filter/DP/DP_FP_TP_all_thresholds_haploid_boxplot.pdf",
    width = 6, height = 12)

ggarrange(
  DP_8_plot_annoated, DP_PARs_plot_annoated, 
  DP_X_nonPARs_plot_annoated, DP_Y_plot_annoated,
  DP_M_plot_annoated,
  ncol = 1,
  nrow = 5
)

dev.off()



pdf("../performance_metrics/EUR/males/per_filter/DP/DP_FP_TP_all_thresholds_haploid_boxplot_2.pdf",
    width = 10, height = 8)

ggarrange(
  DP_8_plot_annoated, DP_PARs_plot_annoated, 
  DP_X_nonPARs_plot_annoated, DP_Y_plot_annoated,
  DP_M_plot_annoated,
  ncol = 2,
  nrow = 3
)

dev.off()


pdf("../performance_metrics/EUR/males/per_filter/DP/DP_FP_TP_all_thresholds_haploid_boxplot_3.pdf",
    width = 10, height = 6)

ggarrange(
  DP_8_plot_annoated, DP_PARs_plot_annoated, 
  DP_X_nonPARs_plot_annoated, DP_Y_plot_annoated,
  ncol = 2,
  nrow = 2
)

dev.off()


#-----------------------------------------------------------------------------#
# DP INFO
#-----------------------------------------------------------------------------#
DP_INFO_thresholds <- c("DP67and201", "DP23and69")

########
# Chr8 #
########
DP_INFO_data_all_thresholds_8 <- c()
for (val in DP_INFO_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/males/per_filter/DP_INFO/ALL_males_chr8_autos_diploid_DP_INFO_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chr8"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  DP_INFO_data_all_thresholds_8 <- rbind(DP_INFO_data_all_thresholds_8, dati)
}


DP_INFO_data_all_thresholds_8_reordered <- DP_INFO_data_all_thresholds_8 %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=DP_INFO_thresholds)) 

DP_INFO_8_TP <- ggplot(DP_INFO_data_all_thresholds_8_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("DP Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 1, linetype="dotted", 
             color = "red")

DP_INFO_8_FP <- ggplot(DP_INFO_data_all_thresholds_8_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("DP Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 1, linetype="dotted", 
             color = "red")
#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")


DP_INFO_8_plot <- ggarrange(
  DP_INFO_8_TP, DP_INFO_8_FP,
  ncol = 2,
  nrow = 1#,
  #labels = c("A", "B")
)
DP_INFO_8_plot_annoated <- annotate_figure(DP_INFO_8_plot, top = text_grob("Chromosome 8", face = "bold", size = 14))
#print(DP_INFO_8_plot_annoated)



########################
# Chromosome X nonPARs # 
########################
# TODO test out what this would look like with diploid called data
DP_INFO_data_all_thresholds_X_nonPARs <- c()
for (val in DP_INFO_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/males/per_filter/DP_INFO/ALL_males_chrX_nonPARs_haploid_DP_INFO_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chrXnonPARs"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  DP_INFO_data_all_thresholds_X_nonPARs <- rbind(DP_INFO_data_all_thresholds_X_nonPARs, dati)
}

# get means
mean(DP_INFO_data_all_thresholds_X_nonPARs[ which(DP_INFO_data_all_thresholds_X_nonPARs$Threshold=='DP67and201'), 2])
mean(DP_INFO_data_all_thresholds_X_nonPARs[ which(DP_INFO_data_all_thresholds_X_nonPARs$Threshold=='DP67and201'), 3])
mean(DP_INFO_data_all_thresholds_X_nonPARs[ which(DP_INFO_data_all_thresholds_X_nonPARs$Threshold=='DP67and201'), 4])

mean(DP_INFO_data_all_thresholds_X_nonPARs[ which(DP_INFO_data_all_thresholds_X_nonPARs$Threshold=='DP23and69'), 2])
mean(DP_INFO_data_all_thresholds_X_nonPARs[ which(DP_INFO_data_all_thresholds_X_nonPARs$Threshold=='DP23and69'), 3])
mean(DP_INFO_data_all_thresholds_X_nonPARs[ which(DP_INFO_data_all_thresholds_X_nonPARs$Threshold=='DP23and69'), 4])

DP_INFO_data_all_thresholds_X_nonPARs_reordered <- DP_INFO_data_all_thresholds_X_nonPARs %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=DP_INFO_thresholds)) 

DP_INFO_X_nonPARs_TP <- ggplot(DP_INFO_data_all_thresholds_X_nonPARs_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("DP Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 1, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")

DP_INFO_X_nonPARs_FP <- ggplot(DP_INFO_data_all_thresholds_X_nonPARs_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("DP Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 1, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")


DP_INFO_X_nonPARs_plot <- ggarrange(
  DP_INFO_X_nonPARs_TP, DP_INFO_X_nonPARs_FP,
  ncol = 2,
  nrow = 1#,
  #labels = c("E", "F")
)
DP_INFO_X_nonPARs_plot_annoated <- annotate_figure(DP_INFO_X_nonPARs_plot, top = text_grob("Chromosome X non-PARs", face = "bold", size = 14))
#print(DP_INFO_X_nonPARs_plot_annoated)

################
# Y chromosome # 
################
DP_INFO_data_all_thresholds_Y <- c()
for (val in DP_INFO_thresholds) {
  #print(paste("../performance_metrics/EUR/males/per_filter/DP_INFO/ALL_males_chrY_haploid_DP_INFO_",
  #            val, "_golden_vs_called_performance_metrics.txt", sep = ""))
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/males/per_filter/DP_INFO/ALL_males_chrY_haploid_DP_INFO_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chrY"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  DP_INFO_data_all_thresholds_Y <- rbind(DP_INFO_data_all_thresholds_Y, dati)
}

# get means
mean(DP_INFO_data_all_thresholds_Y[ which(DP_INFO_data_all_thresholds_Y$Threshold=='DP67and201'), 2])
mean(DP_INFO_data_all_thresholds_Y[ which(DP_INFO_data_all_thresholds_Y$Threshold=='DP67and201'), 3])
mean(DP_INFO_data_all_thresholds_Y[ which(DP_INFO_data_all_thresholds_Y$Threshold=='DP67and201'), 4])

mean(DP_INFO_data_all_thresholds_Y[ which(DP_INFO_data_all_thresholds_Y$Threshold=='DP23and69'), 2])
mean(DP_INFO_data_all_thresholds_Y[ which(DP_INFO_data_all_thresholds_Y$Threshold=='DP23and69'), 3])
mean(DP_INFO_data_all_thresholds_Y[ which(DP_INFO_data_all_thresholds_Y$Threshold=='DP23and69'), 4])

DP_INFO_data_all_thresholds_Y_reordered <- DP_INFO_data_all_thresholds_Y %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=DP_INFO_thresholds)) 

DP_INFO_Y_TP <- ggplot(DP_INFO_data_all_thresholds_Y_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("DP Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 1, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")

DP_INFO_Y_FP <- ggplot(DP_INFO_data_all_thresholds_Y_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("DP Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 1, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")


DP_INFO_Y_plot <- ggarrange(
  DP_INFO_Y_TP, DP_INFO_Y_FP,
  ncol = 2,
  nrow = 1#,
  #labels = c("G", "H")
)
DP_INFO_Y_plot_annoated <- annotate_figure(DP_INFO_Y_plot, top = text_grob("Chromosome Y", face = "bold", size = 14))
#print(DP_INFO_Y_plot_annoated)




########
# PLOT #
########
pdf("../performance_metrics/EUR/males/per_filter/DP_INFO/DP_INFO_FP_TP_all_thresholds_haploid_boxplot.pdf",
    width = 6, height = 11)

ggarrange(
  DP_INFO_8_plot_annoated, 
  DP_INFO_X_nonPARs_plot_annoated, DP_INFO_Y_plot_annoated,
  nrow = 3,
  ncol = 1
)

dev.off()




#----------------------------------------------------------#
# NEW - prep for one large plot (potential main text figure)
#----------------------------------------------------------#
DP_X_nonPARs_plot <- ggarrange(
  DP_X_nonPARs_TP, DP_X_nonPARs_FP,
  ncol = 2,
  nrow = 1
)
DP_X_nonPARs_plot_annoated <- annotate_figure(DP_X_nonPARs_plot, top = text_grob("Chromosome X non-PARs", face = "bold", size = 14))

DP_Y_plot <- ggarrange(
  DP_Y_TP, DP_Y_FP,
  ncol = 2,
  nrow = 1
)
DP_Y_plot_annoated <- annotate_figure(DP_Y_plot, top = text_grob("Chromosome Y", face = "bold", size = 14))


DP_plot_males <- ggarrange(
  DP_X_nonPARs_plot_annoated, DP_Y_plot_annoated,
  ncol = 2,
  nrow = 1
)


#-----------------------------------------------------------------------------#
# AN 
#-----------------------------------------------------------------------------#
AN_thresholds <- c("1", "2", "3", "4", "5", "10", "15", "20")

########
# Chr8 #
########
AN_data_all_thresholds_8 <- c()
for (val in AN_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/males/per_filter/AN/ALL_males_chr8_autos_diploid_AN_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chr8"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  AN_data_all_thresholds_8 <- rbind(AN_data_all_thresholds_8, dati)
}


AN_data_all_thresholds_8_reordered <- AN_data_all_thresholds_8 %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=AN_thresholds)) 

AN_8_TP <- ggplot(AN_data_all_thresholds_8_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("AN Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 8, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")

AN_8_FP <- ggplot(AN_data_all_thresholds_8_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("AN Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 8, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")


AN_8_plot <- ggarrange(
  AN_8_TP, AN_8_FP,
  ncol = 2,
  nrow = 1#,
  #labels = c("A", "B")
)
AN_8_plot_annoated <- annotate_figure(AN_8_plot, top = text_grob("Chromosome 8", face = "bold", size = 14))
#print(AN_8_plot_annoated)


########
# PARs # 
########
AN_data_all_thresholds_PARs <- c()
for (val in AN_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/males/per_filter/AN/ALL_males_chrX_PARs_diploid_AN_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chrXPARs"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  AN_data_all_thresholds_PARs <- rbind(AN_data_all_thresholds_PARs, dati)
}

AN_data_all_thresholds_PARs_reordered <- AN_data_all_thresholds_PARs %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=AN_thresholds)) 

AN_PARs_TP <- ggplot(AN_data_all_thresholds_PARs_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("AN Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 8, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")

AN_PARs_FP <- ggplot(AN_data_all_thresholds_PARs_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("AN Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 8, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")


AN_PARs_plot <- ggarrange(
  AN_PARs_TP, AN_PARs_FP,
  ncol = 2,
  nrow = 1,
  labels = c("C", "D")
)
AN_PARs_plot_annoated <- annotate_figure(AN_PARs_plot, top = text_grob("PARs", face = "bold", size = 14))
#print(AN_PARs_plot_annoated)


########################
# Chromosome X nonPARs # 
########################
# TODO test out what this would look like with diploid called data
AN_data_all_thresholds_X_nonPARs <- c()
for (val in AN_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/males/per_filter/AN/ALL_males_chrX_nonPARs_haploid_AN_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chrXnonPARs"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  AN_data_all_thresholds_X_nonPARs <- rbind(AN_data_all_thresholds_X_nonPARs, dati)
}


AN_data_all_thresholds_X_nonPARs_reordered <- AN_data_all_thresholds_X_nonPARs %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=AN_thresholds)) 

AN_X_nonPARs_TP <- ggplot(AN_data_all_thresholds_X_nonPARs_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("AN Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 8, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")

AN_X_nonPARs_FP <- ggplot(AN_data_all_thresholds_X_nonPARs_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("AN Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 8, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")


AN_X_nonPARs_plot <- ggarrange(
  AN_X_nonPARs_TP, AN_X_nonPARs_FP,
  ncol = 2,
  nrow = 1#,
  #labels = c("E", "F")
)
AN_X_nonPARs_plot_annoated <- annotate_figure(AN_X_nonPARs_plot, top = text_grob("Chromosome X non-PARs", face = "bold", size = 14))
#print(AN_X_nonPARs_plot_annoated)

################
# Y chromosome # 
################
AN_data_all_thresholds_Y <- c()
for (val in AN_thresholds) {
  #print(paste("../performance_metrics/EUR/males/per_filter/AN/ALL_males_chrY_haploid_AN_",
  #            val, "_golden_vs_called_performance_metrics.txt", sep = ""))
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/males/per_filter/AN/ALL_males_chrY_haploid_AN_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chrY"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  AN_data_all_thresholds_Y <- rbind(AN_data_all_thresholds_Y, dati)
}

AN_data_all_thresholds_Y_reordered <- AN_data_all_thresholds_Y %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=AN_thresholds)) 

AN_Y_TP <- ggplot(AN_data_all_thresholds_Y_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("AN Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 8, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")

AN_Y_FP <- ggplot(AN_data_all_thresholds_Y_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("AN Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 8, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")


AN_Y_plot <- ggarrange(
  AN_Y_TP, AN_Y_FP,
  ncol = 2,
  nrow = 1#,
  #labels = c("G", "H")
)
AN_Y_plot_annoated <- annotate_figure(AN_Y_plot, top = text_grob("Chromosome Y", face = "bold", size = 14))
#print(AN_Y_plot_annoated)


#########
# mtDNA # 
#########
AN_data_all_thresholds_M <- c()
for (val in AN_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/males/per_filter/AN/ALL_males_chrM_haploid_AN_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chrM"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  AN_data_all_thresholds_M <- rbind(AN_data_all_thresholds_M, dati)
}


AN_data_all_thresholds_M_reordered <- AN_data_all_thresholds_M %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=AN_thresholds)) 

AN_M_TP <- ggplot(AN_data_all_thresholds_M_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("AN Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 8, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")

AN_M_FP <- ggplot(AN_data_all_thresholds_M_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(
    height = .1, width = .2, 
    #position=position_jitter(0.2), 
    color = "#440154FF"
  ) +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("AN Threshold") +
  scale_y_continuous(limits = c(0,15)) +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 8, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")


AN_M_plot <- ggarrange(
  AN_M_TP, AN_M_FP,
  ncol = 2,
  nrow = 1,
  labels = c("I", "J")
)
AN_M_plot_annoated <- annotate_figure(AN_M_plot, top = text_grob("mtDNA", face = "bold", size = 14))
#print(AN_M_plot_annoated)


########
# PLOT #
########
pdf("../performance_metrics/EUR/males/per_filter/AN/AN_FP_TP_all_thresholds_haploid_boxplot.pdf",
    width = 6, height = 12)

ggarrange(
  AN_8_plot_annoated, AN_PARs_plot_annoated, 
  AN_X_nonPARs_plot_annoated, AN_Y_plot_annoated,
  AN_M_plot_annoated,
  ncol = 1,
  nrow = 5
)

dev.off()



pdf("../performance_metrics/EUR/males/per_filter/AN/AN_FP_TP_all_thresholds_haploid_boxplot_2.pdf",
    width = 10, height = 8)

ggarrange(
  AN_8_plot_annoated, AN_PARs_plot_annoated, 
  AN_X_nonPARs_plot_annoated, AN_Y_plot_annoated,
  AN_M_plot_annoated,
  ncol = 2,
  nrow = 3
)

dev.off()

pdf("../performance_metrics/EUR/males/per_filter/AN/AN_FP_TP_all_thresholds_haploid_boxplot_3.pdf",
    width = 6, height = 11)

ggarrange(
  AN_8_plot_annoated, 
  AN_X_nonPARs_plot_annoated, AN_Y_plot_annoated,
  ncol = 1,
  nrow = 3
)

dev.off()

#----------------------------------------------------------#
# NEW - prep for one large plot (potential main text figure)
#----------------------------------------------------------#
AN_X_nonPARs_plot <- ggarrange(
  AN_X_nonPARs_TP, AN_X_nonPARs_FP,
  ncol = 2,
  nrow = 1
)
AN_X_nonPARs_plot_annoated <- annotate_figure(AN_X_nonPARs_plot, top = text_grob("Chromosome X non-PARs", face = "bold", size = 14))

AN_Y_plot <- ggarrange(
  AN_Y_TP, AN_Y_FP,
  ncol = 2,
  nrow = 1
)
AN_Y_plot_annoated <- annotate_figure(AN_Y_plot, top = text_grob("Chromosome Y", face = "bold", size = 14))


AN_plot_males <- ggarrange(
  AN_X_nonPARs_plot_annoated, AN_Y_plot_annoated,
  ncol = 2,
  nrow = 1
)

#ggarrange(
#  QD_plot_males, MQ_plot_males, DP_plot_males, AN_plot_males,
#  ncol = 1,
#  nrow = 4
#)

###############################################################################
#
# FEMALE DATA #
#
###############################################################################
#-----------------------------------------------------------------------------#
# QD #
#-----------------------------------------------------------------------------#
# Plot the different thresholds
QD_thresholds <- c("1.0", "1.5","2.0", "12.0", "16.0", "20.0", "28.0")

########
# Chr8 #
########
QD_data_all_thresholds_8 <- c()
for (val in QD_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/females/per_filter/QD/ALL_females_chr8_autos_diploid_QD_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chr8"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  QD_data_all_thresholds_8 <- rbind(QD_data_all_thresholds_8, dati)
}


QD_data_all_thresholds_8_reordered <- QD_data_all_thresholds_8 %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=QD_thresholds)) 

QD_8_TP <- ggplot(QD_data_all_thresholds_8_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("QD Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")

QD_8_FP <- ggplot(QD_data_all_thresholds_8_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("QD Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")


QD_8_plot <- ggarrange(
  QD_8_TP, QD_8_FP,
  ncol = 2,
  nrow = 1#,
  #labels = c("A", "B")
)
QD_8_plot_annoated <- annotate_figure(QD_8_plot, top = text_grob("Chromosome 8", face = "bold", size = 14))
#print(QD_8_plot_annoated)

# to get means
summarySE(QD_data_all_thresholds_8_reordered, measurevar=c("TP"), groupvars=c("Threshold"))
summarySE(QD_data_all_thresholds_8_reordered, measurevar=c("FP"), groupvars=c("Threshold"))
summarySE(QD_data_all_thresholds_8_reordered, measurevar=c("FN"), groupvars=c("Threshold"))

########
# PARs # 
########
QD_data_all_thresholds_PARs <- c()
for (val in QD_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/females/per_filter/QD/ALL_females_chrX_PARs_diploid_QD_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chrXPARs"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  QD_data_all_thresholds_PARs <- rbind(QD_data_all_thresholds_PARs, dati)
}

QD_data_all_thresholds_PARs_reordered <- QD_data_all_thresholds_PARs %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=QD_thresholds)) 

QD_PARs_TP <- ggplot(QD_data_all_thresholds_PARs_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("QD Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")

QD_PARs_FP <- ggplot(QD_data_all_thresholds_PARs_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("QD Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")


QD_PARs_plot <- ggarrange(
  QD_PARs_TP, QD_PARs_FP,
  ncol = 2,
  nrow = 1,
  labels = c("C", "D")
)
QD_PARs_plot_annoated <- annotate_figure(QD_PARs_plot, top = text_grob("PARs", face = "bold", size = 14))
#print(QD_PARs_plot_annoated)


########################
# Chromosome X nonPARs # 
########################
# TODO test out what this would look like with diploid called data
QD_data_all_thresholds_X_nonPARs <- c()
for (val in QD_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/females/per_filter/QD/ALL_females_chrX_nonPARs_diploid_QD_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chrXnonPARs"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  QD_data_all_thresholds_X_nonPARs <- rbind(QD_data_all_thresholds_X_nonPARs, dati)
}


QD_data_all_thresholds_X_nonPARs_reordered <- QD_data_all_thresholds_X_nonPARs %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=QD_thresholds)) 

QD_X_nonPARs_TP <- ggplot(QD_data_all_thresholds_X_nonPARs_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("QD Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")

QD_X_nonPARs_FP <- ggplot(QD_data_all_thresholds_X_nonPARs_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("QD Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")


QD_X_nonPARs_plot <- ggarrange(
  QD_X_nonPARs_TP, QD_X_nonPARs_FP,
  ncol = 2,
  nrow = 1#,
  #labels = c("E", "F")
)
QD_X_nonPARs_plot_annoated <- annotate_figure(QD_X_nonPARs_plot, top = text_grob("Chromosome X non-PARs", face = "bold", size = 14))
#print(QD_X_nonPARs_plot_annoated)


#########
# mtDNA # 
#########
QD_data_all_thresholds_M <- c()
for (val in QD_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/females/per_filter/QD/ALL_females_chrM_haploid_QD_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chrM"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  QD_data_all_thresholds_M <- rbind(QD_data_all_thresholds_M, dati)
}


QD_data_all_thresholds_M_reordered <- QD_data_all_thresholds_M %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=QD_thresholds)) 

QD_M_TP <- ggplot(QD_data_all_thresholds_M_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("QD Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")

QD_M_FP <- ggplot(QD_data_all_thresholds_M_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(
    height = .1, width = .2, 
    #position=position_jitter(0.2), 
    color = "#440154FF"
  ) +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("QD Threshold") +
  scale_y_continuous(limits = c(0,15)) +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")


QD_M_plot <- ggarrange(
  QD_M_TP, QD_M_FP,
  ncol = 2,
  nrow = 1,
  labels = c("I", "J")
)
QD_M_plot_annoated <- annotate_figure(QD_M_plot, top = text_grob("mtDNA", face = "bold", size = 14))
#print(QD_M_plot_annoated)

summarySE(QD_data_all_thresholds_M_reordered, measurevar=c("TP"), groupvars=c("Threshold"))
summarySE(QD_data_all_thresholds_M_reordered, measurevar=c("FP"), groupvars=c("Threshold"))
summarySE(QD_data_all_thresholds_M_reordered, measurevar=c("FN"), groupvars=c("Threshold"))


########
# PLOT #
########
pdf("../performance_metrics/EUR/females/per_filter/QD/QD_FP_TP_all_thresholds_haploid_boxplot.pdf",
    width = 6, height = 10)

ggarrange(
  QD_8_plot_annoated, QD_PARs_plot_annoated, 
  QD_X_nonPARs_plot_annoated, QD_M_plot_annoated,
  ncol = 1,
  nrow = 4
)

dev.off()



pdf("../performance_metrics/EUR/females/per_filter/QD/QD_FP_TP_all_thresholds_haploid_boxplot_2.pdf",
    width = 10, height = 6)

ggarrange(
  QD_8_plot_annoated, QD_PARs_plot_annoated, 
  QD_X_nonPARs_plot_annoated, QD_M_plot_annoated,
  ncol = 2,
  nrow = 2
)

dev.off()

pdf("../performance_metrics/EUR/females/per_filter/QD/QD_FP_TP_all_thresholds_diploid_boxplot.pdf",
    width = 6, height = 7)

ggarrange(
  QD_8_plot_annoated, 
  QD_X_nonPARs_plot_annoated,
  nrow = 2,
  ncol = 1
)

dev.off()


#----------------------------------------------------------#
# NEW - prep for one large plot (potential main text figure)
#----------------------------------------------------------#
QD_X_nonPARs_plot <- ggarrange(
  QD_X_nonPARs_TP, QD_X_nonPARs_FP,
  ncol = 2,
  nrow = 1
)
QD_plot_females <- annotate_figure(QD_X_nonPARs_plot, top = text_grob("Chromosome X non-PARs", face = "bold", size = 14))


QD_plot_both <- ggarrange(
  QD_plot_females, QD_plot_males,
  ncol = 2,
  nrow = 1,
  widths = c(.5,1)
)


#-----------------------------------------------------------------------------#
# MQ #
#-----------------------------------------------------------------------------#

# Plot the different thresholds
MQ_thresholds <- c("20.0", "30.0", "40.0", "50.0", "60.0")

########
# Chr8 #
########
MQ_data_all_thresholds_8 <- c()
for (val in MQ_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/females/per_filter/MQ/ALL_females_chr8_autos_diploid_MQ_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chr8"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  MQ_data_all_thresholds_8 <- rbind(MQ_data_all_thresholds_8, dati)
}


MQ_data_all_thresholds_8_reordered <- MQ_data_all_thresholds_8 %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=MQ_thresholds)) 

MQ_8_TP <- ggplot(MQ_data_all_thresholds_8_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("MQ Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")

MQ_8_FP <- ggplot(MQ_data_all_thresholds_8_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("MQ Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")


MQ_8_plot <- ggarrange(
  MQ_8_TP, MQ_8_FP,
  ncol = 2,
  nrow = 1,
  labels = c("A", "B")
)
MQ_8_plot_annoated <- annotate_figure(MQ_8_plot, top = text_grob("Chromosome 8", face = "bold", size = 14))
#print(MQ_8_plot_annoated)


########
# PARs # 
########
MQ_data_all_thresholds_PARs <- c()
for (val in MQ_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/females/per_filter/MQ/ALL_females_chrX_PARs_diploid_MQ_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chrXPARs"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  MQ_data_all_thresholds_PARs <- rbind(MQ_data_all_thresholds_PARs, dati)
}

MQ_data_all_thresholds_PARs_reordered <- MQ_data_all_thresholds_PARs %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=MQ_thresholds)) 

MQ_PARs_TP <- ggplot(MQ_data_all_thresholds_PARs_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("MQ Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")

MQ_PARs_FP <- ggplot(MQ_data_all_thresholds_PARs_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("MQ Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")


MQ_PARs_plot <- ggarrange(
  MQ_PARs_TP, MQ_PARs_FP,
  ncol = 2,
  nrow = 1,
  labels = c("C", "D")
)
MQ_PARs_plot_annoated <- annotate_figure(MQ_PARs_plot, top = text_grob("PARs", face = "bold", size = 14))
#print(MQ_PARs_plot_annoated)


########################
# Chromosome X nonPARs # 
########################
# TODO test out what this would look like with diploid called data
MQ_data_all_thresholds_X_nonPARs <- c()
for (val in MQ_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/females/per_filter/MQ/ALL_females_chrX_nonPARs_diploid_MQ_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chrXnonPARs"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  MQ_data_all_thresholds_X_nonPARs <- rbind(MQ_data_all_thresholds_X_nonPARs, dati)
}


MQ_data_all_thresholds_X_nonPARs_reordered <- MQ_data_all_thresholds_X_nonPARs %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=MQ_thresholds)) 

MQ_X_nonPARs_TP <- ggplot(MQ_data_all_thresholds_X_nonPARs_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("MQ Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")

MQ_X_nonPARs_FP <- ggplot(MQ_data_all_thresholds_X_nonPARs_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("MQ Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")


MQ_X_nonPARs_plot <- ggarrange(
  MQ_X_nonPARs_TP, MQ_X_nonPARs_FP,
  ncol = 2,
  nrow = 1,
  labels = c("E", "F")
)
MQ_X_nonPARs_plot_annoated <- annotate_figure(MQ_X_nonPARs_plot, top = text_grob("Chromosome X non-PARs", face = "bold", size = 14))
#print(MQ_X_nonPARs_plot_annoated)



#########
# mtDNA # 
#########
MQ_data_all_thresholds_M <- c()
for (val in MQ_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/females/per_filter/MQ/ALL_females_chrM_haploid_MQ_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chrM"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  MQ_data_all_thresholds_M <- rbind(MQ_data_all_thresholds_M, dati)
}


MQ_data_all_thresholds_M_reordered <- MQ_data_all_thresholds_M %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=MQ_thresholds)) 

MQ_M_TP <- ggplot(MQ_data_all_thresholds_M_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("MQ Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")

MQ_M_FP <- ggplot(MQ_data_all_thresholds_M_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(
    height = .1, width = .2, 
    #position=position_jitter(0.2), 
    color = "#440154FF"
  ) +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("MQ Threshold") +
  scale_y_continuous(limits = c(0,15)) +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 3, linetype="dotted", 
             color = "red")


MQ_M_plot <- ggarrange(
  MQ_M_TP, MQ_M_FP,
  ncol = 2,
  nrow = 1,
  labels = c("I", "J")
)
MQ_M_plot_annoated <- annotate_figure(MQ_M_plot, top = text_grob("mtDNA", face = "bold", size = 14))
#print(MQ_M_plot_annoated)

summarySE(MQ_data_all_thresholds_M_reordered, measurevar=c("TP"), groupvars=c("Threshold"))
summarySE(MQ_data_all_thresholds_M_reordered, measurevar=c("FP"), groupvars=c("Threshold"))
summarySE(MQ_data_all_thresholds_M_reordered, measurevar=c("FN"), groupvars=c("Threshold"))


########
# PLOT #
########
pdf("../performance_metrics/EUR/females/per_filter/MQ/MQ_FP_TP_all_thresholds_haploid_boxplot.pdf",
    width = 6, height = 10)

ggarrange(
  MQ_8_plot_annoated, MQ_PARs_plot_annoated, 
  MQ_X_nonPARs_plot_annoated,
  MQ_M_plot_annoated,
  ncol = 1,
  nrow = 4
)

dev.off()



pdf("../performance_metrics/EUR/females/per_filter/MQ/MQ_FP_TP_all_thresholds_haploid_boxplot_2.pdf",
    width = 10, height = 6)

ggarrange(
  MQ_8_plot_annoated, MQ_PARs_plot_annoated, 
  MQ_X_nonPARs_plot_annoated,
  MQ_M_plot_annoated,
  ncol = 2,
  nrow = 2
)

dev.off()


#----------------------------------------------------------#
# NEW - prep for one large plot (potential main text figure)
#----------------------------------------------------------#
MQ_X_nonPARs_plot <- ggarrange(
  MQ_X_nonPARs_TP, MQ_X_nonPARs_FP,
  ncol = 2,
  nrow = 1
)
MQ_plot_females <- annotate_figure(MQ_X_nonPARs_plot, top = text_grob("Chromosome X non-PARs", face = "bold", size = 14))


MQ_plot_both <- ggarrange(
  MQ_plot_females, MQ_plot_males,
  ncol = 2,
  nrow = 1,
  widths = c(.5,1)
)


#-----------------------------------------------------------------------------#
# DP 
#-----------------------------------------------------------------------------#
DP_thresholds <- c("1", "2", "3", "4", "5", "10", "15", "20")
DP_thresholds <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "15", "20")

########
# Chr8 #
########
DP_data_all_thresholds_8 <- c()
for (val in DP_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/females/per_filter/DP/ALL_females_chr8_autos_diploid_DP_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chr8"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  DP_data_all_thresholds_8 <- rbind(DP_data_all_thresholds_8, dati)
}


DP_data_all_thresholds_8_reordered <- DP_data_all_thresholds_8 %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=DP_thresholds)) 

DP_8_TP <- ggplot(DP_data_all_thresholds_8_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("DP Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 6, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")

DP_8_FP <- ggplot(DP_data_all_thresholds_8_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("DP Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 6, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")


DP_8_plot <- ggarrange(
  DP_8_TP, DP_8_FP,
  ncol = 2,
  nrow = 1,
  labels = c("A", "B")
)
DP_8_plot_annoated <- annotate_figure(DP_8_plot, top = text_grob("Chromosome 8", face = "bold", size = 14))
#print(DP_8_plot_annoated)


########
# PARs # 
########
DP_data_all_thresholds_PARs <- c()
for (val in DP_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/females/per_filter/DP/ALL_females_chrX_PARs_diploid_DP_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chrXPARs"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  DP_data_all_thresholds_PARs <- rbind(DP_data_all_thresholds_PARs, dati)
}

DP_data_all_thresholds_PARs_reordered <- DP_data_all_thresholds_PARs %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=DP_thresholds)) 

DP_PARs_TP <- ggplot(DP_data_all_thresholds_PARs_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("DP Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 6, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")

DP_PARs_FP <- ggplot(DP_data_all_thresholds_PARs_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("DP Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 6, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")


DP_PARs_plot <- ggarrange(
  DP_PARs_TP, DP_PARs_FP,
  ncol = 2,
  nrow = 1,
  labels = c("C", "D")
)
DP_PARs_plot_annoated <- annotate_figure(DP_PARs_plot, top = text_grob("PARs", face = "bold", size = 14))
#print(DP_PARs_plot_annoated)


########################
# Chromosome X nonPARs # 
########################
# TODO test out what this would look like with diploid called data
DP_data_all_thresholds_X_nonPARs <- c()
for (val in DP_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/females/per_filter/DP/ALL_females_chrX_nonPARs_diploid_DP_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chrXnonPARs"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  DP_data_all_thresholds_X_nonPARs <- rbind(DP_data_all_thresholds_X_nonPARs, dati)
}


DP_data_all_thresholds_X_nonPARs_reordered <- DP_data_all_thresholds_X_nonPARs %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=DP_thresholds)) 

DP_X_nonPARs_TP <- ggplot(DP_data_all_thresholds_X_nonPARs_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("DP Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 6, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")

DP_X_nonPARs_FP <- ggplot(DP_data_all_thresholds_X_nonPARs_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("DP Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 6, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")


DP_X_nonPARs_plot <- ggarrange(
  DP_X_nonPARs_TP, DP_X_nonPARs_FP,
  ncol = 2,
  nrow = 1,
  labels = c("E", "F")
)
DP_X_nonPARs_plot_annoated <- annotate_figure(DP_X_nonPARs_plot, top = text_grob("Chromosome X non-PARs", face = "bold", size = 14))
#print(DP_X_nonPARs_plot_annoated)



#########
# mtDNA # 
#########
DP_data_all_thresholds_M <- c()
for (val in DP_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/females/per_filter/DP/ALL_females_chrM_haploid_DP_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chrM"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  DP_data_all_thresholds_M <- rbind(DP_data_all_thresholds_M, dati)
}


DP_data_all_thresholds_M_reordered <- DP_data_all_thresholds_M %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=DP_thresholds)) 

DP_M_TP <- ggplot(DP_data_all_thresholds_M_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("DP Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 6, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")

DP_M_FP <- ggplot(DP_data_all_thresholds_M_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(
    height = .1, width = .2, 
    #position=position_jitter(0.2), 
    color = "#440154FF"
  ) +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("DP Threshold") +
  scale_y_continuous(limits = c(0,15)) +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 6, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")


DP_M_plot <- ggarrange(
  DP_M_TP, DP_M_FP,
  ncol = 2,
  nrow = 1,
  labels = c("I", "J")
)
DP_M_plot_annoated <- annotate_figure(DP_M_plot, top = text_grob("mtDNA", face = "bold", size = 14))
#print(DP_M_plot_annoated)

summarySE(DP_data_all_thresholds_M_reordered, measurevar=c("TP"), groupvars=c("Threshold"))
summarySE(DP_data_all_thresholds_M_reordered, measurevar=c("FP"), groupvars=c("Threshold"))
summarySE(DP_data_all_thresholds_M_reordered, measurevar=c("FN"), groupvars=c("Threshold"))


########
# PLOT #
########
pdf("../performance_metrics/EUR/females/per_filter/DP/DP_FP_TP_all_thresholds_haploid_boxplot.pdf",
    width = 6, height = 10)

ggarrange(
  DP_8_plot_annoated, DP_PARs_plot_annoated, 
  DP_X_nonPARs_plot_annoated,
  DP_M_plot_annoated,
  ncol = 1,
  nrow = 4
)

dev.off()



pdf("../performance_metrics/EUR/females/per_filter/DP/DP_FP_TP_all_thresholds_haploid_boxplot_2.pdf",
    width = 10, height = 6)

ggarrange(
  DP_8_plot_annoated, DP_PARs_plot_annoated, 
  DP_X_nonPARs_plot_annoated,
  DP_M_plot_annoated,
  ncol = 2,
  nrow = 2
)

dev.off()

pdf("../performance_metrics/EUR/females/per_filter/DP/DP_FP_TP_all_thresholds_haploid_boxplot_3.pdf",
    width = 10, height = 6)

ggarrange(
  DP_8_plot_annoated, DP_PARs_plot_annoated, 
  DP_X_nonPARs_plot_annoated,
  ncol = 2,
  nrow = 2
)

dev.off()


#----------------------------------------------------------#
# NEW - prep for one large plot (potential main text figure)
#----------------------------------------------------------#
DP_X_nonPARs_plot <- ggarrange(
  DP_X_nonPARs_TP, DP_X_nonPARs_FP,
  ncol = 2,
  nrow = 1
)
DP_plot_females <- annotate_figure(DP_X_nonPARs_plot, top = text_grob("Chromosome X non-PARs", face = "bold", size = 14))


DP_plot_both <- ggarrange(
  DP_plot_females, DP_plot_males,
  ncol = 2,
  nrow = 1,
  widths = c(.5,1)
)



#-----------------------------------------------------------------------------#
# AN 
#-----------------------------------------------------------------------------#
AN_thresholds <- c("1", "2", "3", "4", "5", "10", "15", "20")

########
# Chr8 #
########
AN_data_all_thresholds_8 <- c()
for (val in AN_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/females/per_filter/AN/ALL_females_chr8_autos_diploid_AN_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chr8"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  AN_data_all_thresholds_8 <- rbind(AN_data_all_thresholds_8, dati)
}


AN_data_all_thresholds_8_reordered <- AN_data_all_thresholds_8 %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=AN_thresholds)) 

AN_8_TP <- ggplot(AN_data_all_thresholds_8_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("AN Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 8, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")

AN_8_FP <- ggplot(AN_data_all_thresholds_8_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("AN Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 8, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")


AN_8_plot <- ggarrange(
  AN_8_TP, AN_8_FP,
  ncol = 2,
  nrow = 1#,
  #labels = c("A", "B")
)
AN_8_plot_annoated <- annotate_figure(AN_8_plot, top = text_grob("Chromosome 8", face = "bold", size = 14))
#print(AN_8_plot_annoated)


########
# PARs # 
########
AN_data_all_thresholds_PARs <- c()
for (val in AN_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/females/per_filter/AN/ALL_females_chrX_PARs_diploid_AN_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chrXPARs"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  AN_data_all_thresholds_PARs <- rbind(AN_data_all_thresholds_PARs, dati)
}

AN_data_all_thresholds_PARs_reordered <- AN_data_all_thresholds_PARs %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=AN_thresholds)) 

AN_PARs_TP <- ggplot(AN_data_all_thresholds_PARs_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("AN Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 8, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")

AN_PARs_FP <- ggplot(AN_data_all_thresholds_PARs_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("AN Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 8, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")


AN_PARs_plot <- ggarrange(
  AN_PARs_TP, AN_PARs_FP,
  ncol = 2,
  nrow = 1,
  labels = c("C", "D")
)
AN_PARs_plot_annoated <- annotate_figure(AN_PARs_plot, top = text_grob("PARs", face = "bold", size = 14))
#print(AN_PARs_plot_annoated)


########################
# Chromosome X nonPARs # 
########################
# TODO test out what this would look like with diploid called data
AN_data_all_thresholds_X_nonPARs <- c()
for (val in AN_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/females/per_filter/AN/ALL_females_chrX_nonPARs_diploid_AN_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chrXnonPARs"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  AN_data_all_thresholds_X_nonPARs <- rbind(AN_data_all_thresholds_X_nonPARs, dati)
}


AN_data_all_thresholds_X_nonPARs_reordered <- AN_data_all_thresholds_X_nonPARs %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=AN_thresholds)) 

AN_X_nonPARs_TP <- ggplot(AN_data_all_thresholds_X_nonPARs_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("AN Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 8, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")

AN_X_nonPARs_FP <- ggplot(AN_data_all_thresholds_X_nonPARs_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("AN Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 8, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")


AN_X_nonPARs_plot <- ggarrange(
  AN_X_nonPARs_TP, AN_X_nonPARs_FP,
  ncol = 2,
  nrow = 1#,
  #labels = c("E", "F")
)
AN_X_nonPARs_plot_annoated <- annotate_figure(AN_X_nonPARs_plot, top = text_grob("Chromosome X non-PARs", face = "bold", size = 14))
#print(AN_X_nonPARs_plot_annoated)



#########
# mtDNA # 
#########
AN_data_all_thresholds_M <- c()
for (val in AN_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/females/per_filter/AN/ALL_females_chrM_haploid_AN_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chrM"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  AN_data_all_thresholds_M <- rbind(AN_data_all_thresholds_M, dati)
}


AN_data_all_thresholds_M_reordered <- AN_data_all_thresholds_M %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=AN_thresholds)) 

AN_M_TP <- ggplot(AN_data_all_thresholds_M_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("AN Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 8, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")

AN_M_FP <- ggplot(AN_data_all_thresholds_M_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(
    height = .1, width = .2, 
    #position=position_jitter(0.2), 
    color = "#440154FF"
  ) +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("AN Threshold") +
  scale_y_continuous(limits = c(0,15)) +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 8, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")


AN_M_plot <- ggarrange(
  AN_M_TP, AN_M_FP,
  ncol = 2,
  nrow = 1,
  labels = c("I", "J")
)
AN_M_plot_annoated <- annotate_figure(AN_M_plot, top = text_grob("mtDNA", face = "bold", size = 14))
#print(AN_M_plot_annoated)

summarySE(AN_data_all_thresholds_M_reordered, measurevar=c("TP"), groupvars=c("Threshold"))
summarySE(AN_data_all_thresholds_M_reordered, measurevar=c("FP"), groupvars=c("Threshold"))
summarySE(AN_data_all_thresholds_M_reordered, measurevar=c("FN"), groupvars=c("Threshold"))


########
# PLOT #
########
pdf("../performance_metrics/EUR/females/per_filter/AN/AN_FP_TP_all_thresholds_haploid_boxplot.pdf",
    width = 6, height = 10)

ggarrange(
  AN_8_plot_annoated, AN_PARs_plot_annoated, 
  AN_X_nonPARs_plot_annoated,
  AN_M_plot_annoated,
  ncol = 1,
  nrow = 4
)

dev.off()



pdf("../performance_metrics/EUR/females/per_filter/AN/AN_FP_TP_all_thresholds_haploid_boxplot_2.pdf",
    width = 10, height = 6)

ggarrange(
  AN_8_plot_annoated, AN_PARs_plot_annoated, 
  AN_X_nonPARs_plot_annoated,
  AN_M_plot_annoated,
  ncol = 2,
  nrow = 2
)

dev.off()

pdf("../performance_metrics/EUR/females/per_filter/AN/AN_FP_TP_all_thresholds_diploid_boxplot.pdf",
    width = 6, height = 7)

ggarrange(
  AN_8_plot_annoated, 
  AN_X_nonPARs_plot_annoated,
  nrow = 2,
  ncol = 1
)

dev.off()


#----------------------------------------------------------#
# NEW - prep for one large plot (potential main text figure)
#----------------------------------------------------------#
AN_X_nonPARs_plot <- ggarrange(
  AN_X_nonPARs_TP, AN_X_nonPARs_FP,
  ncol = 2,
  nrow = 1
)
AN_plot_females <- annotate_figure(AN_X_nonPARs_plot, top = text_grob("Chromosome X non-PARs", face = "bold", size = 14))


AN_plot_both <- ggarrange(
  AN_plot_females, AN_plot_males,
  ncol = 2,
  nrow = 1,
  widths = c(.5,1)
)



#-----------------------------------------------------------------------------#
# DP INFO
#-----------------------------------------------------------------------------#
DP_INFO_thresholds <- c("DP67and201", "DP23and69")

########
# Chr8 #
########
DP_INFO_data_all_thresholds_8 <- c()
for (val in DP_INFO_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/females/per_filter/DP_INFO/ALL_females_chr8_autos_diploid_DP_INFO_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chr8"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  DP_INFO_data_all_thresholds_8 <- rbind(DP_INFO_data_all_thresholds_8, dati)
}


DP_INFO_data_all_thresholds_8_reordered <- DP_INFO_data_all_thresholds_8 %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=DP_INFO_thresholds)) 

DP_INFO_8_TP <- ggplot(DP_INFO_data_all_thresholds_8_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("DP Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 1, linetype="dotted", 
             color = "red")

DP_INFO_8_FP <- ggplot(DP_INFO_data_all_thresholds_8_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("DP Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 1, linetype="dotted", 
             color = "red")
#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")


DP_INFO_8_plot <- ggarrange(
  DP_INFO_8_TP, DP_INFO_8_FP,
  ncol = 2,
  nrow = 1#,
  #labels = c("A", "B")
)
DP_INFO_8_plot_annoated <- annotate_figure(DP_INFO_8_plot, top = text_grob("Chromosome 8", face = "bold", size = 14))
#print(DP_INFO_8_plot_annoated)



########################
# Chromosome X nonPARs # 
########################
# TODO test out what this would look like with diploid called data
DP_INFO_data_all_thresholds_X_nonPARs <- c()
for (val in DP_INFO_thresholds) {
  # read in file
  dati <- read.table(paste("../performance_metrics/EUR/females/per_filter/DP_INFO/ALL_females_chrX_nonPARs_diploid_DP_INFO_",
                           val, "_golden_vs_called_performance_metrics.txt", sep = ""), 
                     header = T)
  # add chromosome number column
  dati$Chromosome <- "chrXnonPARs"
  # add threshold value
  dati$Threshold <- val
  # add to larger DF
  DP_INFO_data_all_thresholds_X_nonPARs <- rbind(DP_INFO_data_all_thresholds_X_nonPARs, dati)
}

# get means
mean(DP_INFO_data_all_thresholds_X_nonPARs[ which(DP_INFO_data_all_thresholds_X_nonPARs$Threshold=='DP67and201'), 2])
mean(DP_INFO_data_all_thresholds_X_nonPARs[ which(DP_INFO_data_all_thresholds_X_nonPARs$Threshold=='DP67and201'), 3])
mean(DP_INFO_data_all_thresholds_X_nonPARs[ which(DP_INFO_data_all_thresholds_X_nonPARs$Threshold=='DP67and201'), 4])

mean(DP_INFO_data_all_thresholds_X_nonPARs[ which(DP_INFO_data_all_thresholds_X_nonPARs$Threshold=='DP23and69'), 2])
mean(DP_INFO_data_all_thresholds_X_nonPARs[ which(DP_INFO_data_all_thresholds_X_nonPARs$Threshold=='DP23and69'), 3])
mean(DP_INFO_data_all_thresholds_X_nonPARs[ which(DP_INFO_data_all_thresholds_X_nonPARs$Threshold=='DP23and69'), 4])

DP_INFO_data_all_thresholds_X_nonPARs_reordered <- DP_INFO_data_all_thresholds_X_nonPARs %>%
  arrange(TP) %>%
  mutate(Threshold = factor(Threshold, 
                            levels=DP_INFO_thresholds)) 

DP_INFO_X_nonPARs_TP <- ggplot(DP_INFO_data_all_thresholds_X_nonPARs_reordered, aes(x=Threshold, y=TP)) + 
  geom_boxplot(color = "#21908CFF") +
  geom_jitter(position=position_jitter(0.2), color = "#21908CFF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("True Positives") +
  xlab("DP Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#21908CFF", face="bold")) +
  geom_vline(xintercept = 1, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")

DP_INFO_X_nonPARs_FP <- ggplot(DP_INFO_data_all_thresholds_X_nonPARs_reordered, aes(x=Threshold, y=FP)) + 
  geom_boxplot(color = "#440154FF") +
  geom_jitter(position=position_jitter(0.2), color = "#440154FF") +
  theme_classic() +
  #ggtitle("Chromosome 8") + 
  ylab("False Positives") +
  xlab("DP Threshold") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1),
    axis.title.y = element_text(color = "#440154FF", face="bold")) +
  geom_vline(xintercept = 1, linetype="dotted", 
             color = "red")#+
#geom_vline(xintercept = 3, linetype="dotted", 
#color = "red")


DP_INFO_X_nonPARs_plot <- ggarrange(
  DP_INFO_X_nonPARs_TP, DP_INFO_X_nonPARs_FP,
  ncol = 2,
  nrow = 1#,
  #labels = c("E", "F")
)
DP_INFO_X_nonPARs_plot_annoated <- annotate_figure(DP_INFO_X_nonPARs_plot, top = text_grob("Chromosome X non-PARs", face = "bold", size = 14))
#print(DP_INFO_X_nonPARs_plot_annoated)




########
# PLOT #
########
pdf("../performance_metrics/EUR/females/per_filter/DP_INFO/DP_INFO_FP_TP_all_thresholds_diploid_boxplot.pdf",
    width = 6, height = 7)

ggarrange(
  DP_INFO_8_plot_annoated, 
  DP_INFO_X_nonPARs_plot_annoated,
  nrow = 2,
  ncol = 1
)

dev.off()




#----------------------------------------------------------#
# NEW - prep for one large plot (potential main text figure)
#----------------------------------------------------------#
DP_X_nonPARs_plot <- ggarrange(
  DP_X_nonPARs_TP, DP_X_nonPARs_FP,
  ncol = 2,
  nrow = 1
)
DP_X_nonPARs_plot_annoated <- annotate_figure(DP_X_nonPARs_plot, top = text_grob("Chromosome X non-PARs", face = "bold", size = 14))

DP_Y_plot <- ggarrange(
  DP_Y_TP, DP_Y_FP,
  ncol = 2,
  nrow = 1
)
DP_Y_plot_annoated <- annotate_figure(DP_Y_plot, top = text_grob("Chromosome Y", face = "bold", size = 14))


DP_plot_females <- ggarrange(
  DP_X_nonPARs_plot_annoated, DP_Y_plot_annoated,
  ncol = 2,
  nrow = 1
)



# PLOT ALL #
pdf("../plots/test_filters_on_one.pdf", 
    width = 11)
ggarrange(
  QD_plot_both, MQ_plot_both, DP_plot_both, AN_plot_both,
  ncol = 1,
  nrow = 4
)
dev.off()
