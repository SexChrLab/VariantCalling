# Potential Figures - Figure 1
# This will be a multi panel plot with the difference in FNs between default 
# and SCC alignment


#----------------#
# Load libraries # 
#----------------#
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggpubr)
require(scales)


#--------------#
# Read in data #
#--------------#
# variables to subset
myvars <- c("V1", "V2", "V3", "average")

# Per site data #
metric = "TP"

# a. Females chrX: non PARs + PARs 

females_chrX_nonPARs_metric_scc <-
  read.table(
    paste("../performance_metrics/EUR/windows/scc/females/all_females_chrX_nonPARs_scc_diploid_golden_vs_called_", 
    #paste("../performance_metrics/EUR/windows/scc_gencode/females/all_females_chrX_nonPARs_scc_diploid_golden_vs_called_", 
          metric, "_counts_50kb_windows.txt", sep = ""),
    header = F)

females_chrX_nonPARs_metric_scc$average <- rowMeans(females_chrX_nonPARs_metric_scc[,4:13])

females_chrX_nonPARs_metric_scc_sub <- females_chrX_nonPARs_metric_scc[myvars]
colnames(females_chrX_nonPARs_metric_scc_sub) <- c("chr", "start", "end", "SCC")

females_chrX_nonPARs_metric_default <-
  read.table(
    paste("../performance_metrics/EUR/windows/default/females/all_females_chrX_nonPARs_default_diploid_golden_vs_called_", 
          metric, "_counts_50kb_windows.txt", sep = ""),
    header = F)

females_chrX_nonPARs_metric_default$average <- rowMeans(females_chrX_nonPARs_metric_default[,4:13])

females_chrX_nonPARs_metric_default_sub <- females_chrX_nonPARs_metric_default[myvars]
colnames(females_chrX_nonPARs_metric_default_sub) <- c("chr", "start", "end", "Default")


females_chrX_nonPARs_metric_merge <- merge(females_chrX_nonPARs_metric_default_sub, females_chrX_nonPARs_metric_scc_sub, 
                                           by = c("chr", "start", "end"))
females_chrX_nonPARs_metric_merge$FNdifference <- females_chrX_nonPARs_metric_merge$SCC - females_chrX_nonPARs_metric_merge$Default

# remove PAR entries from data frame
females_chrX_nonPARs_metric_merge_edit <-females_chrX_nonPARs_metric_merge[!(females_chrX_nonPARs_metric_merge$start<=2781479 | females_chrX_nonPARs_metric_merge$start>=155701383),]

females_chrX_PARs_metric_scc <-
  read.table(
    paste("../performance_metrics/EUR/windows/scc/females/all_females_chrX_PARs_scc_diploid_golden_vs_called_", 
    #paste("../performance_metrics/EUR/windows/scc_gencode/females/all_females_chrX_PARs_scc_diploid_golden_vs_called_", 
          metric, "_counts_50kb_windows.txt", sep = ""),
    header = F)

females_chrX_PARs_metric_scc$average <- rowMeans(females_chrX_PARs_metric_scc[,4:13])

females_chrX_PARs_metric_scc_sub <- females_chrX_PARs_metric_scc[myvars]
colnames(females_chrX_PARs_metric_scc_sub) <- c("chr", "start", "end", "SCC")

females_chrX_PARs_metric_default <-
  read.table(
    paste("../performance_metrics/EUR/windows/default/females/all_females_chrX_PARs_default_diploid_golden_vs_called_", 
          metric, "_counts_50kb_windows.txt", sep = ""),
    header = F)

females_chrX_PARs_metric_default$average <- rowMeans(females_chrX_PARs_metric_default[,4:13])

females_chrX_PARs_metric_default_sub <- females_chrX_PARs_metric_default[myvars]
colnames(females_chrX_PARs_metric_default_sub) <- c("chr", "start", "end", "Default")


females_chrX_PARs_metric_merge <- merge(females_chrX_PARs_metric_default_sub, females_chrX_PARs_metric_scc_sub, 
                                        by = c("chr", "start", "end"))
females_chrX_PARs_metric_merge$FNdifference <- females_chrX_PARs_metric_merge$SCC - females_chrX_PARs_metric_merge$Default

# remove nonPAR entries from data frame
females_chrX_PARs_metric_merge_edit <-females_chrX_PARs_metric_merge[(females_chrX_PARs_metric_merge$start<=2781479 | females_chrX_PARs_metric_merge$start>=155701383),]

# b. Males chrX: non PARs + PARs 

chrX_nonPARs_metric_scc <-
  read.table(
    paste("../performance_metrics/EUR/windows/scc/males/all_males_chrX_nonPARs_scc_haploid_golden_vs_called_", 
    #paste("../performance_metrics/EUR/windows/scc_gencode/males/all_males_chrX_nonPARs_scc_haploid_golden_vs_called_", 
          metric, "_counts_50kb_windows.txt", sep = ""),
    header = F)

chrX_nonPARs_metric_scc$average <- rowMeans(chrX_nonPARs_metric_scc[,4:13])

chrX_nonPARs_metric_scc_sub <- chrX_nonPARs_metric_scc[myvars]
colnames(chrX_nonPARs_metric_scc_sub) <- c("chr", "start", "end", "SCC")

chrX_nonPARs_metric_default <-
  read.table(
    paste("../performance_metrics/EUR/windows/default/males/all_males_chrX_nonPARs_default_haploid_golden_vs_called_", 
          metric, "_counts_50kb_windows.txt", sep = ""),
    header = F)

chrX_nonPARs_metric_default$average <- rowMeans(chrX_nonPARs_metric_default[,4:13])

chrX_nonPARs_metric_default_sub <- chrX_nonPARs_metric_default[myvars]
colnames(chrX_nonPARs_metric_default_sub) <- c("chr", "start", "end", "Default")


chrX_nonPARs_metric_merge <- merge(chrX_nonPARs_metric_default_sub, chrX_nonPARs_metric_scc_sub, 
                                   by = c("chr", "start", "end"))
chrX_nonPARs_metric_merge$FNdifference <- chrX_nonPARs_metric_merge$SCC - chrX_nonPARs_metric_merge$Default

# remove PAR entries from data frame
chrX_nonPARs_metric_merge_edit <-chrX_nonPARs_metric_merge[!(chrX_nonPARs_metric_merge$start<=2781479 | chrX_nonPARs_metric_merge$start>=155701383),]

chrX_PARs_metric_scc <-
  read.table(
    paste("../performance_metrics/EUR/windows/scc/males/all_males_chrX_PARs_scc_diploid_golden_vs_called_", 
    #paste("../performance_metrics/EUR/windows/scc_gencode/males/all_males_chrX_PARs_scc_diploid_golden_vs_called_", 
          metric, "_counts_50kb_windows.txt", sep = ""),
    header = F)

chrX_PARs_metric_scc$average <- rowMeans(chrX_PARs_metric_scc[,4:13])

chrX_PARs_metric_scc_sub <- chrX_PARs_metric_scc[myvars]
colnames(chrX_PARs_metric_scc_sub) <- c("chr", "start", "end", "SCC")

chrX_PARs_metric_default <-
  read.table(
    paste("../performance_metrics/EUR/windows/default/males/all_males_chrX_PARs_default_diploid_golden_vs_called_", 
          metric, "_counts_50kb_windows.txt", sep = ""),
    header = F)

chrX_PARs_metric_default$average <- rowMeans(chrX_PARs_metric_default[,4:13])

chrX_PARs_metric_default_sub <- chrX_PARs_metric_default[myvars]
colnames(chrX_PARs_metric_default_sub) <- c("chr", "start", "end", "Default")


chrX_PARs_metric_merge <- merge(chrX_PARs_metric_default_sub, chrX_PARs_metric_scc_sub, 
                                by = c("chr", "start", "end"))
chrX_PARs_metric_merge$FNdifference <- chrX_PARs_metric_merge$SCC - chrX_PARs_metric_merge$Default

# remove nonPAR entries from data frame
chrX_PARs_metric_merge_edit <-chrX_PARs_metric_merge[(chrX_PARs_metric_merge$start<=2781479 | chrX_PARs_metric_merge$start>=155701383),]

# b. Males chrY 
chrY_metric_scc <-
  read.table(
    paste("../performance_metrics/EUR/windows/scc/males/all_males_chrY_scc_haploid_golden_vs_called_",
    #paste("../performance_metrics/EUR/windows/scc_gencode/males/all_males_chrY_scc_haploid_golden_vs_called_",
          metric, "_counts_50kb_windows.txt", sep = ""),
    header = F)

chrY_metric_scc$average <- rowMeans(chrY_metric_scc[,4:13])

chrY_metric_scc_sub <- chrY_metric_scc[myvars]
colnames(chrY_metric_scc_sub) <- c("chr", "start", "end", "SCC")

chrY_metric_default <-
  read.table(
    paste("../performance_metrics/EUR/windows/default/males/all_males_chrY_default_haploid_golden_vs_called_",
          metric, "_counts_50kb_windows.txt", sep =  ""),
    header = F)

chrY_metric_default$average <- rowMeans(chrY_metric_default[,4:13])

chrY_metric_default_sub <- chrY_metric_default[myvars]
colnames(chrY_metric_default_sub) <- c("chr", "start", "end", "Default")


chrY_metric_merge <- merge(chrY_metric_default_sub, chrY_metric_scc_sub, 
                                by = c("chr", "start", "end"))
chrY_metric_merge$FNdifference <- chrY_metric_merge$SCC - chrY_metric_merge$Default


# remove PAR entries from data frame
chrY_metric_merge_merge_edit <-chrY_metric_merge[!(chrY_metric_merge$start<=2781479 | chrY_metric_merge$start>=56887903),]

#--------------#
# Plot figures # 
#--------------#
# Y chromosome # 
yaxisvaluemax_Y <- (max(chrY_metric_merge$FNdifference)) + 1
yaxisvaluemin_Y <- (min(chrY_metric_merge$FNdifference)) - 1 

male_chrY <- ggplot(data=chrY_metric_merge_merge_edit, aes(x=start, y=FNdifference)) +
  geom_point(size = 0.1) +
  #geom_line() +
  #scale_x_continuous(limits = c(2781479,56887903), labels = comma) +
  scale_x_continuous(limits = c(0,57227415), labels = comma) +
  scale_y_continuous(limits = c(yaxisvaluemin_Y,yaxisvaluemax_Y)) +
  theme_classic() +
  # Amplicons
  annotate("rect", xmin=6235111, xmax=6532906, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=7574481, xmax=10266944, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=13984473, xmax=14058733, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=15875093, xmax=15905214, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=16159795, xmax=16425965, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=17456266, xmax=18870334, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=21305953, xmax=26415435, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  # XTRs
  annotate("rect", xmin=3050044, xmax=6235111, ymin=-Inf, ymax=Inf, alpha=0.5, fill="darkmagenta") +
  annotate("rect", xmin=6532906, xmax=6748713, ymin=-Inf, ymax=Inf, alpha=0.5, fill="darkmagenta") +
  # PARs
  annotate("rect", xmin=10001, xmax=2781479, ymin=-Inf, ymax=Inf, alpha=0.5, fill="darkolivegreen4") +
  annotate("rect", xmin=56887903, xmax=57217415, ymin=-Inf, ymax=Inf, alpha=0.5, fill="darkolivegreen4") +
  xlab("Chromosome Y (GRCh38)") +
  ylab(paste("Difference in ", metric, "s", sep = "")) +
  #ggtitle("Simulated males (XY) - Default vs SCC alignment") +
  ggtitle("SCC vs Default alignment") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


yaxisvaluemax_Y_s_d <- (max(c(chrY_metric_merge$Default, chrY_metric_merge$SCC))) 
yaxisvaluemin_Y_s_d <- (min(c(chrY_metric_merge$Default, chrY_metric_merge$SCC))) 

male_chrY_default <- ggplot(data=chrY_metric_merge_merge_edit, aes(x=start, y=Default)) +
  geom_point(size = 0.1) +
  #geom_line() +
  #scale_x_continuous(limits = c(2781479,56887903), labels = comma) +
  scale_x_continuous(limits = c(0,57227415), labels = comma) +
  scale_y_continuous(limits = c(yaxisvaluemin_Y_s_d,yaxisvaluemax_Y_s_d)) +
  theme_classic() +
  # Amplicons
  annotate("rect", xmin=6235111, xmax=6532906, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=7574481, xmax=10266944, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=13984473, xmax=14058733, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=15875093, xmax=15905214, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=16159795, xmax=16425965, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=17456266, xmax=18870334, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=21305953, xmax=26415435, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  # XTRs
  annotate("rect", xmin=3050044, xmax=6235111, ymin=-Inf, ymax=Inf, alpha=0.5, fill="darkmagenta") +
  annotate("rect", xmin=6532906, xmax=6748713, ymin=-Inf, ymax=Inf, alpha=0.5, fill="darkmagenta") +
  # PARs
  annotate("rect", xmin=10001, xmax=2781479, ymin=-Inf, ymax=Inf, alpha=0.5, fill="darkolivegreen4") +
  annotate("rect", xmin=56887903, xmax=57217415, ymin=-Inf, ymax=Inf, alpha=0.5, fill="darkolivegreen4") +
  xlab("Chromosome Y (GRCh38)") +
  ylab(paste(metric, "s", sep = "")) +
  #ggtitle("Simulated males (XY) - Default vs SCC alignment") +
  ggtitle("Default alignment") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


male_chrY_scc <- ggplot(data=chrY_metric_merge_merge_edit, aes(x=start, y=SCC)) +
  geom_point(size = 0.1) +
  #geom_line() +
  #scale_x_continuous(limits = c(2781479,56887903), labels = comma) +
  scale_x_continuous(limits = c(0,57227415), labels = comma) +
  scale_y_continuous(limits = c(yaxisvaluemin_Y_s_d,yaxisvaluemax_Y_s_d)) +
  theme_classic() +
  # Amplicons
  annotate("rect", xmin=6235111, xmax=6532906, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=7574481, xmax=10266944, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=13984473, xmax=14058733, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=15875093, xmax=15905214, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=16159795, xmax=16425965, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=17456266, xmax=18870334, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=21305953, xmax=26415435, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  # XTRs
  annotate("rect", xmin=3050044, xmax=6235111, ymin=-Inf, ymax=Inf, alpha=0.5, fill="darkmagenta") +
  annotate("rect", xmin=6532906, xmax=6748713, ymin=-Inf, ymax=Inf, alpha=0.5, fill="darkmagenta") +
  # PARs
  annotate("rect", xmin=10001, xmax=2781479, ymin=-Inf, ymax=Inf, alpha=0.5, fill="darkolivegreen4") +
  annotate("rect", xmin=56887903, xmax=57217415, ymin=-Inf, ymax=Inf, alpha=0.5, fill="darkolivegreen4") +
  xlab("Chromosome Y (GRCh38)") +
  ylab(paste(metric, "s", sep = "")) +
  #ggtitle("Simulated males (XY) - Default vs SCC alignment") +
  ggtitle("SCC alignment") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))



# X chromosome #
# Females
female_yaxisvaluemax_X <- (max(c(females_chrX_nonPARs_metric_merge_edit$FNdifference, females_chrX_PARs_metric_merge_edit$FNdifference))) + 1
female_yaxisvaluemin_X <- (min(c(females_chrX_nonPARs_metric_merge_edit$FNdifference, females_chrX_PARs_metric_merge_edit$FNdifference))) - 1 

female_chrX <- ggplot(data=females_chrX_nonPARs_metric_merge_edit, aes(x=start, y=FNdifference)) +
  geom_point(size = 0.1) +
  geom_point(data = females_chrX_PARs_metric_merge_edit, color = "black", size = 0.1) +
  scale_x_continuous(limits = c(0,156040895), labels = comma) +
  scale_y_continuous(limits = c(female_yaxisvaluemin_X,female_yaxisvaluemax_X)) +
  theme_classic() +
  xlab("Chromosome X (GRCh38)") +
  ylab(paste("Difference in ", metric, "s", sep = "")) +
  #ggtitle("Simulated females (XX) - Default vs SCC alignment") +
  ggtitle("SCC vs Default alignment") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  # PARs
  annotate("rect", xmin=10001, xmax=2781479, ymin=female_yaxisvaluemin_X, ymax=female_yaxisvaluemax_X, alpha=0.5, fill="darkolivegreen4") +
  annotate("rect", xmin=155701383, xmax=156030895, ymin=female_yaxisvaluemin_X, ymax=female_yaxisvaluemax_X, alpha=0.5, fill="darkolivegreen4") +
  # Amplicons
  annotate("rect", xmin=48343310, xmax=48434604, ymin=female_yaxisvaluemin_X, ymax=female_yaxisvaluemax_X, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=52032464, xmax=52223402, ymin=female_yaxisvaluemin_X, ymax=female_yaxisvaluemax_X, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=55437684, xmax=55547739, ymin=female_yaxisvaluemin_X, ymax=female_yaxisvaluemax_X, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=71674267, xmax=71835832, ymin=female_yaxisvaluemin_X, ymax=female_yaxisvaluemax_X, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=72721314, xmax=73105236, ymin=female_yaxisvaluemin_X, ymax=female_yaxisvaluemax_X, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=101563740, xmax=101648990, ymin=female_yaxisvaluemin_X, ymax=female_yaxisvaluemax_X, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=102180805, xmax=102519463, ymin=female_yaxisvaluemin_X, ymax=female_yaxisvaluemax_X, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=103940531, xmax=104117650, ymin=female_yaxisvaluemin_X, ymax=female_yaxisvaluemax_X, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=49157254, xmax=49162986, ymin=female_yaxisvaluemin_X, ymax=female_yaxisvaluemax_X, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=52511000, xmax=52511492, ymin=female_yaxisvaluemin_X, ymax=female_yaxisvaluemax_X, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=52489200, xmax=52501394, ymin=female_yaxisvaluemin_X, ymax=female_yaxisvaluemax_X, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=52520509, xmax=52520537, ymin=female_yaxisvaluemin_X, ymax=female_yaxisvaluemax_X, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=139870917, xmax=139871081, ymin=female_yaxisvaluemin_X, ymax=female_yaxisvaluemax_X, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=63185573, xmax=63193371, ymin=female_yaxisvaluemin_X, ymax=female_yaxisvaluemax_X, alpha=0.5, fill="steelblue") +
  # XTRs
  annotate("rect", xmin=89140830, xmax=93428068, ymin=female_yaxisvaluemin_X, ymax=female_yaxisvaluemax_X, alpha=0.5, fill="darkmagenta") #+

# Default and SCC separately 
female_yaxisvaluemax_X_d_s <- (max(c(females_chrX_nonPARs_metric_merge_edit$Default, females_chrX_PARs_metric_merge_edit$SCC))) 
female_yaxisvaluemin_X_d_s <- (min(c(females_chrX_nonPARs_metric_merge_edit$Default, females_chrX_PARs_metric_merge_edit$SCC)))  

female_chrX_default <- ggplot(data=females_chrX_nonPARs_metric_merge_edit, aes(x=start, y=Default)) +
  geom_point(size = 0.1) +
  geom_point(data = females_chrX_PARs_metric_merge_edit, color = "black", size = 0.1) +
  scale_x_continuous(limits = c(0,156040895), labels = comma) +
  scale_y_continuous(limits = c(female_yaxisvaluemin_X_d_s,female_yaxisvaluemax_X_d_s)) +
  theme_classic() +
  xlab("Chromosome X (GRCh38)") +
  ylab(paste(metric, "s", sep = "")) +
  #ggtitle("Simulated females (XX) - Default alignment") +
  ggtitle("Default alignment") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  # PARs
  annotate("rect", xmin=10001, xmax=2781479, ymin=-Inf, ymax=Inf, alpha=0.5, fill="darkolivegreen4") +
  annotate("rect", xmin=155701383, xmax=156030895, ymin=-Inf, ymax=Inf, alpha=0.5, fill="darkolivegreen4") +
  # Amplicons
  annotate("rect", xmin=48343310, xmax=48434604, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=52032464, xmax=52223402, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=55437684, xmax=55547739, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=71674267, xmax=71835832, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=72721314, xmax=73105236, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=101563740, xmax=101648990, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=102180805, xmax=102519463, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=103940531, xmax=104117650, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=49157254, xmax=49162986, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=52511000, xmax=52511492, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=52489200, xmax=52501394, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=52520509, xmax=52520537, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=139870917, xmax=139871081, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=63185573, xmax=63193371, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  # XTRs
  annotate("rect", xmin=89140830, xmax=93428068, ymin=-Inf, ymax=Inf, alpha=0.5, fill="darkmagenta") #+

female_chrX_scc <- ggplot(data=females_chrX_nonPARs_metric_merge_edit, aes(x=start, y=SCC)) +
  geom_point(size = 0.1) +
  geom_point(data = females_chrX_PARs_metric_merge_edit, color = "black", size = 0.1) +
  scale_x_continuous(limits = c(0,156040895), labels = comma) +
  scale_y_continuous(limits = c(female_yaxisvaluemin_X_d_s,female_yaxisvaluemax_X_d_s)) +
  theme_classic() +
  xlab("Chromosome X (GRCh38)") +
  ylab(paste(metric, "s", sep = "")) +
  #ggtitle("Simulated females (XX) - SCC alignment") +
  ggtitle("SCC alignment") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  # PARs
  annotate("rect", xmin=10001, xmax=2781479, ymin=-Inf, ymax=Inf, alpha=0.5, fill="darkolivegreen4") +
  annotate("rect", xmin=155701383, xmax=156030895, ymin=-Inf, ymax=Inf, alpha=0.5, fill="darkolivegreen4") +
  # Amplicons
  annotate("rect", xmin=48343310, xmax=48434604, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=52032464, xmax=52223402, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=55437684, xmax=55547739, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=71674267, xmax=71835832, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=72721314, xmax=73105236, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=101563740, xmax=101648990, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=102180805, xmax=102519463, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=103940531, xmax=104117650, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=49157254, xmax=49162986, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=52511000, xmax=52511492, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=52489200, xmax=52501394, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=52520509, xmax=52520537, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=139870917, xmax=139871081, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=63185573, xmax=63193371, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  # XTRs
  annotate("rect", xmin=89140830, xmax=93428068, ymin=-Inf, ymax=Inf, alpha=0.5, fill="darkmagenta") #+


# Males
male_yaxisvaluemax_X <- (max(c(chrX_nonPARs_metric_merge_edit$FNdifference, chrX_PARs_metric_merge_edit$FNdifference))) + 1
male_yaxisvaluemin_X <- (min(c(chrX_nonPARs_metric_merge_edit$FNdifference, chrX_PARs_metric_merge_edit$FNdifference))) - 1 

male_chrX <- ggplot(data=chrX_nonPARs_metric_merge_edit, aes(x=start, y=FNdifference)) +
  geom_point(size = 0.1) +
  geom_point(data = chrX_PARs_metric_merge_edit, color = "black", size = 0.1) +
  scale_x_continuous(limits = c(0,156040895), labels = comma) +
  scale_y_continuous(limits = c(male_yaxisvaluemin_X,male_yaxisvaluemax_X)) +
  theme_classic() +
  xlab("Chromosome X (GRCh38)") +
  ylab(paste("Difference in ", metric, "s", sep = "")) +
  #ggtitle("Simulated males (XX) - Default vs SCC alignment") +
  ggtitle("SCC vs Default alignment") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  # PARs
  annotate("rect", xmin=10001, xmax=2781479, ymin=male_yaxisvaluemin_X, ymax=male_yaxisvaluemax_X, alpha=0.5, fill="darkolivegreen4") +
  annotate("rect", xmin=155701383, xmax=156030895, ymin=male_yaxisvaluemin_X, ymax=male_yaxisvaluemax_X, alpha=0.5, fill="darkolivegreen4") +
  # Amplicons
  annotate("rect", xmin=48343310, xmax=48434604, ymin=male_yaxisvaluemin_X, ymax=male_yaxisvaluemax_X, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=52032464, xmax=52223402, ymin=male_yaxisvaluemin_X, ymax=male_yaxisvaluemax_X, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=55437684, xmax=55547739, ymin=male_yaxisvaluemin_X, ymax=male_yaxisvaluemax_X, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=71674267, xmax=71835832, ymin=male_yaxisvaluemin_X, ymax=male_yaxisvaluemax_X, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=72721314, xmax=73105236, ymin=male_yaxisvaluemin_X, ymax=male_yaxisvaluemax_X, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=101563740, xmax=101648990, ymin=male_yaxisvaluemin_X, ymax=male_yaxisvaluemax_X, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=102180805, xmax=102519463, ymin=male_yaxisvaluemin_X, ymax=male_yaxisvaluemax_X, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=103940531, xmax=104117650, ymin=male_yaxisvaluemin_X, ymax=male_yaxisvaluemax_X, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=49157254, xmax=49162986, ymin=male_yaxisvaluemin_X, ymax=male_yaxisvaluemax_X, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=52511000, xmax=52511492, ymin=male_yaxisvaluemin_X, ymax=male_yaxisvaluemax_X, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=52489200, xmax=52501394, ymin=male_yaxisvaluemin_X, ymax=male_yaxisvaluemax_X, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=52520509, xmax=52520537, ymin=male_yaxisvaluemin_X, ymax=male_yaxisvaluemax_X, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=139870917, xmax=139871081, ymin=male_yaxisvaluemin_X, ymax=male_yaxisvaluemax_X, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=63185573, xmax=63193371, ymin=male_yaxisvaluemin_X, ymax=male_yaxisvaluemax_X, alpha=0.5, fill="steelblue") +
  # XTRs
  annotate("rect", xmin=89140830, xmax=93428068, ymin=male_yaxisvaluemin_X, ymax=male_yaxisvaluemax_X, alpha=0.5, fill="darkmagenta") #+

# Default and SCC separately 
male_yaxisvaluemax_X_d_s <- (max(c(chrX_nonPARs_metric_merge_edit$Default, chrX_PARs_metric_merge_edit$SCC))) 
male_yaxisvaluemin_X_d_s <- (min(c(chrX_nonPARs_metric_merge_edit$Default, chrX_PARs_metric_merge_edit$SCC)))  

male_chrX_default <- ggplot(data=chrX_nonPARs_metric_merge_edit, aes(x=start, y=Default)) +
  geom_point(size = 0.1) +
  geom_point(data = chrX_PARs_metric_merge_edit, color = "black", size = 0.1) +
  scale_x_continuous(limits = c(0,156040895), labels = comma) +
  scale_y_continuous(limits = c(male_yaxisvaluemin_X_d_s,male_yaxisvaluemax_X_d_s)) +
  theme_classic() +
  xlab("Chromosome X (GRCh38)") +
  ylab(paste(metric, "s", sep = "")) +
  #ggtitle("Simulated males (XX) - Default alignment") +
  ggtitle("Default alignment") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  # PARs
  annotate("rect", xmin=10001, xmax=2781479, ymin=-Inf, ymax=Inf, alpha=0.5, fill="darkolivegreen4") +
  annotate("rect", xmin=155701383, xmax=156030895, ymin=-Inf, ymax=Inf, alpha=0.5, fill="darkolivegreen4") +
  # Amplicons
  annotate("rect", xmin=48343310, xmax=48434604, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=52032464, xmax=52223402, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=55437684, xmax=55547739, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=71674267, xmax=71835832, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=72721314, xmax=73105236, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=101563740, xmax=101648990, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=102180805, xmax=102519463, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=103940531, xmax=104117650, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=49157254, xmax=49162986, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=52511000, xmax=52511492, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=52489200, xmax=52501394, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=52520509, xmax=52520537, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=139870917, xmax=139871081, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=63185573, xmax=63193371, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  # XTRs
  annotate("rect", xmin=89140830, xmax=93428068, ymin=-Inf, ymax=Inf, alpha=0.5, fill="darkmagenta") #+

male_chrX_scc <- ggplot(data=chrX_nonPARs_metric_merge_edit, aes(x=start, y=SCC)) +
  geom_point(size = 0.1) +
  geom_point(data = chrX_PARs_metric_merge_edit, color = "black", size = 0.1) +
  scale_x_continuous(limits = c(0,156040895), labels = comma) +
  scale_y_continuous(limits = c(male_yaxisvaluemin_X_d_s,male_yaxisvaluemax_X_d_s)) +
  theme_classic() +
  xlab("Chromosome X (GRCh38)") +
  ylab(paste(metric, "s", sep = "")) +
  #ggtitle("Simulated males (XX) - SCC alignment") +
  ggtitle("SCC alignment") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  # PARs
  annotate("rect", xmin=10001, xmax=2781479, ymin=-Inf, ymax=Inf, alpha=0.5, fill="darkolivegreen4") +
  annotate("rect", xmin=155701383, xmax=156030895, ymin=-Inf, ymax=Inf, alpha=0.5, fill="darkolivegreen4") +
  # Amplicons
  annotate("rect", xmin=48343310, xmax=48434604, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=52032464, xmax=52223402, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=55437684, xmax=55547739, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=71674267, xmax=71835832, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=72721314, xmax=73105236, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=101563740, xmax=101648990, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=102180805, xmax=102519463, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=103940531, xmax=104117650, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=49157254, xmax=49162986, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=52511000, xmax=52511492, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=52489200, xmax=52501394, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=52520509, xmax=52520537, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=139870917, xmax=139871081, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  annotate("rect", xmin=63185573, xmax=63193371, ymin=-Inf, ymax=Inf, alpha=0.5, fill="steelblue") +
  # XTRs
  annotate("rect", xmin=89140830, xmax=93428068, ymin=-Inf, ymax=Inf, alpha=0.5, fill="darkmagenta") #+


# Plot chromosomes
df <- data.frame()
yplotmin <- 0
yplotmax <- 2.2
chrmax <- 2
chrX <- ggplot(df) + geom_point() + 
  xlim(0,156040895) + ylim(c(yplotmin,yplotmax)) +
  theme_classic(base_line_size = 0) +
  annotate("rect", xmin=89140830, xmax=93428068, ymin=yplotmin, ymax=chrmax, fill="darkmagenta") +
  
  annotate("rect", xmin=10001, xmax=2781479, ymin=yplotmin, ymax=chrmax, fill="darkolivegreen4") +
  annotate("rect", xmin=155701383, xmax=156030895, ymin=yplotmin, ymax=chrmax, fill="darkolivegreen4") +
  
  annotate("rect", xmin=48343310, xmax=48434604, ymin=yplotmin, ymax=chrmax, fill="steelblue") +
  annotate("rect", xmin=52032464, xmax=52223402, ymin=yplotmin, ymax=chrmax, fill="steelblue") +
  annotate("rect", xmin=55437684, xmax=55547739, ymin=yplotmin, ymax=chrmax, fill="steelblue") +
  annotate("rect", xmin=71674267, xmax=71835832, ymin=yplotmin, ymax=chrmax, fill="steelblue") +
  annotate("rect", xmin=72721314, xmax=73105236, ymin=yplotmin, ymax=chrmax, fill="steelblue") +
  annotate("rect", xmin=101563740, xmax=101648990, ymin=yplotmin, ymax=chrmax, fill="steelblue") +
  annotate("rect", xmin=102180805, xmax=102519463, ymin=yplotmin, ymax=chrmax, fill="steelblue") +
  annotate("rect", xmin=103940531, xmax=104117650, ymin=yplotmin, ymax=chrmax, fill="steelblue") +
  annotate("rect", xmin=49157254, xmax=49162986, ymin=yplotmin, ymax=chrmax, fill="steelblue") +
  annotate("rect", xmin=52511000, xmax=52511492, ymin=yplotmin, ymax=chrmax, fill="steelblue") +
  annotate("rect", xmin=52489200, xmax=52501394, ymin=yplotmin, ymax=chrmax, fill="steelblue") +
  annotate("rect", xmin=52520509, xmax=52520537, ymin=yplotmin, ymax=chrmax, fill="steelblue") +
  annotate("rect", xmin=139870917, xmax=139871081, ymin=yplotmin, ymax=chrmax, fill="steelblue") +
  annotate("rect", xmin=63185573, xmax=63193371, ymin=yplotmin, ymax=chrmax, fill="steelblue") +
  
  annotate("rect", xmin=0, xmax=156040895, ymin=yplotmin, ymax=chrmax, alpha=0.3, color = "black", fill = "grey") +

  #annotate("text", x=1500000, y=yplotmax, label = "PAR1", color = "darkolivegreen4", fontface =2, size = 3) +
  #annotate("text", x=155701383, y=yplotmax, label = "PAR2", color = "darkolivegreen4", fontface =2, size = 3) +
  #annotate("text", x=91284449, y=yplotmax, label = "XTR", color = "darkmagenta", fontface =2, size = 3) +
  #annotate("text", x=71674267, y=yplotmax, label = "Ampliconic", color = "steelblue", fontface =2) +
  ggtitle("Chromosome X") +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) 


chrY <- ggplot(df) + geom_point() + 
  xlim(0,57227415) + ylim(c(yplotmin,yplotmax)) +
  theme_classic(base_line_size = 0) +
  # Amplicons
  annotate("rect", xmin=6235111, xmax=6532906, ymin=yplotmin, ymax=chrmax, fill="steelblue") +
  annotate("rect", xmin=7574481, xmax=10266944, ymin=yplotmin, ymax=chrmax, fill="steelblue") +
  annotate("rect", xmin=13984473, xmax=14058733, ymin=yplotmin, ymax=chrmax, fill="steelblue") +
  annotate("rect", xmin=15875093, xmax=15905214, ymin=yplotmin, ymax=chrmax, fill="steelblue") +
  annotate("rect", xmin=16159795, xmax=16425965, ymin=yplotmin, ymax=chrmax, fill="steelblue") +
  annotate("rect", xmin=17456266, xmax=18870334, ymin=yplotmin, ymax=chrmax, fill="steelblue") +
  annotate("rect", xmin=21305953, xmax=26415435, ymin=yplotmin, ymax=chrmax, fill="steelblue") +
  # XTRs
  annotate("rect", xmin=3050044, xmax=6235111, ymin=yplotmin, ymax=chrmax, fill="darkmagenta") +
  annotate("rect", xmin=6532906, xmax=6748713, ymin=yplotmin, ymax=chrmax, fill="darkmagenta") +
  # PARs
  annotate("rect", xmin=10001, xmax=2781479, ymin=yplotmin, ymax=chrmax, fill="darkolivegreen4") +
  annotate("rect", xmin=56887903, xmax=57217415, ymin=yplotmin, ymax=chrmax, fill="darkolivegreen4") +
  

  annotate("rect", xmin=0, xmax=57217415, ymin=yplotmin, ymax=chrmax, alpha=0.3, color = "black", fill = "grey") +
  
  #annotate("text", x=0, y=yplotmax, label = "PAR1", color = "darkolivegreen4", fontface =2) +
  #annotate("text", x=56887903, y=yplotmax, label = "PAR2", color = "darkolivegreen4", fontface =2) +
  #annotate("text", x=4899379, y=yplotmax, label = "XTR", color = "darkmagenta", fontface =2) +
  #annotate("text", x=71674267, y=yplotmax, label = "Ampliconic", color = "steelblue", fontface =2) +
  ggtitle("Chromosome Y") +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) 
emptyplot <- ggplot(df) + theme_classic(base_line_size = 0) + theme(panel.grid = element_blank(),
                                                  axis.title = element_blank(),
                                                  axis.text = element_blank(),
                                                  axis.ticks = element_blank())
pdf("../plots/test_Figure_01.pdf",
    width = 9, height = 13)
ggarrange(
  female_chrX, chrX, male_chrX,
  chrY, male_chrY,
  ncol = 1,
  nrow = 5,
  align = c("v"),
  labels = c("A", "B", "C", "D", "E"),
  heights = c(1,.4,1,.4,1)
)

dev.off()


plotX <- ggarrange(
  female_chrX, chrX, male_chrX,
  ncol = 1,
  nrow = 3,
  align = c("v"),
  labels = c("A", "B", "C"),
  heights = c(1,.4,1)
)

plotY <- ggarrange(
  chrY, male_chrY, emptyplot,
  ncol = 1,
  nrow = 3,
  align = c("v"),
  labels = c("D", "E"),
  heights = c(.4,1,1)
)

pdf("../plots/test_Figure_01_rearranged.pdf",
    width = 16, height = 9)
ggarrange(
  plotX, plotY,
  ncol = 2,
  nrow = 1,
  align = c("v")
)
dev.off()




female_X_plot <- ggarrange(
  chrX, female_chrX_default, female_chrX_scc, female_chrX, 
  ncol = 1,
  nrow = 4,
  align = c("v"),
  #labels = c("A", "B", "C", "D"),
  heights = c(.4,1,1,1)
)

female_X_plot_annoated <- annotate_figure(female_X_plot, top = text_grob("Simulated females (XX)", face = "bold", size = 14))
#print(female_X_plot_annoated)


male_X_plot <- ggarrange(
  chrX, male_chrX_default, male_chrX_scc, male_chrX, 
  ncol = 1,
  nrow = 4,
  align = c("v"),
  #labels = c("E", "F", "G", "H"),
  heights = c(.4,1,1,1)
)

male_X_plot_annoated <- annotate_figure(male_X_plot, top = text_grob("Simulated males (XY)", face = "bold", size = 14))
#print(male_X_plot_annoated)

pdf("../plots/test_Figure_01_chrX_TPs_SCC_Default.pdf",
    width = 10, height = 7)
ggarrange(
  female_X_plot_annoated, male_X_plot_annoated, 
  ncol = 2,
  nrow = 1,
  align = c("v"),
  labels = c("A", "B")
)
dev.off()


male_Y_plot <- ggarrange(
  chrY, male_chrY_default, male_chrY_scc, male_chrY, 
  ncol = 1,
  nrow = 4,
  align = c("v"),
  #labels = c("E", "F", "G", "H"),
  heights = c(.4,1,1,1)
)

male_Y_plot_annoated <- annotate_figure(male_Y_plot, top = text_grob("Simulated males (XY)", face = "bold", size = 14))

pdf("../plots/supplemental_Figure_chrY_TPs_SCC_Default.pdf",
    width = 5, height = 7)
print(male_Y_plot_annoated)
dev.off()


pdf("../plots/test_Figure_01_chrX_TPs_SCC_Default_chrY_NEW.pdf",
    width = 10, height = 7)

ggarrange(
  female_X_plot, male_X_plot, male_Y_plot,
  ncol = 3,
  nrow = 1,
  align = c("v"),
  labels = c("A", "B", "C")
)
dev.off()

ggarrange(
  female_X_plot_annoated, male_X_plot_annoated, male_Y_plot_annoated,
  ncol = 3,
  nrow = 1,
  align = c("v"),
  labels = c("A", "B", "C")
)
