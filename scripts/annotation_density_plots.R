# Angela Oill
# 2020-04-28
# Plot density curves for different subsets of sites

# 1. Load libraries # 
library(ggplot2)
library(ggpubr)
library(reshape2)


# 2. Command Line Variables #
args = commandArgs(trailingOnly=TRUE)
vqsr.uniq.fn = args[1] # file name and path to annotations file for sites uniquely retained post VQSR
hardf.uniq.fn = args[2] # file name and path to annotations file for sites uniquely retained post hard filtering
overlap.fn = args[3] # file name and path to annotations file for sites that overlap between VQSR and hard filtering
out.fn.1 = args[4] # output file for DP
out.fn.2 = args[5] # output file for QD


# 3. Read in data #
vqsr.uniq <- read.table(vqsr.uniq.fn, sep = "\t", header = T)
hardf.uniq <- read.table(hardf.uniq.fn, sep = "\t", header = T)
overlap <- read.table(overlap.fn, sep = "\t", header = T)


# 4. Plot #
#---------------------
# DP and Normalized DP 
#---------------------
png(out.fn.1, width=13, height=5, units = "in", res = 300)

# Add DP for all 3 in one data frame
dp.dat.merged <- cbind(vqsr.uniq$DP, hardf.uniq$DP, overlap$DP)
colnames(dp.dat.merged) <- c("VQSR", "Hard Filter", "Overlap")
dat.for.plot.DP <- melt(dp.dat.merged)
p1 <- ggplot(dat.for.plot.DP,aes(x=value, fill=Var2)) + geom_density(alpha=0.25) +
  theme_bw() +
  theme(axis.text.x = element_text(size=16, colour="black"), 
        axis.text.y = element_text(size=16, colour = "black"), 
        axis.title= element_text(size=16), 
        axis.line=element_line(colour = "black"), 
        axis.ticks = element_line(colour = "black", size=1.5), 
        axis.ticks.length = unit(0.3, 'cm')) + 
  guides(fill=guide_legend(title="Key"))

# Calculate DP normalized by AN
vqsr.uniq$DP_normed <- vqsr.uniq$DP/(vqsr.uniq$AN/2)
hardf.uniq$DP_normed <- hardf.uniq$DP/(hardf.uniq$AN/2)
overlap$DP_normed <- overlap$DP/(overlap$AN/2)
# Combine data from all 3
dp.norm.dat.merged <- cbind(vqsr.uniq$DP_normed, hardf.uniq$DP_normed, overlap$DP_normed)
colnames(dp.norm.dat.merged) <- c("VQSR", "Hard Filter", "Overlap")
dat.for.plot.DP.norm <- melt(dp.norm.dat.merged)
p2 <- ggplot(dat.for.plot.DP.norm,aes(x=value, fill=Var2)) + geom_density(alpha=0.25) +
  theme_bw() +
  theme(axis.text.x = element_text(size=16, colour="black"), 
        axis.text.y = element_text(size=16, colour = "black"), 
        axis.title= element_text(size=16), 
        axis.line=element_line(colour = "black"), 
        axis.ticks = element_line(colour = "black", size=1.5), 
        axis.ticks.length = unit(0.3, 'cm')) + 
  guides(fill=guide_legend(title="Key"))

ggarrange(p1, p2, ncol=2, common.legend = T)
dev.off()

#----
# QD
#----
png(out.fn.2, width=8, height=5, units = "in", res = 300)
qd.dat.merged <- cbind(vqsr.uniq$QD, hardf.uniq$QD, overlap$QD)
colnames(qd.dat.merged) <- c("VQSR", "Hard Filter", "Overlap")
dat.for.plot.QD<- melt(qd.dat.merged)
ggplot(dat.for.plot.QD,aes(x=value, fill=Var2)) + geom_density(alpha=0.25) +
  theme_bw() +
  theme(axis.text.x = element_text(size=16, colour="black"), 
        axis.text.y = element_text(size=16, colour = "black"), 
        axis.title= element_text(size=16), 
        axis.line=element_line(colour = "black"), 
        axis.ticks = element_line(colour = "black", size=1.5), 
        axis.ticks.length = unit(0.3, 'cm')) + 
  guides(fill=guide_legend(title="Key")) 
dev.off()