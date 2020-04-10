# Angela Oill 
# 2020-04-10
# multiple SFSs on 1 plot to compare no filtering, vqsr and hard filtering 
# USAGE: Rscript sfs_plotting.R sfs_raw_file.txt sfs_vqsr_file.txt sfs_hard_filter.txt output.SFSs.png

# Load libraries #
library(ggplot2)
library(viridis)


# Command Line Variables #
args = commandArgs(trailingOnly=TRUE)
raw.dat.fn = args[1]
vqsr.dat.fn = args[2]
hard.filter.dat.fn = args[3]
out.fn = args[4]


# Load data #
raw.dat <- read.table(raw.dat.fn, header = T)
vqsr.dat <- read.table(vqsr.dat.fn, header = T)
hard.filter.dat <- read.table(hard.filter.dat.fn, header = T)


# Set up input data frame for plotting #
# Add column to each dat file with the data type
raw.dat$type <- 'Raw'
vqsr.dat$type <- 'VQSR'
hard.filter.dat$type <- 'Hard Filter'
# Merge all data together
full.dat <- rbind(raw.dat, vqsr.dat, hard.filter.dat)


# Plot #
# Order data by type - I want it in order from left to right raw, vqsr, hard filter
full.dat$type <- factor(full.dat$type, levels = c("Raw", "VQSR", "Hard Filter"))
# Save as png
png(out.fn, width = 11, height = 8, res = 300, units = "in")
# Plot
ggplot(full.dat, aes(x = af_bin, y = count, fill = type)) +
  theme(legend.title=element_blank(), text = element_text(size=14), plot.title = element_text(hjust = 0.5, face="bold")) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  labs(title="Site frequency spectrum", x="Count", y="Number of variants") +
  scale_x_continuous(breaks = full.dat$af_bin) + scale_fill_viridis(discrete=TRUE)
dev.off()
