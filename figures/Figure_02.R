# Angela Oill
# 2022-03-11
# Plotting variants in XTR with different alignment approaches

##################
# Load Libraries #
##################
library(ggplot2)
library(ggpubr)


################
# Read in data #
################
scc_XTR <- read.table("males_scc_XTR_variant_counts.txt", header = T)
sim_XTR <- read.table("males_simulated_XTR_variant_counts.txt", header = T)
colnames(sim_XTR) <- c("SampleID","Num_sim_variants","sim_Category")

X_Y_XTR_dat <- cbind(scc_XTR, sim_XTR$Num_sim_variants, sim_XTR$sim_Category)
colnames(X_Y_XTR_dat) <- c("SampleID", "Num_variants", "Category", "Num_sim_variants","sim_Category")
X_Y_XTR_dat$called_to_sim <- X_Y_XTR_dat$Num_variants/X_Y_XTR_dat$Num_sim_variants


# SCC + XTR masked data
YXTRmasked_XTR <- read.table("males_scc_YXTRmased_XTR_variant_counts.txt", header = T)
sim_XTR_XY_sum <- read.table("males_simulated_XTR_variant_counts_XYsummed.txt", header = T)
colnames(sim_XTR_XY_sum) <- c("SampleID","Num_sim_variants","sim_Category")

YXTRmasked_XTR_dat <- cbind(YXTRmasked_XTR, sim_XTR_XY_sum$Num_sim_variants, sim_XTR_XY_sum$sim_Category)
colnames(YXTRmasked_XTR_dat) <- c("SampleID", "Num_variants", "Category", "Num_sim_variants","sim_Category")
YXTRmasked_XTR_dat$called_to_sim <- YXTRmasked_XTR_dat$Num_variants/YXTRmasked_XTR_dat$Num_sim_variants

all_dat_merged <- rbind(X_Y_XTR_dat,
                        YXTRmasked_XTR_dat)

##########
# Plot 1 #
##########
# This is just number of called variants in XTR per chromosome and alignment strategy
ggplot(all_dat_merged, aes(x=Category, y=Num_variants)) + 
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  theme_classic() +
  ylab("# Called Variants in XTR") +
  scale_x_discrete(labels=c("SCC_chrX_XTR" = "X chromosome\n(SCC)", 
                            "SCC_chrY_XTR" = "Y chromosome\n(SCC)", 
                            "SCC_YXTRmasked_chrX_XTR" = "X chromosome\n(SCC + Y-XTR masked)")
  ) +
  ggtitle("Simulated Males (XY)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1)) 

##########
# Plot 2 #
##########
# This is the proportion of called variants over the number of simulated variants
# in XTR per chromosome and alignment strategy
pdf("prop_called_to_sim_XTR_SCC_vs_SCCwithYXTRmasking.pdf",
    width = 5,
    height = 5)
ggplot(all_dat_merged, aes(x=Category, y=called_to_sim)) + 
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  theme_classic() +
  ylab("# Called Variants in XTR/\n# Simulated Variants in XTR") +
  scale_x_discrete(labels=c("SCC_chrX_XTR" = "X chromosome\n(SCC)", 
                            "SCC_chrY_XTR" = "Y chromosome\n(SCC)", 
                            "SCC_YXTRmasked_chrX_XTR" = "X chromosome\n(SCC + Y-XTR masked)")
  ) +
  ggtitle("Simulated Males (XY)") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold"),
    axis.text.x = element_text(angle = 45, hjust=1)) +
  geom_hline(yintercept=1, 
             linetype="dashed", color = "red") 
dev.off()


