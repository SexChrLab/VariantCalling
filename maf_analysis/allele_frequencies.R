#=============================================================================#
# Angela Oill
# 2023-01-29
# M/F allele frequency plots
#=============================================================================#

library(ggplot2)
library(cowplot)
library(ggpubr)
library(dplyr)
library(stringr)
library(scales)

###############################################################################
###############################################################################
# SCC ----
###############################################################################
###############################################################################

# Read in data
# remove header in a text editor
XnonPARs_M <- read.table("males_SCC_XnonPARs_frq_fix.frq",
                         header = F, sep = "\t")
PARs_M <- read.table("males_SCC_PARs_frq_fix.frq",
                     header = F, sep = "\t")
X_M <- rbind(PARs_M, XnonPARs_M)
colnames(X_M) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "FREQ_A1", "FREQ_A2")
# Get minor allele
X_M$FREQ_A1 <- gsub("^.*:", "", X_M$FREQ_A1)
X_M$FREQ_A2 <- gsub("^.*:", "", X_M$FREQ_A2)

#X_M <- X_M %>%
#  filter(FREQ_A1 != "-nan") %>%
#  filter(FREQ_A2 != "-nan") 

X_M <- X_M %>%
  mutate(
    FREQ_MINOR = case_when(
      X_M$FREQ_A1 > X_M$FREQ_A2 ~ X_M$FREQ_A2,
      X_M$FREQ_A1 < X_M$FREQ_A2 ~ X_M$FREQ_A1,
      TRUE                      ~ X_M$FREQ_A1
    )
  )
X_M$SEX <- "Males"


XnonPARs_F <- read.table("females_SCC_XnonPARs_frq_fix.frq",
                         header = F, sep = "\t")
PARs_F <- read.table("females_SCC_PARs_frq_fix.frq",
                         header = F, sep = "\t")
X_F <- rbind(PARs_F, XnonPARs_F)
colnames(X_F) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "FREQ_A1", "FREQ_A2")
# Get minor allele
X_F$FREQ_A1 <- gsub("^.*:", "", X_F$FREQ_A1)
X_F$FREQ_A2 <- gsub("^.*:", "", X_F$FREQ_A2)

X_F <- X_F %>%
  mutate(
    FREQ_FINOR = case_when(
      X_F$FREQ_A1 > X_F$FREQ_A2 ~ X_F$FREQ_A2,
      X_F$FREQ_A1 < X_F$FREQ_A2 ~ X_F$FREQ_A1,
      TRUE                      ~ X_F$FREQ_A1
    )
  )

X_F$SEX <- "Females"

colnames(X_M) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "FREQ_A1", "FREQ_A2", "FREQ_MINOR", "SEX")
colnames(X_F) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "FREQ_A1", "FREQ_A2", "FREQ_MINOR", "SEX")

X_dat <- rbind(X_F, X_M)

# Scatterplot
ggplot(X_dat, aes(x=POS, y=as.numeric(FREQ_MINOR))) + 
  geom_point() +
  facet_wrap(~SEX,  ncol=1) +
  theme_classic()

X_dat$FREQ_MINOR <- as.numeric(X_dat$FREQ_MINOR)

X_dat <- X_dat %>%
  mutate(
    THRESHOLD = case_when(
      as.numeric(X_dat$POS) <= 2781479 ~ "PAR1",
      as.numeric(X_dat$POS) >= 155701383 ~ "PAR2",
      as.numeric(X_dat$POS) >= 89140830 & as.numeric(X_dat$POS) <= 93428068 ~ "XTR",
      TRUE                      ~ "Non-PAR"
    )
  )


pdf("SCC_F_M_MAF.pdf",
    width = 12,
    height = 10)
ggplot(X_dat, aes(x=POS, y=as.numeric(FREQ_MINOR), color=THRESHOLD)) + 
  geom_point(size = 0.2) +
  facet_wrap(~SEX,  ncol=1) +
  scale_x_continuous(labels = function(x) format(x, scientific = FALSE)) +
  scale_color_manual(values = c("black", "grey", "grey", "grey")) +
  xlab("Chromosome X") + ylab("MAF") +
  theme_classic() + 
  theme(legend.position = "none")
dev.off()
#X_dat_filter <- X_dat %>%
#  filter(FREQ_MINOR >= 0.05)

#ggplot(X_dat_filter, aes(x=POS, y=as.numeric(FREQ_MINOR))) + 
#  geom_point() +
#  facet_wrap(~SEX,  ncol=1) +
#  theme_classic()

X_dat_diff <- as.data.frame(as.numeric(X_F$FREQ_MINOR) - as.numeric(X_M$FREQ_MINOR))
X_dat_diff$POS <- X_F$POS
colnames(X_dat_diff) <- c("FREQ_MINOR_DIFF", "POS")

ggplot(X_dat_diff, aes(x=POS, y=as.numeric(FREQ_MINOR_DIFF))) + 
  geom_point() +
  scale_x_continuous(labels = function(x) format(x, scientific = FALSE)) +
  theme_classic() 

X_dat_diff <- X_dat_diff %>%
  mutate(
    THRESHOLD = case_when(
      as.numeric(X_dat_diff$POS) <= 2781479 ~ "PAR1",
      as.numeric(X_dat_diff$POS) >= 155701383 ~ "PAR2",
      as.numeric(X_dat_diff$POS) >= 89140830 & as.numeric(X_dat_diff$POS) <= 93428068 ~ "XTR",
      TRUE                      ~ "Non-PAR"
  )
  )

pdf("SCC_F_M_diff.pdf",
    width = 12,
    height = 5)
ggplot(X_dat_diff, aes(x=POS, y=as.numeric(FREQ_MINOR_DIFF), color=THRESHOLD)) + 
  geom_point(size = 0.2) +
  scale_x_continuous(labels = function(x) format(x, scientific = FALSE)) +
  scale_color_manual(values = c("black", "grey", "grey", "grey")) +
  xlab("Chromosome X") + ylab("Female - Male MAF") +
  geom_hline(yintercept=0, color = "red") +
  theme_classic() + 
  theme(legend.position = "none")
dev.off()


###############################################################################
###############################################################################
# DEFAULT ----
###############################################################################
###############################################################################

# Read in data
# remove header in a text editor
XnonPARs_M_default <- read.table("males_default_XnonPARs_frq_fix.frq",
                                 header = F, sep = "\t")
PARs_M_default <- read.table("males_default_PARs_frq_fix.frq",
                             header = F, sep = "\t")
X_M_default <- rbind(PARs_M_default, XnonPARs_M_default)
colnames(X_M_default) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "FREQ_A1", "FREQ_A2")
# Get minor allele
X_M_default$FREQ_A1 <- gsub("^.*:", "", X_M_default$FREQ_A1)
X_M_default$FREQ_A2 <- gsub("^.*:", "", X_M_default$FREQ_A2)

#X_M_default <- X_M_default %>%
#  filter(FREQ_A1 != "-nan") %>%
#  filter(FREQ_A2 != "-nan") 

X_M_default <- X_M_default %>%
  mutate(
    FREQ_MINOR = case_when(
      X_M_default$FREQ_A1 > X_M_default$FREQ_A2 ~ X_M_default$FREQ_A2,
      X_M_default$FREQ_A1 < X_M_default$FREQ_A2 ~ X_M_default$FREQ_A1,
      TRUE                      ~ X_M_default$FREQ_A1
    )
  )
X_M_default$SEX <- "Males"


XnonPARs_F_default <- read.table("females_default_XnonPARs_frq_fix.frq",
                                 header = F, sep = "\t")
PARs_F_default <- read.table("females_default_PARs_frq_fix.frq",
                             header = F, sep = "\t")
X_F_default <- rbind(PARs_F_default, XnonPARs_F_default)
colnames(X_F_default) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "FREQ_A1", "FREQ_A2")
# Get minor allele
X_F_default$FREQ_A1 <- gsub("^.*:", "", X_F_default$FREQ_A1)
X_F_default$FREQ_A2 <- gsub("^.*:", "", X_F_default$FREQ_A2)

X_F_default <- X_F_default %>%
  mutate(
    FREQ_FINOR = case_when(
      X_F_default$FREQ_A1 > X_F_default$FREQ_A2 ~ X_F_default$FREQ_A2,
      X_F_default$FREQ_A1 < X_F_default$FREQ_A2 ~ X_F_default$FREQ_A1,
      TRUE                      ~ X_F_default$FREQ_A1
    )
  )

X_F_default$SEX <- "Females"

colnames(X_M_default) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "FREQ_A1", "FREQ_A2", "FREQ_MINOR", "SEX")
colnames(X_F_default) <- c("CHROM", "POS", "N_ALLELES", "N_CHR", "FREQ_A1", "FREQ_A2", "FREQ_MINOR", "SEX")

X_dat_default <- rbind(X_F_default, X_M_default)

# Scatterplot
ggplot(X_dat_default, aes(x=POS, y=as.numeric(FREQ_MINOR))) + 
  geom_point() +
  facet_wrap(~SEX,  ncol=1) +
  theme_classic()

X_dat_default$FREQ_MINOR <- as.numeric(X_dat_default$FREQ_MINOR)

X_dat_default <- X_dat_default %>%
  mutate(
    THRESHOLD = case_when(
      as.numeric(X_dat_default$POS) <= 2781479 ~ "PAR1",
      as.numeric(X_dat_default$POS) >= 155701383 ~ "PAR2",
      as.numeric(X_dat_default$POS) >= 89140830 & as.numeric(X_dat_default$POS) <= 93428068 ~ "XTR",
      TRUE                      ~ "Non-PAR"
    )
  )


pdf("default_F_M_MAF.pdf",
    width = 12,
    height = 10)
ggplot(X_dat_default, aes(x=POS, y=as.numeric(FREQ_MINOR), color=THRESHOLD)) + 
  geom_point(size = 0.2) +
  facet_wrap(~SEX,  ncol=1) +
  scale_x_continuous(labels = function(x) format(x, scientific = FALSE)) +
  scale_color_manual(values = c("black", "grey", "grey", "grey")) +
  xlab("Chromosome X") + ylab("MAF") +
  theme_classic() + 
  theme(legend.position = "none")
dev.off()
#X_dat_default_filter <- X_dat_default %>%
#  filter(FREQ_MINOR >= 0.05)

#ggplot(X_dat_default_filter, aes(x=POS, y=as.numeric(FREQ_MINOR))) + 
#  geom_point() +
#  facet_wrap(~SEX,  ncol=1) +
#  theme_classic()

X_dat_default_diff <- as.data.frame(as.numeric(X_F_default$FREQ_MINOR) - as.numeric(X_M_default$FREQ_MINOR))
X_dat_default_diff$POS <- X_F_default$POS
colnames(X_dat_default_diff) <- c("FREQ_MINOR_DIFF", "POS")

ggplot(X_dat_default_diff, aes(x=POS, y=as.numeric(FREQ_MINOR_DIFF))) + 
  geom_point() +
  scale_x_continuous(labels = function(x) format(x, scientific = FALSE)) +
  theme_classic() 

X_dat_default_diff <- X_dat_default_diff %>%
  mutate(
    THRESHOLD = case_when(
      as.numeric(X_dat_default_diff$POS) <= 2781479 ~ "PAR1",
      as.numeric(X_dat_default_diff$POS) >= 155701383 ~ "PAR2",
      as.numeric(X_dat_default_diff$POS) >= 89140830 & as.numeric(X_dat_default_diff$POS) <= 93428068 ~ "XTR",
      TRUE                      ~ "Non-PAR"
    )
  )

pdf("default_F_M_diff.pdf",
    width = 12,
    height = 5)
ggplot(X_dat_default_diff, aes(x=POS, y=as.numeric(FREQ_MINOR_DIFF), color=THRESHOLD)) + 
  geom_point(size = 0.2) +
  scale_x_continuous(labels = function(x) format(x, scientific = FALSE)) +
  scale_color_manual(values = c("black", "grey", "grey", "grey")) +
  xlab("Chromosome X") + ylab("Female - Male MAF") +
  geom_hline(yintercept=0, color = "red") +
  theme_classic() + 
  theme(legend.position = "none")
dev.off()


###############################################################################
###############################################################################
# SCC/DEFAULT DIFFERENCE ----
###############################################################################
###############################################################################
#---------#
# FEMALES #
#---------#
# Difference in allele frequencies in females (default vs SCC)
X_F_default$ALIGNMENT <- "Default"
X_F$ALIGNMENT <- "SCC"

X_F_alignments_merged <- merge(X_F, X_F_default, by = "POS", all = TRUE)
#X_F_alignments_merged[is.na(X_F_alignments_merged)] <- 0
# Might want to keep NAs?
X_F_alignments_merged$DIFF <- as.numeric(X_F_alignments_merged$FREQ_MINOR.x) - 
  as.numeric(X_F_alignments_merged$FREQ_MINOR.y)

# Or make them 0
X_F_alignments_merged_NAs0 <- X_F_alignments_merged
X_F_alignments_merged_NAs0[is.na(X_F_alignments_merged_NAs0)] <- 0
X_F_alignments_merged_NAs0$DIFF <- as.numeric(X_F_alignments_merged_NAs0$FREQ_MINOR.x) - 
  as.numeric(X_F_alignments_merged_NAs0$FREQ_MINOR.y)



X_F_alignments_merged <- X_F_alignments_merged %>%
  mutate(
    THRESHOLD = case_when(
      as.numeric(X_F_alignments_merged$POS) <= 2781479 ~ "PAR1",
      as.numeric(X_F_alignments_merged$POS) >= 155701383 ~ "PAR2",
      as.numeric(X_F_alignments_merged$POS) >= 89140830 & as.numeric(X_F_alignments_merged$POS) <= 93428068 ~ "XTR",
      TRUE                      ~ "Non-PAR"
    )
  )


pdf("scc_default_F_diff_noNAs.pdf",
    width = 12,
    height = 5)
ggplot(X_F_alignments_merged, aes(x=POS, y=as.numeric(DIFF), color=THRESHOLD)) + 
  geom_point(size = 0.2) +
  scale_x_continuous(labels = function(x) format(x, scientific = FALSE)) +
  scale_color_manual(values = c("black", "grey", "grey", "grey")) +
  xlab("Chromosome X") + ylab("SCC - Default") +
  geom_hline(yintercept=0, color = "red") +
  theme_classic() + 
  ggtitle("Female (XX) individuals") +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()



X_F_alignments_merged_NAs0 <- X_F_alignments_merged_NAs0 %>%
  mutate(
    THRESHOLD = case_when(
      as.numeric(X_F_alignments_merged_NAs0$POS) <= 2781479 ~ "PAR1",
      as.numeric(X_F_alignments_merged_NAs0$POS) >= 155701383 ~ "PAR2",
      as.numeric(X_F_alignments_merged_NAs0$POS) >= 89140830 & as.numeric(X_F_alignments_merged_NAs0$POS) <= 93428068 ~ "XTR",
      TRUE                      ~ "Non-PAR"
    )
  )

pdf("scc_default_F_diff_NAs0.pdf",
    width = 12,
    height = 5)
ggplot(X_F_alignments_merged_NAs0, aes(x=POS, y=as.numeric(DIFF), color=THRESHOLD)) + 
  geom_point(size = 0.2) +
  scale_x_continuous(labels = function(x) format(x, scientific = FALSE)) +
  scale_color_manual(values = c("black", "grey", "grey", "grey")) +
  xlab("Chromosome X") + ylab("SCC - Default") +
  geom_hline(yintercept=0, color = "red") +
  theme_classic() + 
  ggtitle("Female (XX) individuals") +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()



#-------#
# MALES #
#-------#
# Difference in allele frequencies in females (default vs SCC)
X_M_default$ALIGNMENT <- "Default"
X_M$ALIGNMENT <- "SCC"

X_M_alignments_merged <- merge(X_M, X_M_default, by = "POS", all = TRUE)
#X_M_alignments_merged[is.na(X_M_alignments_merged)] <- 0
# Might want to keep NAs?
X_M_alignments_merged$DIFF <- as.numeric(X_M_alignments_merged$FREQ_MINOR.x) - 
  as.numeric(X_M_alignments_merged$FREQ_MINOR.y)

# Or make them 0
X_M_alignments_merged_NAs0 <- X_M_alignments_merged
X_M_alignments_merged_NAs0[is.na(X_M_alignments_merged_NAs0)] <- 0
X_M_alignments_merged_NAs0$DIFF <- as.numeric(X_M_alignments_merged_NAs0$FREQ_MINOR.x) - 
  as.numeric(X_M_alignments_merged_NAs0$FREQ_MINOR.y)



X_M_alignments_merged <- X_M_alignments_merged %>%
  mutate(
    THRESHOLD = case_when(
      as.numeric(X_M_alignments_merged$POS) <= 2781479 ~ "PAR1",
      as.numeric(X_M_alignments_merged$POS) >= 155701383 ~ "PAR2",
      as.numeric(X_M_alignments_merged$POS) >= 89140830 & as.numeric(X_M_alignments_merged$POS) <= 93428068 ~ "XTR",
      TRUE                      ~ "Non-PAR"
    )
  )


pdf("scc_default_M_diff_noNAs.pdf",
    width = 12,
    height = 5)
ggplot(X_M_alignments_merged, aes(x=POS, y=as.numeric(DIFF), color=THRESHOLD)) + 
  geom_point(size = 0.2) +
  scale_x_continuous(labels = function(x) format(x, scientific = FALSE)) +
  scale_color_manual(values = c("black", "grey", "grey", "grey")) +
  xlab("Chromosome X") + ylab("SCC - Default") +
  geom_hline(yintercept=0, color = "red") +
  theme_classic() + 
  ggtitle("Male (XX) individuals") +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()



X_M_alignments_merged_NAs0 <- X_M_alignments_merged_NAs0 %>%
  mutate(
    THRESHOLD = case_when(
      as.numeric(X_M_alignments_merged_NAs0$POS) <= 2781479 ~ "PAR1",
      as.numeric(X_M_alignments_merged_NAs0$POS) >= 155701383 ~ "PAR2",
      as.numeric(X_M_alignments_merged_NAs0$POS) >= 89140830 & as.numeric(X_M_alignments_merged_NAs0$POS) <= 93428068 ~ "XTR",
      TRUE                      ~ "Non-PAR"
    )
  )

pdf("scc_default_M_diff_NAs0.pdf",
    width = 12,
    height = 5)
ggplot(X_M_alignments_merged_NAs0, aes(x=POS, y=as.numeric(DIFF), color=THRESHOLD)) + 
  geom_point(size = 0.2) +
  scale_x_continuous(labels = function(x) format(x, scientific = FALSE)) +
  scale_color_manual(values = c("black", "grey", "grey", "grey")) +
  xlab("Chromosome X") + ylab("SCC - Default") +
  geom_hline(yintercept=0, color = "red") +
  theme_classic() + 
  ggtitle("Male (XX) individuals") +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()


#-----------------------------------------------------------------#
# Difference in the sex difference in allele freq (default vs SCC)
#-----------------------------------------------------------------#
# X_dat_default_diff
# X_dat_diff
X_dat_default_diff$ALIGNMENT <- "Default"
X_dat_diff$ALIGNMENT <- "SCC"

X_diff_alignments_merged <- merge(X_dat_diff, X_dat_default_diff, by = "POS", all = TRUE)
X_diff_alignments_merged$DIFF <- as.numeric(X_diff_alignments_merged$FREQ_MINOR_DIFF.x) - 
  as.numeric(X_diff_alignments_merged$FREQ_MINOR_DIFF.y)

# Or make them 0
X_diff_alignments_merged_NAs0 <- X_diff_alignments_merged
X_diff_alignments_merged_NAs0[is.na(X_diff_alignments_merged_NAs0)] <- 0
X_diff_alignments_merged_NAs0$DIFF <- as.numeric(X_diff_alignments_merged_NAs0$FREQ_MINOR_DIFF.x) - 
  as.numeric(X_diff_alignments_merged_NAs0$FREQ_MINOR_DIFF.y)



X_diff_alignments_merged <- X_diff_alignments_merged %>%
  mutate(
    THRESHOLD = case_when(
      as.numeric(X_diff_alignments_merged$POS) <= 2781479 ~ "PAR1",
      as.numeric(X_diff_alignments_merged$POS) >= 155701383 ~ "PAR2",
      as.numeric(X_diff_alignments_merged$POS) >= 89140830 & as.numeric(X_diff_alignments_merged$POS) <= 93428068 ~ "XTR",
      TRUE                      ~ "Non-PAR"
    )
  )

pdf("scc_default_maf_diff_diff_noNAs.pdf",
    width = 12,
    height = 5)
ggplot(X_diff_alignments_merged, aes(x=POS, y=as.numeric(DIFF), color=THRESHOLD)) + 
  geom_point(size = 0.2) +
  scale_x_continuous(labels = function(x) format(x, scientific = FALSE)) +
  scale_color_manual(values = c("black", "grey", "grey", "grey")) +
  xlab("Chromosome X") + ylab("Difference in Female - Male MAF between alignments") +
  geom_hline(yintercept=0, color = "red") +
  theme_classic() + 
  ggtitle("SCC - Default") +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()


X_diff_alignments_merged_NAs0 <- X_diff_alignments_merged_NAs0 %>%
  mutate(
    THRESHOLD = case_when(
      as.numeric(X_diff_alignments_merged_NAs0$POS) <= 2781479 ~ "PAR1",
      as.numeric(X_diff_alignments_merged_NAs0$POS) >= 155701383 ~ "PAR2",
      as.numeric(X_diff_alignments_merged_NAs0$POS) >= 89140830 & as.numeric(X_diff_alignments_merged_NAs0$POS) <= 93428068 ~ "XTR",
      TRUE                      ~ "Non-PAR"
    )
  )


pdf("scc_default_maf_diff_diff_NAs0.pdf",
    width = 12,
    height = 5)
ggplot(X_diff_alignments_merged_NAs0, aes(x=POS, y=as.numeric(DIFF), color=THRESHOLD)) + 
  geom_point(size = 0.2) +
  scale_x_continuous(labels = function(x) format(x, scientific = FALSE)) +
  scale_color_manual(values = c("black", "grey", "grey", "grey")) +
  xlab("Chromosome X") + ylab("Difference in Female - Male MAF between alignments") +
  geom_hline(yintercept=0, color = "red") +
  theme_classic() + 
  ggtitle("SCC - Default") +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()
