# Angela Oill
# 2020-04-16
# Plot pi results for variant calling project
# USAGE:


# Load libraries #
library(ggplot2)
library(viridis)
require(scales)


# Read in all data
results.main.dir <- "Dropbox (ASU)/00_Lab_Research_Projects/Kenya/00_Main_Projects/variant_calling/output/pop_gen_results/"

# All called variants #
# Raw files and add columns for chr and filtering strategy
chr.8.f.raw <- read.table(paste(results.main.dir, 
                                "pre_filter/pi/autosomes/chr8/chr8.gatk.called.raw.females.pi.per.site.txt", 
                                sep = ""),
                          header = T)
chr.8.f.raw$chromosome <- 'Chr 8 (Females)'
chr.8.f.raw$type <- 'Raw'

chr.8.m.raw <- read.table(paste(results.main.dir, 
                                "pre_filter/pi/autosomes/chr8/chr8.gatk.called.raw.males.pi.per.site.txt", 
                                sep = ""),
                          header = T)
chr.8.m.raw$chromosome <- 'Chr 8 (Males)'
chr.8.m.raw$type <- 'Raw'

chr.X.f.raw <- read.table(paste(results.main.dir, 
                                "pre_filter/pi/sex_chromosomes/chrX/chrX.gatk.called.raw.females.noPARs.noXTR.pi.per.site.txt", 
                                sep = ""),
                          header = T)
chr.X.f.raw$chromosome <- 'X chr (Females)'
chr.X.f.raw$type <- 'Raw'


chr.X.m.hap.raw <- read.table(paste(results.main.dir, 
                                "pre_filter/pi/sex_chromosomes/chrX/chrX.gatk.called.raw.haploid.noPARs.noXTR.pi.per.site.txt", 
                                sep = ""),
                          header = T)
chr.X.m.hap.raw$chromosome <- 'X chr (Males haploid)'
chr.X.m.hap.raw$type <- 'Raw'

chr.X.m.dip.raw <- read.table(paste(results.main.dir, 
                                    "pre_filter/pi/sex_chromosomes/chrX/chrX.gatk.called.raw.diploid.noPARs.noXTR.pi.per.site.txt", 
                                    sep = ""),
                              header = T)
chr.X.m.dip.raw$chromosome <- 'X chr (Males diploid)'
chr.X.m.dip.raw$type <- 'Raw'


# VQSR filtered files 
chr.8.f.vqsr <- read.table(paste(results.main.dir, 
                                "vqsr/pi/autosomes/chr8/chr8.gatk.called.females.vqsr.sv.pi.per.site.txt", 
                                sep = ""),
                          header = T)
chr.8.f.vqsr$chromosome <- 'Chr 8 (Females)'
chr.8.f.vqsr$type <- 'VQSR'

chr.8.m.vqsr <- read.table(paste(results.main.dir, 
                                "vqsr/pi/autosomes/chr8/chr8.gatk.called.males.vqsr.sv.pi.per.site.txt", 
                                sep = ""),
                          header = T)
chr.8.m.vqsr$chromosome <- 'Chr 8 (Males)'
chr.8.m.vqsr$type <- 'VQSR'

chr.X.f.vqsr <- read.table(paste(results.main.dir, 
                                "vqsr/pi/sex_chromosomes/chrX/chrX.gatk.called.females.vqsr.sv.noPARs.noXTR.pi.per.site.txt", 
                                sep = ""),
                          header = T)
chr.X.f.vqsr$chromosome <- 'X chr (Females)'
chr.X.f.vqsr$type <- 'VQSR'

chr.X.m.hap.vqsr <- read.table(paste(results.main.dir, 
                                    "vqsr/pi/sex_chromosomes/chrX/chrX.gatk.called.haploid.vqsr.sv.noPARs.noXTR.pi.per.site.txt", 
                                    sep = ""),
                              header = T)
chr.X.m.hap.vqsr$chromosome <- 'X chr (Males haploid)'
chr.X.m.hap.vqsr$type <- 'VQSR'

chr.X.m.dip.vqsr <- read.table(paste(results.main.dir, 
                                    "vqsr/pi/sex_chromosomes/chrX/chrX.gatk.called.diploid.vqsr.sv.noPARs.noXTR.pi.per.site.txt", 
                                    sep = ""),
                              header = T)
chr.X.m.dip.vqsr$chromosome <- 'X chr (Males diploid)'
chr.X.m.dip.vqsr$type <- 'VQSR'


# Hard filtered files 
chr.8.f.hardf <- read.table(paste(results.main.dir, 
                                 "hard_filter/pi/autosomes/chr8/chr8.gatk.called.females.snps.gatkHardFilter.sv.pi.per.site.txt", 
                                 sep = ""),
                           header = T)
chr.8.f.hardf$chromosome <- 'Chr 8 (Females)'
chr.8.f.hardf$type <- 'Hard Filter'

chr.8.m.hardf <- read.table(paste(results.main.dir, 
                                 "hard_filter/pi/autosomes/chr8/chr8.gatk.called.males.snps.gatkHardFilter.sv.pi.per.site.txt", 
                                 sep = ""),
                           header = T)
chr.8.m.hardf$chromosome <- 'Chr 8 (Males)'
chr.8.m.hardf$type <- 'Hard Filter'

chr.X.f.hardf <- read.table(paste(results.main.dir, 
                                 "hard_filter/pi/sex_chromosomes/chrX/chrX.gatk.called.females.snps.gatkHardFilter.sv.noPARs.noXTR.pi.per.site.txt", 
                                 sep = ""),
                           header = T)
chr.X.f.hardf$chromosome <- 'X chr (Females)'
chr.X.f.hardf$type <- 'Hard Filter'

chr.X.m.hap.hardf <- read.table(paste(results.main.dir, 
                                     "hard_filter/pi/sex_chromosomes/chrX/chrX.gatk.called.haploid.snps.gatkHardFilter.sv.noPARs.noXTR.pi.per.site.txt", 
                                     sep = ""),
                               header = T)
chr.X.m.hap.hardf$chromosome <- 'X chr (Males haploid)'
chr.X.m.hap.hardf$type <- 'Hard Filter'

chr.X.m.dip.hardf <- read.table(paste(results.main.dir, 
                                     "hard_filter/pi/sex_chromosomes/chrX/chrX.gatk.called.diploid.snps.gatkHardFilter.sv.noPARs.noXTR.pi.per.site.txt", 
                                     sep = ""),
                               header = T)
chr.X.m.dip.hardf$chromosome <- 'X chr (Males diploid)'
chr.X.m.dip.hardf$type <- 'Hard Filter'

# Merge all results into 1 data frame
full.dat <- rbind(chr.8.f.raw, chr.8.m.raw, chr.X.f.raw, chr.X.m.dip.raw, chr.X.m.hap.raw,
                  chr.8.f.vqsr, chr.8.m.vqsr, chr.X.f.vqsr, chr.X.m.dip.vqsr, chr.X.m.hap.vqsr,
                  chr.8.f.hardf, chr.8.m.hardf, chr.X.f.hardf, chr.X.m.dip.hardf, chr.X.m.hap.hardf)

# Plot # 
# Order data by type - I want it in order from left to right raw, vqsr, hard filter
full.dat$type <- factor(full.dat$type, levels = c("Raw", "VQSR", "Hard Filter"))

# Plot
png(paste(results.main.dir, "merged/pi/pi.all.called.variants.noPARs.noXTR.results.png", sep = ""), width = 11, height = 8, res = 300, units = "in")

# To make y labels nice looking (not scientific notation)
point <- format_format(big.mark = " ", decimal.mark = ".", scientific = FALSE)

ggplot(full.dat, aes(x = chromosome, y = pi_per_site, fill = type)) +
  theme(legend.title=element_blank(), text = element_text(size=14), plot.title = element_text(hjust = 0.5, face="bold")) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  labs(title="Pi comparisons - all called variants", x="Chromosome", y="Pi") +
  scale_fill_viridis(discrete=TRUE) +
  scale_y_continuous(labels = point)
dev.off()

#------------------------------------------------------------------------------
# X/A - all called sites
# RAW #
# Females no filtering X/A - A is chr 8
x.a.f.raw <- as.data.frame(chr.X.f.raw$pi_per_site / chr.8.f.raw$pi_per_site)
x.a.f.raw.col.n <- colnames(x.a.f.raw)
names(x.a.f.raw)[names(x.a.f.raw) == x.a.f.raw.col.n] <- "X.A"
x.a.f.raw$chromosome <- 'Females'
x.a.f.raw$type <- 'Raw'

# Males diploid no filtering X/A - A is chr 8
x.a.m.dip.raw <- as.data.frame(chr.X.m.dip.raw$pi_per_site / chr.8.m.raw$pi_per_site)
x.a.m.dip.raw.col.n <- colnames(x.a.m.dip.raw)
names(x.a.m.dip.raw)[names(x.a.m.dip.raw) == x.a.m.dip.raw.col.n] <- "X.A"
x.a.m.dip.raw$chromosome <- 'Males (diploid called X)'
x.a.m.dip.raw$type <- 'Raw'

# Males haploid no filtering X/A - A is chr 8
x.a.m.hap.raw <- as.data.frame(chr.X.m.hap.raw$pi_per_site / chr.8.m.raw$pi_per_site)
x.a.m.hap.raw.col.n <- colnames(x.a.m.hap.raw)
names(x.a.m.hap.raw)[names(x.a.m.hap.raw) == x.a.m.hap.raw.col.n] <- "X.A"
x.a.m.hap.raw$chromosome <- 'Males (haploid called X)'
x.a.m.hap.raw$type <- 'Raw'

# VQSR #
# Females vqsr X/A - A is chr 8
x.a.f.vqsr <- as.data.frame(chr.X.f.vqsr$pi_per_site / chr.8.f.vqsr$pi_per_site)
x.a.f.vqsr.col.n <- colnames(x.a.f.vqsr)
names(x.a.f.vqsr)[names(x.a.f.vqsr) == x.a.f.vqsr.col.n] <- "X.A"
x.a.f.vqsr$chromosome <- 'Females'
x.a.f.vqsr$type <- 'VQSR'

# Males diploid vqsr X/A - A is chr 8
x.a.m.dip.vqsr <- as.data.frame(chr.X.m.dip.vqsr$pi_per_site / chr.8.m.vqsr$pi_per_site)
x.a.m.dip.vqsr.col.n <- colnames(x.a.m.dip.vqsr)
names(x.a.m.dip.vqsr)[names(x.a.m.dip.vqsr) == x.a.m.dip.vqsr.col.n] <- "X.A"
x.a.m.dip.vqsr$chromosome <- 'Males (diploid called X)'
x.a.m.dip.vqsr$type <- 'VQSR'

# Males haploid vqsr X/A - A is chr 8
x.a.m.hap.vqsr <- as.data.frame(chr.X.m.hap.vqsr$pi_per_site / chr.8.m.vqsr$pi_per_site)
x.a.m.hap.vqsr.col.n <- colnames(x.a.m.hap.vqsr)
names(x.a.m.hap.vqsr)[names(x.a.m.hap.vqsr) == x.a.m.hap.vqsr.col.n] <- "X.A"
x.a.m.hap.vqsr$chromosome <- 'Males (haploid called X)'
x.a.m.hap.vqsr$type <- 'VQSR'

# HARD FILTER #
# Females Hard Filter X/A - A is chr 8
x.a.f.hardf <- as.data.frame(chr.X.f.hardf$pi_per_site / chr.8.f.hardf$pi_per_site)
x.a.f.hardf.col.n <- colnames(x.a.f.hardf)
names(x.a.f.hardf)[names(x.a.f.hardf) == x.a.f.hardf.col.n] <- "X.A"
x.a.f.hardf$chromosome <- 'Females'
x.a.f.hardf$type <- 'Hard Filter'

# Males diploid Hard Filter X/A - A is chr 8
x.a.m.dip.hardf <- as.data.frame(chr.X.m.dip.hardf$pi_per_site / chr.8.m.hardf$pi_per_site)
x.a.m.dip.hardf.col.n <- colnames(x.a.m.dip.hardf)
names(x.a.m.dip.hardf)[names(x.a.m.dip.hardf) == x.a.m.dip.hardf.col.n] <- "X.A"
x.a.m.dip.hardf$chromosome <- 'Males (diploid called X)'
x.a.m.dip.hardf$type <- 'Hard Filter'

# Males haploid Hard Filter X/A - A is chr 8
x.a.m.hap.hardf <- as.data.frame(chr.X.m.hap.hardf$pi_per_site / chr.8.m.hardf$pi_per_site)
x.a.m.hap.hardf.col.n <- colnames(x.a.m.hap.hardf)
names(x.a.m.hap.hardf)[names(x.a.m.hap.hardf) == x.a.m.hap.hardf.col.n] <- "X.A"
x.a.m.hap.hardf$chromosome <- 'Males (haploid called X)'
x.a.m.hap.hardf$type <- 'Hard Filter'

x.a.full.dat <- rbind(x.a.f.raw, x.a.m.dip.raw, x.a.m.hap.raw,
                      x.a.f.vqsr, x.a.m.dip.vqsr, x.a.m.hap.vqsr,
                      x.a.f.hardf, x.a.m.dip.hardf, x.a.m.hap.hardf)
# Order data by type - I want it in order from left to right raw, vqsr, hard filter
x.a.full.dat$type <- factor(x.a.full.dat$type, levels = c("Raw", "VQSR", "Hard Filter"))

# To make y labels nice looking (not scientific notation)
point <- format_format(big.mark = " ", decimal.mark = ".", scientific = FALSE)

# Plot
png(paste(results.main.dir, "merged/pi/X.A.all.called.variants.noPARs.noXTR.results.png", sep = ""), width = 11, height = 8, res = 300, units = "in")

ggplot(x.a.full.dat, aes(x = chromosome, y = X.A, fill = type)) +
  theme(legend.title=element_blank(), text = element_text(size=14), plot.title = element_text(hjust = 0.5, face="bold")) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  labs(title="X/A - all called variants", x="", y="Pi") +
  scale_fill_viridis(discrete=TRUE) +
  scale_y_continuous(labels = point)

dev.off()

#------------------------------------------------------------------------------
# Variants restricted to the array #
# Raw files and add columns for chr and filtering strategy
chr.8.f.raw <- read.table(paste(results.main.dir, 
                                "pre_filter/pi/autosomes/chr8/chr8.gatk.called.raw.females.array.sites.pi.per.site.txt", 
                                sep = ""),
                          header = T)
chr.8.f.raw$chromosome <- 'Chr 8 (Females)'
chr.8.f.raw$type <- 'Raw'

chr.8.m.raw <- read.table(paste(results.main.dir, 
                                "pre_filter/pi/autosomes/chr8/chr8.gatk.called.raw.males.array.sites.pi.per.site.txt", 
                                sep = ""),
                          header = T)
chr.8.m.raw$chromosome <- 'Chr 8 (Males)'
chr.8.m.raw$type <- 'Raw'

chr.X.f.raw <- read.table(paste(results.main.dir, 
                                "pre_filter/pi/sex_chromosomes/chrX/chrX.gatk.called.raw.females.array.sites.noPARs.noXTR.pi.per.site.txt", 
                                sep = ""),
                          header = T)
chr.X.f.raw$chromosome <- 'X chr (Females)'
chr.X.f.raw$type <- 'Raw'


chr.X.m.hap.raw <- read.table(paste(results.main.dir, 
                                    "pre_filter/pi/sex_chromosomes/chrX/chrX.gatk.called.raw.haploid.array.sites.noPARs.noXTR.pi.per.site.txt", 
                                    sep = ""),
                              header = T)
chr.X.m.hap.raw$chromosome <- 'X chr (Males haploid)'
chr.X.m.hap.raw$type <- 'Raw'

chr.X.m.dip.raw <- read.table(paste(results.main.dir, 
                                    "pre_filter/pi/sex_chromosomes/chrX/chrX.gatk.called.raw.diploid.array.sites.noPARs.noXTR.pi.per.site.txt", 
                                    sep = ""),
                              header = T)
chr.X.m.dip.raw$chromosome <- 'X chr (Males diploid)'
chr.X.m.dip.raw$type <- 'Raw'


# VQSR filtered files 
chr.8.f.vqsr <- read.table(paste(results.main.dir, 
                                 "vqsr/pi/autosomes/chr8/chr8.gatk.called.females.vqsr.sv.array.sites.pi.per.site.txt", 
                                 sep = ""),
                           header = T)
chr.8.f.vqsr$chromosome <- 'Chr 8 (Females)'
chr.8.f.vqsr$type <- 'VQSR'

chr.8.m.vqsr <- read.table(paste(results.main.dir, 
                                 "vqsr/pi/autosomes/chr8/chr8.gatk.called.males.vqsr.sv.array.sites.pi.per.site.txt", 
                                 sep = ""),
                           header = T)
chr.8.m.vqsr$chromosome <- 'Chr 8 (Males)'
chr.8.m.vqsr$type <- 'VQSR'

chr.X.f.vqsr <- read.table(paste(results.main.dir, 
                                 "vqsr/pi/sex_chromosomes/chrX/chrX.gatk.called.females.vqsr.sv.array.sites.noPARs.noXTR.pi.per.site.txt", 
                                 sep = ""),
                           header = T)
chr.X.f.vqsr$chromosome <- 'X chr (Females)'
chr.X.f.vqsr$type <- 'VQSR'

chr.X.m.hap.vqsr <- read.table(paste(results.main.dir, 
                                     "vqsr/pi/sex_chromosomes/chrX/chrX.gatk.called.haploid.vqsr.sv.array.sites.noPARs.noXTR.pi.per.site.txt", 
                                     sep = ""),
                               header = T)
chr.X.m.hap.vqsr$chromosome <- 'X chr (Males haploid)'
chr.X.m.hap.vqsr$type <- 'VQSR'

chr.X.m.dip.vqsr <- read.table(paste(results.main.dir, 
                                     "vqsr/pi/sex_chromosomes/chrX/chrX.gatk.called.diploid.vqsr.sv.array.sites.noPARs.noXTR.pi.per.site.txt", 
                                     sep = ""),
                               header = T)
chr.X.m.dip.vqsr$chromosome <- 'X chr (Males diploid)'
chr.X.m.dip.vqsr$type <- 'VQSR'


# Hard filtered files 
chr.8.f.hardf <- read.table(paste(results.main.dir, 
                                  "hard_filter/pi/autosomes/chr8/chr8.gatk.called.females.snps.gatkHardFilter.sv.array.sites.pi.per.site.txt", 
                                  sep = ""),
                            header = T)
chr.8.f.hardf$chromosome <- 'Chr 8 (Females)'
chr.8.f.hardf$type <- 'Hard Filter'

chr.8.m.hardf <- read.table(paste(results.main.dir, 
                                  "hard_filter/pi/autosomes/chr8/chr8.gatk.called.males.snps.gatkHardFilter.sv.array.sites.pi.per.site.txt", 
                                  sep = ""),
                            header = T)
chr.8.m.hardf$chromosome <- 'Chr 8 (Males)'
chr.8.m.hardf$type <- 'Hard Filter'

chr.X.f.hardf <- read.table(paste(results.main.dir, 
                                  "hard_filter/pi/sex_chromosomes/chrX/chrX.gatk.called.females.snps.gatkHardFilter.sv.array.sites.noPARs.noXTR.pi.per.site.txt", 
                                  sep = ""),
                            header = T)
chr.X.f.hardf$chromosome <- 'X chr (Females)'
chr.X.f.hardf$type <- 'Hard Filter'

chr.X.m.hap.hardf <- read.table(paste(results.main.dir, 
                                      "hard_filter/pi/sex_chromosomes/chrX/chrX.gatk.called.haploid.snps.gatkHardFilter.sv.array.sites.noPARs.noXTR.pi.per.site.txt", 
                                      sep = ""),
                                header = T)
chr.X.m.hap.hardf$chromosome <- 'X chr (Males haploid)'
chr.X.m.hap.hardf$type <- 'Hard Filter'

chr.X.m.dip.hardf <- read.table(paste(results.main.dir, 
                                      "hard_filter/pi/sex_chromosomes/chrX/chrX.gatk.called.diploid.snps.gatkHardFilter.sv.array.sites.noPARs.noXTR.pi.per.site.txt", 
                                      sep = ""),
                                header = T)
chr.X.m.dip.hardf$chromosome <- 'X chr (Males diploid)'
chr.X.m.dip.hardf$type <- 'Hard Filter'

# Merge all results into 1 data frame
full.dat <- rbind(chr.8.f.raw, chr.8.m.raw, chr.X.f.raw, chr.X.m.dip.raw, chr.X.m.hap.raw,
                  chr.8.f.vqsr, chr.8.m.vqsr, chr.X.f.vqsr, chr.X.m.dip.vqsr, chr.X.m.hap.vqsr,
                  chr.8.f.hardf, chr.8.m.hardf, chr.X.f.hardf, chr.X.m.dip.hardf, chr.X.m.hap.hardf)

# Plot # 
# Order data by type - I want it in order from left to right raw, vqsr, hard filter
full.dat$type <- factor(full.dat$type, levels = c("Raw", "VQSR", "Hard Filter"))

# Plot
png(paste(results.main.dir, "merged/pi/pi.array.sites.noPARs.noXTR.results.png", sep = ""), width = 11, height = 8, res = 300, units = "in")

# To make y labels nice looking (not scientific notation)
point <- format_format(big.mark = " ", decimal.mark = ".", scientific = FALSE)

ggplot(full.dat, aes(x = chromosome, y = pi_per_site, fill = type)) +
  theme(legend.title=element_blank(), text = element_text(size=14), plot.title = element_text(hjust = 0.5, face="bold")) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  labs(title="Pi comparisons - array sites", x="Chromosome", y="Pi") +
  scale_fill_viridis(discrete=TRUE) +
  scale_y_continuous(labels = point)

dev.off()

#------------------------------------------------------------------------------
# X/A - array sites
# RAW #
# Females no filtering X/A - A is chr 8
x.a.f.raw <- as.data.frame(chr.X.f.raw$pi_per_site / chr.8.f.raw$pi_per_site)
x.a.f.raw.col.n <- colnames(x.a.f.raw)
names(x.a.f.raw)[names(x.a.f.raw) == x.a.f.raw.col.n] <- "X.A"
x.a.f.raw$chromosome <- 'Females'
x.a.f.raw$type <- 'Raw'

# Males diploid no filtering X/A - A is chr 8
x.a.m.dip.raw <- as.data.frame(chr.X.m.dip.raw$pi_per_site / chr.8.m.raw$pi_per_site)
x.a.m.dip.raw.col.n <- colnames(x.a.m.dip.raw)
names(x.a.m.dip.raw)[names(x.a.m.dip.raw) == x.a.m.dip.raw.col.n] <- "X.A"
x.a.m.dip.raw$chromosome <- 'Males (diploid called X)'
x.a.m.dip.raw$type <- 'Raw'

# Males haploid no filtering X/A - A is chr 8
x.a.m.hap.raw <- as.data.frame(chr.X.m.hap.raw$pi_per_site / chr.8.m.raw$pi_per_site)
x.a.m.hap.raw.col.n <- colnames(x.a.m.hap.raw)
names(x.a.m.hap.raw)[names(x.a.m.hap.raw) == x.a.m.hap.raw.col.n] <- "X.A"
x.a.m.hap.raw$chromosome <- 'Males (haploid called X)'
x.a.m.hap.raw$type <- 'Raw'

# VQSR #
# Females vqsr X/A - A is chr 8
x.a.f.vqsr <- as.data.frame(chr.X.f.vqsr$pi_per_site / chr.8.f.vqsr$pi_per_site)
x.a.f.vqsr.col.n <- colnames(x.a.f.vqsr)
names(x.a.f.vqsr)[names(x.a.f.vqsr) == x.a.f.vqsr.col.n] <- "X.A"
x.a.f.vqsr$chromosome <- 'Females'
x.a.f.vqsr$type <- 'VQSR'

# Males diploid vqsr X/A - A is chr 8
x.a.m.dip.vqsr <- as.data.frame(chr.X.m.dip.vqsr$pi_per_site / chr.8.m.vqsr$pi_per_site)
x.a.m.dip.vqsr.col.n <- colnames(x.a.m.dip.vqsr)
names(x.a.m.dip.vqsr)[names(x.a.m.dip.vqsr) == x.a.m.dip.vqsr.col.n] <- "X.A"
x.a.m.dip.vqsr$chromosome <- 'Males (diploid called X)'
x.a.m.dip.vqsr$type <- 'VQSR'

# Males haploid vqsr X/A - A is chr 8
x.a.m.hap.vqsr <- as.data.frame(chr.X.m.hap.vqsr$pi_per_site / chr.8.m.vqsr$pi_per_site)
x.a.m.hap.vqsr.col.n <- colnames(x.a.m.hap.vqsr)
names(x.a.m.hap.vqsr)[names(x.a.m.hap.vqsr) == x.a.m.hap.vqsr.col.n] <- "X.A"
x.a.m.hap.vqsr$chromosome <- 'Males (haploid called X)'
x.a.m.hap.vqsr$type <- 'VQSR'

# HARD FILTER #
# Females Hard Filter X/A - A is chr 8
x.a.f.hardf <- as.data.frame(chr.X.f.hardf$pi_per_site / chr.8.f.hardf$pi_per_site)
x.a.f.hardf.col.n <- colnames(x.a.f.hardf)
names(x.a.f.hardf)[names(x.a.f.hardf) == x.a.f.hardf.col.n] <- "X.A"
x.a.f.hardf$chromosome <- 'Females'
x.a.f.hardf$type <- 'Hard Filter'

# Males diploid Hard Filter X/A - A is chr 8
x.a.m.dip.hardf <- as.data.frame(chr.X.m.dip.hardf$pi_per_site / chr.8.m.hardf$pi_per_site)
x.a.m.dip.hardf.col.n <- colnames(x.a.m.dip.hardf)
names(x.a.m.dip.hardf)[names(x.a.m.dip.hardf) == x.a.m.dip.hardf.col.n] <- "X.A"
x.a.m.dip.hardf$chromosome <- 'Males (diploid called X)'
x.a.m.dip.hardf$type <- 'Hard Filter'

# Males haploid Hard Filter X/A - A is chr 8
x.a.m.hap.hardf <- as.data.frame(chr.X.m.hap.hardf$pi_per_site / chr.8.m.hardf$pi_per_site)
x.a.m.hap.hardf.col.n <- colnames(x.a.m.hap.hardf)
names(x.a.m.hap.hardf)[names(x.a.m.hap.hardf) == x.a.m.hap.hardf.col.n] <- "X.A"
x.a.m.hap.hardf$chromosome <- 'Males (haploid called X)'
x.a.m.hap.hardf$type <- 'Hard Filter'

x.a.full.dat <- rbind(x.a.f.raw, x.a.m.dip.raw, x.a.m.hap.raw,
                      x.a.f.vqsr, x.a.m.dip.vqsr, x.a.m.hap.vqsr,
                      x.a.f.hardf, x.a.m.dip.hardf, x.a.m.hap.hardf)
# Order data by type - I want it in order from left to right raw, vqsr, hard filter
x.a.full.dat$type <- factor(x.a.full.dat$type, levels = c("Raw", "VQSR", "Hard Filter"))

# To make y labels nice looking (not scientific notation)
point <- format_format(big.mark = " ", decimal.mark = ".", scientific = FALSE)

# Plot
png(paste(results.main.dir, "merged/pi/X.A.array.sites.noPARs.noXTR.results.png", sep = ""), width = 11, height = 8, res = 300, units = "in")

ggplot(x.a.full.dat, aes(x = chromosome, y = X.A, fill = type)) +
  theme(legend.title=element_blank(), text = element_text(size=14), plot.title = element_text(hjust = 0.5, face="bold")) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  labs(title="X/A - array sites", x="", y="Pi") +
  scale_fill_viridis(discrete=TRUE) +
  scale_y_continuous(labels = point)

dev.off()
