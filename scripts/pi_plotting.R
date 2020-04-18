# Angela Oill
# 2020-04-16
# Plot pi results for variant calling project
# USAGE:


# Load libraries #
library(ggplot2)
library(viridis)
require(scales)


# Read in all data
# I need to read in a bunch of data 
# chr 8 F, chr 8 M, chr X F, chr X M dip, chr X M hap - Raw
# chr 8 F, chr 8 M, chr X F, chr X M dip, chr X M hap - VQSR
# chr 8 F, chr 8 M, chr X F, chr X M dip, chr X M hap - Hard filter
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
                                "pre_filter/pi/sex_chromosomes/chrX/chrX.gatk.called.raw.females.pi.per.site.txt", 
                                sep = ""),
                          header = T)
chr.X.f.raw$chromosome <- 'X chr (Females)'
chr.X.f.raw$type <- 'Raw'


chr.X.m.hap.raw <- read.table(paste(results.main.dir, 
                                "pre_filter/pi/sex_chromosomes/chrX/chrX.gatk.called.raw.haploid.pi.per.site.txt", 
                                sep = ""),
                          header = T)
chr.X.m.hap.raw$chromosome <- 'X chr (Males haploid)'
chr.X.m.hap.raw$type <- 'Raw'

chr.X.m.dip.raw <- read.table(paste(results.main.dir, 
                                    "pre_filter/pi/sex_chromosomes/chrX/chrX.gatk.called.raw.diploid.pi.per.site.txt", 
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
                                "vqsr/pi/sex_chromosomes/chrX/chrX.gatk.called.females.vqsr.sv.pi.per.site.txt", 
                                sep = ""),
                          header = T)
chr.X.f.vqsr$chromosome <- 'X chr (Females)'
chr.X.f.vqsr$type <- 'VQSR'

chr.X.m.hap.vqsr <- read.table(paste(results.main.dir, 
                                    "vqsr/pi/sex_chromosomes/chrX/chrX.gatk.called.haploid.vqsr.sv.pi.per.site.txt", 
                                    sep = ""),
                              header = T)
chr.X.m.hap.vqsr$chromosome <- 'X chr (Males haploid)'
chr.X.m.hap.vqsr$type <- 'VQSR'

chr.X.m.dip.vqsr <- read.table(paste(results.main.dir, 
                                    "vqsr/pi/sex_chromosomes/chrX/chrX.gatk.called.diploid.vqsr.sv.pi.per.site.txt", 
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
                                 "hard_filter/pi/sex_chromosomes/chrX/chrX.gatk.called.females.snps.gatkHardFilter.sv.pi.per.site.txt", 
                                 sep = ""),
                           header = T)
chr.X.f.hardf$chromosome <- 'X chr (Females)'
chr.X.f.hardf$type <- 'Hard Filter'

chr.X.m.hap.hardf <- read.table(paste(results.main.dir, 
                                     "hard_filter/pi/sex_chromosomes/chrX/chrX.gatk.called.haploid.snps.gatkHardFilter.sv.pi.per.site.txt", 
                                     sep = ""),
                               header = T)
chr.X.m.hap.hardf$chromosome <- 'X chr (Males haploid)'
chr.X.m.hap.hardf$type <- 'Hard Filter'

chr.X.m.dip.hardf <- read.table(paste(results.main.dir, 
                                     "hard_filter/pi/sex_chromosomes/chrX/chrX.gatk.called.diploid.snps.gatkHardFilter.sv.pi.per.site.txt", 
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
png(paste(results.main.dir, "merged/pi/pi.all.called.variants.results.png", sep = ""), width = 11, height = 8, res = 300, units = "in")

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
                                "pre_filter/pi/sex_chromosomes/chrX/chrX.gatk.called.raw.females.array.sites.pi.per.site.txt", 
                                sep = ""),
                          header = T)
chr.X.f.raw$chromosome <- 'X chr (Females)'
chr.X.f.raw$type <- 'Raw'


chr.X.m.hap.raw <- read.table(paste(results.main.dir, 
                                    "pre_filter/pi/sex_chromosomes/chrX/chrX.gatk.called.raw.haploid.array.sites.pi.per.site.txt", 
                                    sep = ""),
                              header = T)
chr.X.m.hap.raw$chromosome <- 'X chr (Males haploid)'
chr.X.m.hap.raw$type <- 'Raw'

chr.X.m.dip.raw <- read.table(paste(results.main.dir, 
                                    "pre_filter/pi/sex_chromosomes/chrX/chrX.gatk.called.raw.diploid.array.sites.pi.per.site.txt", 
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
                                 "vqsr/pi/sex_chromosomes/chrX/chrX.gatk.called.females.vqsr.sv.array.sites.pi.per.site.txt", 
                                 sep = ""),
                           header = T)
chr.X.f.vqsr$chromosome <- 'X chr (Females)'
chr.X.f.vqsr$type <- 'VQSR'

chr.X.m.hap.vqsr <- read.table(paste(results.main.dir, 
                                     "vqsr/pi/sex_chromosomes/chrX/chrX.gatk.called.haploid.vqsr.sv.array.sites.pi.per.site.txt", 
                                     sep = ""),
                               header = T)
chr.X.m.hap.vqsr$chromosome <- 'X chr (Males haploid)'
chr.X.m.hap.vqsr$type <- 'VQSR'

chr.X.m.dip.vqsr <- read.table(paste(results.main.dir, 
                                     "vqsr/pi/sex_chromosomes/chrX/chrX.gatk.called.diploid.vqsr.sv.array.sites.pi.per.site.txt", 
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
                                  "hard_filter/pi/sex_chromosomes/chrX/chrX.gatk.called.females.snps.gatkHardFilter.sv.array.sites.pi.per.site.txt", 
                                  sep = ""),
                            header = T)
chr.X.f.hardf$chromosome <- 'X chr (Females)'
chr.X.f.hardf$type <- 'Hard Filter'

chr.X.m.hap.hardf <- read.table(paste(results.main.dir, 
                                      "hard_filter/pi/sex_chromosomes/chrX/chrX.gatk.called.haploid.snps.gatkHardFilter.sv.array.sites.pi.per.site.txt", 
                                      sep = ""),
                                header = T)
chr.X.m.hap.hardf$chromosome <- 'X chr (Males haploid)'
chr.X.m.hap.hardf$type <- 'Hard Filter'

chr.X.m.dip.hardf <- read.table(paste(results.main.dir, 
                                      "hard_filter/pi/sex_chromosomes/chrX/chrX.gatk.called.diploid.snps.gatkHardFilter.sv.array.sites.pi.per.site.txt", 
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
png(paste(results.main.dir, "merged/pi/pi.array.sites.results.png", sep = ""), width = 11, height = 8, res = 300, units = "in")

# To make y labels nice looking (not scientific notation)
point <- format_format(big.mark = " ", decimal.mark = ".", scientific = FALSE)

ggplot(full.dat, aes(x = chromosome, y = pi_per_site, fill = type)) +
  theme(legend.title=element_blank(), text = element_text(size=14), plot.title = element_text(hjust = 0.5, face="bold")) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  labs(title="Pi comparisons - array sites", x="Chromosome", y="Pi") +
  scale_fill_viridis(discrete=TRUE) +
  scale_y_continuous(labels = point)

dev.off()