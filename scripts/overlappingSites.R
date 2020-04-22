# Angela Oill
# 2020-04-21
# Make venn diagrams to show overlapping sites and unique sites for each
# filtering strategy.
# USAGE: Rscript overlappingSites.R <chr.pos.vqsr.txt> <chr.pos.hard.filter.txt> <output.png>


# 1. Load libraries #
library(VennDiagram)
library(viridis)


# 2. Command Line Variables #
args = commandArgs(trailingOnly=TRUE)
vqsr.dir.fn = args[1] # file name and path to chr.pos file for vqsr filtering strategy
hardf.dir.fn = args[2] # file name and path to chr.pos file for hard filtering strategy
out.dir.fn = args[3] # file name and path to output png file. 


# 3. Read in data # 
#vqsr <- unlist(read.table("chr8.gatk.called.females.vqsr.sv.chr.pos.txt"))
vqsr <- unlist(read.table(vqsr.dir.fn))
#hardf <- unlist(read.table("chr8.gatk.called.females.snps.gatkHardFilter.sv.chr.pos.txt"))
hardf <- unlist(read.table(hardf.dir.fn))


# 4. Plot #
# color based on virdis
vcolors <- viridis(3) #first entry is raw
vcolors.plot <- vcolors[2:3]
venn.diagram(
  x = list(vqsr, hardf),
  category.names = c("VQSR" , "Hard Filter"),
  filename = out.dir.fn,
  output=TRUE,
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = vcolors.plot
)
