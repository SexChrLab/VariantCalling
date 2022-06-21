#!/bin/bash
#SBATCH --job-name=metrics_per_window_20samples # Job name
#SBATCH -o slurm.metrics_per_window_20samples.out
#SBATCH -e slurm.metrics_per_window_20samples.err
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mail-user=amtarave@asu.edu # send-to address

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/scripts

# Command (for autosomes):
# cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/EUR/males/autos
# cat *_chr8_autos_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_samples_chr8_autos_golden_vs_called_performance_metrics.txt

################################################################################
# MALES - SCC #
################################################################################
#-----------------------------#
# DIPLOID CALLING - Autosomes #
#-----------------------------#
ancestry="EUR"
sex="males"
ploidy="diploid"

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/autos

for i in {1..22}
do
	echo "chr${i}"

  cat *_chr${i}_autos_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chr${i}_${ploidy}_golden_vs_called_performance_metrics.txt
done

#------------------------------------#
# DIPLOID CALLING - Sex chromosomes  #
#------------------------------------#

# PARs
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/chrX_PARs
cat *_chrX_PARs_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrX_PARs_${ploidy}_golden_vs_called_performance_metrics.txt

# ChrX non-PARs
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/chrX_nonPARs
cat *_chrX_nonPARs_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrX_nonPARs_${ploidy}_golden_vs_called_performance_metrics.txt

# ChrY non-PARs
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/chrY
cat *_chrY_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrY_${ploidy}_golden_vs_called_performance_metrics.txt

# ChrM
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/chrM
cat *_chrM_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrM_${ploidy}_golden_vs_called_performance_metrics.txt

#-----------------#
# HAPLOID CALLING #
#-----------------#
ploidy="haploid"

# ChrX non-PARs
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/${ploidy}/compare_VCFs/${ancestry}/${sex}/chrX_nonPARs
cat *_chrX_nonPARs_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrX_nonPARs_${ploidy}_golden_vs_called_performance_metrics.txt

# ChrY non-PARs
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/${ploidy}/compare_VCFs/${ancestry}/${sex}/chrY
cat *_chrY_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrY_${ploidy}_golden_vs_called_performance_metrics.txt

# ChrM
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/${ploidy}/compare_VCFs/${ancestry}/${sex}/chrM
cat *_chrM_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrM_${ploidy}_golden_vs_called_performance_metrics.txt




################################################################################
# MALES - DEFAULT #
################################################################################
#-----------------------------#
# DIPLOID CALLING - Autosomes #
#-----------------------------#
ancestry="EUR"
sex="males"
ploidy="diploid"

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/default/autos

for i in {1..22}
do
	echo "chr${i}"

  cat *_chr${i}_autos_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chr${i}_${ploidy}_default_golden_vs_called_performance_metrics.txt
done


#------------------------------------#
# DIPLOID CALLING - Sex chromosomes  #
#------------------------------------#

# PARs
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/default/chrX_PARs
cat *_chrX_PARs_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrX_PARs_${ploidy}_default_golden_vs_called_performance_metrics.txt

# ChrX non-PARs
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/default/chrX_nonPARs
cat *_chrX_nonPARs_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrX_nonPARs_${ploidy}_default_golden_vs_called_performance_metrics.txt

# ChrY non-PARs
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/default/chrY
cat *_chrY_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrY_${ploidy}_default_golden_vs_called_performance_metrics.txt

# ChrM
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/default/chrM
cat *_chrM_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrM_${ploidy}_default_golden_vs_called_performance_metrics.txt

#-----------------#
# HAPLOID CALLING #
#-----------------#
ploidy="haploid"

# ChrX non-PARs
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/${ploidy}/compare_VCFs/${ancestry}/${sex}/default/chrX_nonPARs
cat *_chrX_nonPARs_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrX_nonPARs_${ploidy}_default_golden_vs_called_performance_metrics.txt

# ChrY non-PARs
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/${ploidy}/compare_VCFs/${ancestry}/${sex}/default/chrY
cat *_chrY_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrY_${ploidy}_default_golden_vs_called_performance_metrics.txt

# ChrM
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/${ploidy}/compare_VCFs/${ancestry}/${sex}/default/chrM
cat *_chrM_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrM_${ploidy}_default_golden_vs_called_performance_metrics.txt



################################################################################
# FEMALES #
################################################################################

#-----------------------------#
# DIPLOID CALLING - Autosomes #
#-----------------------------#
ancestry="EUR"
sex="females"
ploidy="diploid"

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/autos

for i in {1..22}
do
	echo "chr${i}"

  cat *_chr${i}_autos_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chr${i}_${ploidy}_golden_vs_called_performance_metrics.txt
done

#------------------------------------#
# DIPLOID CALLING - Sex chromosomes  #
#------------------------------------#

# PARs
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/default/chrX_PARs
cat *_chrX_PARs_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrX_PARs_${ploidy}_golden_vs_called_performance_metrics.txt

# ChrX non-PARs
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/default/chrX_nonPARs
cat *_chrX_nonPARs_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrX_nonPARs_${ploidy}_golden_vs_called_performance_metrics.txt


# ChrM
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/default/chrM
cat *_chrM_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrM_${ploidy}_golden_vs_called_performance_metrics.txt


#-----------------#
# HAPLOID CALLING #
#-----------------#
ploidy="haploid"

# ChrM
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/${ploidy}/compare_VCFs/${ancestry}/${sex}/chrM
cat *_chrM_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrM_${ploidy}_golden_vs_called_performance_metrics.txt



################################################################################
# FEMALES - DEFAULT #
################################################################################
#-----------------------------#
# DIPLOID CALLING - Autosomes #
#-----------------------------#
ancestry="EUR"
sex="females"
ploidy="diploid"

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/default/autos

for i in {1..22}
do
	echo "chr${i}"

  cat *_chr${i}_autos_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chr${i}_${ploidy}_default_golden_vs_called_performance_metrics.txt
done


#------------------------------------#
# DIPLOID CALLING - Sex chromosomes  #
#------------------------------------#

# PARs
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/default/chrX_PARs
cat *_chrX_PARs_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrX_PARs_${ploidy}_default_golden_vs_called_performance_metrics.txt

# ChrX non-PARs
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/default/chrX_nonPARs
cat *_chrX_nonPARs_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrX_nonPARs_${ploidy}_default_golden_vs_called_performance_metrics.txt


# ChrM
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/default/chrM
cat *_chrM_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrM_${ploidy}_default_golden_vs_called_performance_metrics.txt


#-----------------#
# HAPLOID CALLING #
#-----------------#
ploidy="haploid"

# ChrM
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/${ploidy}/compare_VCFs/${ancestry}/${sex}/default/chrM
cat *_chrM_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrM_${ploidy}_default_golden_vs_called_performance_metrics.txt





################################################################################
################################################################################
################################################################################
# 20 SAMPLES #
################################################################################
################################################################################
################################################################################


################################################################################
# MALES - SCC - 20 samples #
################################################################################
#-----------------------------#
# DIPLOID CALLING - Autosomes #
#-----------------------------#
ancestry="EUR"
sex="males"
ploidy="diploid"

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs_20_samples/${ancestry}/${sex}/autos

#for i in {1..22}
for i in {8..8}
do
	echo "chr${i}"

  cat *_chr${i}_autos_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chr${i}_${ploidy}_golden_vs_called_performance_metrics.txt
done


#------------------------------------#
# DIPLOID CALLING - Sex chromosomes  #
#------------------------------------#

# PARs
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs_20_samples/${ancestry}/${sex}/chrX_PARs
cat *_chrX_PARs_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrX_PARs_${ploidy}_golden_vs_called_performance_metrics.txt

# ChrX non-PARs
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs_20_samples/${ancestry}/${sex}/chrX_nonPARs
cat *_chrX_nonPARs_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrX_nonPARs_${ploidy}_golden_vs_called_performance_metrics.txt

# ChrY non-PARs
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs_20_samples/${ancestry}/${sex}/chrY
cat *_chrY_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrY_${ploidy}_golden_vs_called_performance_metrics.txt

# ChrM
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs_20_samples/${ancestry}/${sex}/chrM
cat *_chrM_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrM_${ploidy}_golden_vs_called_performance_metrics.txt

#-----------------#
# HAPLOID CALLING #
#-----------------#
ploidy="haploid"

# ChrX non-PARs
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/${ploidy}/compare_VCFs_20_samples/${ancestry}/${sex}/chrX_nonPARs
cat *_chrX_nonPARs_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrX_nonPARs_${ploidy}_golden_vs_called_performance_metrics.txt

# ChrY non-PARs
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/${ploidy}/compare_VCFs_20_samples/${ancestry}/${sex}/chrY
cat *_chrY_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrY_${ploidy}_golden_vs_called_performance_metrics.txt

# ChrM
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/${ploidy}/compare_VCFs_20_samples/${ancestry}/${sex}/chrM
cat *_chrM_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrM_${ploidy}_golden_vs_called_performance_metrics.txt



################################################################################
# FEMALES - SCC - 20 samples #
################################################################################
#-----------------------------#
# DIPLOID CALLING - Autosomes #
#-----------------------------#
ancestry="EUR"
sex="females"
ploidy="diploid"

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs_20_samples/${ancestry}/${sex}/autos

#for i in {1..22}
for i in {8..8}
do
	echo "chr${i}"

  cat *_chr${i}_autos_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chr${i}_${ploidy}_golden_vs_called_performance_metrics.txt
done


#------------------------------------#
# DIPLOID CALLING - Sex chromosomes  #
#------------------------------------#
# PARs
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs_20_samples/${ancestry}/${sex}/chrX_PARs
cat *_chrX_PARs_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrX_PARs_${ploidy}_golden_vs_called_performance_metrics.txt

# ChrX non-PARs
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs_20_samples/${ancestry}/${sex}/chrX_nonPARs
cat *_chrX_nonPARs_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrX_nonPARs_${ploidy}_golden_vs_called_performance_metrics.txt

# ChrM
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs_20_samples/${ancestry}/${sex}/chrM
cat *_chrM_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ${ancestry}_${sex}_chrM_${ploidy}_golden_vs_called_performance_metrics.txt

#-----------------#
# HAPLOID CALLING #
#-----------------#
# havent run haploid on mtDNA for females yet



# DEFAULT hasnt been run yet for 20 males and females
