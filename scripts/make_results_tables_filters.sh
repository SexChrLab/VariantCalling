#!/bin/bash
#SBATCH --job-name=metrics_per_window_20samples # Job name
#SBATCH -o slurm.metrics_per_window_20samples.out
#SBATCH -e slurm.metrics_per_window_20samples.err
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mail-user=amtarave@asu.edu # send-to address

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/scripts


################################################################################
# MALES #
################################################################################
#-----------------#
# DIPLOID CALLING #
#-----------------#
# QD
filter="QD"
ancestry="EUR"
sex="males"
#for i in 1.0 1.5 2.0 12.0 16.0 20.0 28.0
for i in 1.0 1.5
do
	echo "${i}"
  # Chr8 #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/autos/
  cat *_chr8_autos_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chr8_autos_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX non-PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_nonPARs/
  cat *_chrX_nonPARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrX_nonPARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_PARs/
  cat *_chrX_PARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrX_PARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrY #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrY/
  cat *_chrY_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrY_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrM_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done

# QUAL
filter="QUAL"
ancestry="EUR"
sex="males"
for i in 30.0
do
	echo "${i}"
  # Chr8 #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/autos/
  cat *_chr8_autos_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chr8_autos_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX non-PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_nonPARs/
  cat *_chrX_nonPARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrX_nonPARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_PARs/
  cat *_chrX_PARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrX_PARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrY #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrY/
  cat *_chrY_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrY_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrM_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done

# SOR
filter="SOR"
ancestry="EUR"
sex="males"
for i in 3.0
do
	echo "${i}"
  # Chr8 #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/autos/
  cat *_chr8_autos_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chr8_autos_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX non-PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_nonPARs/
  cat *_chrX_nonPARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrX_nonPARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_PARs/
  cat *_chrX_PARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrX_PARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrY #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrY/
  cat *_chrY_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrY_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrM_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done


# FS
filter="FS"
ancestry="EUR"
sex="males"
for i in 60.0
do
	echo "${i}"
  # Chr8 #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/autos/
  cat *_chr8_autos_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chr8_autos_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX non-PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_nonPARs/
  cat *_chrX_nonPARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrX_nonPARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_PARs/
  cat *_chrX_PARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrX_PARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrY #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrY/
  cat *_chrY_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrY_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrM_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done


# MQ
filter="MQ"
ancestry="EUR"
sex="males"
for i in 20.0 30.0 40.0 50.0 60.0
do
	echo "${i}"
  # Chr8 #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/autos/
  cat *_chr8_autos_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chr8_autos_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX non-PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_nonPARs/
  cat *_chrX_nonPARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrX_nonPARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_PARs/
  cat *_chrX_PARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrX_PARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrY #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrY/
  cat *_chrY_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrY_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrM_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done


# MQRankSum
filter="MQRankSum"
ancestry="EUR"
sex="males"
for i in -12.5
do
	echo "${i}"
  # Chr8 #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/autos/
  cat *_chr8_autos_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chr8_autos_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX non-PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_nonPARs/
  cat *_chrX_nonPARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrX_nonPARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_PARs/
  cat *_chrX_PARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrX_PARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrY #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrY/
  cat *_chrY_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrY_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrM_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done


# ReadPosRankSum
filter="ReadPosRankSum"
ancestry="EUR"
sex="males"
for i in -8.0
do
	echo "${i}"
  # Chr8 #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/autos/
  cat *_chr8_autos_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chr8_autos_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX non-PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_nonPARs/
  cat *_chrX_nonPARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrX_nonPARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_PARs/
  cat *_chrX_PARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrX_PARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrY #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrY/
  cat *_chrY_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrY_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrM_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done

# AN
filter="AN"
ancestry="EUR"
sex="males"
#for i in 20 15 10 5 4 3 2 1
for i in 4 3 2 1
do
	echo "${i}"
  # Chr8 #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/autos/
  cat *_chr8_autos_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chr8_autos_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX non-PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_nonPARs/
  cat *_chrX_nonPARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrX_nonPARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_PARs/
  cat *_chrX_PARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrX_PARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrY #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrY/
  cat *_chrY_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrY_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrM_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done

# DP
filter="DP"
ancestry="EUR"
sex="males"
#for i in 20 15 10 5 4 3 2 1
for i in 6 7 8 9
do
	echo "${i}"
  # Chr8 #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/autos/
  cat *_chr8_autos_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chr8_autos_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX non-PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_nonPARs/
  cat *_chrX_nonPARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrX_nonPARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_PARs/
  cat *_chrX_PARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrX_PARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrY #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrY/
  cat *_chrY_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrY_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrM_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done

# DP INFO
filter="DP_INFO"
ancestry="EUR"
sex="males"
for i in DP67and201 DP23and69
do
	echo "${i}"
  # Chr8 #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/autos/
  cat *_chr8_autos_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chr8_autos_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX non-PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_nonPARs/
  cat *_chrX_nonPARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrX_nonPARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_PARs/
  cat *_chrX_PARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrX_PARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrY #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrY/
  cat *_chrY_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrY_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrM #
  #cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  #cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrM_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done

#-----------------#
# HAPLOID CALLING #
#-----------------#
# QD
filter="QD"
ancestry="EUR"
sex="males"
for i in 1.0 1.5 2.0 12.0 16.0 20.0 28.0
do
	echo "${i}"

  # ChrX non-PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_nonPARs/
  cat *_chrX_nonPARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrX_nonPARs_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrY #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrY/
  cat *_chrY_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrY_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrM_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done

# QUAL
filter="QUAL"
ancestry="EUR"
sex="males"
for i in 30.0
do
	echo "${i}"

  # ChrX non-PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_nonPARs/
  cat *_chrX_nonPARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrX_nonPARs_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrY #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrY/
  cat *_chrY_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrY_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrM_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done

# SOR
filter="SOR"
ancestry="EUR"
sex="males"
for i in 3.0
do
	echo "${i}"

  # ChrX non-PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_nonPARs/
  cat *_chrX_nonPARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrX_nonPARs_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrY #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrY/
  cat *_chrY_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrY_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrM_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done


# FS
filter="FS"
ancestry="EUR"
sex="males"
for i in 60.0
do
	echo "${i}"

  # ChrX non-PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_nonPARs/
  cat *_chrX_nonPARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrX_nonPARs_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrY #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrY/
  cat *_chrY_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrY_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrM_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done

# MQ
filter="MQ"
ancestry="EUR"
sex="males"
for i in 20.0 30.0 40.0 50.0 60.0
do
	echo "${i}"

  # ChrX non-PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_nonPARs/
  cat *_chrX_nonPARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrX_nonPARs_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrY #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrY/
  cat *_chrY_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrY_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrM_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done

# MQRankSum
filter="MQRankSum"
ancestry="EUR"
sex="males"
for i in -12.5
do
	echo "${i}"

  # ChrX non-PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_nonPARs/
  cat *_chrX_nonPARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrX_nonPARs_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrY #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrY/
  cat *_chrY_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrY_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrM_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done

# ReadPosRankSum
filter="ReadPosRankSum"
ancestry="EUR"
sex="males"
for i in -8.0
do
	echo "${i}"

  # ChrX non-PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_nonPARs/
  cat *_chrX_nonPARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrX_nonPARs_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrY #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrY/
  cat *_chrY_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrY_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrM_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done

# AN
filter="AN"
ancestry="EUR"
sex="males"
#for i in 20 15 10 5 4 3 2 1
for i in 4 3 2 1
do
	echo "${i}"

  # ChrX non-PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_nonPARs/
  cat *_chrX_nonPARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrX_nonPARs_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrY #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrY/
  cat *_chrY_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrY_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrM_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done

# DP
filter="DP"
ancestry="EUR"
sex="males"
#for i in 20 15 10 5 4 3 2 1
for i in 6 7 8 9
do
	echo "${i}"

  # ChrX non-PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_nonPARs/
  cat *_chrX_nonPARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrX_nonPARs_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrY #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrY/
  cat *_chrY_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrY_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrM_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done



# DP
filter="DP_INFO"
ancestry="EUR"
sex="males"
for i in DP67and201 DP23and69
do
	echo "${i}"

  # ChrX non-PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_nonPARs/
  cat *_chrX_nonPARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrX_nonPARs_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrY #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrY/
  cat *_chrY_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrY_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrM #
  #cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  #cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_males_chrM_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done


################################################################################
# FEMALES #
################################################################################

#-----------------#
# DIPLOID CALLING #
#-----------------#
# QD
filter="QD"
ancestry="EUR"
sex="females"
for i in 1.0 1.5 2.0 12.0 16.0 20.0 28.0
do
	echo "${i}"
  # Chr8 #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/autos/
  cat *_chr8_autos_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chr8_autos_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX non-PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_nonPARs/
  cat *_chrX_nonPARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrX_nonPARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_PARs/
  cat *_chrX_PARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrX_PARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt


  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrM_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done

# QUAL
filter="QUAL"
ancestry="EUR"
sex="females"
for i in 30.0
do
	echo "${i}"
  # Chr8 #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/autos/
  cat *_chr8_autos_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chr8_autos_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX non-PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_nonPARs/
  cat *_chrX_nonPARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrX_nonPARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_PARs/
  cat *_chrX_PARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrX_PARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrM_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done

# SOR
filter="SOR"
ancestry="EUR"
sex="females"
for i in 3.0
do
	echo "${i}"
  # Chr8 #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/autos/
  cat *_chr8_autos_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chr8_autos_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX non-PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_nonPARs/
  cat *_chrX_nonPARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrX_nonPARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_PARs/
  cat *_chrX_PARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrX_PARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrM_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done


# FS
filter="FS"
ancestry="EUR"
sex="females"
for i in 60.0
do
	echo "${i}"
  # Chr8 #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/autos/
  cat *_chr8_autos_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chr8_autos_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX non-PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_nonPARs/
  cat *_chrX_nonPARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrX_nonPARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_PARs/
  cat *_chrX_PARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrX_PARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt


  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrM_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done


# MQ
filter="MQ"
ancestry="EUR"
sex="females"
for i in 20.0 30.0 40.0 50.0 60.0
do
	echo "${i}"
  # Chr8 #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/autos/
  cat *_chr8_autos_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chr8_autos_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX non-PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_nonPARs/
  cat *_chrX_nonPARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrX_nonPARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_PARs/
  cat *_chrX_PARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrX_PARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt


  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrM_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done


# MQRankSum
filter="MQRankSum"
ancestry="EUR"
sex="females"
for i in -12.5
do
	echo "${i}"
  # Chr8 #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/autos/
  cat *_chr8_autos_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chr8_autos_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX non-PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_nonPARs/
  cat *_chrX_nonPARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrX_nonPARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_PARs/
  cat *_chrX_PARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrX_PARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrM_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done


# ReadPosRankSum
filter="ReadPosRankSum"
ancestry="EUR"
sex="females"
for i in -8.0
do
	echo "${i}"
  # Chr8 #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/autos/
  cat *_chr8_autos_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chr8_autos_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX non-PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_nonPARs/
  cat *_chrX_nonPARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrX_nonPARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_PARs/
  cat *_chrX_PARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrX_PARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt


  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrM_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done

# AN
filter="AN"
ancestry="EUR"
sex="females"
#for i in 20 15 10 5 4 3 2 1
for i in 4 3 2 1
do
	echo "${i}"
  # Chr8 #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/autos/
  cat *_chr8_autos_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chr8_autos_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX non-PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_nonPARs/
  cat *_chrX_nonPARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrX_nonPARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_PARs/
  cat *_chrX_PARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrX_PARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt


  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrM_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done

# DP
filter="DP"
ancestry="EUR"
sex="females"
#for i in 20 15 10 5 4 3 2 1
for i in 6 7 8 9
do
	echo "${i}"
  # Chr8 #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/autos/
  cat *_chr8_autos_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chr8_autos_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX non-PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_nonPARs/
  cat *_chrX_nonPARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrX_nonPARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_PARs/
  cat *_chrX_PARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrX_PARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrM_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done


filter="DP_INFO"
ancestry="EUR"
sex="females"
for i in DP67and201 DP23and69
do
	echo "${i}"
  # Chr8 #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/autos/
  cat *_chr8_autos_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chr8_autos_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX non-PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_nonPARs/
  cat *_chrX_nonPARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrX_nonPARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrX PARs #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrX_PARs/
  cat *_chrX_PARs_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrX_PARs_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

  # ChrM #
  #cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  #cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrM_diploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done

#-----------------#
# HAPLOID CALLING #
#-----------------#
# QD
filter="QD"
ancestry="EUR"
sex="females"
for i in 1.0 1.5 2.0 12.0 16.0 20.0 28.0
do
  echo "${i}"
  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrM_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done

# QUAL
filter="QUAL"
ancestry="EUR"
sex="females"
for i in 30.0
do
	echo "${i}"

  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrM_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done

# SOR
filter="SOR"
ancestry="EUR"
sex="females"
for i in 3.0
do
	echo "${i}"

  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrM_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done


# FS
filter="FS"
ancestry="EUR"
sex="females"
for i in 60.0
do
	echo "${i}"

  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrM_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done


# MQ
filter="MQ"
ancestry="EUR"
sex="females"
for i in 20.0 30.0 40.0 50.0 60.0
do
	echo "${i}"

  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrM_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done


# MQRankSum
filter="MQRankSum"
ancestry="EUR"
sex="females"
for i in -12.5
do
	echo "${i}"

  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrM_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done


# ReadPosRankSum
filter="ReadPosRankSum"
ancestry="EUR"
sex="females"
for i in -8.0
do
	echo "${i}"

  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrM_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done

# AN
filter="AN"
ancestry="EUR"
sex="females"
for i in 20 15 10 5 4 3 2 1
do
	echo "${i}"

  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrM_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done

# DP
filter="DP"
ancestry="EUR"
sex="females"
for i in 20 15 10 5 4 3 2 1
do
	echo "${i}"

  # ChrM #
  cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/${ancestry}/${sex}/${filter}/chrM/
  cat *_chrM_${filter}_${i}_golden_vs_called_performance_metrics.txt | sort | uniq | sort -r > ALL_${sex}_chrM_haploid_${filter}_${i}_golden_vs_called_performance_metrics.txt

done
