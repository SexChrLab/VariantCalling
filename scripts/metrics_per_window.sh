#!/bin/bash
#SBATCH --job-name=metrics_per_window # Job name
#SBATCH -o slurm.metrics_per_window.out
#SBATCH -e slurm.metrics_per_window.err
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mail-user=amtarave@asu.edu # send-to address

cd /scratch/amtarave/variant_calling_simulations_project/window_analysis/

conda activate variant_calling_simulations_project

# First make window file
bedtools makewindows -g chrY_genome_fix.txt -w 50000 > chrY_genome_50kb_windows.bed # can change this later
bedtools makewindows -g chrX_genome_fix.txt -w 50000 > chrX_genome_50kb_windows.bed
bedtools makewindows -g chr8_genome_fix.txt -w 50000 > chr8_genome_50kb_windows.bed

# TO DO - try this on another autosomes
# bedtools makewindows -g chr15_genome_fix.txt -w 50000 > chr15_genome_50kb_windows.bed

# Then loop through and get counts per 50kb window
###############
# MALES - SCC #
###############
cd /scratch/amtarave/variant_calling_simulations_project/window_analysis/scc/males/

for i in NA06984 NA06986 NA06994 NA07048 NA07051 NA07347 NA07357 NA10851 NA11829 NA11831
do
  echo $i
  for t in FP FN TP
  do
    echo $t
    # Y
    # make a bed file out of one of the FN results file
    #awk '{ print $1"\t"$2-1"\t"$2 }' /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/EUR/males/chrY/${i}_chrY_golden_vs_called_${t}_pos.txt > /data/CEM/wilsonlab/projects/variant_calling_simulations_project/window_analysis/scc/males/${i}_chrY_scc_haploid_golden_vs_called_${t}_pos.bed
    # get counts across Y
    #bedtools coverage -a /scratch/amtarave/variant_calling_simulations_project/window_analysis/chrY_genome_50kb_windows.bed -b /data/CEM/wilsonlab/projects/variant_calling_simulations_project/window_analysis/scc/males/${i}_chrY_scc_haploid_golden_vs_called_${t}_pos.bed -counts > ${i}_chrY_scc_haploid_golden_vs_called_${t}_counts_50kb_windows.txt


    #awk '{ print $1"\t"$2-1"\t"$2 }' /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/EUR/males/chrY/${i}_chrY_golden_vs_called_${t}_pos.txt > /data/CEM/wilsonlab/projects/variant_calling_simulations_project/window_analysis/scc/males/${i}_chrY_scc_diploid_golden_vs_called_${t}_pos.bed
    # get counts across Y
    #bedtools coverage -a /scratch/amtarave/variant_calling_simulations_project/window_analysis/chrY_genome_50kb_windows.bed -b /data/CEM/wilsonlab/projects/variant_calling_simulations_project/window_analysis/scc/males/${i}_chrY_scc_diploid_golden_vs_called_${t}_pos.bed -counts > ${i}_chrY_scc_diploid_golden_vs_called_${t}_counts_50kb_windows.txt


    # X nonPARs
    #awk '{ print $1"\t"$2-1"\t"$2 }' /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/EUR/males/chrX_nonPARs/${i}_chrX_nonPARs_golden_vs_called_${t}_pos.txt > /data/CEM/wilsonlab/projects/variant_calling_simulations_project/window_analysis/scc/males/${i}_chrX_nonPARs_scc_haploid_golden_vs_called_${t}_pos.bed
    # get FN counts across Y
    #bedtools coverage -a /scratch/amtarave/variant_calling_simulations_project/window_analysis/chrX_genome_50kb_windows.bed -b /data/CEM/wilsonlab/projects/variant_calling_simulations_project/window_analysis/scc/males/${i}_chrX_nonPARs_scc_haploid_golden_vs_called_${t}_pos.bed -counts > ${i}_chrX_nonPARs_scc_haploid_golden_vs_called_${t}_counts_50kb_windows.txt

    #awk '{ print $1"\t"$2-1"\t"$2 }' /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/EUR/males/chrX_nonPARs/${i}_chrX_nonPARs_golden_vs_called_${t}_pos.txt > /data/CEM/wilsonlab/projects/variant_calling_simulations_project/window_analysis/scc/males/${i}_chrX_nonPARs_scc_diploid_golden_vs_called_${t}_pos.bed
    # get counts across Y
    #bedtools coverage -a /scratch/amtarave/variant_calling_simulations_project/window_analysis/chrX_genome_50kb_windows.bed -b /data/CEM/wilsonlab/projects/variant_calling_simulations_project/window_analysis/scc/males/${i}_chrX_nonPARs_scc_diploid_golden_vs_called_${t}_pos.bed -counts > ${i}_chrX_nonPARs_scc_diploid_golden_vs_called_${t}_counts_50kb_windows.txt

    # X PARs
    #awk '{ print $1"\t"$2-1"\t"$2 }' /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/EUR/males/chrX_PARs/${i}_chrX_PARs_golden_vs_called_${t}_pos.txt > /data/CEM/wilsonlab/projects/variant_calling_simulations_project/window_analysis/scc/males/${i}_chrX_PARs_scc_diploid_golden_vs_called_${t}_pos.bed
    # get counts across Y
    #bedtools coverage -a /scratch/amtarave/variant_calling_simulations_project/window_analysis/chrX_genome_50kb_windows.bed -b /data/CEM/wilsonlab/projects/variant_calling_simulations_project/window_analysis/scc/males/${i}_chrX_PARs_scc_diploid_golden_vs_called_${t}_pos.bed -counts > ${i}_chrX_PARs_scc_diploid_golden_vs_called_${t}_counts_50kb_windows.txt

    # chr8
    # make bed files
    #awk '{ print $1"\t"$2-1"\t"$2 }' /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/EUR/males/autos/${i}_chr8_autos_golden_vs_called_${t}_pos.txt > /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/EUR/males/autos/${i}_chr8_autos_scc_diploid_golden_vs_called_${t}_pos.bed
    # then run bedtools coverage
    #bedtools coverage -a /scratch/amtarave/variant_calling_simulations_project/window_analysis/chr8_genome_50kb_windows.bed -b /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/EUR/males/autos/${i}_chr8_autos_scc_diploid_golden_vs_called_${t}_pos.bed -counts > ${i}_chr8_autos_scc_diploid_golden_vs_called_${t}_counts_50kb_windows.txt

  done
  # merge results across metrics per sample

  # chrY - haploid
  paste ${i}_chrY_scc_haploid_golden_vs_called_TP_counts_50kb_windows.txt ${i}_chrY_scc_haploid_golden_vs_called_FP_counts_50kb_windows.txt ${i}_chrY_scc_haploid_golden_vs_called_FN_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12 }' > ${i}_chrY_scc_haploid_golden_vs_called_all_metrics_counts_50kb_windows.txt

  # chrY - diploid
  paste ${i}_chrY_scc_diploid_golden_vs_called_TP_counts_50kb_windows.txt ${i}_chrY_scc_diploid_golden_vs_called_FP_counts_50kb_windows.txt ${i}_chrY_scc_diploid_golden_vs_called_FN_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12 }' > ${i}_chrY_scc_diploid_golden_vs_called_all_metrics_counts_50kb_windows.txt

  # chrX nonPARs - haploid
  paste ${i}_chrX_nonPARs_scc_haploid_golden_vs_called_TP_counts_50kb_windows.txt ${i}_chrX_nonPARs_scc_haploid_golden_vs_called_FP_counts_50kb_windows.txt ${i}_chrX_nonPARs_scc_haploid_golden_vs_called_FN_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12 }' > ${i}_chrX_nonPARs_scc_haploid_golden_vs_called_all_metrics_counts_50kb_windows.txt

  # chrX nonPARs - diploid
  paste ${i}_chrX_nonPARs_scc_diploid_golden_vs_called_TP_counts_50kb_windows.txt ${i}_chrX_nonPARs_scc_diploid_golden_vs_called_FP_counts_50kb_windows.txt ${i}_chrX_nonPARs_scc_diploid_golden_vs_called_FN_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12 }' > ${i}_chrX_nonPARs_scc_diploid_golden_vs_called_all_metrics_counts_50kb_windows.txt

  # chrX PARs
  paste ${i}_chrX_PARs_scc_diploid_golden_vs_called_TP_counts_50kb_windows.txt ${i}_chrX_PARs_scc_diploid_golden_vs_called_FP_counts_50kb_windows.txt ${i}_chrX_PARs_scc_diploid_golden_vs_called_FN_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12 }' > ${i}_chrX_PARs_scc_diploid_golden_vs_called_all_metrics_counts_50kb_windows.txt

done

# Then after, merege the results across samples
for t in FP FN TP
do
  paste *_chrY_scc_haploid_golden_vs_called_${t}_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12,$16,$20,$24,$28,$32,$36,$40 }' > all_males_chrY_scc_haploid_golden_vs_called_${t}_counts_50kb_windows.txt

  paste *_chrY_scc_diploid_golden_vs_called_${t}_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12,$16,$20,$24,$28,$32,$36,$40 }' > all_males_chrY_scc_diploid_golden_vs_called_${t}_counts_50kb_windows.txt

  paste *_chrX_nonPARs_scc_haploid_golden_vs_called_${t}_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12,$16,$20,$24,$28,$32,$36,$40 }' > all_males_chrX_nonPARs_scc_haploid_golden_vs_called_${t}_counts_50kb_windows.txt

  paste *_chrX_nonPARs_scc_diploid_golden_vs_called_${t}_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12,$16,$20,$24,$28,$32,$36,$40 }' > all_males_chrX_nonPARs_scc_diploid_golden_vs_called_${t}_counts_50kb_windows.txt

  paste *_chrX_PARs_scc_diploid_golden_vs_called_${t}_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12,$16,$20,$24,$28,$32,$36,$40 }' > all_males_chrX_PARs_scc_diploid_golden_vs_called_${t}_counts_50kb_windows.txt

  paste *_chr8_autos_scc_diploid_golden_vs_called_${t}_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12,$16,$20,$24,$28,$32,$36,$40 }' > all_males_chr8_autos_scc_diploid_golden_vs_called_${t}_counts_50kb_windows.txt
done

#################
# FEMALES - SCC #
#################
cd /scratch/amtarave/variant_calling_simulations_project/window_analysis/scc/females/

for i in NA06985 NA06989 NA07000 NA07037 NA07056 NA10847 NA11830 NA11832 NA11840 NA11892
do
  echo $i
  for t in FP FN TP
  do
    echo $t
    # X nonPARs
    #awk '{ print $1"\t"$2-1"\t"$2 }' /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/EUR/females/chrX_nonPARs/${i}_chrX_nonPARs_golden_vs_called_${t}_pos.txt > /data/CEM/wilsonlab/projects/variant_calling_simulations_project/window_analysis/scc/females/${i}_chrX_nonPARs_scc_diploid_golden_vs_called_${t}_pos.bed
    # get counts across chr
    #bedtools coverage -a /scratch/amtarave/variant_calling_simulations_project/window_analysis/chrX_genome_50kb_windows.bed -b /data/CEM/wilsonlab/projects/variant_calling_simulations_project/window_analysis/scc/females/${i}_chrX_nonPARs_scc_diploid_golden_vs_called_${t}_pos.bed -counts > ${i}_chrX_nonPARs_scc_diploid_golden_vs_called_${t}_counts_50kb_windows.txt

    # X PARs
    #awk '{ print $1"\t"$2-1"\t"$2 }' /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/EUR/females/chrX_PARs/${i}_chrX_PARs_golden_vs_called_${t}_pos.txt > /data/CEM/wilsonlab/projects/variant_calling_simulations_project/window_analysis/scc/females/${i}_chrX_PARs_scc_diploid_golden_vs_called_${t}_pos.bed
    # get counts across chr
    #bedtools coverage -a /scratch/amtarave/variant_calling_simulations_project/window_analysis/chrX_genome_50kb_windows.bed -b /data/CEM/wilsonlab/projects/variant_calling_simulations_project/window_analysis/scc/females/${i}_chrX_PARs_scc_diploid_golden_vs_called_${t}_pos.bed -counts > ${i}_chrX_PARs_scc_diploid_golden_vs_called_${t}_counts_50kb_windows.txt

    # chr8
    # make bed files
    #awk '{ print $1"\t"$2-1"\t"$2 }' /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/EUR/females/autos/${i}_chr8_autos_golden_vs_called_${t}_pos.txt > /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/EUR/females/autos/${i}_chr8_autos_scc_diploid_golden_vs_called_${t}_pos.bed
    # then run bedtools coverage
    #bedtools coverage -a /scratch/amtarave/variant_calling_simulations_project/window_analysis/chr8_genome_50kb_windows.bed -b /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/EUR/females/autos/${i}_chr8_autos_scc_diploid_golden_vs_called_${t}_pos.bed -counts > ${i}_chr8_autos_scc_diploid_golden_vs_called_${t}_counts_50kb_windows.txt

  done
  # chrX nonPARs - diploid
  paste ${i}_chrX_nonPARs_scc_diploid_golden_vs_called_TP_counts_50kb_windows.txt ${i}_chrX_nonPARs_scc_diploid_golden_vs_called_FP_counts_50kb_windows.txt ${i}_chrX_nonPARs_scc_diploid_golden_vs_called_FN_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12 }' > ${i}_chrX_nonPARs_scc_diploid_golden_vs_called_all_metrics_counts_50kb_windows.txt

  # chrX PARs
  paste ${i}_chrX_PARs_scc_diploid_golden_vs_called_TP_counts_50kb_windows.txt ${i}_chrX_PARs_scc_diploid_golden_vs_called_FP_counts_50kb_windows.txt ${i}_chrX_PARs_scc_diploid_golden_vs_called_FN_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12 }' > ${i}_chrX_PARs_scc_diploid_golden_vs_called_all_metrics_counts_50kb_windows.txt

done

for t in FP FN TP
do

  paste *_chrX_nonPARs_scc_diploid_golden_vs_called_${t}_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12,$16,$20,$24,$28,$32,$36,$40 }' > all_females_chrX_nonPARs_scc_diploid_golden_vs_called_${t}_counts_50kb_windows.txt

  paste *_chrX_PARs_scc_diploid_golden_vs_called_${t}_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12,$16,$20,$24,$28,$32,$36,$40 }' > all_females_chrX_PARs_scc_diploid_golden_vs_called_${t}_counts_50kb_windows.txt

  paste *_chr8_autos_scc_diploid_golden_vs_called_${t}_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12,$16,$20,$24,$28,$32,$36,$40 }' > all_females_chr8_autos_scc_diploid_golden_vs_called_${t}_counts_50kb_windows.txt
done

###################
# MALES - DEFAULT #
###################
cd /scratch/amtarave/variant_calling_simulations_project/window_analysis/default/males/

for i in NA06984 NA06986 NA06994 NA07048 NA07051 NA07347 NA07357 NA10851 NA11829 NA11831
do
  echo $i
  for t in FP FN TP
  do
    echo $t # TO DO
    # Y
    # make a bed file out of one of the FN results file
    awk '{ print $1"\t"$2-1"\t"$2 }' /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/EUR/males/default/chrY/${i}_chrY_golden_vs_called_${t}_pos.txt > /data/CEM/wilsonlab/projects/variant_calling_simulations_project/window_analysis/default/males/${i}_chrY_default_haploid_golden_vs_called_${t}_pos.bed
    # get counts across Y
    bedtools coverage -a /scratch/amtarave/variant_calling_simulations_project/window_analysis/chrY_genome_50kb_windows.bed -b /data/CEM/wilsonlab/projects/variant_calling_simulations_project/window_analysis/default/males/${i}_chrY_default_haploid_golden_vs_called_${t}_pos.bed -counts > ${i}_chrY_default_haploid_golden_vs_called_${t}_counts_50kb_windows.txt


    awk '{ print $1"\t"$2-1"\t"$2 }' /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/EUR/males/default/chrY/${i}_chrY_golden_vs_called_${t}_pos.txt > /data/CEM/wilsonlab/projects/variant_calling_simulations_project/window_analysis/default/males/${i}_chrY_default_diploid_golden_vs_called_${t}_pos.bed
    # get counts across Y
    bedtools coverage -a /scratch/amtarave/variant_calling_simulations_project/window_analysis/chrY_genome_50kb_windows.bed -b /data/CEM/wilsonlab/projects/variant_calling_simulations_project/window_analysis/default/males/${i}_chrY_default_diploid_golden_vs_called_${t}_pos.bed -counts > ${i}_chrY_default_diploid_golden_vs_called_${t}_counts_50kb_windows.txt


    # X nonPARs
    awk '{ print $1"\t"$2-1"\t"$2 }' /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/compare_VCFs/EUR/males/default/chrX_nonPARs/${i}_chrX_nonPARs_golden_vs_called_${t}_pos.txt > /data/CEM/wilsonlab/projects/variant_calling_simulations_project/window_analysis/default/males/${i}_chrX_nonPARs_default_haploid_golden_vs_called_${t}_pos.bed
    # get FN counts across Y
    bedtools coverage -a /scratch/amtarave/variant_calling_simulations_project/window_analysis/chrX_genome_50kb_windows.bed -b /data/CEM/wilsonlab/projects/variant_calling_simulations_project/window_analysis/default/males/${i}_chrX_nonPARs_default_haploid_golden_vs_called_${t}_pos.bed -counts > ${i}_chrX_nonPARs_default_haploid_golden_vs_called_${t}_counts_50kb_windows.txt

    awk '{ print $1"\t"$2-1"\t"$2 }' /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/EUR/males/default/chrX_nonPARs/${i}_chrX_nonPARs_golden_vs_called_${t}_pos.txt > /data/CEM/wilsonlab/projects/variant_calling_simulations_project/window_analysis/default/males/${i}_chrX_nonPARs_default_diploid_golden_vs_called_${t}_pos.bed
    # get counts across Y
    bedtools coverage -a /scratch/amtarave/variant_calling_simulations_project/window_analysis/chrX_genome_50kb_windows.bed -b /data/CEM/wilsonlab/projects/variant_calling_simulations_project/window_analysis/default/males/${i}_chrX_nonPARs_default_diploid_golden_vs_called_${t}_pos.bed -counts > ${i}_chrX_nonPARs_default_diploid_golden_vs_called_${t}_counts_50kb_windows.txt

    # X PARs
    awk '{ print $1"\t"$2-1"\t"$2 }' /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/EUR/males/default/chrX_PARs/${i}_chrX_PARs_golden_vs_called_${t}_pos.txt > /data/CEM/wilsonlab/projects/variant_calling_simulations_project/window_analysis/default/males/${i}_chrX_PARs_default_diploid_golden_vs_called_${t}_pos.bed
    # get counts across Y
    bedtools coverage -a /scratch/amtarave/variant_calling_simulations_project/window_analysis/chrX_genome_50kb_windows.bed -b /data/CEM/wilsonlab/projects/variant_calling_simulations_project/window_analysis/default/males/${i}_chrX_PARs_default_diploid_golden_vs_called_${t}_pos.bed -counts > ${i}_chrX_PARs_default_diploid_golden_vs_called_${t}_counts_50kb_windows.txt

    # chr8
    # make bed files
    awk '{ print $1"\t"$2-1"\t"$2 }' /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/EUR/males/default/autos/${i}_chr8_autos_golden_vs_called_${t}_pos.txt > /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/EUR/males/default/autos/${i}_chr8_autos_default_diploid_golden_vs_called_${t}_pos.bed
    # then run bedtools coverage
    bedtools coverage -a /scratch/amtarave/variant_calling_simulations_project/window_analysis/chr8_genome_50kb_windows.bed -b /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/EUR/males/default/autos/${i}_chr8_autos_default_diploid_golden_vs_called_${t}_pos.bed -counts > ${i}_chr8_autos_default_diploid_golden_vs_called_${t}_counts_50kb_windows.txt

  done
  # chrY - haploid
  paste ${i}_chrY_default_haploid_golden_vs_called_TP_counts_50kb_windows.txt ${i}_chrY_default_haploid_golden_vs_called_FP_counts_50kb_windows.txt ${i}_chrY_default_haploid_golden_vs_called_FN_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12 }' > ${i}_chrY_default_haploid_golden_vs_called_all_metrics_counts_50kb_windows.txt

  # chrY - diploid
  paste ${i}_chrY_default_diploid_golden_vs_called_TP_counts_50kb_windows.txt ${i}_chrY_default_diploid_golden_vs_called_FP_counts_50kb_windows.txt ${i}_chrY_default_diploid_golden_vs_called_FN_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12 }' > ${i}_chrY_default_diploid_golden_vs_called_all_metrics_counts_50kb_windows.txt

  # chrX nonPARs - haploid
  paste ${i}_chrX_nonPARs_default_haploid_golden_vs_called_TP_counts_50kb_windows.txt ${i}_chrX_nonPARs_default_haploid_golden_vs_called_FP_counts_50kb_windows.txt ${i}_chrX_nonPARs_default_haploid_golden_vs_called_FN_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12 }' > ${i}_chrX_nonPARs_default_haploid_golden_vs_called_all_metrics_counts_50kb_windows.txt

  # chrX nonPARs - diploid
  paste ${i}_chrX_nonPARs_default_diploid_golden_vs_called_TP_counts_50kb_windows.txt ${i}_chrX_nonPARs_default_diploid_golden_vs_called_FP_counts_50kb_windows.txt ${i}_chrX_nonPARs_default_diploid_golden_vs_called_FN_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12 }' > ${i}_chrX_nonPARs_default_diploid_golden_vs_called_all_metrics_counts_50kb_windows.txt

  # chrX PARs
  paste ${i}_chrX_PARs_default_diploid_golden_vs_called_TP_counts_50kb_windows.txt ${i}_chrX_PARs_default_diploid_golden_vs_called_FP_counts_50kb_windows.txt ${i}_chrX_PARs_default_diploid_golden_vs_called_FN_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12 }' > ${i}_chrX_PARs_default_diploid_golden_vs_called_all_metrics_counts_50kb_windows.txt

done

for t in FP FN TP
do
  paste *_chrY_default_haploid_golden_vs_called_${t}_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12,$16,$20,$24,$28,$32,$36,$40 }' > all_males_chrY_default_haploid_golden_vs_called_${t}_counts_50kb_windows.txt

  paste *_chrY_default_diploid_golden_vs_called_${t}_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12,$16,$20,$24,$28,$32,$36,$40 }' > all_males_chrY_default_diploid_golden_vs_called_${t}_counts_50kb_windows.txt

  paste *_chrX_nonPARs_default_haploid_golden_vs_called_${t}_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12,$16,$20,$24,$28,$32,$36,$40 }' > all_males_chrX_nonPARs_default_haploid_golden_vs_called_${t}_counts_50kb_windows.txt

  paste *_chrX_nonPARs_default_diploid_golden_vs_called_${t}_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12,$16,$20,$24,$28,$32,$36,$40 }' > all_males_chrX_nonPARs_default_diploid_golden_vs_called_${t}_counts_50kb_windows.txt

  paste *_chrX_PARs_default_diploid_golden_vs_called_${t}_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12,$16,$20,$24,$28,$32,$36,$40 }' > all_males_chrX_PARs_default_diploid_golden_vs_called_${t}_counts_50kb_windows.txt

  paste *_chr8_autos_default_diploid_golden_vs_called_${t}_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12,$16,$20,$24,$28,$32,$36,$40 }' > all_males_chr8_autos_default_diploid_golden_vs_called_${t}_counts_50kb_windows.txt
done

#####################
# FEMALES - DEFAULT #
#####################
cd /scratch/amtarave/variant_calling_simulations_project/window_analysis/default/females/

for i in NA06985 NA06989 NA07000 NA07037 NA07056 NA10847 NA11830 NA11832 NA11840 NA11892
do
  echo $i
  for t in FP FN TP
  do
    echo $t
    # X nonPARs
    awk '{ print $1"\t"$2-1"\t"$2 }' /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/EUR/females/default/chrX_nonPARs/${i}_chrX_nonPARs_golden_vs_called_${t}_pos.txt > /data/CEM/wilsonlab/projects/variant_calling_simulations_project/window_analysis/default/females/${i}_chrX_nonPARs_default_diploid_golden_vs_called_${t}_pos.bed
    # get counts across chr
    bedtools coverage -a /scratch/amtarave/variant_calling_simulations_project/window_analysis/chrX_genome_50kb_windows.bed -b /data/CEM/wilsonlab/projects/variant_calling_simulations_project/window_analysis/default/females/${i}_chrX_nonPARs_default_diploid_golden_vs_called_${t}_pos.bed -counts > ${i}_chrX_nonPARs_default_diploid_golden_vs_called_${t}_counts_50kb_windows.txt

    # X PARs
    awk '{ print $1"\t"$2-1"\t"$2 }' /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/EUR/females/default/chrX_PARs/${i}_chrX_PARs_golden_vs_called_${t}_pos.txt > /data/CEM/wilsonlab/projects/variant_calling_simulations_project/window_analysis/default/females/${i}_chrX_PARs_default_diploid_golden_vs_called_${t}_pos.bed
    # get counts across chr
    bedtools coverage -a /scratch/amtarave/variant_calling_simulations_project/window_analysis/chrX_genome_50kb_windows.bed -b /data/CEM/wilsonlab/projects/variant_calling_simulations_project/window_analysis/default/females/${i}_chrX_PARs_default_diploid_golden_vs_called_${t}_pos.bed -counts > ${i}_chrX_PARs_default_diploid_golden_vs_called_${t}_counts_50kb_windows.txt

    # chr8
    # make bed files
    awk '{ print $1"\t"$2-1"\t"$2 }' /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/EUR/females/default/autos/${i}_chr8_autos_golden_vs_called_${t}_pos.txt > /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/EUR/females/default/autos/${i}_chr8_autos_default_diploid_golden_vs_called_${t}_pos.bed
    # then run bedtools coverage
    bedtools coverage -a /scratch/amtarave/variant_calling_simulations_project/window_analysis/chr8_genome_50kb_windows.bed -b /data/CEM/wilsonlab/projects/variant_calling_simulations_project/compare_VCFs/EUR/females/default/autos/${i}_chr8_autos_default_diploid_golden_vs_called_${t}_pos.bed -counts > ${i}_chr8_autos_default_diploid_golden_vs_called_${t}_counts_50kb_windows.txt

  done
  # chrX nonPARs - diploid
  paste ${i}_chrX_nonPARs_default_diploid_golden_vs_called_TP_counts_50kb_windows.txt ${i}_chrX_nonPARs_default_diploid_golden_vs_called_FP_counts_50kb_windows.txt ${i}_chrX_nonPARs_default_diploid_golden_vs_called_FN_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12 }' > ${i}_chrX_nonPARs_default_diploid_golden_vs_called_all_metrics_counts_50kb_windows.txt

  # chrX PARs
  paste ${i}_chrX_PARs_default_diploid_golden_vs_called_TP_counts_50kb_windows.txt ${i}_chrX_PARs_default_diploid_golden_vs_called_FP_counts_50kb_windows.txt ${i}_chrX_PARs_default_diploid_golden_vs_called_FN_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12 }' > ${i}_chrX_PARs_default_diploid_golden_vs_called_all_metrics_counts_50kb_windows.txt

done

for t in FP FN TP
do

  paste *_chrX_nonPARs_default_diploid_golden_vs_called_${t}_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12,$16,$20,$24,$28,$32,$36,$40 }' > all_females_chrX_nonPARs_default_diploid_golden_vs_called_${t}_counts_50kb_windows.txt

  paste *_chrX_PARs_default_diploid_golden_vs_called_${t}_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12,$16,$20,$24,$28,$32,$36,$40 }' > all_females_chrX_PARs_default_diploid_golden_vs_called_${t}_counts_50kb_windows.txt

  paste *_chr8_autos_default_diploid_golden_vs_called_${t}_counts_50kb_windows.txt | awk '{ print $1,$2,$3,$4,$8,$12,$16,$20,$24,$28,$32,$36,$40 }' > all_females_chr8_autos_default_diploid_golden_vs_called_${t}_counts_50kb_windows.txt
done

# Now plot results in R. Average across samples.
