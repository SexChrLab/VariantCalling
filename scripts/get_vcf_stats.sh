#!/bin/bash
#SBATCH --job-name=get_vcf_stats # Job name
#SBATCH -o slurm.get_vcf_stats.out
#SBATCH -e slurm.get_vcf_stats.err
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mail-user=amtarave@asu.edu # send-to address

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/scripts

conda actiavate variant_calling_simulations_project
# example:
# python3 /home/amtarave/packages/vcfhelper/extract_stats_from_vcf.py QD QUAL SOR FS MQ MQRankSum ReadPosRankSum AN DP \
# --vcf /data/CEM/wilsonlab/projects/variant_calling_simulations_project/hard_filtered_vcfs/EUR/males/all/chrX_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_PARs.vcf.gz \
# --outfile test_males_chrX_PARs_vcf_stats.txt


################################################################################
# MALES - SCC #
################################################################################
#-----------------------------#
# DIPLOID CALLING - Autosomes #
#-----------------------------#
ancestry="EUR"
sex="males"

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/hard_filtered_vcfs/${ancestry}/${sex}/all/

for i in {1..22}
do
	echo "chr${i}"
  python3 /home/amtarave/packages/vcfhelper/extract_stats_from_vcf.py QD QUAL SOR FS MQ MQRankSum ReadPosRankSum AN DP --vcf chr${i}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz --outfile chr${i}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_stats.txt
done

#------------------------------------#
# DIPLOID CALLING - Sex chromosomes  #
#------------------------------------#

# PARs
python3 /home/amtarave/packages/vcfhelper/extract_stats_from_vcf.py QD QUAL SOR FS MQ MQRankSum ReadPosRankSum AN DP --vcf chrX_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_PARs.vcf.gz --outfile chrX_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_PARs_stats.txt


#-----------------#
# HAPLOID CALLING #
#-----------------#
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/hard_filtered_vcfs/${ancestry}/${sex}/all/

# ChrX non-PARs
python3 /home/amtarave/packages/vcfhelper/extract_stats_from_vcf.py QD QUAL SOR FS MQ MQRankSum ReadPosRankSum AN DP --vcf chrX_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_gatkHardFilter_filtered_nonPARs.vcf.gz --outfile chrX_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_gatkHardFilter_filtered_nonPARs_stats.txt

# ChrY non-PARs
python3 /home/amtarave/packages/vcfhelper/extract_stats_from_vcf.py QD QUAL SOR FS MQ MQRankSum ReadPosRankSum AN DP --vcf chrY_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_gatkHardFilter_filtered.vcf.gz --outfile chrY_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_gatkHardFilter_filtered_stats.txt


################################################################################
# MALES - DEFAULT #
################################################################################
#-----------------------------#
# DIPLOID CALLING - Autosomes #
#-----------------------------#
ancestry="EUR"
sex="males"

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/hard_filtered_vcfs/${ancestry}/${sex}/default/all/

for i in {1..22}
do
	echo "chr${i}"
  python3 /home/amtarave/packages/vcfhelper/extract_stats_from_vcf.py QD QUAL SOR FS MQ MQRankSum ReadPosRankSum AN DP --vcf chr${i}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz --outfile chr${i}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_stats.txt
done

#------------------------------------#
# DIPLOID CALLING - Sex chromosomes  #
#------------------------------------#

# PARs
python3 /home/amtarave/packages/vcfhelper/extract_stats_from_vcf.py QD QUAL SOR FS MQ MQRankSum ReadPosRankSum AN DP --vcf chrX_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_PARs.vcf.gz --outfile chrX_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_PARs_stats.txt


#-----------------#
# HAPLOID CALLING #
#-----------------#
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/hard_filtered_vcfs/${ancestry}/${sex}/default/all/

# ChrX non-PARs
python3 /home/amtarave/packages/vcfhelper/extract_stats_from_vcf.py QD QUAL SOR FS MQ MQRankSum ReadPosRankSum AN DP --vcf chrX_GRCh38_default_gatk_haploid_called_SNPs_gatkHardFilter_filtered_nonPARs.vcf.gz --outfile chrX_GRCh38_default_gatk_haploid_called_SNPs_gatkHardFilter_filtered_nonPARs_stats.txt

# ChrY non-PARs
python3 /home/amtarave/packages/vcfhelper/extract_stats_from_vcf.py QD QUAL SOR FS MQ MQRankSum ReadPosRankSum AN DP --vcf chrY_GRCh38_default_gatk_haploid_called_SNPs_gatkHardFilter_filtered.vcf.gz --outfile chrY_GRCh38_default_gatk_haploid_called_SNPs_gatkHardFilter_filtered_stats.txt



################################################################################
# FEMALES - SCC #
################################################################################
#-----------------------------#
# DIPLOID CALLING - Autosomes #
#-----------------------------#
ancestry="EUR"
sex="females"

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/hard_filtered_vcfs/${ancestry}/${sex}/all/

for i in {1..22}
do
	echo "chr${i}"
  python3 /home/amtarave/packages/vcfhelper/extract_stats_from_vcf.py QD QUAL SOR FS MQ MQRankSum ReadPosRankSum AN DP --vcf chr${i}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz --outfile chr${i}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_stats.txt
done

#------------------------------------#
# DIPLOID CALLING - Sex chromosomes  #
#------------------------------------#

# PARs
python3 /home/amtarave/packages/vcfhelper/extract_stats_from_vcf.py QD QUAL SOR FS MQ MQRankSum ReadPosRankSum AN DP --vcf chrX_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_PARs.vcf.gz --outfile chrX_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_PARs_stats.txt

# ChrX non-PARs
python3 /home/amtarave/packages/vcfhelper/extract_stats_from_vcf.py QD QUAL SOR FS MQ MQRankSum ReadPosRankSum AN DP --vcf chrX_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs.vcf.gz --outfile chrX_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs_stats.txt


################################################################################
# FEMALES - DEFAULT #
################################################################################
#-----------------------------#
# DIPLOID CALLING - Autosomes #
#-----------------------------#
ancestry="EUR"
sex="females"

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/hard_filtered_vcfs/${ancestry}/${sex}/default/all/

for i in {1..22}
do
	echo "chr${i}"
  python3 /home/amtarave/packages/vcfhelper/extract_stats_from_vcf.py QD QUAL SOR FS MQ MQRankSum ReadPosRankSum AN DP --vcf chr${i}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz --outfile chr${i}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_stats.txt
done

#------------------------------------#
# DIPLOID CALLING - Sex chromosomes  #
#------------------------------------#

# PARs
python3 /home/amtarave/packages/vcfhelper/extract_stats_from_vcf.py QD QUAL SOR FS MQ MQRankSum ReadPosRankSum AN DP --vcf chrX_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_PARs.vcf.gz --outfile chrX_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_PARs_stats.txt

# ChrX non-PARs
python3 /home/amtarave/packages/vcfhelper/extract_stats_from_vcf.py QD QUAL SOR FS MQ MQRankSum ReadPosRankSum AN DP --vcf chrX_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs.vcf.gz --outfile chrX_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs_stats.txt




################################################################################
################################################################################
# Repeat but just with the raw snp (unfiltered VCFs)
################################################################################
################################################################################

################################################################################
# MALES - SCC #
################################################################################
#-----------------------------#
# DIPLOID CALLING - Autosomes #
#-----------------------------#
ancestry="EUR"
sex="males"

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/hard_filtered_vcfs/${ancestry}/${sex}/all/

for i in {1..22}
do
	echo "chr${i}"
  python3 /home/amtarave/packages/vcfhelper/extract_stats_from_vcf.py QD QUAL SOR FS MQ MQRankSum ReadPosRankSum AN DP --vcf chr${i}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_gatkHardFilter.vcf.gz --outfile chr${i}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_gatkHardFilter_stats.txt
done

#------------------------------------#
# DIPLOID CALLING - Sex chromosomes  #
#------------------------------------#

# PARs
# extract PARs from unfiltered vcf first
gatk SelectVariants -R /scratch/amtarave/variant_calling_simulations_project/reference_genomes/GRCh38_full_analysis_set_plus_decoy_hla_YPARmasked_XX.fa -V chrX_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_gatkHardFilter.vcf.gz -L /home/amtarave/projects/variant_calling/lists/XPARs.interval.list -O chrX_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_gatkHardFilter_PARs.vcf.gz

python3 /home/amtarave/packages/vcfhelper/extract_stats_from_vcf.py QD QUAL SOR FS MQ MQRankSum ReadPosRankSum AN DP --vcf chrX_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_gatkHardFilter_PARs.vcf.gz --outfile chrX_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_gatkHardFilter_PARs_stats.txt


#-----------------#
# HAPLOID CALLING #
#-----------------#
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/hard_filtered_vcfs/${ancestry}/${sex}/all/

# ChrX non-PARs
# extract non-PARs from unfiltered vcf first
gatk SelectVariants -R /scratch/amtarave/variant_calling_simulations_project/reference_genomes/GRCh38_full_analysis_set_plus_decoy_hla_YPARmasked_XX.fa -V chrX_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs_gatkHardFilter.vcf.gz -XL /home/amtarave/projects/variant_calling/lists/XPARs.interval.list -O chrX_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs_gatkHardFilter_nonPARs.vcf.gz

python3 /home/amtarave/packages/vcfhelper/extract_stats_from_vcf.py QD QUAL SOR FS MQ MQRankSum ReadPosRankSum AN DP --vcf chrX_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs_gatkHardFilter_nonPARs.vcf.gz --outfile chrX_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs_gatkHardFilter_nonPARs_stats.txt

# ChrY non-PARs
python3 /home/amtarave/packages/vcfhelper/extract_stats_from_vcf.py QD QUAL SOR FS MQ MQRankSum ReadPosRankSum AN DP --vcf chrY_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs_gatkHardFilter.vcf.gz --outfile chrY_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs_gatkHardFilter_stats.txt


################################################################################
# MALES - DEFAULT #
################################################################################
#-----------------------------#
# DIPLOID CALLING - Autosomes #
#-----------------------------#
ancestry="EUR"
sex="males"

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/hard_filtered_vcfs/${ancestry}/${sex}/default/all/

for i in {1..22}
do
	echo "chr${i}"
  python3 /home/amtarave/packages/vcfhelper/extract_stats_from_vcf.py QD QUAL SOR FS MQ MQRankSum ReadPosRankSum AN DP --vcf chr${i}_GRCh38_default_gatk_diploid_called_raw_SNPs_gatkHardFilter.vcf.gz --outfile chr${i}_GRCh38_default_gatk_diploid_called_raw_SNPs_gatkHardFilter_stats.txt
done

#------------------------------------#
# DIPLOID CALLING - Sex chromosomes  #
#------------------------------------#

# PARs
# extract PARs from unfiltered vcf first
gatk SelectVariants -R /scratch/amtarave/variant_calling_simulations_project/reference_genomes/GRCh38_full_analysis_set_plus_decoy_hla_DEFAULT_reformat.fa -V chrX_GRCh38_default_gatk_diploid_called_raw_SNPs_gatkHardFilter.vcf.gz -L /home/amtarave/projects/variant_calling/lists/XPARs.interval.list -O chrX_GRCh38_default_gatk_diploid_called_raw_SNPs_gatkHardFilter_PARs.vcf.gz

python3 /home/amtarave/packages/vcfhelper/extract_stats_from_vcf.py QD QUAL SOR FS MQ MQRankSum ReadPosRankSum AN DP --vcf chrX_GRCh38_default_gatk_diploid_called_raw_SNPs_gatkHardFilter_PARs.vcf.gz --outfile chrX_GRCh38_default_gatk_diploid_called_raw_SNPs_gatkHardFilter_PARs_stats.txt


#-----------------#
# HAPLOID CALLING #
#-----------------#
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/hard_filtered_vcfs/${ancestry}/${sex}/default/all/

# ChrX non-PARs
# extract non-PARs from unfiltered vcf first
gatk SelectVariants -R /scratch/amtarave/variant_calling_simulations_project/reference_genomes/GRCh38_full_analysis_set_plus_decoy_hla_DEFAULT_reformat.fa -V chrX_GRCh38_default_gatk_haploid_called_raw_SNPs_gatkHardFilter.vcf.gz -XL /home/amtarave/projects/variant_calling/lists/XPARs.interval.list -O chrX_GRCh38_default_gatk_haploid_called_raw_SNPs_gatkHardFilter_nonPARs.vcf.gz

python3 /home/amtarave/packages/vcfhelper/extract_stats_from_vcf.py QD QUAL SOR FS MQ MQRankSum ReadPosRankSum AN DP --vcf chrX_GRCh38_default_gatk_haploid_called_raw_SNPs_gatkHardFilter_nonPARs.vcf.gz --outfile chrX_GRCh38_default_gatk_haploid_called_raw_SNPs_gatkHardFilter_nonPARs_stats.txt

# ChrY non-PARs
python3 /home/amtarave/packages/vcfhelper/extract_stats_from_vcf.py QD QUAL SOR FS MQ MQRankSum ReadPosRankSum AN DP --vcf chrY_GRCh38_default_gatk_haploid_called_raw_SNPs_gatkHardFilter.vcf.gz --outfile chrY_GRCh38_default_gatk_haploid_called_raw_SNPs_gatkHardFilter_stats.txt



################################################################################
# FEMALES - SCC #
################################################################################
#-----------------------------#
# DIPLOID CALLING - Autosomes #
#-----------------------------#
ancestry="EUR"
sex="females"

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/hard_filtered_vcfs/${ancestry}/${sex}/all/

for i in {1..22}
do
	echo "chr${i}"
  python3 /home/amtarave/packages/vcfhelper/extract_stats_from_vcf.py QD QUAL SOR FS MQ MQRankSum ReadPosRankSum AN DP --vcf chr${i}_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs_gatkHardFilter.vcf.gz --outfile chr${i}_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs_gatkHardFilter_stats.txt
done

#------------------------------------#
# DIPLOID CALLING - Sex chromosomes  #
#------------------------------------#

# PARs
# extract PARs from unfiltered vcf first
gatk SelectVariants -R /scratch/amtarave/variant_calling_simulations_project/reference_genomes/GRCh38_full_analysis_set_plus_decoy_hla_Ymasked_XX.fa -V chrX_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs_gatkHardFilter.vcf.gz -L /home/amtarave/projects/variant_calling/lists/XPARs.interval.list -O chrX_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs_gatkHardFilter_PARs.vcf.gz

python3 /home/amtarave/packages/vcfhelper/extract_stats_from_vcf.py QD QUAL SOR FS MQ MQRankSum ReadPosRankSum AN DP --vcf chrX_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs_gatkHardFilter_PARs.vcf.gz --outfile chrX_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs_gatkHardFilter_PARs_stats.txt

# ChrX non-PARs
# extract PARs from unfiltered vcf first
gatk SelectVariants -R /scratch/amtarave/variant_calling_simulations_project/reference_genomes/GRCh38_full_analysis_set_plus_decoy_hla_Ymasked_XX.fa -V chrX_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs_gatkHardFilter.vcf.gz -XL /home/amtarave/projects/variant_calling/lists/XPARs.interval.list -O chrX_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs_gatkHardFilter_nonPARs.vcf.gz

python3 /home/amtarave/packages/vcfhelper/extract_stats_from_vcf.py QD QUAL SOR FS MQ MQRankSum ReadPosRankSum AN DP --vcf chrX_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs_gatkHardFilter_nonPARs.vcf.gz --outfile chrX_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs_gatkHardFilter_nonPARs_stats.txt


################################################################################
# FEMALES - DEFAULT #
################################################################################
#-----------------------------#
# DIPLOID CALLING - Autosomes #
#-----------------------------#
ancestry="EUR"
sex="females"

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/hard_filtered_vcfs/${ancestry}/${sex}/default/all/

for i in {1..22}
do
	echo "chr${i}"
  python3 /home/amtarave/packages/vcfhelper/extract_stats_from_vcf.py QD QUAL SOR FS MQ MQRankSum ReadPosRankSum AN DP --vcf chr${i}_GRCh38_default_gatk_diploid_called_raw_SNPs_gatkHardFilter.vcf.gz --outfile chr${i}_GRCh38_default_gatk_diploid_called_raw_SNPs_gatkHardFilter_stats.txt
done

#------------------------------------#
# DIPLOID CALLING - Sex chromosomes  #
#------------------------------------#

# PARs
# extract PARs from unfiltered vcf first
gatk SelectVariants -R /scratch/amtarave/variant_calling_simulations_project/reference_genomes/GRCh38_full_analysis_set_plus_decoy_hla_DEFAULT_reformat.fa -V chrX_GRCh38_default_gatk_diploid_called_raw_SNPs_gatkHardFilter.vcf.gz -L /home/amtarave/projects/variant_calling/lists/XPARs.interval.list -O chrX_GRCh38_default_gatk_diploid_called_raw_SNPs_gatkHardFilter_PARs.vcf.gz

python3 /home/amtarave/packages/vcfhelper/extract_stats_from_vcf.py QD QUAL SOR FS MQ MQRankSum ReadPosRankSum AN DP --vcf chrX_GRCh38_default_gatk_diploid_called_raw_SNPs_gatkHardFilter_PARs.vcf.gz --outfile chrX_GRCh38_default_gatk_diploid_called_raw_SNPs_gatkHardFilter_PARs_stats.txt

# ChrX non-PARs
# extract non-PARs from unfiltered vcf first
gatk SelectVariants -R /scratch/amtarave/variant_calling_simulations_project/reference_genomes/GRCh38_full_analysis_set_plus_decoy_hla_DEFAULT_reformat.fa -V chrX_GRCh38_default_gatk_diploid_called_raw_SNPs_gatkHardFilter.vcf.gz -XL /home/amtarave/projects/variant_calling/lists/XPARs.interval.list -O chrX_GRCh38_default_gatk_diploid_called_raw_SNPs_gatkHardFilter_nonPARs.vcf.gz

python3 /home/amtarave/packages/vcfhelper/extract_stats_from_vcf.py QD QUAL SOR FS MQ MQRankSum ReadPosRankSum AN DP --vcf chrX_GRCh38_default_gatk_diploid_called_raw_SNPs_gatkHardFilter_nonPARs.vcf.gz --outfile chrX_GRCh38_default_gatk_diploid_called_raw_SNPs_gatkHardFilter_nonPARs_stats.txt
