#!/bin/bash
#SBATCH --job-name=get_high_quality_sites # Job name
#SBATCH -o slurm.get_high_quality_sites.out
#SBATCH -e slurm.get_high_quality_sites.err
#SBATCH -n 2
#SBATCH -t 7-00:00:00
#SBATCH --mail-user=amtarave@asu.edu # send-to address

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/scripts

conda actiavate variant_calling_simulations_project

# Need to get bed file for all autosomes and X for default and SCC

# Autosome filters
# AN: 10 # IT NEEDS TO BE CHANGED TO 20
# DP: 67 - 201

# X filters:
# AN: 5 # IT NEEDS TO BE CHANGED TO 10
# DP: 23 - 69

# Thought...should I have ran this on the filtered set?

################################################################################
# MALES - SCC #
################################################################################
#-----------------------------#
# DIPLOID CALLING - Autosomes #
#-----------------------------#
ancestry="EUR"
sex="males"

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/joint_called_vcfs/${ancestry}/${sex}/

for i in {1..22}
do
	echo "chr${i}"
  # first filter for AN
  bcftools filter -e'INFO/AN<20' chr${i}_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites.vcf.gz -Oz -o chr${i}_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20.vcf.gz
  tabix -p vcf chr${i}_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20.vcf.gz
  # then filter for DP
  bcftools filter -e'INFO/DP<67 || INFO/DP>201' chr${i}_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20.vcf.gz -Oz -o chr${i}_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201.vcf.gz
  tabix -p vcf chr${i}_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201.vcf.gz

  # then make bed file? (or can just intersect the vcf with the neutral bed file)
done


#------------------------------#
# HAPLOID CALLING - X non-PARs #
#------------------------------#
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/joint_called_vcfs/${ancestry}/${sex}/

# first remove PARs
gatk SelectVariants -R /scratch/amtarave/variant_calling_simulations_project/reference_genomes/GRCh38_full_analysis_set_plus_decoy_hla_YPARmasked_XX.fa -V chrX_GRCh38_YPARsMasked_gatk_haploid_called_raw_allSites.vcf.gz -XL /home/amtarave/projects/variant_calling/lists/XPARs.interval.list -O chrX_GRCh38_YPARsMasked_gatk_haploid_called_raw_allSites_nonPARs.vcf.gz

bcftools filter -e'INFO/AN<10' chrX_GRCh38_YPARsMasked_gatk_haploid_called_raw_allSites_nonPARs.vcf.gz -Oz -o chrX_GRCh38_YPARsMasked_gatk_haploid_called_raw_allSites_nonPARs_AN10.vcf.gz
tabix -p vcf chrX_GRCh38_YPARsMasked_gatk_haploid_called_raw_allSites_nonPARs_AN10.vcf.gz
# then filter for DP
bcftools filter -e'INFO/DP<23 || INFO/DP>69' chrX_GRCh38_YPARsMasked_gatk_haploid_called_raw_allSites_nonPARs_AN10.vcf.gz -Oz -o chrX_GRCh38_YPARsMasked_gatk_haploid_called_raw_allSites_nonPARs_AN10_DP23and69.vcf.gz
tabix -p vcf chrX_GRCh38_YPARsMasked_gatk_haploid_called_raw_allSites_nonPARs_AN10_DP23and69.vcf.gz



################################################################################
# MALES - DEFAULT #
################################################################################
#-----------------------------#
# DIPLOID CALLING - Autosomes #
#-----------------------------#
ancestry="EUR"
sex="males"

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/joint_called_vcfs/${ancestry}/${sex}/default/

for i in {1..22}
do
  echo "chr${i}"
  # first filter for AN
  bcftools filter -e'INFO/AN<20' chr${i}_GRCh38_default_gatk_diploid_called_raw_allSites.vcf.gz -Oz -o chr${i}_GRCh38_default_gatk_diploid_called_raw_allSites_AN20.vcf.gz
  tabix -p vcf chr${i}_GRCh38_default_gatk_diploid_called_raw_allSites_AN20.vcf.gz
  # then filter for DP
  bcftools filter -e'INFO/DP<67 || INFO/DP>201' chr${i}_GRCh38_default_gatk_diploid_called_raw_allSites_AN20.vcf.gz -Oz -o chr${i}_GRCh38_default_gatk_diploid_called_raw_allSites_AN20_DP67and201.vcf.gz
  tabix -p vcf chr${i}_GRCh38_default_gatk_diploid_called_raw_allSites_AN20_DP67and201.vcf.gz

  # then make bed file? (or can just intersect the vcf with the neutral bed file)
done


#------------------------------#
# HAPLOID CALLING - X non-PARs #
#------------------------------#
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/joint_called_vcfs/${ancestry}/${sex}/default/

# first remove PARs
gatk SelectVariants -R /scratch/amtarave/variant_calling_simulations_project/reference_genomes/GRCh38_full_analysis_set_plus_decoy_hla_DEFAULT_reformat.fa -V chrX_GRCh38_default_gatk_haploid_called_raw_allSites.vcf.gz -XL /home/amtarave/projects/variant_calling/lists/XPARs.interval.list -O chrX_GRCh38_default_gatk_haploid_called_raw_allSites_nonPARs.vcf.gz

bcftools filter -e'INFO/AN<10' chrX_GRCh38_default_gatk_haploid_called_raw_allSites_nonPARs.vcf.gz -Oz -o chrX_GRCh38_default_gatk_haploid_called_raw_allSites_nonPARs_AN10.vcf.gz
tabix -p vcf chrX_GRCh38_default_gatk_haploid_called_raw_allSites_nonPARs_AN10.vcf.gz
# then filter for DP
bcftools filter -e'INFO/DP<23 || INFO/DP>69' chrX_GRCh38_default_gatk_haploid_called_raw_allSites_nonPARs_AN10.vcf.gz -Oz -o chrX_GRCh38_default_gatk_haploid_called_raw_allSites_nonPARs_AN10_DP23and69.vcf.gz
tabix -p vcf chrX_GRCh38_default_gatk_haploid_called_raw_allSites_nonPARs_AN10_DP23and69.vcf.gz


################################################################################
# FEMALES - SCC #
################################################################################
#-----------------------------#
# DIPLOID CALLING - Autosomes #
#-----------------------------#
ancestry="EUR"
sex="females"

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/joint_called_vcfs/${ancestry}/${sex}/

for i in {1..22}
do
	echo "chr${i}"
  # first filter for AN
  bcftools filter -e'INFO/AN<20' chr${i}_GRCh38_YHardMasked_gatk_diploid_called_raw_allSites.vcf.gz -Oz -o chr${i}_GRCh38_YHardMasked_gatk_diploid_called_raw_allSites_AN20.vcf.gz
  tabix -p vcf chr${i}_GRCh38_YHardMasked_gatk_diploid_called_raw_allSites_AN20.vcf.gz
  # then filter for DP
  bcftools filter -e'INFO/DP<67 || INFO/DP>201' chr${i}_GRCh38_YHardMasked_gatk_diploid_called_raw_allSites_AN20.vcf.gz -Oz -o chr${i}_GRCh38_YHardMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201.vcf.gz
  tabix -p vcf chr${i}_GRCh38_YHardMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201.vcf.gz

  # then make bed file? (or can just intersect the vcf with the neutral bed file)
done


#------------#
# X non-PARs #
#------------#
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/joint_called_vcfs/${ancestry}/${sex}/

# first remove PARs
gatk SelectVariants -R /scratch/amtarave/variant_calling_simulations_project/reference_genomes/GRCh38_full_analysis_set_plus_decoy_hla_Ymasked_XX.fa -V chrX_GRCh38_YHardMasked_gatk_diploid_called_raw_allSites.vcf.gz -XL /home/amtarave/projects/variant_calling/lists/XPARs.interval.list -O chrX_GRCh38_YHardMasked_gatk_diploid_called_raw_allSites_nonPARs.vcf.gz

bcftools filter -e'INFO/AN<20' chrX_GRCh38_YHardMasked_gatk_diploid_called_raw_allSites_nonPARs.vcf.gz -Oz -o chrX_GRCh38_YHardMasked_gatk_diploid_called_raw_allSites_nonPARs_AN20.vcf.gz
tabix -p vcf chrX_GRCh38_YHardMasked_gatk_diploid_called_raw_allSites_nonPARs_AN20.vcf.gz
# then filter for DP
bcftools filter -e'INFO/DP<67 || INFO/DP>201' chrX_GRCh38_YHardMasked_gatk_diploid_called_raw_allSites_nonPARs_AN20.vcf.gz -Oz -o chrX_GRCh38_YHardMasked_gatk_diploid_called_raw_allSites_nonPARs_AN20_DP67and201.vcf.gz
tabix -p vcf chrX_GRCh38_YHardMasked_gatk_diploid_called_raw_allSites_nonPARs_AN20_DP67and201.vcf.gz



################################################################################
# FEMALES - DEFAULT #
################################################################################
#-----------------------------#
# DIPLOID CALLING - Autosomes #
#-----------------------------#
ancestry="EUR"
sex="females"

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/joint_called_vcfs/${ancestry}/${sex}/default/

for i in {1..22}
do
  echo "chr${i}"
  # first filter for AN
  bcftools filter -e'INFO/AN<20' chr${i}_GRCh38_default_gatk_diploid_called_raw_allSites.vcf.gz -Oz -o chr${i}_GRCh38_default_gatk_diploid_called_raw_allSites_AN20.vcf.gz
  tabix -p vcf chr${i}_GRCh38_default_gatk_diploid_called_raw_allSites_AN20.vcf.gz
  # then filter for DP
  bcftools filter -e'INFO/DP<67 || INFO/DP>201' chr${i}_GRCh38_default_gatk_diploid_called_raw_allSites_AN20.vcf.gz -Oz -o chr${i}_GRCh38_default_gatk_diploid_called_raw_allSites_AN20_DP67and201.vcf.gz
  tabix -p vcf chr${i}_GRCh38_default_gatk_diploid_called_raw_allSites_AN20_DP67and201.vcf.gz

  # then make bed file? (or can just intersect the vcf with the neutral bed file)
done


#------------#
# X non-PARs #
#------------#
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/joint_called_vcfs/${ancestry}/${sex}/default/

# first remove PARs
gatk SelectVariants -R /scratch/amtarave/variant_calling_simulations_project/reference_genomes/GRCh38_full_analysis_set_plus_decoy_hla_DEFAULT_reformat.fa -V chrX_GRCh38_default_gatk_diploid_called_raw_allSites.vcf.gz -XL /home/amtarave/projects/variant_calling/lists/XPARs.interval.list -O chrX_GRCh38_default_gatk_diploid_called_raw_allSites_nonPARs.vcf.gz

bcftools filter -e'INFO/AN<10' chrX_GRCh38_default_gatk_diploid_called_raw_allSites_nonPARs.vcf.gz -Oz -o chrX_GRCh38_default_gatk_diploid_called_raw_allSites_nonPARs_AN20.vcf.gz
tabix -p vcf chrX_GRCh38_default_gatk_diploid_called_raw_allSites_nonPARs_AN20.vcf.gz
# then filter for DP
bcftools filter -e'INFO/DP<67 || INFO/DP>201' chrX_GRCh38_default_gatk_diploid_called_raw_allSites_nonPARs_AN20.vcf.gz -Oz -o chrX_GRCh38_default_gatk_diploid_called_raw_allSites_nonPARs_AN20_DP67and201.vcf.gz
tabix -p vcf chrX_GRCh38_default_gatk_diploid_called_raw_allSites_nonPARs_AN20_DP67and201.vcf.gz
