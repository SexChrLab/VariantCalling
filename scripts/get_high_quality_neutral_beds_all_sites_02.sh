#!/bin/bash
#SBATCH --job-name=get_high_quality_neutral_beds # Job name
#SBATCH -o slurm.get_high_quality_neutral_beds.out
#SBATCH -e slurm.get_high_quality_neutral_beds.err
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem-per-cpu=16000
#SBATCH --mail-user=amtarave@asu.edu # send-to address

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/scripts

conda activate variant_calling_simulations_project

# For each of the VCFs generated in get_high_quality_sites.sh run bedtools
# intersect with the putatively neutral bed file to get bed file of neutral sites
# that are also high quality.

# Path to neutral sites bed file:
neutralbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/UCSC_genome_browser/for_diversity_calculation/all.chromosomes.putatively.neutral.bed"
neutralbedpath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/UCSC_genome_browser/for_diversity_calculation/"
#grep chrX ${neutralbed} > ${neutralbedpath}chrX_neutral.bed
#grep -v chrY ${neutralbed} | grep -v chrX > ${neutralbedpath}autosomes_neutral.bed

# Path to output bed files will be where the vcf files are

################################################################################
# MALES - SCC #
################################################################################
#-----------------------------#
# DIPLOID CALLING - Autosomes #
#-----------------------------#
ancestry="EUR"
sex="males"
ploidy="diploid"

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/joint_called_vcfs/${ancestry}/${sex}/

for i in {1..22}
do
	echo "chr${i}"

  bedtools merge -i chr${i}_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201.vcf.gz > chr${i}_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge.bed

  bedtools intersect -a ${neutralbed} -b chr${i}_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge.bed > chr${i}_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed

done

# merge all autosomes together: vcfs and bed files
cat chr*_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed | sort -k1,1 -k2,2n > autosomes_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed

# note to self --- bed filr does not appear to be sorted (chr 9 is at the bottom of the file)
#cat chr1_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed chr2_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed chr3_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed chr4_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed chr5_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed chr6_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed chr7_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed chr8_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed chr9_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed chr10_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed chr11_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed chr12_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed chr13_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed chr14_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed chr15_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed chr16_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed chr17_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed chr18_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed chr19_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed chr20_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed chr21_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed chr22_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed > autosomes_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual_sorted.bed
#------------------------------#
# HAPLOID CALLING - X non-PARs #
#------------------------------#
ploidy="haploid"

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/joint_called_vcfs/${ancestry}/${sex}/

bedtools merge -i chrX_GRCh38_YPARsMasked_gatk_haploid_called_raw_allSites_nonPARs_AN10_DP23and69.vcf.gz > chrX_GRCh38_YPARsMasked_gatk_haploid_called_raw_allSites_nonPARs_AN10_DP23and69_merge.bed

bedtools intersect -a ${neutralbed} -b chrX_GRCh38_YPARsMasked_gatk_haploid_called_raw_allSites_nonPARs_AN10_DP23and69_merge.bed > chrX_GRCh38_YPARsMasked_gatk_haploid_called_raw_allSites_nonPARs_AN10_DP23and69_merge_neutral_highqual.bed


################################################################################
# MALES - DEFAULT #
################################################################################
#-----------------------------#
# DIPLOID CALLING - Autosomes #
#-----------------------------#
ancestry="EUR"
sex="males"
ploidy="diploid"

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/joint_called_vcfs/${ancestry}/${sex}/default/

#for i in {1..11}
#do
#	#echo "${i}"
#	bedtools merge -i chr${i}_GRCh38_default_gatk_diploid_called_raw_allSites_AN20_DP67and201.vcf.gz > chr${i}_GRCh38_default_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge.bed
#	bedtools intersect -a ${neutralbed} -b chr${i}_GRCh38_default_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge.bed > chr${i}_GRCh38_default_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed
#done
#for i in {21..22}; do echo "${i}"; bedtools merge -i chr${i}_GRCh38_default_gatk_diploid_called_raw_allSites_AN20_DP67and201.vcf.gz > chr${i}_GRCh38_default_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge.bed; bedtools intersect -a ${neutralbed} -b chr${i}_GRCh38_default_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge.bed > chr${i}_GRCh38_default_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed; done

for i in {1..22}
do
	echo "chr${i}"

	bedtools merge -i chr${i}_GRCh38_default_gatk_diploid_called_raw_allSites_AN20_DP67and201.vcf.gz > chr${i}_GRCh38_default_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge.bed
	bedtools intersect -a ${neutralbed} -b chr${i}_GRCh38_default_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge.bed > chr${i}_GRCh38_default_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed

done

# merge all autosomes together: vcfs and bed files
cat chr*_GRCh38_default_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed | sort -k1,1 -k2,2n > autosomes_GRCh38_default_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed


#------------------------------#
# HAPLOID CALLING - X non-PARs #
#------------------------------#
ploidy="haploid"

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/joint_called_vcfs/${ancestry}/${sex}/default/

bedtools merge -i chrX_GRCh38_default_gatk_haploid_called_raw_allSites_nonPARs_AN10_DP23and69.vcf.gz > chrX_GRCh38_default_gatk_haploid_called_raw_allSites_nonPARs_AN10_DP23and69_merge.bed

bedtools intersect -a ${neutralbed} -b chrX_GRCh38_default_gatk_haploid_called_raw_allSites_nonPARs_AN10_DP23and69_merge.bed > chrX_GRCh38_default_gatk_haploid_called_raw_allSites_nonPARs_AN10_DP23and69_merge_neutral_highqual.bed



################################################################################
# FEMALES - SCC #
################################################################################
#-----------------------------#
# DIPLOID CALLING - Autosomes #
#-----------------------------#
ancestry="EUR"
sex="females"
ploidy="diploid"

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/joint_called_vcfs/${ancestry}/${sex}/

for i in {1..22}
do
	echo "chr${i}"

  bedtools merge -i chr${i}_GRCh38_YHardMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201.vcf.gz > chr${i}_GRCh38_YHardMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge.bed

  bedtools intersect -a ${neutralbed} -b chr${i}_GRCh38_YHardMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge.bed > chr${i}_GRCh38_YHardMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed

done

# merge all autosomes together: vcfs and bed files
cat chr*_GRCh38_YHardMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed | sort -k1,1 -k2,2n > autosomes_GRCh38_YHardMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed


#------------#
# X non-PARs #
#------------#
bedtools merge -i chrX_GRCh38_YHardMasked_gatk_diploid_called_raw_allSites_nonPARs_AN20_DP67and201.vcf.gz > chrX_GRCh38_YHardMasked_gatk_diploid_called_raw_allSites_nonPARs_AN20_DP67and201_merge.bed

bedtools intersect -a ${neutralbed} -b chrX_GRCh38_YHardMasked_gatk_diploid_called_raw_allSites_nonPARs_AN20_DP67and201_merge.bed > chrX_GRCh38_YHardMasked_gatk_diploid_called_raw_allSites_nonPARs_AN20_DP67and201_merge_neutral_highqual.bed


################################################################################
# FEMALES - DEFAULT #
################################################################################
#-----------------------------#
# DIPLOID CALLING - Autosomes #
#-----------------------------#
ancestry="EUR"
sex="females"
ploidy="diploid"

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/joint_called_vcfs/${ancestry}/${sex}/default/

for i in {1..22}
do
	echo "chr${i}"

  bedtools merge -i chr${i}_GRCh38_default_gatk_diploid_called_raw_allSites_AN20_DP67and201.vcf.gz > chr${i}_GRCh38_default_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge.bed

  bedtools intersect -a ${neutralbed} -b chr${i}_GRCh38_default_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge.bed > chr${i}_GRCh38_default_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed

done

# merge all autosomes together: vcfs and bed files
cat chr*_GRCh38_default_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed | sort -k1,1 -k2,2n > autosomes_GRCh38_default_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed


#------------#
# X non-PARs #
#------------#
bedtools merge -i chrX_GRCh38_default_gatk_diploid_called_raw_allSites_nonPARs_AN20_DP67and201.vcf.gz > chrX_GRCh38_default_gatk_diploid_called_raw_allSites_nonPARs_AN20_DP67and201_merge.bed

bedtools intersect -a ${neutralbed} -b chrX_GRCh38_default_gatk_diploid_called_raw_allSites_nonPARs_AN20_DP67and201_merge.bed > chrX_GRCh38_default_gatk_diploid_called_raw_allSites_nonPARs_AN20_DP67and201_merge_neutral_highqual.bed
