#!/bin/bash
#SBATCH --job-name=calculate_pi # Job name
#SBATCH -o slurm.calculate_pi.out
#SBATCH -e slurm.calculate_pi.err
#SBATCH -n 2
#SBATCH --mem-per-cpu=16000
#SBATCH -t 48:00:00
#SBATCH --mail-user=amtarave@asu.edu # send-to address

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/scripts

conda activate variant_calling_simulations_project

# Only calculating pi in neutral regions (not also high quality regions because
# I am still figuring out how to do this)

# Path to neutral sites bed file:
neutralbed="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/UCSC_genome_browser/for_diversity_calculation/all.chromosomes.putatively.neutral.bed"
neutralbedpath="/data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/UCSC_genome_browser/for_diversity_calculation/"
#grep chrX ${neutralbed} > ${neutralbedpath}chrX_neutral.bed
#grep -v chrY ${neutralbed} | grep -v chrX > ${neutralbedpath}autosomes_neutral.bed


################################################################################
# MALES - SCC #
################################################################################
#-----------------------------#
# DIPLOID CALLING - Autosomes #
#-----------------------------#
ancestry="EUR"
sex="males"
ploidy="diploid"

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/hard_filtered_vcfs/${ancestry}/${sex}/all/

# merge all autosomes together: vcfs and bed files
bcftools concat chr1_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr2_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr3_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr4_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr5_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr6_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr7_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr8_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr9_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr10_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr11_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr12_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr13_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr14_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr15_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr16_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr17_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr18_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr19_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr20_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr21_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr22_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz -Oz -o autosomes_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz


# Then run diversity script on merged vcf
python3 /home/amtarave/packages/popgen_tools/popgen_tools.py --vcf_file autosomes_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz --pi --ploidy ${ploidy} --pi_target --target_bed /data/CEM/wilsonlab/projects/variant_calling_simulations_project/joint_called_vcfs/${ancestry}/${sex}/autosomes_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed --pi_target_out autosomes_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_merge_neutral_highqual_pi_target_out_02.txt --pi_target_per_site_out autosomes_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_merge_neutral_highqual_pi_target_per_site_out_02.txt

python3 /home/amtarave/packages/popgen_tools/popgen_tools.py --vcf_file autosomes_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz --pi --ploidy ${ploidy} --pi_target --target_bed /data/CEM/wilsonlab/projects/variant_calling_simulations_project/joint_called_vcfs/${ancestry}/${sex}/autosomes_GRCh38_YPARsMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual_sorted.bed --pi_target_out autosomes_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_merge_neutral_highqual_pi_target_out_sorted.txt --pi_target_per_site_out autosomes_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_merge_neutral_highqual_pi_target_per_site_out_sorted.txt

#------------------------------#
# HAPLOID CALLING - X non-PARs #
#------------------------------#
ploidy="haploid"

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/hard_filtered_vcfs/${ancestry}/${sex}/all/


python3 /home/amtarave/packages/popgen_tools/popgen_tools.py --vcf_file chrX_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_gatkHardFilter_filtered_nonPARs.vcf.gz --pi --ploidy ${ploidy} --pi_target --target_bed /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/joint_called_vcfs/${ancestry}/${sex}/chrX_GRCh38_YPARsMasked_gatk_haploid_called_raw_allSites_nonPARs_AN10_DP23and69_merge_neutral_highqual.bed --pi_target_out chrX_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_gatkHardFilter_filtered_nonPARs_merge_neutral_highqual_pi_target_out_02.txt --pi_target_per_site_out chrX_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_gatkHardFilter_filtered_nonPARs_merge_neutral_highqual_pi_target_per_site_out_02.txt

##############
# Remove XTR #
##############
bedtools subtract -a chrX_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_gatkHardFilter_filtered_nonPARs.vcf.gz -b /home/amtarave/projects/variant_calling/lists/X_XTR.interval_fix.bed | cut -f1,2 > chrX_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_gatkHardFilter_filtered_nonPARs_nonXTRpositions.txt

bcftools view -R chrX_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_gatkHardFilter_filtered_nonPARs_nonXTRpositions.txt chrX_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_gatkHardFilter_filtered_nonPARs.vcf.gz -Oz -o chrX_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_gatkHardFilter_filtered_nonPARs_nonXTR.vcf.gz


python3 /home/amtarave/packages/popgen_tools/popgen_tools.py --vcf_file chrX_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_gatkHardFilter_filtered_nonPARs_nonXTR.vcf.gz --pi --ploidy ${ploidy} --pi_target --target_bed /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/joint_called_vcfs/${ancestry}/${sex}/chrX_GRCh38_YPARsMasked_gatk_haploid_called_raw_allSites_nonPARs_AN10_DP23and69_merge_neutral_highqual.bed --pi_target_out chrX_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_gatkHardFilter_filtered_nonPARs_nonXTR_merge_neutral_highqual_pi_target_out_02.txt --pi_target_per_site_out chrX_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_gatkHardFilter_filtered_nonPARs_nonXTR_merge_neutral_highqual_pi_target_per_site_out_02.txt


################################################################################
# MALES - DEFAULT #
################################################################################
#-----------------------------#
# DIPLOID CALLING - Autosomes #
#-----------------------------#
ancestry="EUR"
sex="males"
ploidy="diploid"

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/hard_filtered_vcfs/${ancestry}/${sex}/default/all/

# merge all autosomes together: vcfs and bed files
bcftools concat chr1_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr2_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr3_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr4_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr5_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr6_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr7_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr8_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr9_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr10_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr11_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr12_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr13_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr14_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr15_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr16_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr17_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr18_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr19_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr20_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr21_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr22_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz -Oz -o autosomes_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz


# Then run diversity script on merged vcf
python3 /home/amtarave/packages/popgen_tools/popgen_tools.py --vcf_file autosomes_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz --pi --ploidy ${ploidy} --pi_target --target_bed /data/CEM/wilsonlab/projects/variant_calling_simulations_project/joint_called_vcfs/${ancestry}/${sex}/default/autosomes_GRCh38_default_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed --pi_target_out autosomes_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_merge_neutral_highqual_pi_target_out_02.txt --pi_target_per_site_out autosomes_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_merge_neutral_highqual_pi_target_per_site_out_02.txt

#------------------------------#
# HAPLOID CALLING - X non-PARs #
#------------------------------#
ploidy="haploid"

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/hard_filtered_vcfs/${ancestry}/${sex}/default/all/


python3 /home/amtarave/packages/popgen_tools/popgen_tools.py --vcf_file chrX_GRCh38_default_gatk_haploid_called_SNPs_gatkHardFilter_filtered_nonPARs.vcf.gz --pi --ploidy ${ploidy} --pi_target --target_bed /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/joint_called_vcfs/${ancestry}/${sex}/default/chrX_GRCh38_default_gatk_haploid_called_raw_allSites_nonPARs_AN10_DP23and69_merge_neutral_highqual.bed --pi_target_out chrX_GRCh38_default_gatk_haploid_called_SNPs_gatkHardFilter_filtered_nonPARs_merge_neutral_highqual_pi_target_out_02.txt --pi_target_per_site_out chrX_GRCh38_default_gatk_haploid_called_SNPs_gatkHardFilter_filtered_nonPARs_merge_neutral_highqual_pi_target_per_site_out_02.txt

##############
# Remove XTR #
##############
bedtools subtract -a chrX_GRCh38_default_gatk_haploid_called_SNPs_gatkHardFilter_filtered_nonPARs.vcf.gz -b /home/amtarave/projects/variant_calling/lists/X_XTR.interval_fix.bed | cut -f1,2 > chrX_GRCh38_default_gatk_haploid_called_SNPs_gatkHardFilter_filtered_nonPARs_nonXTRpositions.txt

bcftools view -R chrX_GRCh38_default_gatk_haploid_called_SNPs_gatkHardFilter_filtered_nonPARs_nonXTRpositions.txt chrX_GRCh38_default_gatk_haploid_called_SNPs_gatkHardFilter_filtered_nonPARs.vcf.gz -Oz -o chrX_GRCh38_default_gatk_haploid_called_SNPs_gatkHardFilter_filtered_nonPARs_nonXTR.vcf.gz


python3 /home/amtarave/packages/popgen_tools/popgen_tools.py --vcf_file chrX_GRCh38_default_gatk_haploid_called_SNPs_gatkHardFilter_filtered_nonPARs_nonXTR.vcf.gz --pi --ploidy ${ploidy} --pi_target --target_bed /data/CEM/wilsonlab/projects/variant_calling_simulations_project/haploid/joint_called_vcfs/${ancestry}/${sex}/default/chrX_GRCh38_default_gatk_haploid_called_raw_allSites_nonPARs_AN10_DP23and69_merge_neutral_highqual.bed --pi_target_out chrX_GRCh38_default_gatk_haploid_called_SNPs_gatkHardFilter_filtered_nonPARs_nonXTR_merge_neutral_highqual_pi_target_out_02.txt --pi_target_per_site_out chrX_GRCh38_default_gatk_haploid_called_SNPs_gatkHardFilter_filtered_nonPARs_nonXTR_merge_neutral_highqual_pi_target_per_site_out_02.txt


################################################################################
# FEMALES - SCC #
################################################################################
#-----------------------------#
# DIPLOID CALLING - Autosomes #
#-----------------------------#
ancestry="EUR"
sex="females"
ploidy="diploid"

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/hard_filtered_vcfs/${ancestry}/${sex}/all/

# merge all autosomes together: vcfs and bed files
bcftools concat chr1_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr2_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr3_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr4_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr5_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr6_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr7_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr8_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr9_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr10_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr11_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr12_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr13_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr14_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr15_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr16_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr17_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr18_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr19_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr20_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr21_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr22_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz -Oz -o autosomes_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz


# Then run diversity script on merged vcf
python3 /home/amtarave/packages/popgen_tools/popgen_tools.py --vcf_file autosomes_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz --pi --ploidy ${ploidy} --pi_target --target_bed /data/CEM/wilsonlab/projects/variant_calling_simulations_project/joint_called_vcfs/${ancestry}/${sex}/autosomes_GRCh38_YHardMasked_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed --pi_target_out autosomes_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_merge_neutral_highqual_pi_target_out_02.txt --pi_target_per_site_out autosomes_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_merge_neutral_highqual_pi_target_per_site_out_02.txt

#grep chr8 ${neutralbedpath}autosomes_neutral.bed > ${neutralbedpath}chr8_neutral.bed
#python3 /home/amtarave/packages/popgen_tools/popgen_tools.py --vcf_file chr8_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz --pi --ploidy ${ploidy} --pi_target --target_bed ${neutralbedpath}chr8_neutral.bed --pi_target_out chr8_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_merge_neutral_highqual_pi_target_out_02.txt --pi_target_per_site_out chr8_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_merge_neutral_highqual_pi_target_per_site_out_02.txt


#-------------#
# X non-PARs #
#-------------#
python3 /home/amtarave/packages/popgen_tools/popgen_tools.py --vcf_file chrX_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs.vcf.gz --pi --ploidy ${ploidy} --pi_target --target_bed /data/CEM/wilsonlab/projects/variant_calling_simulations_project/joint_called_vcfs/${ancestry}/${sex}/chrX_GRCh38_YHardMasked_gatk_diploid_called_raw_allSites_nonPARs_AN20_DP67and201_merge_neutral_highqual.bed --pi_target_out chrX_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs_merge_neutral_highqual_pi_target_out_AN20_DP67and201.txt --pi_target_per_site_out chrX_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs_merge_neutral_highqual_pi_target_per_site_out_AN20_DP67and201.txt

##############
# Remove XTR #
##############
bedtools subtract -a chrX_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs.vcf.gz -b /home/amtarave/projects/variant_calling/lists/X_XTR.interval_fix.bed | cut -f1,2 > chrX_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs_nonXTRpositions.txt

bcftools view -R chrX_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs_nonXTRpositions.txt chrX_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs.vcf.gz -Oz -o chrX_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs_nonXTR.vcf.gz


python3 /home/amtarave/packages/popgen_tools/popgen_tools.py --vcf_file chrX_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs_nonXTR.vcf.gz --pi --ploidy ${ploidy} --pi_target --target_bed /data/CEM/wilsonlab/projects/variant_calling_simulations_project/joint_called_vcfs/${ancestry}/${sex}/chrX_GRCh38_YHardMasked_gatk_diploid_called_raw_allSites_nonPARs_AN20_DP67and201_merge_neutral_highqual.bed --pi_target_out chrX_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs_nonXTR_merge_neutral_highqual_pi_target_out_AN20_DP67and201.txt --pi_target_per_site_out chrX_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs_nonXTR_merge_neutral_highqual_pi_target_per_site_out_AN20_DP67and201.txt


################################################################################
# FEMALES - DEFAULT #
################################################################################
#-----------------------------#
# DIPLOID CALLING - Autosomes #
#-----------------------------#
ancestry="EUR"
sex="females"
ploidy="diploid"

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/hard_filtered_vcfs/${ancestry}/${sex}/default/all/

# merge all autosomes together: vcfs and bed files
bcftools concat chr1_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr2_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr3_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr4_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr5_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr6_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr7_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr8_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr9_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr10_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr11_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr12_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr13_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr14_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr15_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr16_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr17_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr18_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr19_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr20_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr21_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz chr22_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz -Oz -o autosomes_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz


# Then run diversity script on merged vcf
python3 /home/amtarave/packages/popgen_tools/popgen_tools.py --vcf_file autosomes_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz --pi --ploidy ${ploidy} --pi_target --target_bed /data/CEM/wilsonlab/projects/variant_calling_simulations_project/joint_called_vcfs/${ancestry}/${sex}/default/autosomes_GRCh38_default_gatk_diploid_called_raw_allSites_AN20_DP67and201_merge_neutral_highqual.bed --pi_target_out autosomes_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_merge_neutral_highqual_pi_target_out_02.txt --pi_target_per_site_out autosomes_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_merge_neutral_highqual_pi_target_per_site_out_02.txt

#-------------#
# X non-PARs #
#-------------#
python3 /home/amtarave/packages/popgen_tools/popgen_tools.py --vcf_file chrX_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs.vcf.gz --pi --ploidy ${ploidy} --pi_target --target_bed /data/CEM/wilsonlab/projects/variant_calling_simulations_project/joint_called_vcfs/${ancestry}/${sex}/default/chrX_GRCh38_default_gatk_diploid_called_raw_allSites_nonPARs_AN20_DP67and201_merge_neutral_highqual.bed --pi_target_out chrX_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs_merge_neutral_highqual_pi_target_out_02.txt --pi_target_per_site_out chrX_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs_merge_neutral_highqual_pi_target_per_site_out_02.txt


##############
# Remove XTR #
##############
bedtools subtract -a chrX_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs.vcf.gz -b /home/amtarave/projects/variant_calling/lists/X_XTR.interval_fix.bed | cut -f1,2 > chrX_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs_nonXTRpositions.txt

bcftools view -R chrX_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs_nonXTRpositions.txt chrX_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs.vcf.gz -Oz -o chrX_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs_nonXTR.vcf.gz


python3 /home/amtarave/packages/popgen_tools/popgen_tools.py --vcf_file chrX_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs_nonXTR.vcf.gz --pi --ploidy ${ploidy} --pi_target --target_bed /data/CEM/wilsonlab/projects/variant_calling_simulations_project/joint_called_vcfs/${ancestry}/${sex}/default/chrX_GRCh38_default_gatk_diploid_called_raw_allSites_nonPARs_AN20_DP67and201_merge_neutral_highqual.bed --pi_target_out chrX_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs_nonXTR_merge_neutral_highqual_pi_target_out_02.txt --pi_target_per_site_out chrX_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs_nonXTR_merge_neutral_highqual_pi_target_per_site_out_02.txt
