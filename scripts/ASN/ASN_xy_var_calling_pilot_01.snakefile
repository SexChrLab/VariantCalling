import os

# Environment: variant_calling_simulations_project

configfile: "ASN_xy_var_calling_pilot_01.config.json"


rule all:
    input:
        expand(os.path.join(config["proj_path"], "refs/GRCh38_SCC/GRCh38_full_analysis_set_plus_decoy_hla_YPARmasked_{chrms}.fa"), chrms=config["chromosomes"]),
        expand(os.path.join(config["proj_path"], "refs/GRCh38_SCC/GRCh38_full_analysis_set_plus_decoy_hla_YPARmasked_{chrms}.fa.fai"), chrms=config["chromosomes"]),

        expand(os.path.join(config["proj_path"], "vcfs/males/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_reformat.vcf"), chrms=config["chromosomes_m"], sample=config["males"]),

        expand(os.path.join(config["proj_path"], "vcfs/males/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrM.recalibrated_variants_{sample}_biallelic_SNPs_reformat.vcf"), sample=config["males"]),

        expand(os.path.join(config["proj_path"], "vcfs/males/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_PARs_reformat.vcf"), chrms=config["chrX"], sample=config["males"]),
        expand(os.path.join(config["proj_path"], "vcfs/males/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_nonPARs_reformat.vcf"), chrms=config["chrX"], sample=config["males"]),

        expand(os.path.join(config["proj_path"], "vcfs/females/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_reformat.vcf"), chrms=config["chromosomes_f"], sample=config["females"]),

        expand(os.path.join(config["proj_path"], "vcfs/females/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrM.recalibrated_variants_{sample}_biallelic_SNPs_reformat.vcf"), sample=config["females"]),

        expand(os.path.join(config["proj_path"], "vcfs/females/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_PARs_reformat.vcf"), chrms=config["chrX"], sample=config["females"]),
        expand(os.path.join(config["proj_path"], "vcfs/females/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_nonPARs_reformat.vcf"), chrms=config["chrX"], sample=config["females"])


#------------------------------------------------------------------------------#
# Step: Separate reference genomes by chromosome
#------------------------------------------------------------------------------#
# Extract chromosomes from reference fasta and then index. For simulations.
rule prep_ref_genome:
    input:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"]
        #ref = config["genome_paths"]["Ref_GRCh38_Default"] # use default so we can simulate PARs on both X and Y and then for variant calling, align to SCC ref (so that Y PAR reads align to X (PAR) chr)
    params:
        chrm = "{chrms}"
    output:
        fa = os.path.join(config["proj_path"], "refs/GRCh38_SCC/GRCh38_full_analysis_set_plus_decoy_hla_YPARmasked_{chrms}.fa"),
        fai = os.path.join(config["proj_path"], "refs/GRCh38_SCC/GRCh38_full_analysis_set_plus_decoy_hla_YPARmasked_{chrms}.fa.fai")
    shell:
        """
        samtools faidx {input} {params.chrm} > {output.fa};
        samtools faidx {output.fa}
        """

#------------------------------------------------------------------------------#
# Step: Extract samples from 1000 Genomes VCFs
#------------------------------------------------------------------------------#
# MALES #
# Select biallelic snps only for the specific sample we are analyzing.
rule extract_snps_and_samples_males:
    input:
        vcf = os.path.join(config["1000genomes_vcf_path"], "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        smpl = "{sample}"
    output:
        vcf = temp(os.path.join(config["proj_path"], "vcfs/males/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs.vcf")),
    shell:
        """
        gatk --java-options '-Xmx10g' SelectVariants -R {params.ref} -sn {params.smpl} -V {input.vcf} -O {output.vcf} --select-type-to-include SNP --restrict-alleles-to BIALLELIC
        """


# Select biallelic snps only for the specific sample we are analyzing.
rule extract_snps_and_samples_males_chrM:
    input:
        vcf = os.path.join(config["1000genomes_vcf_path"], "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_others.recalibrated_variants.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        smpl = "{sample}"
    output:
        vcf = temp(os.path.join(config["proj_path"], "vcfs/males/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrM.recalibrated_variants_{sample}_biallelic_SNPs.vcf"))
    shell:
        """
        gatk --java-options '-Xmx10g' SelectVariants -R {params.ref} -sn {params.smpl} -L chrM -V {input.vcf} -O {output.vcf} --select-type-to-include SNP --restrict-alleles-to BIALLELIC
        """

# reformat VCF to replace FILTER column with PASS so that NEAT recognizes variants
rule reformat_vcfs_males_all:
    input:
        vcf = os.path.join(config["proj_path"], "vcfs/males/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs.vcf")
    output:
        vcf = os.path.join(config["proj_path"], "vcfs/males/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_reformat.vcf")
    shell:
        """
        awk '{{ if ($1 ~ "#") print $0; else print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t.\\tPASS\\t.\\t"$9"\\t"$10}}' {input.vcf} > {output.vcf}
        """

'''
# reformat VCF to replace FILTER column with PASS so that NEAT recognizes variants
rule reformat_vcfs_males_chrM:
    input:
        vcf = os.path.join(config["proj_path"], "vcfs/males/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrM.recalibrated_variants_{sample}_biallelic_SNPs.vcf")
    output:
        vcf = os.path.join(config["proj_path"], "vcfs/males/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrM.recalibrated_variants_{sample}_biallelic_SNPs_reformat.vcf"),
    shell:
        """
        awk '{{ if ($1 ~ "#") print $0; else print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t.\\tPASS\\t.\\t"$9"\\t"$10}}' {input.vcf} > {output.vcf}
        """

'''

# In VCF extract PARs and nonPARs and run separately
rule prep_X_vcf_males:
    input:
        vcf = os.path.join(config["proj_path"], "vcfs/males/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs.vcf"),
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        smpl = "{sample}"
    output:
        vcfPARs = temp(os.path.join(config["proj_path"], "vcfs/males/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_PARs.vcf")),
        vcfnonPARs = temp(os.path.join(config["proj_path"], "vcfs/males/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_nonPARs.vcf"))
    shell:
        """
        gatk --java-options '-Xmx10g' SelectVariants -R {params.ref} -sn {params.smpl} -L chrX:10001-2781479 -L chrX:155701383-156030895 -V {input.vcf} -O {output.vcfPARs} --select-type-to-include SNP --restrict-alleles-to BIALLELIC;
        gatk --java-options '-Xmx10g' SelectVariants -R {params.ref} -sn {params.smpl} -XL chrX:10001-2781479 -XL chrX:155701383-156030895 -V {input.vcf} -O {output.vcfnonPARs} --select-type-to-include SNP --restrict-alleles-to BIALLELIC
        """

# reformat VCF to replace FILTER column with PASS so that NEAT recognizes variants
rule reformat_vcfs_X_males:
    input:
        vcfPARs = os.path.join(config["proj_path"], "vcfs/males/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_PARs.vcf"),
        vcfnonPARs = os.path.join(config["proj_path"], "vcfs/males/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_nonPARs.vcf")
    output:
        vcfPARs = os.path.join(config["proj_path"], "vcfs/males/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_PARs_reformat.vcf"),
        vcfnonPARs = os.path.join(config["proj_path"], "vcfs/males/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_nonPARs_reformat.vcf")
    shell:
        """
        awk '{{ if ($1 ~ "#") print $0; else print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t.\\tPASS\\t.\\t"$9"\\t"$10}}' {input.vcfPARs} > {output.vcfPARs};
        awk '{{ if ($1 ~ "#") print $0; else print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t.\\tPASS\\t.\\t"$9"\\t"$10}}' {input.vcfnonPARs} > {output.vcfnonPARs}
        """

################################################################################
# FEMALES #
################################################################################
# Select biallelic snps only for the specific sample we are analyzing.
rule extract_snps_and_samples_females:
    input:
        vcf = os.path.join(config["1000genomes_vcf_path"], "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        smpl = "{sample}"
    output:
        vcf = temp(os.path.join(config["proj_path"], "vcfs/females/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs.vcf")),
    shell:
        """
        gatk --java-options '-Xmx10g' SelectVariants -R {params.ref} -sn {params.smpl} -V {input.vcf} -O {output.vcf} --select-type-to-include SNP --restrict-alleles-to BIALLELIC
        """

rule reformat_vcfs_females_all:
    input:
        vcf = os.path.join(config["proj_path"], "vcfs/females/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs.vcf")
    output:
        vcf = os.path.join(config["proj_path"], "vcfs/females/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_reformat.vcf"),
    shell:
        """
        awk '{{ if ($1 ~ "#") print $0; else print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t.\\tPASS\\t.\\t"$9"\\t"$10}}' {input.vcf} > {output.vcf}
        """


# Select biallelic snps only for the specific sample we are analyzing.
rule extract_snps_and_samples_females_chrM:
    input:
        vcf = os.path.join(config["1000genomes_vcf_path"], "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_others.recalibrated_variants.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        smpl = "{sample}"
    output:
        vcf = temp(os.path.join(config["proj_path"], "vcfs/females/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrM.recalibrated_variants_{sample}_biallelic_SNPs.vcf")),
    shell:
        """
        gatk --java-options '-Xmx10g' SelectVariants -R {params.ref} -sn {params.smpl} -L chrM -V {input.vcf} -O {output.vcf} --select-type-to-include SNP --restrict-alleles-to BIALLELIC
        """

'''
rule reformat_vcfs_females_chrM:
    input:
        vcf = os.path.join(config["proj_path"], "vcfs/females/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrM.recalibrated_variants_{sample}_biallelic_SNPs.vcf")
    output:
        vcf = os.path.join(config["proj_path"], "vcfs/females/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chrM.recalibrated_variants_{sample}_biallelic_SNPs_reformat.vcf"),
    shell:
        """
        awk '{{ if ($1 ~ "#") print $0; else print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t.\\tPASS\\t.\\t"$9"\\t"$10}}' {input.vcf} > {output.vcf}
        """
'''

# In VCF extract PARs and nonPARs and run separately
rule prep_X_vcf_females:
    input:
        vcf = os.path.join(config["proj_path"], "vcfs/females/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs.vcf"),
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        smpl = "{sample}"
    output:
        vcfPARs = temp(os.path.join(config["proj_path"], "vcfs/females/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_PARs.vcf")),
        vcfnonPARs = temp(os.path.join(config["proj_path"], "vcfs/females/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_nonPARs.vcf"))
    shell:
        """
        gatk --java-options '-Xmx10g' SelectVariants -R {params.ref} -sn {params.smpl} -L chrX:10001-2781479 -L chrX:155701383-156030895 -V {input.vcf} -O {output.vcfPARs} --select-type-to-include SNP --restrict-alleles-to BIALLELIC;
        gatk --java-options '-Xmx10g' SelectVariants -R {params.ref} -sn {params.smpl} -XL chrX:10001-2781479 -XL chrX:155701383-156030895 -V {input.vcf} -O {output.vcfnonPARs} --select-type-to-include SNP --restrict-alleles-to BIALLELIC
        """

rule reformat_vcfs_X_females:
    input:
        vcfPARs = os.path.join(config["proj_path"], "vcfs/females/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_PARs.vcf"),
        vcfnonPARs = os.path.join(config["proj_path"], "vcfs/females/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_nonPARs.vcf")
    output:
        vcfPARs = os.path.join(config["proj_path"], "vcfs/females/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_PARs_reformat.vcf"),
        vcfnonPARs = os.path.join(config["proj_path"], "vcfs/females/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_nonPARs_reformat.vcf")
    shell:
        """
        awk '{{ if ($1 ~ "#") print $0; else print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t.\\tPASS\\t.\\t"$9"\\t"$10}}' {input.vcfPARs} > {output.vcfPARs};
        awk '{{ if ($1 ~ "#") print $0; else print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t.\\tPASS\\t.\\t"$9"\\t"$10}}' {input.vcfnonPARs} > {output.vcfnonPARs}
        """
