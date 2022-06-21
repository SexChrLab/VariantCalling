import os

# Environment: NEAT_env

configfile: "AFR_xy_var_calling_pilot_03_haploid.config.json"

# Notes #

        #expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/QD/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_QD_{filter}.vcf.gz"), chrms=config["chromosomes"], filtering_options=config["filtering_options"]),
        #expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/QUAL/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_QUAL_{filter}.vcf.gz"), chrms=config["chromosomes"], filtering_options=config["filtering_options"]),
        #expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/SOR/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_SOR_{filter}.vcf.gz"), chrms=config["chromosomes"], filtering_options=config["filtering_options"]),
        #expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/FS/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_FS_{filter}.vcf.gz"), chrms=config["chromosomes"], filtering_options=config["filtering_options"]),
        #expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/MQ/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_MQ_{filter}.vcf.gz"), chrms=config["chromosomes"], filtering_options=config["filtering_options"]),
        #expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/MQRankSum/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_MQRankSum_{filter}.vcf.gz"), chrms=config["chromosomes"], filtering_options=config["filtering_options"]),
        #expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/ReadPosRankSum/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_ReadPosRankSum_{filter}.vcf.gz"), chrms=config["chromosomes"], filtering_options=config["filtering_options"]),

        #expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/QD/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_QD_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filtering_options=config["filtering_options"]),
        #expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/QUAL/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_QUAL_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filtering_options=config["filtering_options"]),
        #expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/SOR/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_SOR_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filtering_options=config["filtering_options"]),
        #expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/FS/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_FS_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filtering_options=config["filtering_options"]),
        #expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/MQ/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_MQ_nonPARs_{filter}.vcf.gz"), chrms=config["chrX"], filtering_options=config["filtering_options"]),
        #expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/MQRankSum/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_MQRankSum_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filtering_options=config["filtering_options"]),
        #expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/ReadPosRankSum/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_ReadPosRankSum_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filtering_options=config["filtering_options"])
#wildcard_constraints:
#    filterQD= '|'.join([re.escape(x) for x in config["QD_filters"]]),
#    chrmsXnonPARs= '|'.join([re.escape(x) for x in config["chrX"]]),
#    chrms= '|'.join([re.escape(x) for x in config["chromosomes"]])


rule all:
    input:
        expand(os.path.join(config["proj_path"], "haploid/gvcfs/AFR/males/{sample}_GRCh38_YPARsMasked_haploid.g.vcf.gz"), sample=config["males"]),

        expand(os.path.join(config["proj_path"], "haploid/combine_g_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_haploid_combined.g.vcf.gz"), chrms=config["chromosomes"]),
        expand(os.path.join(config["proj_path"], "haploid/joint_called_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_raw.vcf.gz"), chrms=config["chromosomes"]),

        expand(os.path.join(config["proj_path"], "haploid/joint_called_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs.vcf.gz"), chrms=config["chromosomes"]),

        expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/all/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_{filtering_options}_filtered.vcf.gz"), chrms=config["chromosomes"], filtering_options=config["filtering_options"]),

        expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/all/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_{filtering_options}_filtered_nonPARs.vcf.gz"), chrms=config["chrX"], filtering_options=config["filtering_options"]),

        #expand(os.path.join(config["proj_path"], "haploid/bcftools_stats/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_{filtering_options}_filtered_nonPARs.bcftools.stats.PSC.txt"), chrms=config["chrX"], filtering_options=config["filtering_options"]),
        #expand(os.path.join(config["proj_path"], "haploid/bcftools_stats/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_{filtering_options}_filtered.bcftools.stats.PSC.txt"), chrms=config["chrY"], filtering_options=config["filtering_options"]),

        # MALES - DEFAULT MAPPING #
        expand(os.path.join(config["proj_path"], "haploid/gvcfs/AFR/males/default/{sample}_GRCh38_default_haploid.g.vcf.gz"), sample=config["males"]),

        expand(os.path.join(config["proj_path"], "haploid/combine_g_vcfs/AFR/males/default/{chrms}_GRCh38_default_haploid_combined.g.vcf.gz"), chrms=config["chromosomes"]),
        expand(os.path.join(config["proj_path"], "haploid/joint_called_vcfs/AFR/males/default/{chrms}_GRCh38_default_gatk_haploid_called_raw.vcf.gz"), chrms=config["chromosomes"]),

        expand(os.path.join(config["proj_path"], "haploid/joint_called_vcfs/AFR/males/default/{chrms}_GRCh38_default_gatk_haploid_called_raw_SNPs.vcf.gz"), chrms=config["chromosomes"]),

        expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/default/all/{chrms}_GRCh38_default_gatk_haploid_called_SNPs_{filtering_options}_filtered.vcf.gz"), chrms=config["chromosomes"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/default/all/{chrms}_GRCh38_default_gatk_haploid_called_SNPs_{filtering_options}_filtered_nonPARs.vcf.gz"), chrms=config["chrX"], filtering_options=config["filtering_options"]),

'''
        expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/QD/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_QD_{filter}.vcf.gz"), chrms=config["chromosomes"], filter=config["QD_filters"]),
        expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/QUAL/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_QUAL_{filter}.vcf.gz"), chrms=config["chromosomes"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/SOR/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_SOR_{filter}.vcf.gz"), chrms=config["chromosomes"], filter=config["SOR_filters"]),
        expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/FS/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_FS_{filter}.vcf.gz"), chrms=config["chromosomes"], filter=config["FS_filters"]),
        expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/MQ/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_MQ_{filter}.vcf.gz"), chrms=config["chromosomes"], filter=config["MQ_filters"]),
        expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/MQRankSum/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_MQRankSum_{filter}.vcf.gz"), chrms=config["chromosomes"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/ReadPosRankSum/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_ReadPosRankSum_{filter}.vcf.gz"), chrms=config["chromosomes"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/AN/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_AN_{filter}.vcf.gz"), chrms=config["chromosomes"], filter=config["AN_filters"]),
        expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/AN/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_AN_{filter}.vcf.gz.tbi"), chrms=config["chromosomes"], filter=config["AN_filters"]),
        expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/DP/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_DP_{filter}.vcf.gz"), chrms=config["chromosomes"], filter=config["DP_filters"]),
        expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/DP/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_DP_{filter}.vcf.gz.tbi"), chrms=config["chromosomes"], filter=config["DP_filters"]),

        expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/QD/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_QD_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filter=config["QD_filters"]),
        expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/QUAL/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_QUAL_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/SOR/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_SOR_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filter=config["SOR_filters"]),
        expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/FS/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_FS_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filter=config["FS_filters"]),
        expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/MQ/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_MQ_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filter=config["MQ_filters"]),
        expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/MQRankSum/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_MQRankSum_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/ReadPosRankSum/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_ReadPosRankSum_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/AN/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_AN_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filter=config["AN_filters"]),
        expand(os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/DP/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_DP_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filter=config["DP_filters"]),

'''

#------------------------------------------------------------------------------#
# Step: Call variants - haploid mode
#------------------------------------------------------------------------------#
# Just do this for X and Y
# For this I need to specify non PARs for X...then for PARs run as diploid
rule make_gvcfs_males_haploid:
    input:
        bam = os.path.join(config["proj_path"], "bams/AFR/males/{sample}_GRCh38_YPARsMasked_sorted_mkdups.bam"),
        bai = os.path.join(config["proj_path"], "bams/AFR/males/{sample}_GRCh38_YPARsMasked_sorted_mkdups.bam.bai")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        ploidy = 1,
        ilist = config["par_intervals"]
        #chrms = "{chrms}"
    output:
        gvcf = os.path.join(config["proj_path"], "haploid/gvcfs/AFR/males/{sample}_GRCh38_YPARsMasked_haploid.g.vcf.gz")
    shell:
        """
        gatk --java-options '-Xmx4g' HaplotypeCaller -R {params.ref} -L chrX -L chrY -L chrM -XL {params.ilist} -I {input.bam} -ploidy {params.ploidy} -O {output.gvcf} -ERC GVCF
        """

# gatk --java-options '-Xmx4g' HaplotypeCaller -R {params.ref} -I {input.bam} -L {params.chrms} -ploidy {params.ploidy} -O {output.gvcf} -ERC GVCF



#------------------------------------------------------------------------------#
# Step: Combine gvcfs and joint genotype - all together, haploid
#------------------------------------------------------------------------------#
rule CombineGVCFs:
    input:
        gvcf =  expand(os.path.join(config["proj_path"], "haploid/gvcfs/AFR/males/{sample}_GRCh38_YPARsMasked_haploid.g.vcf.gz"), sample=config["males"])
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        gvcfs = expand(("-V " + config["proj_path"] + "/haploid/gvcfs/AFR/males/{sample}_GRCh38_YPARsMasked_haploid.g.vcf.gz"), sample=config["males"]),
        chrms = "{chrms}"
    output:
        os.path.join(config["proj_path"], "haploid/combine_g_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_haploid_combined.g.vcf.gz")
    shell:
        """
        gatk CombineGVCFs -R {params.ref} {params.gvcfs} -L {params.chrms} -O {output}
        """


# this sites: https://pixy.readthedocs.io/en/latest/generating_invar/generating_invar.html
# suggests that you can use a flag called -all-sites
# https://gatk.broadinstitute.org/hc/en-us/articles/4404607598875-GenotypeGVCFs
# --include-non-variant-sites / -all-sites. Include loci found to be non-variant
# after genotyping. boolean  false
rule GenotypeGVCFs:
    input:
        os.path.join(config["proj_path"], "haploid/combine_g_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_haploid_combined.g.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}"
    output:
        os.path.join(config["proj_path"], "haploid/joint_called_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_raw.vcf.gz")
    shell:
        """
        gatk GenotypeGVCFs -R {params.ref} -V {input} -L {params.chrms} -O {output}
        """

#------------------------------------------------------------------------------#
# Step: Extract SNPs
#------------------------------------------------------------------------------#
rule select_SNPs:
    input:
        os.path.join(config["proj_path"], "haploid/joint_called_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_raw.vcf.gz")
    output:
        os.path.join(config["proj_path"], "haploid/joint_called_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}"
    shell:
        """gatk --java-options "-Xmx10g" SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-O {output} """
        """--select-type-to-include SNP """


#------------------------------------------------------------------------------#
# Step: Implement hard filtering on SNPs
#------------------------------------------------------------------------------#
rule hard_filter_variant_filtration_all:
    input:
        vcf = os.path.join(config["proj_path"], "haploid/joint_called_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs.vcf.gz")
    params:
        chrms = "{chrms}",
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        qd = lambda wildcards: config[wildcards.filtering_options]["QD"],
        qual = lambda wildcards: config[wildcards.filtering_options]["QUAL"],
        sor = lambda wildcards: config[wildcards.filtering_options]["SOR"],
        fs = lambda wildcards: config[wildcards.filtering_options]["FS"],
        mq = lambda wildcards: config[wildcards.filtering_options]["MQ"],
        mqranksum = lambda wildcards: config[wildcards.filtering_options]["MQRankSum"],
        readposranksum = lambda wildcards: config[wildcards.filtering_options]["ReadPosRankSum"],
        intermediatevcf = os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/all/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs_{filtering_options}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/all/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_{filtering_options}_filtered.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "QD < {params.qd}" --filter-name "QD{params.qd}" -filter "QUAL < {params.qual}" --filter-name "QUAL{params.qual}" -filter "SOR > {params.sor}" --filter-name "SOR{params.sor}" -filter "FS > {params.fs}" --filter-name "FS{params.fs}" -filter "MQ < {params.mq}" --filter-name "MQ{params.mq}" -filter "MQRankSum < {params.mqranksum}" --filter-name "MQRankSum{params.mqranksum}" -filter "ReadPosRankSum < {params.readposranksum}" --filter-name "ReadPosRankSum{params.readposranksum}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """


#------------------------------------------------------------------------------#
# Step: For X chromosome, extract non PARS (exclude PARs)
#------------------------------------------------------------------------------#
# since I made the edit to make_gvcfs_males_haploid, there shouldnt be PARs in
# this VCF but just in case, remove them
rule extract_X_nonPARs:
    input:
         os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/all/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_{filtering_options}_filtered.vcf.gz")
    output:
         os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/all/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_{filtering_options}_filtered_nonPARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.chrms} -XL {params.ilist} """
        """-O {output} """


#------------------------------------------------------------------------------#
# Step: For Y and X chromosomes (non PARS) get per sample stats
#------------------------------------------------------------------------------#
rule bcftools_PSCs_XnonPARs:
    input:
        os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/all/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_{filtering_options}_filtered_nonPARs.vcf.gz")
    output:
        os.path.join(config["proj_path"], "haploid/bcftools_stats/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_{filtering_options}_filtered_nonPARs.bcftools.stats.PSC.txt")
    shell:
        """
        bcftools stats -s - {input} | grep PSC > {output}
        """

rule bcftools_PSCs_Y:
    input:
        os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/all/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_{filtering_options}_filtered.vcf.gz")
    output:
        os.path.join(config["proj_path"], "haploid/bcftools_stats/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_{filtering_options}_filtered.bcftools.stats.PSC.txt")
    shell:
        """
        bcftools stats -s - {input} | grep PSC > {output}
        """

'''
#------------------------------------------------------------------------------#
# Step: Repeat but for each filter separately
#------------------------------------------------------------------------------#
##############
# Filter: QD #
##############
rule hard_filter_variant_filtration_QD:
    input:
        vcf = os.path.join(config["proj_path"], "haploid/joint_called_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs.vcf.gz")
    params:
        chrms = "{chrms}",
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        qd = "{filter}",
        intermediatevcf = os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/QD/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs_QD_{filter}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/QD/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_QD_{filter}.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "QD < {params.qd}" --filter-name "QD{params.qd}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """

rule extract_X_nonPARs_QD:
    input:
         os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/QD/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_QD_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/QD/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_QD_{filter}_nonPARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.chrms} -XL {params.ilist} """
        """-O {output} """

################
# Filter: QUAL #
################
rule hard_filter_variant_filtration_QUAL:
    input:
        vcf = os.path.join(config["proj_path"], "haploid/joint_called_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs.vcf.gz")
    params:
        chrms = "{chrms}",
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        qual = "{filter}",
        intermediatevcf = os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/QUAL/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs_QUAL_{filter}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/QUAL/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_QUAL_{filter}.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "QUAL < {params.qual}" --filter-name "QUAL{params.qual}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """

rule extract_X_nonPARs_QUAL:
    input:
         os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/QUAL/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_QUAL_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/QUAL/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_QUAL_{filter}_nonPARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.chrms} -XL {params.ilist} """
        """-O {output} """

###############
# Filter: SOR #
###############
rule hard_filter_variant_filtration_SOR:
    input:
        vcf = os.path.join(config["proj_path"], "haploid/joint_called_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs.vcf.gz")
    params:
        chrms = "{chrms}",
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        sor = "{filter}",
        intermediatevcf = os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/SOR/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs_SOR_{filter}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/SOR/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_SOR_{filter}.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "SOR > {params.sor}" --filter-name "SOR{params.sor}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """

rule extract_X_nonPARs_SOR:
    input:
         os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/SOR/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_SOR_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/SOR/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_SOR_{filter}_nonPARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.chrms} -XL {params.ilist} """
        """-O {output} """

##############
# Filter: FS #
##############
rule hard_filter_variant_filtration_FS:
    input:
        vcf = os.path.join(config["proj_path"], "haploid/joint_called_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs.vcf.gz")
    params:
        chrms = "{chrms}",
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        fs = "{filter}",
        intermediatevcf = os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/FS/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs_FS_{filter}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/FS/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_FS_{filter}.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "FS > {params.fs}" --filter-name "FS{params.fs}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """

rule extract_X_nonPARs_FS:
    input:
         os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/FS/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_FS_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/FS/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_FS_{filter}_nonPARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.chrms} -XL {params.ilist} """
        """-O {output} """

##############
# Filter: MQ #
##############
rule hard_filter_variant_filtration_MQ:
    input:
        vcf = os.path.join(config["proj_path"], "haploid/joint_called_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs.vcf.gz")
    params:
        chrms = "{chrms}",
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        mq = "{filter}",
        intermediatevcf = os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/MQ/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs_MQ_{filter}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/MQ/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_MQ_{filter}.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "MQ < {params.mq}" --filter-name "MQ{params.mq}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """

rule extract_X_nonPARs_MQ:
    input:
         os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/MQ/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_MQ_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/MQ/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_MQ_{filter}_nonPARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.chrms} -XL {params.ilist} """
        """-O {output} """

#####################
# Filter: MQRankSum #
#####################
rule hard_filter_variant_filtration_MQRankSum:
    input:
        vcf = os.path.join(config["proj_path"], "haploid/joint_called_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs.vcf.gz")
    params:
        chrms = "{chrms}",
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        mqranksum = "{filter}",
        intermediatevcf = os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/MQRankSum/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs_MQRankSum_{filter}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/MQRankSum/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_MQRankSum_{filter}.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "MQRankSum < {params.mqranksum}" --filter-name "MQRankSum{params.mqranksum}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """

rule extract_X_nonPARs_MQRankSum:
    input:
         os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/MQRankSum/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_MQRankSum_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/MQRankSum/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_MQRankSum_{filter}_nonPARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.chrms} -XL {params.ilist} """
        """-O {output} """

##########################
# Filter: ReadPosRankSum #
##########################
rule hard_filter_variant_filtration_ReadPosRankSum:
    input:
        vcf = os.path.join(config["proj_path"], "haploid/joint_called_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs.vcf.gz")
    params:
        chrms = "{chrms}",
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        readposranksum = "{filter}",
        intermediatevcf = os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/ReadPosRankSum/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs_ReadPosRankSum_{filter}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/ReadPosRankSum/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_ReadPosRankSum_{filter}.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "ReadPosRankSum < {params.readposranksum}" --filter-name "ReadPosRankSum{params.readposranksum}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """

rule extract_X_nonPARs_ReadPosRankSum:
    input:
         os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/ReadPosRankSum/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_ReadPosRankSum_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/ReadPosRankSum/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_ReadPosRankSum_{filter}_nonPARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.chrms} -XL {params.ilist} """
        """-O {output} """


##############
# Filter: AN #
##############
# AN is allele number
rule hard_filter_variant_filtration_AN:
    input:
        vcf = os.path.join(config["proj_path"], "haploid/joint_called_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs.vcf.gz")
    params:
        chrms = "{chrms}",
        an = "{filter}"
    output:
        ovcf = os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/AN/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_AN_{filter}.vcf.gz"),
        ovcfidx = os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/AN/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_AN_{filter}.vcf.gz.tbi")
    shell:
        """
        bcftools filter -e'INFO/AN<{params.an}' {input.vcf} -Oz -o {output.ovcf};
        tabix -p vcf {output.ovcf}
        """

rule extract_X_nonPARs_AN:
    input:
        ovcf = os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/AN/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_AN_{filter}.vcf.gz"),
        ovcfidx = os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/AN/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_AN_{filter}.vcf.gz.tbi")
    output:
        os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/AN/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_AN_{filter}_nonPARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input.ovcf} """
        """-L {params.chrms} -XL {params.ilist} """
        """-O {output} """


##############
# Filter: DP #
##############
# DP is depth at a position for a given sample
# will replace genotype with missing if below a given DP
rule hard_filter_variant_filtration_DP:
    input:
        vcf = os.path.join(config["proj_path"], "haploid/joint_called_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs.vcf.gz")
    params:
        chrms = "{chrms}",
        an = "{filter}"
    output:
        ovcf = os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/DP/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_DP_{filter}.vcf.gz"),
        ovcfidx = os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/DP/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_DP_{filter}.vcf.gz.tbi")
    shell:
        """
        bcftools filter -e'FORMAT/DP<{params.an}' -S . {input.vcf} -Oz -o {output.ovcf};
        tabix -p vcf {output.ovcf}
        """

rule extract_X_nonPARs_DP:
    input:
        ovcf = os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/DP/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_DP_{filter}.vcf.gz"),
        ovcfidx = os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/DP/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_DP_{filter}.vcf.gz.tbi")
    output:
        os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/DP/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_DP_{filter}_nonPARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input.ovcf} """
        """-L {params.chrms} -XL {params.ilist} """
        """-O {output} """

'''

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# DEFAULT #
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Step: Call variants - haploid mode
#------------------------------------------------------------------------------#
# Just do this for X and Y
# For this I need to specify non PARs for X...then for PARs run as diploid
rule default_make_gvcfs_males_haploid:
    input:
        bam = os.path.join(config["proj_path"], "bams/AFR/males/default/{sample}_GRCh38_default_sorted_mkdups.bam"),
        bai = os.path.join(config["proj_path"], "bams/AFR/males/default/{sample}_GRCh38_default_sorted_mkdups.bam.bai")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        ploidy = 1,
        ilist = config["par_intervals"]
        #chrms = "{chrms}"
    output:
        gvcf = os.path.join(config["proj_path"], "haploid/gvcfs/AFR/males/default/{sample}_GRCh38_default_haploid.g.vcf.gz")
    shell:
        """
        gatk --java-options '-Xmx4g' HaplotypeCaller -R {params.ref} -L chrX -L chrY -L chrM -XL {params.ilist} -I {input.bam} -ploidy {params.ploidy} -O {output.gvcf} -ERC GVCF
        """

# gatk --java-options '-Xmx4g' HaplotypeCaller -R {params.ref} -I {input.bam} -L {params.chrms} -ploidy {params.ploidy} -O {output.gvcf} -ERC GVCF



#------------------------------------------------------------------------------#
# Step: Combine gvcfs and joint genotype - all together, haploid
#------------------------------------------------------------------------------#
rule default_CombineGVCFs:
    input:
        gvcf =  expand(os.path.join(config["proj_path"], "haploid/gvcfs/AFR/males/default/{sample}_GRCh38_default_haploid.g.vcf.gz"), sample=config["males"])
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        gvcfs = expand(("-V " + config["proj_path"] + "/haploid/gvcfs/AFR/males/default/{sample}_GRCh38_default_haploid.g.vcf.gz"), sample=config["males"]),
        chrms = "{chrms}"
    output:
        os.path.join(config["proj_path"], "haploid/combine_g_vcfs/AFR/males/default/{chrms}_GRCh38_default_haploid_combined.g.vcf.gz")
    shell:
        """
        gatk CombineGVCFs -R {params.ref} {params.gvcfs} -L {params.chrms} -O {output}
        """


# this sites: https://pixy.readthedocs.io/en/latest/generating_invar/generating_invar.html
# suggests that you can use a flag called -all-sites
# https://gatk.broadinstitute.org/hc/en-us/articles/4404607598875-GenotypeGVCFs
# --include-non-variant-sites / -all-sites. Include loci found to be non-variant
# after genotyping. boolean  false
rule default_GenotypeGVCFs:
    input:
        os.path.join(config["proj_path"], "haploid/combine_g_vcfs/AFR/males/default/{chrms}_GRCh38_default_haploid_combined.g.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        chrms = "{chrms}"
    output:
        os.path.join(config["proj_path"], "haploid/joint_called_vcfs/AFR/males/default/{chrms}_GRCh38_default_gatk_haploid_called_raw.vcf.gz")
    shell:
        """
        gatk GenotypeGVCFs -R {params.ref} -V {input} -L {params.chrms} -O {output}
        """

#------------------------------------------------------------------------------#
# Step: Extract SNPs
#------------------------------------------------------------------------------#
rule default_select_SNPs:
    input:
        os.path.join(config["proj_path"], "haploid/joint_called_vcfs/AFR/males/default/{chrms}_GRCh38_default_gatk_haploid_called_raw.vcf.gz")
    output:
        os.path.join(config["proj_path"], "haploid/joint_called_vcfs/AFR/males/default/{chrms}_GRCh38_default_gatk_haploid_called_raw_SNPs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        chrms = "{chrms}"
    shell:
        """gatk --java-options "-Xmx10g" SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-O {output} """
        """--select-type-to-include SNP """


#------------------------------------------------------------------------------#
# Step: Implement hard filtering on SNPs
#------------------------------------------------------------------------------#
rule default_hard_filter_variant_filtration_all:
    input:
        vcf = os.path.join(config["proj_path"], "haploid/joint_called_vcfs/AFR/males/default/{chrms}_GRCh38_default_gatk_haploid_called_raw_SNPs.vcf.gz")
    params:
        chrms = "{chrms}",
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        qd = lambda wildcards: config[wildcards.filtering_options]["QD"],
        qual = lambda wildcards: config[wildcards.filtering_options]["QUAL"],
        sor = lambda wildcards: config[wildcards.filtering_options]["SOR"],
        fs = lambda wildcards: config[wildcards.filtering_options]["FS"],
        mq = lambda wildcards: config[wildcards.filtering_options]["MQ"],
        mqranksum = lambda wildcards: config[wildcards.filtering_options]["MQRankSum"],
        readposranksum = lambda wildcards: config[wildcards.filtering_options]["ReadPosRankSum"],
        intermediatevcf = os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/default/all/{chrms}_GRCh38_default_gatk_haploid_called_raw_SNPs_{filtering_options}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/default/all/{chrms}_GRCh38_default_gatk_haploid_called_SNPs_{filtering_options}_filtered.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "QD < {params.qd}" --filter-name "QD{params.qd}" -filter "QUAL < {params.qual}" --filter-name "QUAL{params.qual}" -filter "SOR > {params.sor}" --filter-name "SOR{params.sor}" -filter "FS > {params.fs}" --filter-name "FS{params.fs}" -filter "MQ < {params.mq}" --filter-name "MQ{params.mq}" -filter "MQRankSum < {params.mqranksum}" --filter-name "MQRankSum{params.mqranksum}" -filter "ReadPosRankSum < {params.readposranksum}" --filter-name "ReadPosRankSum{params.readposranksum}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """


#------------------------------------------------------------------------------#
# Step: For X chromosome, extract non PARS (exclude PARs)
#------------------------------------------------------------------------------#
# since I made the edit to make_gvcfs_males_haploid, there shouldnt be PARs in
# this VCF but just in case, remove them
rule default_extract_X_nonPARs:
    input:
         os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/default/all/{chrms}_GRCh38_default_gatk_haploid_called_SNPs_{filtering_options}_filtered.vcf.gz")
    output:
         os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/default/all/{chrms}_GRCh38_default_gatk_haploid_called_SNPs_{filtering_options}_filtered_nonPARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.chrms} -XL {params.ilist} """
        """-O {output} """

'''
#------------------------------------------------------------------------------#
# Step: For Y and X chromosomes (non PARS) get per sample stats
#------------------------------------------------------------------------------#
rule default_bcftools_PSCs_XnonPARs:
    input:
        os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/default/all/{chrms}_GRCh38_default_gatk_haploid_called_SNPs_{filtering_options}_filtered_nonPARs.vcf.gz")
    output:
        os.path.join(config["proj_path"], "haploid/bcftools_stats/{chrms}_GRCh38_default_gatk_haploid_called_SNPs_{filtering_options}_filtered_nonPARs.bcftools.stats.PSC.txt")
    shell:
        """
        bcftools stats -s - {input} | grep PSC > {output}
        """

rule default_bcftools_PSCs_Y:
    input:
        os.path.join(config["proj_path"], "haploid/hard_filtered_vcfs/AFR/males/default/all/{chrms}_GRCh38_default_gatk_haploid_called_SNPs_{filtering_options}_filtered.vcf.gz")
    output:
        os.path.join(config["proj_path"], "haploid/bcftools_stats/{chrms}_GRCh38_default_gatk_haploid_called_SNPs_{filtering_options}_filtered.bcftools.stats.PSC.txt")
    shell:
        """
        bcftools stats -s - {input} | grep PSC > {output}
        """
'''
