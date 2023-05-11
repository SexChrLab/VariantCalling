import os

# Environment: variant_calling_simulations_project

configfile: "20230121_xy_var_calling_followups.config.json"

'''
**1) Joint genotype 10 M + 10 F with correct ploidy, SCC alignment, X non PARs
- Calc performance metrics (FP,FN,TP) (TO ADD TO SCRIPT)
**2) Joint genotype 20 M + 20 F with correct ploidy, SCC alignment, X non PARs (also do PARs for later on)
- Calc performance metrics (FP,FN,TP) (TO ADD TO SCRIPT)
3) Make a plot comparing 1 and 2 (TODO)

**4) Joint genotype 20 M + 20 F with correct ploidy, DEFAULT alignment, X non PARs (also do PARs for later on)

5) Calculate allele frequencies across X in females (non PARs and PARs)
- SCC
- DEFAULT
6) Calculate allele frequencies across X in males (non PARs and PARs)
- SCC
- DEFAULT
7) Plot allele frequencies across X
- SCC, females
- SCC, males
- DEFAULT, females
- DEFAULT, males
8) Plot difference in allele frequencies SCC-DEFAULT
- females
- males
9) Plot difference between females and males from 8
'''

wildcard_constraints:
    chrmsXnonPARs= '|'.join([re.escape(x) for x in config["chrX"]]),
    chrmsXPARs= '|'.join([re.escape(x) for x in config["chrX"]])


rule all:
    input:
        expand(os.path.join(config["proj_path"], "20230121_followups/combine_g_vcfs_10M_10F/EUR/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_combined.g.vcf.gz"), chrms=config["chrX"]),
        expand(os.path.join(config["proj_path"], "20230121_followups/joint_called_vcfs_10M_10F/EUR/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_raw.vcf.gz"), chrms=config["chrX"]),
        expand(os.path.join(config["proj_path"], "20230121_followups/joint_called_vcfs_10M_10F/EUR/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_raw_SNPs.vcf.gz"), chrms=config["chrX"]),
        expand(os.path.join(config["proj_path"], "20230121_followups/hard_filtered_vcfs_10M_10F/EUR/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_SNPs_{filtering_options}_filtered.vcf.gz"), chrms=config["chrX"], filtering_options=config["filtering_options"]),

        expand(os.path.join(config["proj_path"], "20230121_followups/combine_g_vcfs_20M_20F/EUR/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_combined.g.vcf.gz"), chrms=config["chrX"]),
        expand(os.path.join(config["proj_path"], "20230121_followups/joint_called_vcfs_20M_20F/EUR/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_raw.vcf.gz"), chrms=config["chrX"]),
        expand(os.path.join(config["proj_path"], "20230121_followups/joint_called_vcfs_20M_20F/EUR/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_raw_SNPs.vcf.gz"), chrms=config["chrX"]),
        expand(os.path.join(config["proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_SNPs_{filtering_options}_filtered.vcf.gz"), chrms=config["chrX"], filtering_options=config["filtering_options"]),

        expand(os.path.join(config["proj_path"], "20230121_followups/combine_g_vcfs_20M_20F/EUR/PARs/{chrms}_nonPARs_GRCh38_males_diploid_females_diploid_combined.g.vcf.gz"), chrms=config["chrX"]),
        expand(os.path.join(config["proj_path"], "20230121_followups/joint_called_vcfs_20M_20F/EUR/PARs/{chrms}_nonPARs_GRCh38_males_diploid_females_diploid_raw.vcf.gz"), chrms=config["chrX"]),
        expand(os.path.join(config["proj_path"], "20230121_followups/joint_called_vcfs_20M_20F/EUR/PARs/{chrms}_nonPARs_GRCh38_males_diploid_females_diploid_raw_SNPs.vcf.gz"), chrms=config["chrX"]),
        expand(os.path.join(config["proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/PARs/{chrms}_nonPARs_GRCh38_males_diploid_females_diploid_SNPs_{filtering_options}_filtered.vcf.gz"), chrms=config["chrX"], filtering_options=config["filtering_options"]),

        expand(os.path.join(config["proj_path"], "20230121_followups/combine_g_vcfs_20M_20F/EUR/default/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_combined.g.vcf.gz"), chrms=config["chrX"]),
        expand(os.path.join(config["proj_path"], "20230121_followups/joint_called_vcfs_20M_20F/EUR/default/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_raw.vcf.gz"), chrms=config["chrX"]),
        expand(os.path.join(config["proj_path"], "20230121_followups/joint_called_vcfs_20M_20F/EUR/default/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_raw_SNPs.vcf.gz"), chrms=config["chrX"]),
        expand(os.path.join(config["proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/default/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_SNPs_{filtering_options}_filtered.vcf.gz"), chrms=config["chrX"], filtering_options=config["filtering_options"]),

        expand(os.path.join(config["proj_path"], "20230121_followups/combine_g_vcfs_20M_20F/EUR/default/PARs/{chrms}_nonPARs_GRCh38_males_diploid_females_diploid_combined.g.vcf.gz"), chrms=config["chrX"]),
        expand(os.path.join(config["proj_path"], "20230121_followups/joint_called_vcfs_20M_20F/EUR/default/PARs/{chrms}_nonPARs_GRCh38_males_diploid_females_diploid_raw.vcf.gz"), chrms=config["chrX"]),
        expand(os.path.join(config["proj_path"], "20230121_followups/joint_called_vcfs_20M_20F/EUR/default/PARs/{chrms}_nonPARs_GRCh38_males_diploid_females_diploid_raw_SNPs.vcf.gz"), chrms=config["chrX"]),
        expand(os.path.join(config["proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/default/PARs/{chrms}_nonPARs_GRCh38_males_diploid_females_diploid_SNPs_{filtering_options}_filtered.vcf.gz"), chrms=config["chrX"], filtering_options=config["filtering_options"])


#------------------------------------------------------------------------------#
# 1. Joint genotype 10 M + 10 F with correct ploidy, SCC alignment, X non PARs #
#------------------------------------------------------------------------------#
rule CombineGVCFs_10M_10F_XnonPARs_SCC:
    input:
        males = expand(os.path.join(config["proj_path"], "haploid/gvcfs/EUR/males/{sample}_GRCh38_YPARsMasked_haploid.g.vcf.gz"), sample=config["males_first10"], chrms=config["chrX"]),
        females = expand(os.path.join(config["proj_path"], "gvcfs/EUR/females/{sample}_{chrms}_GRCh38_YHardMasked_diploid.g.vcf.gz"), sample=config["females_first10"], chrms=config["chrX"])
    params:
        # since calling variants I wasnt sure which ref to use for this step 
        ref = config["genome_paths"]["Ref_GRCh38_Default"], # This shouldn't matter because I already did the alignment with the appropriate ref
        gvcfsM = expand(("-V " + config["proj_path"] + "haploid/gvcfs/EUR/males/{sample}_GRCh38_YPARsMasked_haploid.g.vcf.gz"), sample=config["males_first10"], chrms=config["chrX"]),
        gvcfsF = expand(("-V " + config["proj_path"] + "gvcfs/EUR/females/{sample}_{chrms}_GRCh38_YHardMasked_diploid.g.vcf.gz"), sample=config["females_first10"], chrms=config["chrX"]),
        chrms = "{chrms}"
    output:
        os.path.join(config["proj_path"], "20230121_followups/combine_g_vcfs_10M_10F/EUR/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_combined.g.vcf.gz")
    shell:
        """
        gatk CombineGVCFs -R {params.ref} {params.gvcfsM} {params.gvcfsF} -L {params.chrms} -XL {params.chrms}:10001-2781479 -XL {params.chrms}:155701383-156030895 -O {output}
        """


rule GenotypeGVCFs_10M_10F_XnonPARs_SCC:
    input:
        os.path.join(config["proj_path"], "20230121_followups/combine_g_vcfs_10M_10F/EUR/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_combined.g.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        chrms = "{chrms}"
    output:
        os.path.join(config["proj_path"], "20230121_followups/joint_called_vcfs_10M_10F/EUR/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_raw.vcf.gz")
    shell:
        """
        gatk GenotypeGVCFs -R {params.ref} -V {input} -L {params.chrms} -XL {params.chrms}:10001-2781479 -XL {params.chrms}:155701383-156030895 -O {output}
        """


rule select_SNPs_10M_10F_XnonPARs_SCC:
    input:
        os.path.join(config["proj_path"], "20230121_followups/joint_called_vcfs_10M_10F/EUR/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_raw.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        chrms = "{chrms}"
    output:
        os.path.join(config["proj_path"], "20230121_followups/joint_called_vcfs_10M_10F/EUR/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_raw_SNPs.vcf.gz")
    shell:
        """
        gatk --java-options "-Xmx10g" SelectVariants -R {params.ref} -V {input} -O {output} --select-type-to-include SNP
        """


rule hard_filter_10M_10F_XnonPARs_SCC:
    input:
        vcf = os.path.join(config["proj_path"], "20230121_followups/joint_called_vcfs_10M_10F/EUR/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_raw_SNPs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        chrms = "{chrms}",
        qd = lambda wildcards: config[wildcards.filtering_options]["QD"],
        qual = lambda wildcards: config[wildcards.filtering_options]["QUAL"],
        sor = lambda wildcards: config[wildcards.filtering_options]["SOR"],
        fs = lambda wildcards: config[wildcards.filtering_options]["FS"],
        mq = lambda wildcards: config[wildcards.filtering_options]["MQ"],
        mqranksum = lambda wildcards: config[wildcards.filtering_options]["MQRankSum"],
        readposranksum = lambda wildcards: config[wildcards.filtering_options]["ReadPosRankSum"],
        intermediatevcf = os.path.join(config["proj_path"], "20230121_followups/hard_filtered_vcfs_10M_10F/EUR/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_raw_SNPs_{filtering_options}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "20230121_followups/hard_filtered_vcfs_10M_10F/EUR/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_SNPs_{filtering_options}_filtered.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "QD < {params.qd}" --filter-name "QD{params.qd}" -filter "QUAL < {params.qual}" --filter-name "QUAL{params.qual}" -filter "SOR > {params.sor}" --filter-name "SOR{params.sor}" -filter "FS > {params.fs}" --filter-name "FS{params.fs}" -filter "MQ < {params.mq}" --filter-name "MQ{params.mq}" -filter "MQRankSum < {params.mqranksum}" --filter-name "MQRankSum{params.mqranksum}" -filter "ReadPosRankSum < {params.readposranksum}" --filter-name "ReadPosRankSum{params.readposranksum}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """

# TO DOs HERE - performance metrics
# General steps:
# Separate out each sample into their own vcf
# Run python script compare_VCFs_fix.py on the per sample vcf along with its 
#   corresponding golden vcf, set ploidy to the ploidy of the golden vcf

# might need to separate this rule by males and females too since I might get a 
#   wildcard error in the next rule since it is separated by m and f
rule extract_samples_VCF_XnonPARs_SCC:
    input:
        os.path.join(config["proj_path"], "20230121_followups/hard_filtered_vcfs_10M_10F/EUR/XnonPARs/{chrmsXnonPARs}_nonPARs_GRCh38_males_haploid_females_diploid_SNPs_gatkHardFilter_filtered.vcf.gz")
        #os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs.vcf.gz")
    params:
        smpl = "{sample}", # samples in config has been changed to the 10m and 10f
        vcf = os.path.join(config["scratch_proj_path"], "20230121_followups/hard_filtered_vcfs_10M_10F/EUR/XnonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs")
    output:
        os.path.join(config["scratch_proj_path"], "20230121_followups/hard_filtered_vcfs_10M_10F/EUR/XnonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

# have to run a separate rule for males and females since naming of golden file is different
rule get_metrics_males_XnonPARs_SCC:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "20230121_followups/hard_filtered_vcfs_10M_10F/EUR/XnonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf") # sample will be males_first10
        #called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrms}_nonPARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "20230121_followups/compare_VCFs/EUR/males/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called") 
    output:
        os.path.join(config["scratch_proj_path"], "20230121_followups/compare_VCFs/EUR/males/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule get_metrics_females_XnonPARs_SCC:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "20230121_followups/hard_filtered_vcfs_10M_10F/EUR/XnonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf") # sample will be females_first10
        #called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrms}_nonPARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "20230121_followups/compare_VCFs/EUR/females/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called") 
    output:
        os.path.join(config["scratch_proj_path"], "20230121_followups/compare_VCFs/EUR/females/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# 2. Joint genotype 20 M + 20 F with correct ploidy, SCC alignment, X non PARs #
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
rule CombineGVCFs_20M_20F_XnonPARs_SCC:
    input:
        males = expand(os.path.join(config["proj_path"], "haploid/gvcfs/EUR/males/{sample}_GRCh38_YPARsMasked_haploid.g.vcf.gz"), sample=config["males"], chrms=config["chrX"]), 
        females = expand(os.path.join(config["proj_path"], "gvcfs/EUR/females/{sample}_{chrms}_GRCh38_YHardMasked_diploid.g.vcf.gz"), sample=config["females"], chrms=config["chrX"])
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        gvcfsM = expand(("-V " + config["proj_path"] + "haploid/gvcfs/EUR/males/{sample}_GRCh38_YPARsMasked_haploid.g.vcf.gz"), sample=config["males"], chrms=config["chrX"]),
        gvcfsF = expand(("-V " + config["proj_path"] + "gvcfs/EUR/females/{sample}_{chrms}_GRCh38_YHardMasked_diploid.g.vcf.gz"), sample=config["females"], chrms=config["chrX"]),
        chrms = "{chrms}"
    output:
        os.path.join(config["proj_path"], "20230121_followups/combine_g_vcfs_20M_20F/EUR/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_combined.g.vcf.gz")
    shell:
        """
        gatk CombineGVCFs -R {params.ref} {params.gvcfsM} {params.gvcfsF} -L {params.chrms} -XL {params.chrms}:10001-2781479 -XL {params.chrms}:155701383-156030895 -O {output}
        """


rule GenotypeGVCFs_20M_20F_XnonPARs_SCC:
    input:
        os.path.join(config["proj_path"], "20230121_followups/combine_g_vcfs_20M_20F/EUR/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_combined.g.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        chrms = "{chrms}"
    output:
        os.path.join(config["proj_path"], "20230121_followups/joint_called_vcfs_20M_20F/EUR/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_raw.vcf.gz")
    shell:
        """
        gatk GenotypeGVCFs -R {params.ref} -V {input} -L {params.chrms} -XL {params.chrms}:10001-2781479 -XL {params.chrms}:155701383-156030895 -O {output}
        """


rule select_SNPs_20M_20F_XnonPARs_SCC:
    input:
        os.path.join(config["proj_path"], "20230121_followups/joint_called_vcfs_20M_20F/EUR/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_raw.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        chrms = "{chrms}"
    output:
        os.path.join(config["proj_path"], "20230121_followups/joint_called_vcfs_20M_20F/EUR/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_raw_SNPs.vcf.gz")
    shell:
        """
        gatk --java-options "-Xmx10g" SelectVariants -R {params.ref} -V {input} -O {output} --select-type-to-include SNP
        """


rule hard_filter_20M_20F_XnonPARs_SCC:
    input:
        vcf = os.path.join(config["proj_path"], "20230121_followups/joint_called_vcfs_20M_20F/EUR/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_raw_SNPs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        chrms = "{chrms}",
        qd = lambda wildcards: config[wildcards.filtering_options]["QD"],
        qual = lambda wildcards: config[wildcards.filtering_options]["QUAL"],
        sor = lambda wildcards: config[wildcards.filtering_options]["SOR"],
        fs = lambda wildcards: config[wildcards.filtering_options]["FS"],
        mq = lambda wildcards: config[wildcards.filtering_options]["MQ"],
        mqranksum = lambda wildcards: config[wildcards.filtering_options]["MQRankSum"],
        readposranksum = lambda wildcards: config[wildcards.filtering_options]["ReadPosRankSum"],
        intermediatevcf = os.path.join(config["proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_raw_SNPs_{filtering_options}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_SNPs_{filtering_options}_filtered.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "QD < {params.qd}" --filter-name "QD{params.qd}" -filter "QUAL < {params.qual}" --filter-name "QUAL{params.qual}" -filter "SOR > {params.sor}" --filter-name "SOR{params.sor}" -filter "FS > {params.fs}" --filter-name "FS{params.fs}" -filter "MQ < {params.mq}" --filter-name "MQ{params.mq}" -filter "MQRankSum < {params.mqranksum}" --filter-name "MQRankSum{params.mqranksum}" -filter "ReadPosRankSum < {params.readposranksum}" --filter-name "ReadPosRankSum{params.readposranksum}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """

# TO DOs HERE - allele frequencies
# vcftools --gzvcf /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/hard_filtered_vcfs_20M_20F/EUR/XnonPARs/chrX_nonPARs_GRCh38_males_haploid_females_diploid_SNPs_gatkHardFilter_filtered.vcf.gz --min-alleles 2 --max-alleles 2 --keep /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/females.txt --freq --out /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/allele_frq_vcfs_20M_20F/EUR/XnonPARs/females/females_SCC_XnonPARs_frq
# vcftools --gzvcf /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/hard_filtered_vcfs_20M_20F/EUR/XnonPARs/chrX_nonPARs_GRCh38_males_haploid_females_diploid_SNPs_gatkHardFilter_filtered.vcf.gz --min-alleles 2 --max-alleles 2 --keep /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/males.txt --freq --out /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/allele_frq_vcfs_20M_20F/EUR/XnonPARs/males/males_SCC_XnonPARs_frq

# grep -v "CHROM" /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/allele_frq_vcfs_20M_20F/EUR/XnonPARs/males/males_SCC_XnonPARs_frq.frq > /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/allele_frq_vcfs_20M_20F/EUR/XnonPARs/males/males_SCC_XnonPARs_frq_fix.frq
# grep -v "CHROM" /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/allele_frq_vcfs_20M_20F/EUR/XnonPARs/females/females_SCC_XnonPARs_frq.frq > /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/allele_frq_vcfs_20M_20F/EUR/XnonPARs/females/females_SCC_XnonPARs_frq_fix.frq


#####################################
# TO DOs HERE - performance metrics #
#####################################
rule extract_samples_VCF_XnonPARs_SCC_20MF:
    input:
        os.path.join(config["proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/XnonPARs/{chrmsXnonPARs}_nonPARs_GRCh38_males_haploid_females_diploid_SNPs_gatkHardFilter_filtered.vcf.gz")
        #os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs.vcf.gz")
    params:
        smpl = "{sample}", # samples_20 in config has been changed to the 20m and 20f
        vcf = os.path.join(config["scratch_proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/XnonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs")
    output:
        os.path.join(config["scratch_proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/XnonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

# have to run a separate rule for males and females since naming of golden file is different
rule get_metrics_males_XnonPARs_SCC_20MF:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/XnonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf") # sample will be males_first10
        #called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrms}_nonPARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "20230121_followups/compare_VCFs/EUR/males/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called") 
    output:
        os.path.join(config["scratch_proj_path"], "20230121_followups/compare_VCFs/EUR/males/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule get_metrics_females_XnonPARs_SCC_20MF:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/XnonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf") # sample will be females_first10
        #called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrms}_nonPARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "20230121_followups/compare_VCFs/EUR/females/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called") 
    output:
        os.path.join(config["scratch_proj_path"], "20230121_followups/compare_VCFs/EUR/females/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# 3. Joint genotype 20 M + 20 F with correct ploidy, SCC alignment, X PARs #
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
rule CombineGVCFs_20M_20F_PARs_SCC:
    input:
        males = expand(os.path.join(config["proj_path"], "gvcfs/EUR/males/{sample}_{chrms}_GRCh38_YPARsMasked_diploid.g.vcf.gz"), sample=config["males"], chrms=config["chrX"]), # I forgot to change file name but these are PAR variants
        females = expand(os.path.join(config["proj_path"], "gvcfs/EUR/females/{sample}_{chrms}_GRCh38_YHardMasked_diploid.g.vcf.gz"), sample=config["females"], chrms=config["chrX"])
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        gvcfsM = expand(("-V " + config["proj_path"] + "gvcfs/EUR/males/{sample}_{chrms}_GRCh38_YPARsMasked_diploid.g.vcf.gz"), sample=config["males"], chrms=config["chrX"]),
        gvcfsF = expand(("-V " + config["proj_path"] + "gvcfs/EUR/females/{sample}_{chrms}_GRCh38_YHardMasked_diploid.g.vcf.gz"), sample=config["females"], chrms=config["chrX"]),
        chrms = "{chrms}"
    output:
        os.path.join(config["proj_path"], "20230121_followups/combine_g_vcfs_20M_20F/EUR/PARs/{chrms}_nonPARs_GRCh38_males_diploid_females_diploid_combined.g.vcf.gz")
    shell:
        """
        gatk CombineGVCFs -R {params.ref} {params.gvcfsM} {params.gvcfsF} -L {params.chrms}:10001-2781479 -L {params.chrms}:155701383-156030895 -O {output}
        """


rule GenotypeGVCFs_20M_20F_PARs_SCC:
    input:
        os.path.join(config["proj_path"], "20230121_followups/combine_g_vcfs_20M_20F/EUR/PARs/{chrms}_nonPARs_GRCh38_males_diploid_females_diploid_combined.g.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        chrms = "{chrms}"
    output:
        os.path.join(config["proj_path"], "20230121_followups/joint_called_vcfs_20M_20F/EUR/PARs/{chrms}_nonPARs_GRCh38_males_diploid_females_diploid_raw.vcf.gz")
    shell:
        """
        gatk GenotypeGVCFs -R {params.ref} -V {input} -L {params.chrms}:10001-2781479 -L {params.chrms}:155701383-156030895 -O {output}
        """


rule select_SNPs_20M_20F_PARs_SCC:
    input:
        os.path.join(config["proj_path"], "20230121_followups/joint_called_vcfs_20M_20F/EUR/PARs/{chrms}_nonPARs_GRCh38_males_diploid_females_diploid_raw.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        chrms = "{chrms}"
    output:
        os.path.join(config["proj_path"], "20230121_followups/joint_called_vcfs_20M_20F/EUR/PARs/{chrms}_nonPARs_GRCh38_males_diploid_females_diploid_raw_SNPs.vcf.gz")
    shell:
        """
        gatk --java-options "-Xmx10g" SelectVariants -R {params.ref} -V {input} -O {output} --select-type-to-include SNP
        """


rule hard_filter_20M_20F_PARs_SCC:
    input:
        vcf = os.path.join(config["proj_path"], "20230121_followups/joint_called_vcfs_20M_20F/EUR/PARs/{chrms}_nonPARs_GRCh38_males_diploid_females_diploid_raw_SNPs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        chrms = "{chrms}",
        qd = lambda wildcards: config[wildcards.filtering_options]["QD"],
        qual = lambda wildcards: config[wildcards.filtering_options]["QUAL"],
        sor = lambda wildcards: config[wildcards.filtering_options]["SOR"],
        fs = lambda wildcards: config[wildcards.filtering_options]["FS"],
        mq = lambda wildcards: config[wildcards.filtering_options]["MQ"],
        mqranksum = lambda wildcards: config[wildcards.filtering_options]["MQRankSum"],
        readposranksum = lambda wildcards: config[wildcards.filtering_options]["ReadPosRankSum"],
        intermediatevcf = os.path.join(config["proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/PARs/{chrms}_nonPARs_GRCh38_males_diploid_females_diploid_raw_SNPs_{filtering_options}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/PARs/{chrms}_nonPARs_GRCh38_males_diploid_females_diploid_SNPs_{filtering_options}_filtered.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "QD < {params.qd}" --filter-name "QD{params.qd}" -filter "QUAL < {params.qual}" --filter-name "QUAL{params.qual}" -filter "SOR > {params.sor}" --filter-name "SOR{params.sor}" -filter "FS > {params.fs}" --filter-name "FS{params.fs}" -filter "MQ < {params.mq}" --filter-name "MQ{params.mq}" -filter "MQRankSum < {params.mqranksum}" --filter-name "MQRankSum{params.mqranksum}" -filter "ReadPosRankSum < {params.readposranksum}" --filter-name "ReadPosRankSum{params.readposranksum}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """

# TO DOs HERE - allele frequencies
# vcftools --gzvcf /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/hard_filtered_vcfs_20M_20F/EUR/PARs/chrX_nonPARs_GRCh38_males_diploid_females_diploid_SNPs_gatkHardFilter_filtered.vcf.gz --min-alleles 2 --max-alleles 2 --keep /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/females.txt --freq --out /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/allele_frq_vcfs_20M_20F/EUR/PARs/females/females_SCC_PARs_frq
# vcftools --gzvcf /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/hard_filtered_vcfs_20M_20F/EUR/PARs/chrX_nonPARs_GRCh38_males_diploid_females_diploid_SNPs_gatkHardFilter_filtered.vcf.gz --min-alleles 2 --max-alleles 2 --keep /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/males.txt --freq --out /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/allele_frq_vcfs_20M_20F/EUR/PARs/males/males_SCC_PARs_frq

# grep -v "CHROM" /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/allele_frq_vcfs_20M_20F/EUR/PARs/males/males_SCC_PARs_frq.frq > /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/allele_frq_vcfs_20M_20F/EUR/PARs/males/males_SCC_PARs_frq_fix.frq
# grep -v "CHROM" /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/allele_frq_vcfs_20M_20F/EUR/PARs/females/females_SCC_PARs_frq.frq > /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/allele_frq_vcfs_20M_20F/EUR/PARs/females/females_SCC_PARs_frq_fix.frq


#####################################
# TO DOs HERE - performance metrics
#####################################
rule extract_samples_VCF_PARs_SCC_20MF:
    input:
        os.path.join(config["proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/PARs/{chrmsXPARs}_nonPARs_GRCh38_males_diploid_females_diploid_SNPs_gatkHardFilter_filtered.vcf.gz")
    params:
        smpl = "{sample}", # samples_20 in config has been changed to the 20m and 20f
        vcf = os.path.join(config["scratch_proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/PARs/by_sample/{sample}_{chrmsXPARs}_PARs")
    output:
        os.path.join(config["scratch_proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/PARs/by_sample/{sample}_{chrmsXPARs}_PARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

# have to run a separate rule for males and females since naming of golden file is different
rule get_metrics_males_PARs_SCC_20MF:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/PARs/by_sample/{sample}_{chrmsXPARs}_PARs.recode.vcf") # sample will be males_first10
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "20230121_followups/compare_VCFs/EUR/males/PARs/{sample}_{chrms}_PARs_golden_vs_called") 
    output:
        os.path.join(config["scratch_proj_path"], "20230121_followups/compare_VCFs/EUR/males/PARs/{sample}_{chrms}_PARs_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule get_metrics_females_PARs_SCC_20MF:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/PARs/by_sample/{sample}_{chrmsXPARs}_PARs.recode.vcf") # sample will be females_first10
        #called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrms}_PARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "20230121_followups/compare_VCFs/EUR/females/PARs/{sample}_{chrms}_PARs_golden_vs_called") 
    output:
        os.path.join(config["scratch_proj_path"], "20230121_followups/compare_VCFs/EUR/females/PARs/{sample}_{chrms}_PARs_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """



#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# 4. Joint genotype 20 M + 20 F with correct ploidy, DEFAULT alignment, X non  #
#    PARs                                                                      #
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
rule CombineGVCFs_20M_20F_XnonPARs_DEFAULT:
    input:
        males = expand(os.path.join(config["proj_path"], "haploid/gvcfs/EUR/males/default/{sample}_GRCh38_default_haploid.g.vcf.gz"), sample=config["males"], chrms=config["chrX"]),
        females = expand(os.path.join(config["proj_path"], "gvcfs/EUR/females/default/{sample}_{chrms}_GRCh38_default_diploid.g.vcf.gz"), sample=config["females"], chrms=config["chrX"])
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        gvcfsM = expand(("-V " + config["proj_path"] + "haploid/gvcfs/EUR/males/default/{sample}_GRCh38_default_haploid.g.vcf.gz"), sample=config["males"], chrms=config["chrX"]),
        gvcfsF = expand(("-V " + config["proj_path"] + "gvcfs/EUR/females/default/{sample}_{chrms}_GRCh38_default_diploid.g.vcf.gz"), sample=config["females"], chrms=config["chrX"]),
        chrms = "{chrms}"
    output:
        os.path.join(config["proj_path"], "20230121_followups/combine_g_vcfs_20M_20F/EUR/default/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_combined.g.vcf.gz")
    shell:
        """
        gatk CombineGVCFs -R {params.ref} {params.gvcfsM} {params.gvcfsF} -L {params.chrms} -XL {params.chrms}:10001-2781479 -XL {params.chrms}:155701383-156030895 -O {output}
        """


rule GenotypeGVCFs_20M_20F_XnonPARs_DEFAULT:
    input:
        os.path.join(config["proj_path"], "20230121_followups/combine_g_vcfs_20M_20F/EUR/default/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_combined.g.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        chrms = "{chrms}"
    output:
        os.path.join(config["proj_path"], "20230121_followups/joint_called_vcfs_20M_20F/EUR/default/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_raw.vcf.gz")
    shell:
        """
        gatk GenotypeGVCFs -R {params.ref} -V {input} -L {params.chrms} -XL {params.chrms}:10001-2781479 -XL {params.chrms}:155701383-156030895 -O {output}
        """


rule select_SNPs_20M_20F_XnonPARs_DEFAULT:
    input:
        os.path.join(config["proj_path"], "20230121_followups/joint_called_vcfs_20M_20F/EUR/default/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_raw.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        chrms = "{chrms}"
    output:
        os.path.join(config["proj_path"], "20230121_followups/joint_called_vcfs_20M_20F/EUR/default/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_raw_SNPs.vcf.gz")
    shell:
        """
        gatk --java-options "-Xmx10g" SelectVariants -R {params.ref} -V {input} -O {output} --select-type-to-include SNP
        """


rule hard_filter_20M_20F_XnonPARs_DEFAULT:
    input:
        vcf = os.path.join(config["proj_path"], "20230121_followups/joint_called_vcfs_20M_20F/EUR/default/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_raw_SNPs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        chrms = "{chrms}",
        qd = lambda wildcards: config[wildcards.filtering_options]["QD"],
        qual = lambda wildcards: config[wildcards.filtering_options]["QUAL"],
        sor = lambda wildcards: config[wildcards.filtering_options]["SOR"],
        fs = lambda wildcards: config[wildcards.filtering_options]["FS"],
        mq = lambda wildcards: config[wildcards.filtering_options]["MQ"],
        mqranksum = lambda wildcards: config[wildcards.filtering_options]["MQRankSum"],
        readposranksum = lambda wildcards: config[wildcards.filtering_options]["ReadPosRankSum"],
        intermediatevcf = os.path.join(config["proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/default/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_raw_SNPs_{filtering_options}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/default/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_SNPs_{filtering_options}_filtered.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "QD < {params.qd}" --filter-name "QD{params.qd}" -filter "QUAL < {params.qual}" --filter-name "QUAL{params.qual}" -filter "SOR > {params.sor}" --filter-name "SOR{params.sor}" -filter "FS > {params.fs}" --filter-name "FS{params.fs}" -filter "MQ < {params.mq}" --filter-name "MQ{params.mq}" -filter "MQRankSum < {params.mqranksum}" --filter-name "MQRankSum{params.mqranksum}" -filter "ReadPosRankSum < {params.readposranksum}" --filter-name "ReadPosRankSum{params.readposranksum}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """

# TO DOs HERE - allele frequencies
# vcftools --gzvcf /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/hard_filtered_vcfs_20M_20F/EUR/default/XnonPARs/chrX_nonPARs_GRCh38_males_haploid_females_diploid_SNPs_gatkHardFilter_filtered.vcf.gz --min-alleles 2 --max-alleles 2 --keep /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/females.txt --freq --out /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/allele_frq_vcfs_20M_20F/EUR/default/XnonPARs/females/females_default_XnonPARs_frq
# vcftools --gzvcf /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/hard_filtered_vcfs_20M_20F/EUR/default/XnonPARs/chrX_nonPARs_GRCh38_males_haploid_females_diploid_SNPs_gatkHardFilter_filtered.vcf.gz --min-alleles 2 --max-alleles 2 --keep /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/males.txt --freq --out /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/allele_frq_vcfs_20M_20F/EUR/default/XnonPARs/males/males_default_XnonPARs_frq

# grep -v "CHROM" /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/allele_frq_vcfs_20M_20F/EUR/default/XnonPARs/males/males_default_XnonPARs_frq.frq > /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/allele_frq_vcfs_20M_20F/EUR/default/XnonPARs/males/males_default_XnonPARs_frq_fix.frq
# grep -v "CHROM" /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/allele_frq_vcfs_20M_20F/EUR/default/XnonPARs/females/females_default_XnonPARs_frq.frq > /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/allele_frq_vcfs_20M_20F/EUR/default/XnonPARs/females/females_default_XnonPARs_frq_fix.frq


#####################################
# TO DOs HERE - performance metrics #
#####################################
rule extract_samples_VCF_XnonPARs_DEFAULT_20MF:
    input:
        os.path.join(config["proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/default/XnonPARs/{chrms}_nonPARs_GRCh38_males_haploid_females_diploid_SNPs_gatkHardFilter_filtered.vcf.gz")
    params:
        smpl = "{sample}", # samples_20 in config has been changed to the 20m and 20f
        vcf = os.path.join(config["scratch_proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/default/XnonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs")
    output:
        os.path.join(config["scratch_proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/default/XnonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

# have to run a separate rule for males and females since naming of golden file is different
rule get_metrics_males_XnonPARs_DEFAULT_20MF:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/default/XnonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf") # sample will be males_first10
        #called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrms}_nonPARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "20230121_followups/compare_VCFs/EUR/default/males/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called") 
    output:
        os.path.join(config["scratch_proj_path"], "20230121_followups/compare_VCFs/EUR/default/males/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule get_metrics_females_XnonPARs_DEFAULT_20MF:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/default/XnonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf") # sample will be females_first10
        #called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrms}_nonPARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "20230121_followups/compare_VCFs/EUR/default/females/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called") 
    output:
        os.path.join(config["scratch_proj_path"], "20230121_followups/compare_VCFs/EUR/default/females/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """



#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# 5. Joint genotype 20 M + 20 F with correct ploidy, DEFAULT alignment, X PARs #
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
rule CombineGVCFs_20M_20F_PARs_DEFAULT:
    input:
        males = expand(os.path.join(config["proj_path"], "gvcfs/EUR/males/default/{sample}_{chrms}_GRCh38_default_diploid.g.vcf.gz"), sample=config["males"], chrms=config["chrX"]),
        females = expand(os.path.join(config["proj_path"], "gvcfs/EUR/females/default/{sample}_{chrms}_GRCh38_default_diploid.g.vcf.gz"), sample=config["females"], chrms=config["chrX"])
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        gvcfsM = expand(("-V " + config["proj_path"] + "gvcfs/EUR/males/default/{sample}_{chrms}_GRCh38_default_diploid.g.vcf.gz"), sample=config["males"], chrms=config["chrX"]),
        gvcfsF = expand(("-V " + config["proj_path"] + "gvcfs/EUR/females/default/{sample}_{chrms}_GRCh38_default_diploid.g.vcf.gz"), sample=config["females"], chrms=config["chrX"]),
        chrms = "{chrms}"
    output:
        os.path.join(config["proj_path"], "20230121_followups/combine_g_vcfs_20M_20F/EUR/default/PARs/{chrms}_nonPARs_GRCh38_males_diploid_females_diploid_combined.g.vcf.gz")
    shell:
        """
        gatk CombineGVCFs -R {params.ref} {params.gvcfsM} {params.gvcfsF} -L {params.chrms}:10001-2781479 -L {params.chrms}:155701383-156030895 -O {output}
        """


rule GenotypeGVCFs_20M_20F_PARs_DEFAULT:
    input:
        os.path.join(config["proj_path"], "20230121_followups/combine_g_vcfs_20M_20F/EUR/default/PARs/{chrms}_nonPARs_GRCh38_males_diploid_females_diploid_combined.g.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        chrms = "{chrms}"
    output:
        os.path.join(config["proj_path"], "20230121_followups/joint_called_vcfs_20M_20F/EUR/default/PARs/{chrms}_nonPARs_GRCh38_males_diploid_females_diploid_raw.vcf.gz")
    shell:
        """
        gatk GenotypeGVCFs -R {params.ref} -V {input} -L {params.chrms}:10001-2781479 -L {params.chrms}:155701383-156030895 -O {output}
        """


rule select_SNPs_20M_20F_PARs_DEFAULT:
    input:
        os.path.join(config["proj_path"], "20230121_followups/joint_called_vcfs_20M_20F/EUR/default/PARs/{chrms}_nonPARs_GRCh38_males_diploid_females_diploid_raw.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        chrms = "{chrms}"
    output:
        os.path.join(config["proj_path"], "20230121_followups/joint_called_vcfs_20M_20F/EUR/default/PARs/{chrms}_nonPARs_GRCh38_males_diploid_females_diploid_raw_SNPs.vcf.gz")
    shell:
        """
        gatk --java-options "-Xmx10g" SelectVariants -R {params.ref} -V {input} -O {output} --select-type-to-include SNP
        """


rule hard_filter_20M_20F_PARs_DEFAULT:
    input:
        vcf = os.path.join(config["proj_path"], "20230121_followups/joint_called_vcfs_20M_20F/EUR/default/PARs/{chrms}_nonPARs_GRCh38_males_diploid_females_diploid_raw_SNPs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        chrms = "{chrms}",
        qd = lambda wildcards: config[wildcards.filtering_options]["QD"],
        qual = lambda wildcards: config[wildcards.filtering_options]["QUAL"],
        sor = lambda wildcards: config[wildcards.filtering_options]["SOR"],
        fs = lambda wildcards: config[wildcards.filtering_options]["FS"],
        mq = lambda wildcards: config[wildcards.filtering_options]["MQ"],
        mqranksum = lambda wildcards: config[wildcards.filtering_options]["MQRankSum"],
        readposranksum = lambda wildcards: config[wildcards.filtering_options]["ReadPosRankSum"],
        intermediatevcf = os.path.join(config["proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/default/PARs/{chrms}_nonPARs_GRCh38_males_diploid_females_diploid_raw_SNPs_{filtering_options}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/default/PARs/{chrms}_nonPARs_GRCh38_males_diploid_females_diploid_SNPs_{filtering_options}_filtered.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "QD < {params.qd}" --filter-name "QD{params.qd}" -filter "QUAL < {params.qual}" --filter-name "QUAL{params.qual}" -filter "SOR > {params.sor}" --filter-name "SOR{params.sor}" -filter "FS > {params.fs}" --filter-name "FS{params.fs}" -filter "MQ < {params.mq}" --filter-name "MQ{params.mq}" -filter "MQRankSum < {params.mqranksum}" --filter-name "MQRankSum{params.mqranksum}" -filter "ReadPosRankSum < {params.readposranksum}" --filter-name "ReadPosRankSum{params.readposranksum}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """

# TO DOs HERE - allele frequencies
# vcftools --gzvcf /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/hard_filtered_vcfs_20M_20F/EUR/default/PARs/chrX_nonPARs_GRCh38_males_diploid_females_diploid_SNPs_gatkHardFilter_filtered.vcf.gz --min-alleles 2 --max-alleles 2 --keep /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/females.txt --freq --out /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/allele_frq_vcfs_20M_20F/EUR/default/PARs/females/females_default_PARs_frq
# vcftools --gzvcf /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/hard_filtered_vcfs_20M_20F/EUR/default/PARs/chrX_nonPARs_GRCh38_males_diploid_females_diploid_SNPs_gatkHardFilter_filtered.vcf.gz --min-alleles 2 --max-alleles 2 --keep /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/males.txt --freq --out /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/allele_frq_vcfs_20M_20F/EUR/default/PARs/males/males_default_PARs_frq

# grep -v "CHROM" /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/allele_frq_vcfs_20M_20F/EUR/default/PARs/males/males_default_PARs_frq.frq > /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/allele_frq_vcfs_20M_20F/EUR/default/PARs/males/males_default_PARs_frq_fix.frq
# grep -v "CHROM" /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/allele_frq_vcfs_20M_20F/EUR/default/PARs/females/females_default_PARs_frq.frq > /data/CEM/wilsonlab/projects/variant_calling_simulations_project/20230121_followups/allele_frq_vcfs_20M_20F/EUR/default/PARs/females/females_default_PARs_frq_fix.frq


#####################################
# TO DOs HERE - performance metrics
#####################################
rule extract_samples_VCF_PARs_DEFAULT_20MF:
    input:
        os.path.join(config["proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/default/PARs/{chrms}_nonPARs_GRCh38_males_diploid_females_diploid_SNPs_gatkHardFilter_filtered.vcf.gz")
    params:
        smpl = "{sample}", # samples_20 in config has been changed to the 20m and 20f
        vcf = os.path.join(config["scratch_proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/default/PARs/by_sample/{sample}_{chrmsXPARs}_PARs")
    output:
        os.path.join(config["scratch_proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/default/PARs/by_sample/{sample}_{chrmsXPARs}_PARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

# have to run a separate rule for males and females since naming of golden file is different
rule get_metrics_males_PARs_DEFAULT_20MF:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/default/PARs/by_sample/{sample}_{chrmsXPARs}_PARs.recode.vcf") # sample will be males_first10
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "20230121_followups/compare_VCFs/EUR/default/males/PARs/{sample}_{chrms}_PARs_golden_vs_called") 
    output:
        os.path.join(config["scratch_proj_path"], "20230121_followups/compare_VCFs/EUR/default/males/PARs/{sample}_{chrms}_PARs_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule get_metrics_females_PARs_DEFAULT_20MF:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "20230121_followups/hard_filtered_vcfs_20M_20F/EUR/default/PARs/by_sample/{sample}_{chrmsXPARs}_PARs.recode.vcf") # sample will be females_first10
        #called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrms}_PARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "20230121_followups/compare_VCFs/EUR/default/females/PARs/{sample}_{chrms}_PARs_golden_vs_called") 
    output:
        os.path.join(config["scratch_proj_path"], "20230121_followups/compare_VCFs/EUR/default/females/PARs/{sample}_{chrms}_PARs_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

