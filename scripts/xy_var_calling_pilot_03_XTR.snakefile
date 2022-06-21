import os

# Environment: variant_calling_simulations_project

configfile: "xy_var_calling_pilot_03_XTR.config.json"


wildcard_constraints:
    chrms= '|'.join([re.escape(x) for x in config["chrX"]])

rule all:
    input:
        expand(os.path.join(config["scratch_proj_path"], "XTRmasked/chrY/gvcfs/EUR/males/{sample}_{chrmsY}_GRCh38_YPARsMasked_haploid_nonPAR_minusXTR.g.vcf.gz"), sample=config["males"], chrmsY=config["chrY"]),
        expand(os.path.join(config["scratch_proj_path"], "XTRmasked/chrY/combine_g_vcfs/EUR/males/{chrmsY}_GRCh38_YPARsMasked_haploid_nonPAR_minusXTR_combined.g.vcf.gz"), sample=config["males"], chrmsY=config["chrY"]),
        expand(os.path.join(config["scratch_proj_path"], "XTRmasked/chrY/joint_called_vcfs/EUR/males/{chrmsY}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_raw.vcf.gz"), sample=config["males"], chrmsY=config["chrY"]),
        expand(os.path.join(config["scratch_proj_path"], "XTRmasked/chrY/joint_called_vcfs/EUR/males/{chrmsY}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_raw_SNPs.vcf.gz"), sample=config["males"], chrmsY=config["chrY"]),
        expand(os.path.join(config["scratch_proj_path"], "XTRmasked/chrY/hard_filtered_vcfs/EUR/males/all/{chrmsY}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_SNPs_{filtering_options}_filtered.vcf.gz"), chrmsY=config["chrY"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["scratch_proj_path"], "XTRmasked/chrY/hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrmsY}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_SNPs_{filtering_options}_filtered.recode.vcf"), sample=config["males"], chrmsY=config["chrY"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["scratch_proj_path"], "XTRmasked/chrY/hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrmsY}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_SNPs_{filtering_options}_filtered_ALTvariants.vcf.gz"), sample=config["males"], chrmsY=config["chrY"], filtering_options=config["filtering_options"])

        expand(os.path.join(config["proj_path"], "XTRmasked/bams/EUR/males/{sample}_GRCh38_YPARsMasked_sorted.bam"), sample=config["males"]),
        expand(os.path.join(config["proj_path"], "XTRmasked/bams/EUR/males/{sample}_GRCh38_YPARsMasked_sorted.bam.bai"), sample=config["males"]),

        expand(os.path.join(config["proj_path"], "XTRmasked/bams/EUR/males/{sample}_GRCh38_YPARsMasked_sorted_mkdups.bam"), sample=config["males"]),
        expand(os.path.join(config["proj_path"], "XTRmasked/stats/EUR/males/{sample}_GRCh38_YPARsMasked_sorted_mkdups_metrics.txt"), sample=config["males"]),
        expand(os.path.join(config["proj_path"], "XTRmasked/bams/EUR/males/{sample}_GRCh38_YPARsMasked_sorted_mkdups.bam.bai"), sample=config["males"]),

        expand(os.path.join(config["proj_path"], "XTRmasked/gvcfs/EUR/males/{sample}_{chrms}_GRCh38_YPARsMasked_diploid.g.vcf.gz"), sample=config["males"], chrms=config["chromosomes"]), # for now just 8 X and Y...add mtDNA and also in haploid snakefile

        expand(os.path.join(config["proj_path"], "XTRmasked/combine_g_vcfs/EUR/males/{chrms}_GRCh38_YPARsMasked_diploid_combined.g.vcf.gz"), chrms=config["chromosomes"]),
        expand(os.path.join(config["proj_path"], "XTRmasked/joint_called_vcfs/EUR/males/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_raw.vcf.gz"), chrms=config["chromosomes"]),

        expand(os.path.join(config["proj_path"], "XTRmasked/joint_called_vcfs/EUR/males/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs.vcf.gz"), chrms=config["chromosomes"]),

        expand(os.path.join(config["proj_path"], "XTRmasked/hard_filtered_vcfs/EUR/males/all/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered.vcf.gz"), chrms=config["chromosomes"], filtering_options=config["filtering_options"]),

        expand(os.path.join(config["proj_path"], "XTRmasked/hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered.recode.vcf"), sample=config["males"], chrms=config["chromosomes"], filtering_options=config["filtering_options"]),

        expand(os.path.join(config["proj_path"], "XTRmasked/hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered_ALTvariants.vcf.gz"), sample=config["males"], chrms=config["chromosomes"], filtering_options=config["filtering_options"]),
        #expand(os.path.join(config["proj_path"], "XTRmasked/hard_filtered_vcfs/EUR/males/all/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered_nonPARs.vcf.gz"), chrms=config["chrX"], filtering_options=config["filtering_options"]),
        #expand(os.path.join(config["proj_path"], "XTRmasked/hard_filtered_vcfs/EUR/males/all/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered_PARs.vcf.gz"), chrms=config["chrX"], filtering_options=config["filtering_options"]),

        #expand(os.path.join(config["proj_path"], "bcftools_XTRmasked/stats/EUR/males/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered_nonPARs.bcftools.stats.PSC.txt"), chrms=config["chrX"], filtering_options=config["filtering_options"]),
        #expand(os.path.join(config["proj_path"], "bcftools_XTRmasked/stats/EUR/males/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered.bcftools.stats.PSC.txt"), chrms=config["chrY"], filtering_options=config["filtering_options"]),
        # X non-PARs minus XTR
        expand(os.path.join(config["scratch_proj_path"], "XTRmasked/gvcfs/EUR/males/{sample}_{chrms}_GRCh38_YPARsMasked_haploid_nonPAR_minusXTR.g.vcf.gz"), sample=config["males"], chrms=config["chromosomes"]),
        expand(os.path.join(config["scratch_proj_path"], "XTRmasked/combine_g_vcfs/EUR/males/{chrms}_GRCh38_YPARsMasked_haploid_nonPAR_minusXTR_combined.g.vcf.gz"), sample=config["males"], chrms=config["chromosomes"]),
        expand(os.path.join(config["scratch_proj_path"], "XTRmasked/joint_called_vcfs/EUR/males/{chrms}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_raw.vcf.gz"), sample=config["males"], chrms=config["chromosomes"]),
        expand(os.path.join(config["scratch_proj_path"], "XTRmasked/joint_called_vcfs/EUR/males/{chrms}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_raw_SNPs.vcf.gz"), sample=config["males"], chrms=config["chromosomes"]),
        expand(os.path.join(config["scratch_proj_path"], "XTRmasked/hard_filtered_vcfs/EUR/males/all/{chrms}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_SNPs_{filtering_options}_filtered.vcf.gz"), chrms=config["chromosomes"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["scratch_proj_path"], "XTRmasked/hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrms}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_SNPs_{filtering_options}_filtered.recode.vcf"), sample=config["males"], chrms=config["chromosomes"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["scratch_proj_path"], "XTRmasked/hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrms}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_SNPs_{filtering_options}_filtered_ALTvariants.vcf.gz"), sample=config["males"], chrms=config["chromosomes"], filtering_options=config["filtering_options"]),



#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#                                     MALES                                    #
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Step: Alignment (SCC + Y-XTR masked)
#------------------------------------------------------------------------------#
rule alignment_SCC_males:
    input:
        fq1 = os.path.join(config["proj_path"], "trimmed_fastqs/EUR/males/{sample}_NEAT_trimmed_read1.fastq.gz"),
        fq2 = os.path.join(config["proj_path"], "trimmed_fastqs/EUR/males/{sample}_NEAT_trimmed_read2.fastq.gz")
    output:
        os.path.join(config["proj_path"], "XTRmasked/bams/EUR/males/{sample}_GRCh38_YPARsMasked_sorted.bam")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        id = lambda wildcards: config[wildcards.sample]["ID"],
        sm = lambda wildcards: config[wildcards.sample]["SM"],
        lb = lambda wildcards: config[wildcards.sample]["LB"],
        pl = lambda wildcards: config[wildcards.sample]["PL"],
        threads = 4
    shell:
        """
        bwa mem -t {params.threads} -R '@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPL:{params.pl}' {params.ref} {input.fq1} {input.fq2} | samtools fixmate -O bam - - | samtools sort -O bam -o {output}
        """


#------------------------------------------------------------------------------#
# Step: Index bams from last step
#------------------------------------------------------------------------------#
rule index_bam_males:
    input:
        os.path.join(config["proj_path"], "XTRmasked/bams/EUR/males/{sample}_GRCh38_YPARsMasked_sorted.bam")
    output:
        os.path.join(config["proj_path"], "XTRmasked/bams/EUR/males/{sample}_GRCh38_YPARsMasked_sorted.bam.bai")
    shell:
        "samtools index {input}"


#------------------------------------------------------------------------------#
# Step: Mark duplicates
#------------------------------------------------------------------------------#
rule MarkDups_males:
    input:
        bam = os.path.join(config["proj_path"], "XTRmasked/bams/EUR/males/{sample}_GRCh38_YPARsMasked_sorted.bam"),
        bai = os.path.join(config["proj_path"], "XTRmasked/bams/EUR/males/{sample}_GRCh38_YPARsMasked_sorted.bam.bai")
    output:
        bam = os.path.join(config["proj_path"], "XTRmasked/bams/EUR/males/{sample}_GRCh38_YPARsMasked_sorted_mkdups.bam"),
        metrics = os.path.join(config["proj_path"], "XTRmasked/stats/EUR/males/{sample}_GRCh38_YPARsMasked_sorted_mkdups_metrics.txt")
    shell:
        """
        picard -Xmx14g MarkDuplicates I={input.bam} O={output.bam} M={output.metrics} VALIDATION_STRINGENCY=LENIENT
        """


#------------------------------------------------------------------------------#
# Step: Index bams from last step
#------------------------------------------------------------------------------#
rule index_MarkDups_bam_males:
    input:
        os.path.join(config["proj_path"], "XTRmasked/bams/EUR/males/{sample}_GRCh38_YPARsMasked_sorted_mkdups.bam")
    output:
        os.path.join(config["proj_path"], "XTRmasked/bams/EUR/males/{sample}_GRCh38_YPARsMasked_sorted_mkdups.bam.bai")
    shell:
        "samtools index {input}"


#------------------------------------------------------------------------------#
# Step: Call variants - diploid mode
#------------------------------------------------------------------------------#
# For this we will just do XTR
rule make_gvcfs_males_diploid:
    input:
        bam = os.path.join(config["proj_path"], "XTRmasked/bams/EUR/males/{sample}_GRCh38_YPARsMasked_sorted_mkdups.bam"),
        bai = os.path.join(config["proj_path"], "XTRmasked/bams/EUR/males/{sample}_GRCh38_YPARsMasked_sorted_mkdups.bam.bai")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        ploidy = 2,
        chrms = "{chrms}",
        ilist = config["XTR_intervals"]
    output:
        gvcf = os.path.join(config["proj_path"], "XTRmasked/gvcfs/EUR/males/{sample}_{chrms}_GRCh38_YPARsMasked_diploid.g.vcf.gz")
        #gvcf = os.path.join(config["proj_path"], "XTRmasked/gvcfs/EUR/males/{sample}_GRCh38_YPARsMasked_diploid.g.vcf.gz")
    shell:
        """
        gatk --java-options '-Xmx4g' HaplotypeCaller -R {params.ref} -I {input.bam} -L {params.ilist} -ploidy {params.ploidy} -O {output.gvcf} -ERC GVCF
        """

# gatk --java-options '-Xmx4g' HaplotypeCaller -R {params.ref} -I {input.bam} -L {params.chrms} -ploidy {params.ploidy} -O {output.gvcf} -ERC GVCF
rule make_gvcfs_males_haploid_nonPAR_minusXTR:
    input:
        bam = os.path.join(config["proj_path"], "XTRmasked/bams/EUR/males/{sample}_GRCh38_YPARsMasked_sorted_mkdups.bam"),
        bai = os.path.join(config["proj_path"], "XTRmasked/bams/EUR/males/{sample}_GRCh38_YPARsMasked_sorted_mkdups.bam.bai")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        ploidy = 1,
        chrms = "{chrms}",
        ilist = config["XTR_intervals"],
        parilist = config["par_intervals"]
    output:
        gvcf = os.path.join(config["scratch_proj_path"], "XTRmasked/gvcfs/EUR/males/{sample}_{chrms}_GRCh38_YPARsMasked_haploid_nonPAR_minusXTR.g.vcf.gz")
    shell:
        """
        gatk --java-options '-Xmx4g' HaplotypeCaller -R {params.ref} -I {input.bam} -L {params.chrms} -XL {params.ilist} -XL {params.parilist} -ploidy {params.ploidy} -O {output.gvcf} -ERC GVCF
        """

#####
# chrY - haploid called
rule make_gvcfs_males_haploid_nonPAR_minusXTR_chrY:
    input:
        bam = os.path.join(config["proj_path"], "XTRmasked/bams/EUR/males/{sample}_GRCh38_YPARsMasked_sorted_mkdups.bam"),
        bai = os.path.join(config["proj_path"], "XTRmasked/bams/EUR/males/{sample}_GRCh38_YPARsMasked_sorted_mkdups.bam.bai")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        ploidy = 1,
        chrms = "{chrmsY}",
        ilist = config["XTR_intervals_chrY"]
    output:
        gvcf = os.path.join(config["scratch_proj_path"], "XTRmasked/chrY/gvcfs/EUR/males/{sample}_{chrmsY}_GRCh38_YPARsMasked_haploid_nonPAR_minusXTR.g.vcf.gz")
    shell:
        """
        gatk --java-options '-Xmx4g' HaplotypeCaller -R {params.ref} -I {input.bam} -L {params.chrms} -XL {params.ilist} -ploidy {params.ploidy} -O {output.gvcf} -ERC GVCF
        """



#------------------------------------------------------------------------------#
# Step: Combine gvcfs and joint genotype - all together, diploid
#------------------------------------------------------------------------------#
rule CombineGVCFs:
    input:
        gvcf =  expand(os.path.join(config["proj_path"], "XTRmasked/gvcfs/EUR/males/{sample}_{chrms}_GRCh38_YPARsMasked_diploid.g.vcf.gz"), sample=config["males"], chrms=config["chromosomes"])
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        gvcfs = expand(("-V " + config["proj_path"] + "/XTRmasked/gvcfs/EUR/males/{sample}_{chrms}_GRCh38_YPARsMasked_diploid.g.vcf.gz"), sample=config["males"], chrms=config["chromosomes"]),
        chrms = "{chrms}",
        ilist = config["XTR_intervals"]
    output:
        os.path.join(config["proj_path"], "XTRmasked/combine_g_vcfs/EUR/males/{chrms}_GRCh38_YPARsMasked_diploid_combined.g.vcf.gz")
    shell:
        """
        gatk CombineGVCFs -R {params.ref} {params.gvcfs} -L {params.ilist} -O {output}
        """

rule CombineGVCFs_nonPAR_minusXTR:
    input:
        gvcf =  expand(os.path.join(config["scratch_proj_path"], "XTRmasked/gvcfs/EUR/males/{sample}_{chrms}_GRCh38_YPARsMasked_haploid_nonPAR_minusXTR.g.vcf.gz"), sample=config["males"], chrms=config["chromosomes"])
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        gvcfs = expand(("-V " + config["scratch_proj_path"] + "/XTRmasked/gvcfs/EUR/males/{sample}_{chrms}_GRCh38_YPARsMasked_haploid_nonPAR_minusXTR.g.vcf.gz"), sample=config["males"], chrms=config["chromosomes"]),
        chrms = "{chrms}",
        ilist = config["XTR_intervals"],
        parilist = config["par_intervals"]
    output:
        os.path.join(config["scratch_proj_path"], "XTRmasked/combine_g_vcfs/EUR/males/{chrms}_GRCh38_YPARsMasked_haploid_nonPAR_minusXTR_combined.g.vcf.gz")
    shell:
        """
        gatk CombineGVCFs -R {params.ref} {params.gvcfs} -L {params.chrms} -XL {params.ilist} -XL {params.parilist} -O {output}
        """


#####
# chrY
rule CombineGVCFs_nonPAR_minusXTR_chrY:
    input:
        gvcf =  expand(os.path.join(config["scratch_proj_path"], "XTRmasked/chrY/gvcfs/EUR/males/{sample}_{chrmsY}_GRCh38_YPARsMasked_haploid_nonPAR_minusXTR.g.vcf.gz"), sample=config["males"], chrmsY=config["chrY"])
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        gvcfs = expand(("-V " + config["scratch_proj_path"] + "XTRmasked/chrY/gvcfs/EUR/males/{sample}_{chrmsY}_GRCh38_YPARsMasked_haploid_nonPAR_minusXTR.g.vcf.gz"), sample=config["males"], chrmsY=config["chrY"]),
        chrms = "{chrmsY}",
        ilist = config["XTR_intervals_chrY"]
    output:
        os.path.join(config["scratch_proj_path"], "XTRmasked/chrY/combine_g_vcfs/EUR/males/{chrmsY}_GRCh38_YPARsMasked_haploid_nonPAR_minusXTR_combined.g.vcf.gz")
    shell:
        """
        gatk CombineGVCFs -R {params.ref} {params.gvcfs} -L {params.chrms} -XL {params.ilist} -O {output}
        """



rule GenotypeGVCFs:
    input:
        os.path.join(config["proj_path"], "XTRmasked/combine_g_vcfs/EUR/males/{chrms}_GRCh38_YPARsMasked_diploid_combined.g.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}",
        ilist = config["XTR_intervals"]
    output:
        os.path.join(config["proj_path"], "XTRmasked/joint_called_vcfs/EUR/males/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_raw.vcf.gz")
    shell:
        """
        gatk GenotypeGVCFs -R {params.ref} -V {input} -L {params.ilist} -O {output}
        """

# PARS are included here!!
rule GenotypeGVCFs_nonPAR_minusXTR:
    input:
        os.path.join(config["scratch_proj_path"], "XTRmasked/combine_g_vcfs/EUR/males/{chrms}_GRCh38_YPARsMasked_haploid_nonPAR_minusXTR_combined.g.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}",
        ilist = config["XTR_intervals"],
        parilist = config["par_intervals"]
    output:
        os.path.join(config["scratch_proj_path"], "XTRmasked/joint_called_vcfs/EUR/males/{chrms}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_raw.vcf.gz")
    shell:
        """
        gatk GenotypeGVCFs -R {params.ref} -V {input} -L {params.chrms} -XL {params.ilist} -XL {params.parilist} -O {output}
        """

#####
# chrY
rule GenotypeGVCFs_nonPAR_minusXTR_chrY:
    input:
        os.path.join(config["scratch_proj_path"], "XTRmasked/chrY/combine_g_vcfs/EUR/males/{chrmsY}_GRCh38_YPARsMasked_haploid_nonPAR_minusXTR_combined.g.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrmsY}",
        ilist = config["XTR_intervals_chrY"]
    output:
        os.path.join(config["scratch_proj_path"], "XTRmasked/chrY/joint_called_vcfs/EUR/males/{chrmsY}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_raw.vcf.gz")
    shell:
        """
        gatk GenotypeGVCFs -R {params.ref} -V {input} -L {params.chrms} -XL {params.ilist} -O {output}
        """


#------------------------------------------------------------------------------#
# Step: Extract SNPs
#------------------------------------------------------------------------------#
rule select_SNPs:
    input:
        os.path.join(config["proj_path"], "XTRmasked/joint_called_vcfs/EUR/males/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_raw.vcf.gz")
    output:
        os.path.join(config["proj_path"], "XTRmasked/joint_called_vcfs/EUR/males/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}"
    shell:
        """gatk --java-options "-Xmx10g" SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-O {output} """
        """--select-type-to-include SNP """

rule select_SNPs_nonPAR_minusXTR:
    input:
        os.path.join(config["scratch_proj_path"], "XTRmasked/joint_called_vcfs/EUR/males/{chrms}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_raw.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "XTRmasked/joint_called_vcfs/EUR/males/{chrms}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_raw_SNPs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}"
    shell:
        """gatk --java-options "-Xmx10g" SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-O {output} """
        """--select-type-to-include SNP """

#####
# chrY
rule select_SNPs_nonPAR_minusXTR_chrY:
    input:
        os.path.join(config["scratch_proj_path"], "XTRmasked/chrY/joint_called_vcfs/EUR/males/{chrmsY}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_raw.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "XTRmasked/chrY/joint_called_vcfs/EUR/males/{chrmsY}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_raw_SNPs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrmsY}"
    shell:
        """gatk --java-options "-Xmx10g" SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-O {output} """
        """--select-type-to-include SNP """

#------------------------------------------------------------------------------#
# Step: Implement hard filtering on SNPs
#------------------------------------------------------------------------------#
# TO DO - Do each hard filter separately.
rule hard_filter_variant_filtration_all:
    input:
        vcf = os.path.join(config["proj_path"], "XTRmasked/joint_called_vcfs/EUR/males/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
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
        intermediatevcf = os.path.join(config["proj_path"], "XTRmasked/hard_filtered_vcfs/EUR/males/all/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_{filtering_options}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "XTRmasked/hard_filtered_vcfs/EUR/males/all/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "QD < {params.qd}" --filter-name "QD{params.qd}" -filter "QUAL < {params.qual}" --filter-name "QUAL{params.qual}" -filter "SOR > {params.sor}" --filter-name "SOR{params.sor}" -filter "FS > {params.fs}" --filter-name "FS{params.fs}" -filter "MQ < {params.mq}" --filter-name "MQ{params.mq}" -filter "MQRankSum < {params.mqranksum}" --filter-name "MQRankSum{params.mqranksum}" -filter "ReadPosRankSum < {params.readposranksum}" --filter-name "ReadPosRankSum{params.readposranksum}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """


rule hard_filter_variant_filtration_all_nonPAR_minusXTR:
    input:
        vcf = os.path.join(config["scratch_proj_path"], "XTRmasked/joint_called_vcfs/EUR/males/{chrms}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_raw_SNPs.vcf.gz")
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
        intermediatevcf = os.path.join(config["scratch_proj_path"], "XTRmasked/hard_filtered_vcfs/EUR/males/all/{chrms}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_raw_SNPs_{filtering_options}.vcf.gz")
    output:
        ovcf = os.path.join(config["scratch_proj_path"], "XTRmasked/hard_filtered_vcfs/EUR/males/all/{chrms}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_SNPs_{filtering_options}_filtered.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "QD < {params.qd}" --filter-name "QD{params.qd}" -filter "QUAL < {params.qual}" --filter-name "QUAL{params.qual}" -filter "SOR > {params.sor}" --filter-name "SOR{params.sor}" -filter "FS > {params.fs}" --filter-name "FS{params.fs}" -filter "MQ < {params.mq}" --filter-name "MQ{params.mq}" -filter "MQRankSum < {params.mqranksum}" --filter-name "MQRankSum{params.mqranksum}" -filter "ReadPosRankSum < {params.readposranksum}" --filter-name "ReadPosRankSum{params.readposranksum}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """

#####
# chrY
rule hard_filter_variant_filtration_all_nonPAR_minusXTR_chrY:
    input:
        vcf = os.path.join(config["scratch_proj_path"], "XTRmasked/chrY/joint_called_vcfs/EUR/males/{chrmsY}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_raw_SNPs.vcf.gz")
    params:
        chrms = "{chrmsY}",
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        qd = lambda wildcards: config[wildcards.filtering_options]["QD"],
        qual = lambda wildcards: config[wildcards.filtering_options]["QUAL"],
        sor = lambda wildcards: config[wildcards.filtering_options]["SOR"],
        fs = lambda wildcards: config[wildcards.filtering_options]["FS"],
        mq = lambda wildcards: config[wildcards.filtering_options]["MQ"],
        mqranksum = lambda wildcards: config[wildcards.filtering_options]["MQRankSum"],
        readposranksum = lambda wildcards: config[wildcards.filtering_options]["ReadPosRankSum"],
        intermediatevcf = os.path.join(config["scratch_proj_path"], "XTRmasked/chrY/hard_filtered_vcfs/EUR/males/all/{chrmsY}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_raw_SNPs_{filtering_options}.vcf.gz")
    output:
        ovcf = os.path.join(config["scratch_proj_path"], "XTRmasked/chrY/hard_filtered_vcfs/EUR/males/all/{chrmsY}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_SNPs_{filtering_options}_filtered.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "QD < {params.qd}" --filter-name "QD{params.qd}" -filter "QUAL < {params.qual}" --filter-name "QUAL{params.qual}" -filter "SOR > {params.sor}" --filter-name "SOR{params.sor}" -filter "FS > {params.fs}" --filter-name "FS{params.fs}" -filter "MQ < {params.mq}" --filter-name "MQ{params.mq}" -filter "MQRankSum < {params.mqranksum}" --filter-name "MQRankSum{params.mqranksum}" -filter "ReadPosRankSum < {params.readposranksum}" --filter-name "ReadPosRankSum{params.readposranksum}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """


#------------------------------------------------------------------------------#
# Step: Extract samples
#------------------------------------------------------------------------------#
rule extract_males_call_VCF_XTR:
    input:
        os.path.join(config["proj_path"], "XTRmasked/hard_filtered_vcfs/EUR/males/all/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["proj_path"], "XTRmasked/hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered")
    output:
        os.path.join(config["proj_path"], "XTRmasked/hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

# Keep only het or homo alt sites
# bcftools view -i 'GT[*]="alt"' input.vcf
rule keep_alt_sites:
    input:
        os.path.join(config["proj_path"], "XTRmasked/hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered.recode.vcf")
    output:
        os.path.join(config["proj_path"], "XTRmasked/hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered_ALTvariants.vcf.gz")
    shell:
        """
        bcftools view -i 'GT[*]="alt"' {input} -Oz -o {output}
        """

# AND Remove PARs at this step
rule extract_males_call_VCF_nonPAR_minusXTR:
    input:
        os.path.join(config["scratch_proj_path"], "XTRmasked/hard_filtered_vcfs/EUR/males/all/{chrms}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_SNPs_{filtering_options}_filtered.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "XTRmasked/hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrms}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_SNPs_{filtering_options}_filtered")
    output:
        os.path.join(config["scratch_proj_path"], "XTRmasked/hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrms}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_SNPs_{filtering_options}_filtered.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

# Keep only het or homo alt sites
# bcftools view -i 'GT[*]="alt"' input.vcf
rule keep_alt_sites_nonPAR_minusXTR:
    input:
        os.path.join(config["scratch_proj_path"], "XTRmasked/hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrms}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_SNPs_{filtering_options}_filtered.recode.vcf")
    output:
        os.path.join(config["scratch_proj_path"], "XTRmasked/hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrms}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_SNPs_{filtering_options}_filtered_ALTvariants.vcf.gz")
    shell:
        """
        bcftools view -i 'GT[*]="alt"' {input} -Oz -o {output}
        """

#####
# chrY
rule extract_males_call_VCF_nonPAR_minusXTR_chrY:
    input:
        os.path.join(config["scratch_proj_path"], "XTRmasked/chrY/hard_filtered_vcfs/EUR/males/all/{chrmsY}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_SNPs_{filtering_options}_filtered.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "XTRmasked/chrY/hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrmsY}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_SNPs_{filtering_options}_filtered")
    output:
        os.path.join(config["scratch_proj_path"], "XTRmasked/chrY/hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrmsY}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_SNPs_{filtering_options}_filtered.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule keep_alt_sites_nonPAR_minusXTR_chrY:
    input:
        os.path.join(config["scratch_proj_path"], "XTRmasked/chrY/hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrmsY}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_SNPs_{filtering_options}_filtered.recode.vcf")
    output:
        os.path.join(config["scratch_proj_path"], "XTRmasked/chrY/hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrmsY}_GRCh38_YPARsMasked_gatk_haploid_nonPAR_minusXTR_called_SNPs_{filtering_options}_filtered_ALTvariants.vcf.gz")
    shell:
        """
        bcftools view -i 'GT[*]="alt"' {input} -Oz -o {output}
        """
