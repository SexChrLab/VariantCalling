import os

# Environment: variant_calling_simulations_project

configfile: "AFR_xy_var_calling_pilot_03.config.json"


rule all:
    input:
        expand(os.path.join(config["proj_path"], "fastqc/AFR/males/{sample}_NEAT_read1_fastqc.html"), sample=config["males"]),
        expand(os.path.join(config["proj_path"], "fastqc/AFR/males/{sample}_NEAT_read2_fastqc.html"), sample=config["males"]),

        os.path.join(config["proj_path"], "multiqc/AFR/males/multiqc_report.html"),

        expand(os.path.join(config["proj_path"], "trimmed_fastqs/AFR/males/{sample}_NEAT_trimmed_read1.fastq.gz"), sample=config["males"]),
        expand(os.path.join(config["proj_path"], "trimmed_fastqs/AFR/males/{sample}_NEAT_trimmed_read2.fastq.gz"), sample=config["males"]),

        #expand(os.path.join(config["proj_path"], "fastqc_trimmed/AFR/males/{sample}_NEAT_trimmed_read1_fastqc.html"), sample=config["males"]),
        #expand(os.path.join(config["proj_path"], "fastqc_trimmed/AFR/males/{sample}_NEAT_trimmed_read2_fastqc.html"), sample=config["males"]),

        #os.path.join(config["proj_path"], "multiqc_trimmed/AFR/males/multiqc_report.html"),

        #expand(os.path.join(config["proj_path"], "bams/AFR/males/{sample}_GRCh38_YPARsMasked_sorted.bam"), sample=config["males"]),
        #expand(os.path.join(config["proj_path"], "bams/AFR/males/{sample}_GRCh38_YPARsMasked_sorted.bam.bai"), sample=config["males"]),


        #expand(os.path.join(config["proj_path"], "bams/AFR/males/{sample}_GRCh38_YPARsMasked_sorted_mkdups.bam"), sample=config["males"]),
        #expand(os.path.join(config["proj_path"], "stats/AFR/males/{sample}_GRCh38_YPARsMasked_sorted_mkdups_metrics.txt"), sample=config["males"]),
        #expand(os.path.join(config["proj_path"], "bams/AFR/males/{sample}_GRCh38_YPARsMasked_sorted_mkdups.bam.bai"), sample=config["males"]),


        #expand(os.path.join(config["proj_path"], "gvcfs/AFR/males/{sample}_{chrms}_GRCh38_YPARsMasked_diploid.g.vcf.gz"), sample=config["males"], chrms=config["chromosomes"]), # for now just 8 X and Y...add mtDNA and also in haploid snakefile

        #expand(os.path.join(config["proj_path"], "combine_g_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_diploid_combined.g.vcf.gz"), chrms=config["chromosomes"]),
        #expand(os.path.join(config["proj_path"], "joint_called_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_raw.vcf.gz"), chrms=config["chromosomes"]),

        #expand(os.path.join(config["proj_path"], "joint_called_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs.vcf.gz"), chrms=config["chromosomes"]),

        #expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/all/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered.vcf.gz"), chrms=config["chromosomes"], filtering_options=config["filtering_options"]),

        #expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/all/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered_nonPARs.vcf.gz"), chrms=config["chrX"], filtering_options=config["filtering_options"]),
        #expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/all/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered_PARs.vcf.gz"), chrms=config["chrX"], filtering_options=config["filtering_options"]),

        #expand(os.path.join(config["proj_path"], "bcftools_stats/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered_nonPARs.bcftools.stats.PSC.txt"), chrms=config["chrX"], filtering_options=config["filtering_options"]),
        #expand(os.path.join(config["proj_path"], "bcftools_stats/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered.bcftools.stats.PSC.txt"), chrms=config["chrY"], filtering_options=config["filtering_options"]),

        ##                         ##
        # MALES - DEFAULT MAPPING #
        expand(os.path.join(config["proj_path"], "bams/AFR/males/default/{sample}_GRCh38_default_sorted.bam"), sample=config["males"]),
        expand(os.path.join(config["proj_path"], "bams/AFR/males/default/{sample}_GRCh38_default_sorted.bam.bai"), sample=config["males"]),

        expand(os.path.join(config["proj_path"], "bams/AFR/males/default/{sample}_GRCh38_default_sorted_mkdups.bam"), sample=config["males"]),
        expand(os.path.join(config["proj_path"], "stats/AFR/males/default/{sample}_GRCh38_default_sorted_mkdups_metrics.txt"), sample=config["males"]),
        expand(os.path.join(config["proj_path"], "bams/AFR/males/default/{sample}_GRCh38_default_sorted_mkdups.bam.bai"), sample=config["males"]),

        expand(os.path.join(config["proj_path"], "gvcfs/AFR/males/default/{sample}_{chrms}_GRCh38_default_diploid.g.vcf.gz"), sample=config["males"], chrms=config["chromosomes"]), # for now just 8 X and Y...add mtDNA and also in haploid snakefile

        expand(os.path.join(config["proj_path"], "combine_g_vcfs/AFR/males/default/{chrms}_GRCh38_default_diploid_combined.g.vcf.gz"), chrms=config["chromosomes"]),
        expand(os.path.join(config["proj_path"], "joint_called_vcfs/AFR/males/default/{chrms}_GRCh38_default_gatk_diploid_called_raw.vcf.gz"), chrms=config["chromosomes"]),

        expand(os.path.join(config["proj_path"], "joint_called_vcfs/AFR/males/default/{chrms}_GRCh38_default_gatk_diploid_called_raw_SNPs.vcf.gz"), chrms=config["chromosomes"]),

        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/default/all/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered.vcf.gz"), chrms=config["chromosomes"], filtering_options=config["filtering_options"]),

        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/default/all/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered_nonPARs.vcf.gz"), chrms=config["chrX"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/default/all/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered_PARs.vcf.gz"), chrms=config["chrX"], filtering_options=config["filtering_options"]),

        #expand(os.path.join(config["proj_path"], "bcftools_stats/AFR/males/default/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered_nonPARs.bcftools.stats.PSC.txt"), chrms=config["chrX"], filtering_options=config["filtering_options"]),
        #expand(os.path.join(config["proj_path"], "bcftools_stats/AFR/males/default/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered.bcftools.stats.PSC.txt"), chrms=config["chrY"], filtering_options=config["filtering_options"]),

        ##                         ##
        # FEMALES - SCC MAPPING #
        expand(os.path.join(config["proj_path"], "fastqc/AFR/females/{sample}_NEAT_read1_fastqc.html"), sample=config["females"]),
        expand(os.path.join(config["proj_path"], "fastqc/AFR/females/{sample}_NEAT_read2_fastqc.html"), sample=config["females"]),

        os.path.join(config["proj_path"], "multiqc/AFR/females/multiqc_report.html"),

        expand(os.path.join(config["proj_path"], "trimmed_fastqs/AFR/females/{sample}_NEAT_trimmed_read1.fastq.gz"), sample=config["females"]),
        expand(os.path.join(config["proj_path"], "trimmed_fastqs/AFR/females/{sample}_NEAT_trimmed_read2.fastq.gz"), sample=config["females"]),

        #expand(os.path.join(config["proj_path"], "fastqc_trimmed/AFR/females/{sample}_NEAT_trimmed_read1_fastqc.html"), sample=config["females"]),
        #expand(os.path.join(config["proj_path"], "fastqc_trimmed/AFR/females/{sample}_NEAT_trimmed_read2_fastqc.html"), sample=config["females"]),

        #os.path.join(config["proj_path"], "multiqc_trimmed/AFR/females/multiqc_report.html"),

        #expand(os.path.join(config["proj_path"], "bams/AFR/females/{sample}_GRCh38_YHardMasked_sorted.bam"), sample=config["females"]),
        #expand(os.path.join(config["proj_path"], "bams/AFR/females/{sample}_GRCh38_YHardMasked_sorted.bam.bai"), sample=config["females"]),

        #expand(os.path.join(config["proj_path"], "bams/AFR/females/{sample}_GRCh38_YHardMasked_sorted_mkdups.bam"), sample=config["females"]),
        #expand(os.path.join(config["proj_path"], "stats/AFR/females/{sample}_GRCh38_YHardMasked_sorted_mkdups_metrics.txt"), sample=config["females"]),
        #expand(os.path.join(config["proj_path"], "bams/AFR/females/{sample}_GRCh38_YHardMasked_sorted_mkdups.bam.bai"), sample=config["females"]),

        #expand(os.path.join(config["proj_path"], "gvcfs/AFR/females/{sample}_{chrms}_GRCh38_YHardMasked_diploid.g.vcf.gz"), sample=config["females"], chrms=config["chromosomes"]), # for now just 8 X and Y

        #expand(os.path.join(config["proj_path"], "combine_g_vcfs/AFR/females/{chrms}_GRCh38_YHardMasked_diploid_combined.g.vcf.gz"), chrms=config["chromosomes"]),
        #expand(os.path.join(config["proj_path"], "joint_called_vcfs/AFR/females/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_raw.vcf.gz"), chrms=config["chromosomes"]),

        #expand(os.path.join(config["proj_path"], "joint_called_vcfs/AFR/females/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs.vcf.gz"), chrms=config["chromosomes"]),

        #expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/all/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered.vcf.gz"), chrms=config["chromosomes"], filtering_options=config["filtering_options"]),

        #expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/all/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered_nonPARs.vcf.gz"), chrms=config["chrX"], filtering_options=config["filtering_options"]),
        #expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/all/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered_PARs.vcf.gz"), chrms=config["chrX"], filtering_options=config["filtering_options"]),

        ##                         ##
        # FEMALES - DEFAULT MAPPING #
        expand(os.path.join(config["proj_path"], "bams/AFR/females/default/{sample}_GRCh38_default_sorted.bam"), sample=config["females"]),
        expand(os.path.join(config["proj_path"], "bams/AFR/females/default/{sample}_GRCh38_default_sorted.bam.bai"), sample=config["females"]),

        expand(os.path.join(config["proj_path"], "bams/AFR/females/default/{sample}_GRCh38_default_sorted_mkdups.bam"), sample=config["females"]),
        expand(os.path.join(config["proj_path"], "stats/AFR/females/default/{sample}_GRCh38_default_sorted_mkdups_metrics.txt"), sample=config["females"]),
        expand(os.path.join(config["proj_path"], "bams/AFR/females/default/{sample}_GRCh38_default_sorted_mkdups.bam.bai"), sample=config["females"]),

        expand(os.path.join(config["proj_path"], "gvcfs/AFR/females/default/{sample}_{chrms}_GRCh38_default_diploid.g.vcf.gz"), sample=config["females"], chrms=config["chromosomes"]), # for now just 8 X and Y

        expand(os.path.join(config["proj_path"], "combine_g_vcfs/AFR/females/default/{chrms}_GRCh38_default_diploid_combined.g.vcf.gz"), chrms=config["chromosomes"]),
        expand(os.path.join(config["proj_path"], "joint_called_vcfs/AFR/females/default/{chrms}_GRCh38_default_gatk_diploid_called_raw.vcf.gz"), chrms=config["chromosomes"]),

        expand(os.path.join(config["proj_path"], "joint_called_vcfs/AFR/females/default/{chrms}_GRCh38_default_gatk_diploid_called_raw_SNPs.vcf.gz"), chrms=config["chromosomes"]),

        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/default/all/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered.vcf.gz"), chrms=config["chromosomes"], filtering_options=config["filtering_options"]),

        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/default/all/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered_nonPARs.vcf.gz"), chrms=config["chrX"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/default/all/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered_PARs.vcf.gz"), chrms=config["chrX"], filtering_options=config["filtering_options"]),

        ###expand(os.path.join(config["proj_path"], "bcftools_stats/AFR/females/default/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered_nonPARs.bcftools.stats.PSC.txt"), chrms=config["chrX"], filtering_options=config["filtering_options"]),
        #expand(os.path.join(config["proj_path"], "bcftools_stats/AFR/females/default/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered.bcftools.stats.PSC.txt"), chrms=config["chrY"], filtering_options=config["filtering_options"]),


'''
# Males SCC filters TO DO later
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/QD/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}.vcf.gz"), chrms=config["chromosomes"], filter=config["QD_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/QUAL/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}.vcf.gz"), chrms=config["chromosomes"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/SOR/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}.vcf.gz"), chrms=config["chromosomes"], filter=config["SOR_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/FS/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}.vcf.gz"), chrms=config["chromosomes"], filter=config["FS_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/MQ/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}.vcf.gz"), chrms=config["chromosomes"], filter=config["MQ_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/MQRankSum/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}.vcf.gz"), chrms=config["chromosomes"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/ReadPosRankSum/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}.vcf.gz"), chrms=config["chromosomes"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/AN/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}.vcf.gz"), chrms=config["chromosomes"], filter=config["AN_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/AN/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}.vcf.gz.tbi"), chrms=config["chromosomes"], filter=config["AN_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/DP/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}.vcf.gz"), chrms=config["chromosomes"], filter=config["DP_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/DP/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}.vcf.gz.tbi"), chrms=config["chromosomes"], filter=config["DP_filters"]),

        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/QD/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filter=config["QD_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/QUAL/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/SOR/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filter=config["SOR_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/FS/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filter=config["FS_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/MQ/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filter=config["MQ_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/MQRankSum/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/ReadPosRankSum/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/AN/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filter=config["AN_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/DP/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filter=config["DP_filters"]),

        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/QD/chrX_PARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}_PARs.vcf.gz"), chrms=config["chrX"], filter=config["QD_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/QUAL/chrX_PARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}_PARs.vcf.gz"), chrms=config["chrX"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/SOR/chrX_PARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}_PARs.vcf.gz"), chrms=config["chrX"], filter=config["SOR_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/FS/chrX_PARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}_PARs.vcf.gz"), chrms=config["chrX"], filter=config["FS_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/MQ/chrX_PARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}_PARs.vcf.gz"), chrms=config["chrX"], filter=config["MQ_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/MQRankSum/chrX_PARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}_PARs.vcf.gz"), chrms=config["chrX"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/ReadPosRankSum/chrX_PARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}_PARs.vcf.gz"), chrms=config["chrX"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/AN/chrX_PARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}_PARs.vcf.gz"), chrms=config["chrX"], filter=config["AN_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/DP/chrX_PARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}_PARs.vcf.gz"), chrms=config["chrX"], filter=config["DP_filters"]),
'''

'''
        ##                         ##
        # MALES - DEFAULT MAPPING #
        expand(os.path.join(config["proj_path"], "bams/AFR/males/default/{sample}_GRCh38_default_sorted.bam"), sample=config["males"]),
        expand(os.path.join(config["proj_path"], "bams/AFR/males/default/{sample}_GRCh38_default_sorted.bam.bai"), sample=config["males"]),

        expand(os.path.join(config["proj_path"], "bams/AFR/males/default/{sample}_GRCh38_default_sorted_mkdups.bam"), sample=config["males"]),
        expand(os.path.join(config["proj_path"], "stats/AFR/males/default/{sample}_GRCh38_default_sorted_mkdups_metrics.txt"), sample=config["males"]),
        expand(os.path.join(config["proj_path"], "bams/AFR/males/default/{sample}_GRCh38_default_sorted_mkdups.bam.bai"), sample=config["males"]),

        expand(os.path.join(config["proj_path"], "gvcfs/AFR/males/default/{sample}_{chrms}_GRCh38_default_diploid.g.vcf.gz"), sample=config["males"], chrms=config["chromosomes"]), # for now just 8 X and Y...add mtDNA and also in haploid snakefile

        expand(os.path.join(config["proj_path"], "combine_g_vcfs/AFR/males/default/{chrms}_GRCh38_default_diploid_combined.g.vcf.gz"), chrms=config["chromosomes"]),
        expand(os.path.join(config["proj_path"], "joint_called_vcfs/AFR/males/default/{chrms}_GRCh38_default_gatk_diploid_called_raw.vcf.gz"), chrms=config["chromosomes"]),

        expand(os.path.join(config["proj_path"], "joint_called_vcfs/AFR/males/default/{chrms}_GRCh38_default_gatk_diploid_called_raw_SNPs.vcf.gz"), chrms=config["chromosomes"]),

        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/default/all/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered.vcf.gz"), chrms=config["chromosomes"], filtering_options=config["filtering_options"]),

        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/default/all/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered_nonPARs.vcf.gz"), chrms=config["chrX"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/default/all/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered_PARs.vcf.gz"), chrms=config["chrX"], filtering_options=config["filtering_options"]),

        expand(os.path.join(config["proj_path"], "bcftools_stats/AFR/males/default/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered_nonPARs.bcftools.stats.PSC.txt"), chrms=config["chrX"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["proj_path"], "bcftools_stats/AFR/males/default/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered.bcftools.stats.PSC.txt"), chrms=config["chrY"], filtering_options=config["filtering_options"]),

        ##                         ##
        # FEMALES - SCC MAPPING #
        #expand(os.path.join(config["proj_path"], "fastqc/AFR/females/{sample}_NEAT_read1_fastqc.html"), sample=config["females"]),
        #expand(os.path.join(config["proj_path"], "fastqc/AFR/females/{sample}_NEAT_read2_fastqc.html"), sample=config["females"]),

        #os.path.join(config["proj_path"], "multiqc/AFR/females/multiqc_report.html"),

        #expand(os.path.join(config["proj_path"], "trimmed_fastqs/AFR/females/{sample}_NEAT_trimmed_read1.fastq.gz"), sample=config["females"]),
        #expand(os.path.join(config["proj_path"], "trimmed_fastqs/AFR/females/{sample}_NEAT_trimmed_read2.fastq.gz"), sample=config["females"]),

        #expand(os.path.join(config["proj_path"], "fastqc_trimmed/AFR/females/{sample}_NEAT_trimmed_read1_fastqc.html"), sample=config["females"]),
        #expand(os.path.join(config["proj_path"], "fastqc_trimmed/AFR/females/{sample}_NEAT_trimmed_read2_fastqc.html"), sample=config["females"]),

        #os.path.join(config["proj_path"], "multiqc_trimmed/AFR/females/multiqc_report.html"),

        expand(os.path.join(config["proj_path"], "bams/AFR/females/{sample}_GRCh38_YHardMasked_sorted.bam"), sample=config["females"]),
        expand(os.path.join(config["proj_path"], "bams/AFR/females/{sample}_GRCh38_YHardMasked_sorted.bam.bai"), sample=config["females"]),

        expand(os.path.join(config["proj_path"], "bams/AFR/females/{sample}_GRCh38_YHardMasked_sorted_mkdups.bam"), sample=config["females"]),
        expand(os.path.join(config["proj_path"], "stats/AFR/females/{sample}_GRCh38_YHardMasked_sorted_mkdups_metrics.txt"), sample=config["females"]),
        expand(os.path.join(config["proj_path"], "bams/AFR/females/{sample}_GRCh38_YHardMasked_sorted_mkdups.bam.bai"), sample=config["females"]),

        expand(os.path.join(config["proj_path"], "gvcfs/AFR/females/{sample}_{chrms}_GRCh38_YHardMasked_diploid.g.vcf.gz"), sample=config["females"], chrms=config["chromosomes"]), # for now just 8 X and Y

        expand(os.path.join(config["proj_path"], "combine_g_vcfs/AFR/females/{chrms}_GRCh38_YHardMasked_diploid_combined.g.vcf.gz"), chrms=config["chromosomes"]),
        expand(os.path.join(config["proj_path"], "joint_called_vcfs/AFR/females/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_raw.vcf.gz"), chrms=config["chromosomes"]),

        expand(os.path.join(config["proj_path"], "joint_called_vcfs/AFR/females/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs.vcf.gz"), chrms=config["chromosomes"]),

        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/all/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered.vcf.gz"), chrms=config["chromosomes"], filtering_options=config["filtering_options"]),

        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/all/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered_nonPARs.vcf.gz"), chrms=config["chrX"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/all/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered_PARs.vcf.gz"), chrms=config["chrX"], filtering_options=config["filtering_options"]),

        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/QD/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}.vcf.gz"), chrms=config["chromosomes"], filter=config["QD_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/QUAL/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}.vcf.gz"), chrms=config["chromosomes"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/SOR/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}.vcf.gz"), chrms=config["chromosomes"], filter=config["SOR_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/FS/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}.vcf.gz"), chrms=config["chromosomes"], filter=config["FS_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/MQ/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}.vcf.gz"), chrms=config["chromosomes"], filter=config["MQ_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/MQRankSum/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}.vcf.gz"), chrms=config["chromosomes"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/ReadPosRankSum/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}.vcf.gz"), chrms=config["chromosomes"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/AN/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}.vcf.gz"), chrms=config["chromosomes"], filter=config["AN_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/AN/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}.vcf.gz.tbi"), chrms=config["chromosomes"], filter=config["AN_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/DP/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}.vcf.gz"), chrms=config["chromosomes"], filter=config["DP_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/DP/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}.vcf.gz.tbi"), chrms=config["chromosomes"], filter=config["DP_filters"]),

        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/QD/chrX_nonPARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filter=config["QD_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/QUAL/chrX_nonPARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/SOR/chrX_nonPARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filter=config["SOR_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/FS/chrX_nonPARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filter=config["FS_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/MQ/chrX_nonPARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filter=config["MQ_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/MQRankSum/chrX_nonPARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/ReadPosRankSum/chrX_nonPARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/AN/chrX_nonPARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filter=config["AN_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/DP/chrX_nonPARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}_nonPARs.vcf.gz"), chrms=config["chrX"], filter=config["DP_filters"]),

        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/QD/chrX_PARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}_PARs.vcf.gz"), chrms=config["chrX"], filter=config["QD_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/QUAL/chrX_PARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}_PARs.vcf.gz"), chrms=config["chrX"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/SOR/chrX_PARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}_PARs.vcf.gz"), chrms=config["chrX"], filter=config["SOR_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/FS/chrX_PARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}_PARs.vcf.gz"), chrms=config["chrX"], filter=config["FS_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/MQ/chrX_PARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}_PARs.vcf.gz"), chrms=config["chrX"], filter=config["MQ_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/MQRankSum/chrX_PARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}_PARs.vcf.gz"), chrms=config["chrX"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/ReadPosRankSum/chrX_PARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}_PARs.vcf.gz"), chrms=config["chrX"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/AN/chrX_PARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}_PARs.vcf.gz"), chrms=config["chrX"], filter=config["AN_filters"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/DP/chrX_PARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}_PARs.vcf.gz"), chrms=config["chrX"], filter=config["DP_filters"]),



        ##                         ##
        # FEMALES - DEFAULT MAPPING #
        expand(os.path.join(config["proj_path"], "bams/AFR/females/default/{sample}_GRCh38_default_sorted.bam"), sample=config["females"]),
        expand(os.path.join(config["proj_path"], "bams/AFR/females/default/{sample}_GRCh38_default_sorted.bam.bai"), sample=config["females"]),

        expand(os.path.join(config["proj_path"], "bams/AFR/females/default/{sample}_GRCh38_default_sorted_mkdups.bam"), sample=config["females"]),
        expand(os.path.join(config["proj_path"], "stats/AFR/females/default/{sample}_GRCh38_default_sorted_mkdups_metrics.txt"), sample=config["females"]),
        expand(os.path.join(config["proj_path"], "bams/AFR/females/default/{sample}_GRCh38_default_sorted_mkdups.bam.bai"), sample=config["females"]),

        expand(os.path.join(config["proj_path"], "gvcfs/AFR/females/default/{sample}_{chrms}_GRCh38_default_diploid.g.vcf.gz"), sample=config["females"], chrms=config["chromosomes"]), # for now just 8 X and Y

        expand(os.path.join(config["proj_path"], "combine_g_vcfs/AFR/females/default/{chrms}_GRCh38_default_diploid_combined.g.vcf.gz"), chrms=config["chromosomes"]),
        expand(os.path.join(config["proj_path"], "joint_called_vcfs/AFR/females/default/{chrms}_GRCh38_default_gatk_diploid_called_raw.vcf.gz"), chrms=config["chromosomes"]),

        expand(os.path.join(config["proj_path"], "joint_called_vcfs/AFR/females/default/{chrms}_GRCh38_default_gatk_diploid_called_raw_SNPs.vcf.gz"), chrms=config["chromosomes"]),

        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/default/all/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered.vcf.gz"), chrms=config["chromosomes"], filtering_options=config["filtering_options"]),

        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/default/all/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered_nonPARs.vcf.gz"), chrms=config["chrX"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/default/all/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered_PARs.vcf.gz"), chrms=config["chrX"], filtering_options=config["filtering_options"]),

        ###expand(os.path.join(config["proj_path"], "bcftools_stats/AFR/females/default/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered_nonPARs.bcftools.stats.PSC.txt"), chrms=config["chrX"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["proj_path"], "bcftools_stats/AFR/females/default/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered.bcftools.stats.PSC.txt"), chrms=config["chrY"], filtering_options=config["filtering_options"]),
'''


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#                                     MALES                                    #
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Step: fastqc and multiQC
#------------------------------------------------------------------------------#
rule fastqc_analysis_read1:
    input:
        fq1 = os.path.join(config["proj_path"], "simulations/AFR/males/{sample}/{sample}_NEAT_read1.fastq.gz"),
    output:
        f1 = os.path.join(config["proj_path"], "fastqc/AFR/males/{sample}_NEAT_read1_fastqc.html"),
    params:
        fastqc = os.path.join(config["proj_path"], "fastqc/AFR/males/")
    shell:
        """
        PERL5LIB="";
        fastqc -o {params.fastqc} {input.fq1}
        """

rule fastqc_analysis_read2:
    input:
        fq2 = os.path.join(config["proj_path"], "simulations/AFR/males/{sample}/{sample}_NEAT_read2.fastq.gz")
    output:
        f2 = os.path.join(config["proj_path"], "fastqc/AFR/males/{sample}_NEAT_read2_fastqc.html")
    params:
        fastqc = os.path.join(config["proj_path"], "fastqc/AFR/males/")
    shell:
        """
        PERL5LIB="";
        fastqc -o {params.fastqc} {input.fq2}
        """


rule multiqc_analysis:
    input:
        f1 = expand(os.path.join(config["proj_path"], "fastqc/AFR/males/{sample}_NEAT_read1_fastqc.html"),
        sample=config["males"]),
        f2 = expand(os.path.join(config["proj_path"], "fastqc/AFR/males/{sample}_NEAT_read2_fastqc.html"),
        sample=config["males"])
    output:
        os.path.join(config["proj_path"], "multiqc/AFR/males/multiqc_report.html")
    params:
        fastqc = os.path.join(config["proj_path"], "fastqc/AFR/males/"),
        multiqc = os.path.join(config["proj_path"], "multiqc/AFR/males/")
    shell:
        "export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 &&"
        "multiqc --interactive -o {params.multiqc} {params.fastqc}"


#------------------------------------------------------------------------------#
# Step: Trimming
#------------------------------------------------------------------------------#
rule trim_bbduk:
    input:
        fq1 = os.path.join(config["proj_path"], "simulations/AFR/males/{sample}/{sample}_NEAT_read1.fastq.gz"),
        fq2 = os.path.join(config["proj_path"], "simulations/AFR/males/{sample}/{sample}_NEAT_read2.fastq.gz")
    output:
        f1 = os.path.join(config["proj_path"], "trimmed_fastqs/AFR/males/{sample}_NEAT_trimmed_read1.fastq.gz"),
        f2 = os.path.join(config["proj_path"], "trimmed_fastqs/AFR/males/{sample}_NEAT_trimmed_read2.fastq.gz")
    shell:
        "bbduk.sh -Xmx3g in1={input.fq1} in2={input.fq2} "
        "out1={output.f1} out2={output.f2} "
        "qtrim=rl trimq=20 minlen=75 maq=20"


#------------------------------------------------------------------------------#
# Step: fastqc and multiQC
#------------------------------------------------------------------------------#
rule fastqc_analysis_trimmed_1:
    input:
        fq1 = os.path.join(config["proj_path"], "trimmed_fastqs/AFR/males/{sample}_NEAT_trimmed_read1.fastq.gz"),
    output:
        f1 = os.path.join(config["proj_path"], "fastqc_trimmed/AFR/males/{sample}_NEAT_trimmed_read1_fastqc.html"),
    params:
        fastqc = os.path.join(config["proj_path"], "fastqc_trimmed/AFR/males/")
    shell:
        """
        PERL5LIB="";
        fastqc -o {params.fastqc} {input.fq1}
        """

rule fastqc_analysis_trimmed_2:
    input:
        fq2 = os.path.join(config["proj_path"], "trimmed_fastqs/AFR/males/{sample}_NEAT_trimmed_read2.fastq.gz")
    output:
        f2 = os.path.join(config["proj_path"], "fastqc_trimmed/AFR/males/{sample}_NEAT_trimmed_read2_fastqc.html")
    params:
        fastqc = os.path.join(config["proj_path"], "fastqc_trimmed/AFR/males/")
    shell:
        """
        PERL5LIB="";
        fastqc -o {params.fastqc} {input.fq2}
        """

rule multiqc_analysis_trimmed:
    input:
        f1 = expand(os.path.join(config["proj_path"], "fastqc_trimmed/AFR/males/{sample}_NEAT_trimmed_read1_fastqc.html"),
        sample=config["males"]),
        f2 = expand(os.path.join(config["proj_path"], "fastqc_trimmed/AFR/males/{sample}_NEAT_trimmed_read2_fastqc.html"),
        sample=config["males"])
    output:
        os.path.join(config["proj_path"], "multiqc_trimmed/AFR/males/multiqc_report.html")
    params:
        fastqc = os.path.join(config["proj_path"], "fastqc_trimmed/AFR/males/"),
        multiqc = os.path.join(config["proj_path"], "multiqc_trimmed/AFR/males/")
    shell:
        "export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 &&"
        "multiqc --interactive -o {params.multiqc} {params.fastqc}"


#------------------------------------------------------------------------------#
# Step: Alignment (SCC)
#------------------------------------------------------------------------------#
rule alignment_SCC_males:
    input:
        fq1 = os.path.join(config["proj_path"], "trimmed_fastqs/AFR/males/{sample}_NEAT_trimmed_read1.fastq.gz"),
        fq2 = os.path.join(config["proj_path"], "trimmed_fastqs/AFR/males/{sample}_NEAT_trimmed_read2.fastq.gz")
    output:
        os.path.join(config["proj_path"], "bams/AFR/males/{sample}_GRCh38_YPARsMasked_sorted.bam")
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
        os.path.join(config["proj_path"], "bams/AFR/males/{sample}_GRCh38_YPARsMasked_sorted.bam")
    output:
        os.path.join(config["proj_path"], "bams/AFR/males/{sample}_GRCh38_YPARsMasked_sorted.bam.bai")
    shell:
        "samtools index {input}"


#------------------------------------------------------------------------------#
# Step: Mark duplicates
#------------------------------------------------------------------------------#
rule MarkDups_males:
    input:
        bam = os.path.join(config["proj_path"], "bams/AFR/males/{sample}_GRCh38_YPARsMasked_sorted.bam"),
        bai = os.path.join(config["proj_path"], "bams/AFR/males/{sample}_GRCh38_YPARsMasked_sorted.bam.bai")
    output:
        bam = os.path.join(config["proj_path"], "bams/AFR/males/{sample}_GRCh38_YPARsMasked_sorted_mkdups.bam"),
        metrics = os.path.join(config["proj_path"], "stats/AFR/males/{sample}_GRCh38_YPARsMasked_sorted_mkdups_metrics.txt")
    shell:
        """
        picard -Xmx14g MarkDuplicates I={input.bam} O={output.bam} M={output.metrics} VALIDATION_STRINGENCY=LENIENT
        """


#------------------------------------------------------------------------------#
# Step: Index bams from last step
#------------------------------------------------------------------------------#
rule index_MarkDups_bam_males:
    input:
        os.path.join(config["proj_path"], "bams/AFR/males/{sample}_GRCh38_YPARsMasked_sorted_mkdups.bam")
    output:
        os.path.join(config["proj_path"], "bams/AFR/males/{sample}_GRCh38_YPARsMasked_sorted_mkdups.bam.bai")
    shell:
        "samtools index {input}"


#------------------------------------------------------------------------------#
# Step: Call variants - diploid mode
#------------------------------------------------------------------------------#
# TO DO: Once scratch is back up:
# - run this step with EMIT_ALL_ACTIVE_SITES so that I can get ref sites that are
#   called. This will help me calculate true negatives

rule make_gvcfs_males_diploid:
    input:
        bam = os.path.join(config["proj_path"], "bams/AFR/males/{sample}_GRCh38_YPARsMasked_sorted_mkdups.bam"),
        bai = os.path.join(config["proj_path"], "bams/AFR/males/{sample}_GRCh38_YPARsMasked_sorted_mkdups.bam.bai")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        ploidy = 2,
        chrms = "{chrms}"
    output:
        gvcf = os.path.join(config["proj_path"], "gvcfs/AFR/males/{sample}_{chrms}_GRCh38_YPARsMasked_diploid.g.vcf.gz")
        #gvcf = os.path.join(config["proj_path"], "gvcfs/AFR/males/{sample}_GRCh38_YPARsMasked_diploid.g.vcf.gz")
    shell:
        """
        gatk --java-options '-Xmx4g' HaplotypeCaller -R {params.ref} -I {input.bam} -L {params.chrms} -ploidy {params.ploidy} -O {output.gvcf} -ERC GVCF
        """

# gatk --java-options '-Xmx4g' HaplotypeCaller -R {params.ref} -I {input.bam} -L {params.chrms} -ploidy {params.ploidy} -O {output.gvcf} -ERC GVCF



#------------------------------------------------------------------------------#
# Step: Combine gvcfs and joint genotype - all together, diploid
#------------------------------------------------------------------------------#
rule CombineGVCFs:
    input:
        gvcf =  expand(os.path.join(config["proj_path"], "gvcfs/AFR/males/{sample}_{chrms}_GRCh38_YPARsMasked_diploid.g.vcf.gz"), sample=config["males"], chrms=config["chromosomes"])
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        gvcfs = expand(("-V " + config["proj_path"] + "/gvcfs/AFR/males/{sample}_{chrms}_GRCh38_YPARsMasked_diploid.g.vcf.gz"), sample=config["males"], chrms=config["chromosomes"]),
        chrms = "{chrms}"
    output:
        os.path.join(config["proj_path"], "combine_g_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_diploid_combined.g.vcf.gz")
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
        os.path.join(config["proj_path"], "combine_g_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_diploid_combined.g.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}"
    output:
        os.path.join(config["proj_path"], "joint_called_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_raw.vcf.gz")
    shell:
        """
        gatk GenotypeGVCFs -R {params.ref} -V {input} -L {params.chrms} -O {output}
        """

#------------------------------------------------------------------------------#
# Step: Extract SNPs
#------------------------------------------------------------------------------#
rule select_SNPs:
    input:
        os.path.join(config["proj_path"], "joint_called_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_raw.vcf.gz")
    output:
        os.path.join(config["proj_path"], "joint_called_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
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
# TO DO - Do each hard filter separately.
rule hard_filter_variant_filtration_all:
    input:
        vcf = os.path.join(config["proj_path"], "joint_called_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
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
        intermediatevcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/all/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_{filtering_options}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/all/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "QD < {params.qd}" --filter-name "QD{params.qd}" -filter "QUAL < {params.qual}" --filter-name "QUAL{params.qual}" -filter "SOR > {params.sor}" --filter-name "SOR{params.sor}" -filter "FS > {params.fs}" --filter-name "FS{params.fs}" -filter "MQ < {params.mq}" --filter-name "MQ{params.mq}" -filter "MQRankSum < {params.mqranksum}" --filter-name "MQRankSum{params.mqranksum}" -filter "ReadPosRankSum < {params.readposranksum}" --filter-name "ReadPosRankSum{params.readposranksum}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """


#------------------------------------------------------------------------------#
# Step: For X chromosome, extract non PARS (exclude PARs)
#------------------------------------------------------------------------------#
rule extract_X_nonPARs:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/all/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered.vcf.gz")
    output:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/all/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered_nonPARs.vcf.gz")
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
# TODO: Step: For X chromosome, extract PARS (we will also analyze PARs)
#------------------------------------------------------------------------------#
rule extract_X_PARs:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/all/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered.vcf.gz")
    output:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/all/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered_PARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.ilist} """
        """-O {output} """


#------------------------------------------------------------------------------#
# Step: For Y and X chromosomes (non PARS) get per sample stats
#------------------------------------------------------------------------------#
rule bcftools_PSCs_XnonPARs:
    input:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/all/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered_nonPARs.vcf.gz")
    output:
        os.path.join(config["proj_path"], "bcftools_stats/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered_nonPARs.bcftools.stats.PSC.txt")
    shell:
        """
        bcftools stats -s - {input} | grep PSC > {output}
        """

rule bcftools_PSCs_Y:
    input:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/all/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered.vcf.gz")
    output:
        os.path.join(config["proj_path"], "bcftools_stats/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered.bcftools.stats.PSC.txt")
    shell:
        """
        bcftools stats -s - {input} | grep PSC > {output}
        """


#------------------------------------------------------------------------------#
# Step: Repeat but for each filter separately
#------------------------------------------------------------------------------#
##############
# Filter: QD #
##############
rule hard_filter_variant_filtration_QD:
    input:
        vcf = os.path.join(config["proj_path"], "joint_called_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    params:
        chrms = "{chrms}",
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        qd = "{filter}",
        intermediatevcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/QD/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_QD_{filter}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/QD/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "QD < {params.qd}" --filter-name "QD{params.qd}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """

rule extract_X_nonPARs_QD:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/QD/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/QD/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}_nonPARs.vcf.gz")
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


rule extract_X_PARs_QD:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/QD/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/QD/chrX_PARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}_PARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.ilist} """
        """-O {output} """


################
# Filter: QUAL #
################
rule hard_filter_variant_filtration_QUAL:
    input:
        vcf = os.path.join(config["proj_path"], "joint_called_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    params:
        chrms = "{chrms}",
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        qual = "{filter}",
        intermediatevcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/QUAL/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_QUAL_{filter}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/QUAL/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "QUAL < {params.qual}" --filter-name "QUAL{params.qual}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """

rule extract_X_nonPARs_QUAL:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/QUAL/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/QUAL/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}_nonPARs.vcf.gz")
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

rule extract_X_PARs_QUAL:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/QUAL/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/QUAL/chrX_PARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}_PARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.ilist} """
        """-O {output} """

###############
# Filter: SOR #
###############
rule hard_filter_variant_filtration_SOR:
    input:
        vcf = os.path.join(config["proj_path"], "joint_called_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    params:
        chrms = "{chrms}",
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        sor = "{filter}",
        intermediatevcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/SOR/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_SOR_{filter}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/SOR/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "SOR > {params.sor}" --filter-name "SOR{params.sor}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """

rule extract_X_nonPARs_SOR:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/SOR/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/SOR/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}_nonPARs.vcf.gz")
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

rule extract_X_PARs_SOR:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/SOR/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/SOR/chrX_PARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}_PARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.ilist} """
        """-O {output} """

##############
# Filter: FS #
##############
rule hard_filter_variant_filtration_FS:
    input:
        vcf = os.path.join(config["proj_path"], "joint_called_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    params:
        chrms = "{chrms}",
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        fs = "{filter}",
        intermediatevcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/FS/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_FS_{filter}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/FS/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "FS > {params.fs}" --filter-name "FS{params.fs}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """

rule extract_X_nonPARs_FS:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/FS/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/FS/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}_nonPARs.vcf.gz")
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

rule extract_X_PARs_FS:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/FS/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/FS/chrX_PARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}_PARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.ilist} """
        """-O {output} """

##############
# Filter: MQ #
##############
rule hard_filter_variant_filtration_MQ:
    input:
        vcf = os.path.join(config["proj_path"], "joint_called_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    params:
        chrms = "{chrms}",
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        mq = "{filter}",
        intermediatevcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/MQ/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_MQ_{filter}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/MQ/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "MQ < {params.mq}" --filter-name "MQ{params.mq}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """

rule extract_X_nonPARs_MQ:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/MQ/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/MQ/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}_nonPARs.vcf.gz")
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

rule extract_X_PARs_MQ:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/MQ/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/MQ/chrX_PARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}_PARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.ilist} """
        """-O {output} """

#####################
# Filter: MQRankSum #
#####################
rule hard_filter_variant_filtration_MQRankSum:
    input:
        vcf = os.path.join(config["proj_path"], "joint_called_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    params:
        chrms = "{chrms}",
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        mqranksum = "{filter}",
        intermediatevcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/MQRankSum/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_MQRankSum_{filter}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/MQRankSum/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "MQRankSum < {params.mqranksum}" --filter-name "MQRankSum{params.mqranksum}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """

rule extract_X_nonPARs_MQRankSum:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/MQRankSum/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/MQRankSum/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}_nonPARs.vcf.gz")
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

rule extract_X_PARs_MQRankSum:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/MQRankSum/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/MQRankSum/chrX_PARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}_PARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.ilist} """
        """-O {output} """

##########################
# Filter: ReadPosRankSum #
##########################
rule hard_filter_variant_filtration_ReadPosRankSum:
    input:
        vcf = os.path.join(config["proj_path"], "joint_called_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    params:
        chrms = "{chrms}",
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        readposranksum = "{filter}",
        intermediatevcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/ReadPosRankSum/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_ReadPosRankSum_{filter}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/ReadPosRankSum/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "ReadPosRankSum < {params.readposranksum}" --filter-name "ReadPosRankSum{params.readposranksum}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """

rule extract_X_nonPARs_ReadPosRankSum:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/ReadPosRankSum/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/ReadPosRankSum/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}_nonPARs.vcf.gz")
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

rule extract_X_PARs_ReadPosRankSum:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/ReadPosRankSum/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/ReadPosRankSum/chrX_PARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}_PARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.ilist} """
        """-O {output} """

##########################
# Filter: AN #
##########################
rule hard_filter_variant_filtration_AN:
    input:
        vcf = os.path.join(config["proj_path"], "joint_called_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    params:
        chrms = "{chrms}",
        AN = "{filter}"
    output:
        ovcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/AN/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}.vcf.gz"),
        ovcfidx = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/AN/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}.vcf.gz.tbi")
    shell:
        """
        bcftools filter -e'INFO/AN<{params.AN}' {input.vcf} -Oz -o {output.ovcf};
        tabix -p vcf {output.ovcf}
        """

rule extract_X_nonPARs_AN:
    input:
         ovcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/AN/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}.vcf.gz"),
         ovcfidx = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/AN/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}.vcf.gz.tbi")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/AN/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}_nonPARs.vcf.gz")
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

rule extract_X_PARs_AN:
    input:
         ovcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/AN/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}.vcf.gz"),
         ovcfidx = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/AN/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}.vcf.gz.tbi")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/AN/chrX_PARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}_PARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input.ovcf} """
        """-L {params.ilist} """
        """-O {output} """

##########################
# Filter: DP #
##########################
rule hard_filter_variant_filtration_DP:
    input:
        vcf = os.path.join(config["proj_path"], "joint_called_vcfs/AFR/males/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    params:
        chrms = "{chrms}",
        DP = "{filter}"
    output:
        ovcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/DP/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}.vcf.gz"),
        ovcfidx = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/DP/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}.vcf.gz.tbi")
    shell:
        """
        bcftools filter -e'FORMAT/DP<{params.DP}' -S . {input.vcf} -Oz -o {output.ovcf};
        tabix -p vcf {output.ovcf}
        """

rule extract_X_nonPARs_DP:
    input:
         ovcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/DP/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}.vcf.gz"),
         ovcfidx = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/DP/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}.vcf.gz.tbi")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/DP/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}_nonPARs.vcf.gz")
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

rule extract_X_PARs_DP:
    input:
         ovcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/DP/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}.vcf.gz"),
         ovcfidx = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/DP/{chrms}/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}.vcf.gz.tbi")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/DP/chrX_PARs/{chrms}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}_PARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input.ovcf} """
        """-L {params.ilist} """
        """-O {output} """



#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#                            MALES (DEFAULT)                                   #
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Step: Alignment (DEFAULT)
#------------------------------------------------------------------------------#
rule default_alignment_males:
    input:
        fq1 = os.path.join(config["proj_path"], "trimmed_fastqs/AFR/males/{sample}_NEAT_trimmed_read1.fastq.gz"),
        fq2 = os.path.join(config["proj_path"], "trimmed_fastqs/AFR/males/{sample}_NEAT_trimmed_read2.fastq.gz")
    output:
        os.path.join(config["proj_path"], "bams/AFR/males/default/{sample}_GRCh38_default_sorted.bam")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
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
rule default_index_bam_males:
    input:
        os.path.join(config["proj_path"], "bams/AFR/males/default/{sample}_GRCh38_default_sorted.bam")
    output:
        os.path.join(config["proj_path"], "bams/AFR/males/default/{sample}_GRCh38_default_sorted.bam.bai")
    shell:
        "samtools index {input}"


#------------------------------------------------------------------------------#
# Step: Mark duplicates
#------------------------------------------------------------------------------#
rule default_MarkDups_males:
    input:
        bam = os.path.join(config["proj_path"], "bams/AFR/males/default/{sample}_GRCh38_default_sorted.bam"),
        bai = os.path.join(config["proj_path"], "bams/AFR/males/default/{sample}_GRCh38_default_sorted.bam.bai")
    output:
        bam = os.path.join(config["proj_path"], "bams/AFR/males/default/{sample}_GRCh38_default_sorted_mkdups.bam"),
        metrics = os.path.join(config["proj_path"], "stats/AFR/males/default/{sample}_GRCh38_default_sorted_mkdups_metrics.txt")
    shell:
        """
        picard -Xmx14g MarkDuplicates I={input.bam} O={output.bam} M={output.metrics} VALIDATION_STRINGENCY=LENIENT
        """


#------------------------------------------------------------------------------#
# Step: Index bams from last step
#------------------------------------------------------------------------------#
rule default_index_MarkDups_bam_males:
    input:
        os.path.join(config["proj_path"], "bams/AFR/males/default/{sample}_GRCh38_default_sorted_mkdups.bam")
    output:
        os.path.join(config["proj_path"], "bams/AFR/males/default/{sample}_GRCh38_default_sorted_mkdups.bam.bai")
    shell:
        "samtools index {input}"


#------------------------------------------------------------------------------#
# Step: Call variants - diploid mode
#------------------------------------------------------------------------------#
# TO DO: Once scratch is back up:
# - run this step with EMIT_ALL_ACTIVE_SITES so that I can get ref sites that are
#   called. This will help me calculate true negatives

rule default_make_gvcfs_males_diploid:
    input:
        bam = os.path.join(config["proj_path"], "bams/AFR/males/default/{sample}_GRCh38_default_sorted_mkdups.bam"),
        bai = os.path.join(config["proj_path"], "bams/AFR/males/default/{sample}_GRCh38_default_sorted_mkdups.bam.bai")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        ploidy = 2,
        chrms = "{chrms}"
    output:
        gvcf = os.path.join(config["proj_path"], "gvcfs/AFR/males/default/{sample}_{chrms}_GRCh38_default_diploid.g.vcf.gz")
        #gvcf = os.path.join(config["proj_path"], "gvcfs/AFR/males/default/{sample}_GRCh38_default_diploid.g.vcf.gz")
    shell:
        """
        gatk --java-options '-Xmx4g' HaplotypeCaller -R {params.ref} -I {input.bam} -L {params.chrms} -ploidy {params.ploidy} -O {output.gvcf} -ERC GVCF
        """

# gatk --java-options '-Xmx4g' HaplotypeCaller -R {params.ref} -I {input.bam} -L {params.chrms} -ploidy {params.ploidy} -O {output.gvcf} -ERC GVCF



#------------------------------------------------------------------------------#
# Step: Combine gvcfs and joint genotype - all together, diploid
#------------------------------------------------------------------------------#
rule default_CombineGVCFs:
    input:
        gvcf =  expand(os.path.join(config["proj_path"], "gvcfs/AFR/males/default/{sample}_{chrms}_GRCh38_default_diploid.g.vcf.gz"), sample=config["males"], chrms=config["chromosomes"])
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        gvcfs = expand(("-V " + config["proj_path"] + "/gvcfs/AFR/males/default/{sample}_{chrms}_GRCh38_default_diploid.g.vcf.gz"), sample=config["males"], chrms=config["chromosomes"]),
        chrms = "{chrms}"
    output:
        os.path.join(config["proj_path"], "combine_g_vcfs/AFR/males/default/{chrms}_GRCh38_default_diploid_combined.g.vcf.gz")
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
        os.path.join(config["proj_path"], "combine_g_vcfs/AFR/males/default/{chrms}_GRCh38_default_diploid_combined.g.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        chrms = "{chrms}"
    output:
        os.path.join(config["proj_path"], "joint_called_vcfs/AFR/males/default/{chrms}_GRCh38_default_gatk_diploid_called_raw.vcf.gz")
    shell:
        """
        gatk GenotypeGVCFs -R {params.ref} -V {input} -L {params.chrms} -O {output}
        """

#------------------------------------------------------------------------------#
# Step: Extract SNPs
#------------------------------------------------------------------------------#
rule default_select_SNPs:
    input:
        os.path.join(config["proj_path"], "joint_called_vcfs/AFR/males/default/{chrms}_GRCh38_default_gatk_diploid_called_raw.vcf.gz")
    output:
        os.path.join(config["proj_path"], "joint_called_vcfs/AFR/males/default/{chrms}_GRCh38_default_gatk_diploid_called_raw_SNPs.vcf.gz")
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
# TO DO - Do each hard filter separately.
rule default_hard_filter_variant_filtration_all:
    input:
        vcf = os.path.join(config["proj_path"], "joint_called_vcfs/AFR/males/default/{chrms}_GRCh38_default_gatk_diploid_called_raw_SNPs.vcf.gz")
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
        intermediatevcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/default/all/{chrms}_GRCh38_default_gatk_diploid_called_raw_SNPs_{filtering_options}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/default/all/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "QD < {params.qd}" --filter-name "QD{params.qd}" -filter "QUAL < {params.qual}" --filter-name "QUAL{params.qual}" -filter "SOR > {params.sor}" --filter-name "SOR{params.sor}" -filter "FS > {params.fs}" --filter-name "FS{params.fs}" -filter "MQ < {params.mq}" --filter-name "MQ{params.mq}" -filter "MQRankSum < {params.mqranksum}" --filter-name "MQRankSum{params.mqranksum}" -filter "ReadPosRankSum < {params.readposranksum}" --filter-name "ReadPosRankSum{params.readposranksum}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """


#------------------------------------------------------------------------------#
# Step: For X chromosome, extract non PARS (exclude PARs)
#------------------------------------------------------------------------------#
rule default_extract_X_nonPARs:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/default/all/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered.vcf.gz")
    output:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/default/all/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered_nonPARs.vcf.gz")
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


#------------------------------------------------------------------------------#
# TODO: Step: For X chromosome, extract PARS (we will also analyze PARs)
#------------------------------------------------------------------------------#
rule default_extract_X_PARs:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/default/all/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered.vcf.gz")
    output:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/default/all/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered_PARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.ilist} """
        """-O {output} """

# I NEED TO CHECK TO SEE IF THERE ARE VARIANTS ON Y PARS? I AM NOT GETTING ANY
# VARIANTS ON X PARS...

#------------------------------------------------------------------------------#
# Step: For Y and X chromosomes (non PARS) get per sample stats
#------------------------------------------------------------------------------#
rule default_bcftools_PSCs_XnonPARs:
    input:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/default/all/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered_nonPARs.vcf.gz")
    output:
        os.path.join(config["proj_path"], "bcftools_stats/AFR/males/default/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered_nonPARs.bcftools.stats.PSC.txt")
    shell:
        """
        bcftools stats -s - {input} | grep PSC > {output}
        """

rule default_bcftools_PSCs_Y:
    input:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/males/default/all/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered.vcf.gz")
    output:
        os.path.join(config["proj_path"], "bcftools_stats/AFR/males/default/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered.bcftools.stats.PSC.txt")
    shell:
        """
        bcftools stats -s - {input} | grep PSC > {output}
        """


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#                                    FEMALES                                   #
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Step: fastqc and multiQC
#------------------------------------------------------------------------------#
rule fastqc_analysis_read1_females:
    input:
        fq1 = os.path.join(config["proj_path"], "simulations/AFR/females/{sample}/{sample}_NEAT_read1.fastq.gz"),
    output:
        f1 = os.path.join(config["proj_path"], "fastqc/AFR/females/{sample}_NEAT_read1_fastqc.html"),
    params:
        fastqc = os.path.join(config["proj_path"], "fastqc/AFR/females/")
    shell:
        """
        PERL5LIB="";
        fastqc -o {params.fastqc} {input.fq1}
        """

rule fastqc_analysis_read2_females:
    input:
        fq2 = os.path.join(config["proj_path"], "simulations/AFR/females/{sample}/{sample}_NEAT_read2.fastq.gz")
    output:
        f2 = os.path.join(config["proj_path"], "fastqc/AFR/females/{sample}_NEAT_read2_fastqc.html")
    params:
        fastqc = os.path.join(config["proj_path"], "fastqc/AFR/females/")
    shell:
        """
        PERL5LIB="";
        fastqc -o {params.fastqc} {input.fq2}
        """


rule multiqc_analysis_females:
    input:
        f1 = expand(os.path.join(config["proj_path"], "fastqc/AFR/females/{sample}_NEAT_read1_fastqc.html"),
        sample=config["females"]),
        f2 = expand(os.path.join(config["proj_path"], "fastqc/AFR/females/{sample}_NEAT_read2_fastqc.html"),
        sample=config["females"])
    output:
        os.path.join(config["proj_path"], "multiqc/AFR/females/multiqc_report.html")
    params:
        fastqc = os.path.join(config["proj_path"], "fastqc/AFR/females/"),
        multiqc = os.path.join(config["proj_path"], "multiqc/AFR/females/")
    shell:
        "export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 &&"
        "multiqc --interactive -o {params.multiqc} {params.fastqc}"


#------------------------------------------------------------------------------#
# Step: Trimming
#------------------------------------------------------------------------------#
rule trim_bbduk_females:
    input:
        fq1 = os.path.join(config["proj_path"], "simulations/AFR/females/{sample}/{sample}_NEAT_read1.fastq.gz"),
        fq2 = os.path.join(config["proj_path"], "simulations/AFR/females/{sample}/{sample}_NEAT_read2.fastq.gz")
    output:
        f1 = os.path.join(config["proj_path"], "trimmed_fastqs/AFR/females/{sample}_NEAT_trimmed_read1.fastq.gz"),
        f2 = os.path.join(config["proj_path"], "trimmed_fastqs/AFR/females/{sample}_NEAT_trimmed_read2.fastq.gz")
    shell:
        "bbduk.sh -Xmx3g in1={input.fq1} in2={input.fq2} "
        "out1={output.f1} out2={output.f2} "
        "qtrim=rl trimq=20 minlen=75 maq=20"


#------------------------------------------------------------------------------#
# Step: fastqc and multiQC
#------------------------------------------------------------------------------#
rule fastqc_analysis_trimmed_1_females:
    input:
        fq1 = os.path.join(config["proj_path"], "trimmed_fastqs/AFR/females/{sample}_NEAT_trimmed_read1.fastq.gz"),
    output:
        f1 = os.path.join(config["proj_path"], "fastqc_trimmed/AFR/females/{sample}_NEAT_trimmed_read1_fastqc.html"),
    params:
        fastqc = os.path.join(config["proj_path"], "fastqc_trimmed/AFR/females/")
    shell:
        """
        PERL5LIB="";
        fastqc -o {params.fastqc} {input.fq1}
        """

rule fastqc_analysis_trimmed_2_females:
    input:
        fq2 = os.path.join(config["proj_path"], "trimmed_fastqs/AFR/females/{sample}_NEAT_trimmed_read2.fastq.gz")
    output:
        f2 = os.path.join(config["proj_path"], "fastqc_trimmed/AFR/females/{sample}_NEAT_trimmed_read2_fastqc.html")
    params:
        fastqc = os.path.join(config["proj_path"], "fastqc_trimmed/AFR/females/")
    shell:
        """
        PERL5LIB="";
        fastqc -o {params.fastqc} {input.fq2}
        """

rule multiqc_analysis_trimmed_females:
    input:
        f1 = expand(os.path.join(config["proj_path"], "fastqc_trimmed/AFR/females/{sample}_NEAT_trimmed_read1_fastqc.html"),
        sample=config["females"]),
        f2 = expand(os.path.join(config["proj_path"], "fastqc_trimmed/AFR/females/{sample}_NEAT_trimmed_read2_fastqc.html"),
        sample=config["females"])
    output:
        os.path.join(config["proj_path"], "multiqc_trimmed/AFR/females/multiqc_report.html")
    params:
        fastqc = os.path.join(config["proj_path"], "fastqc_trimmed/AFR/females/"),
        multiqc = os.path.join(config["proj_path"], "multiqc_trimmed/AFR/females/")
    shell:
        "export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 &&"
        "multiqc --interactive -o {params.multiqc} {params.fastqc}"


#------------------------------------------------------------------------------#
# Step: Alignment (SCC)
#------------------------------------------------------------------------------#
# CHANGE REFERENCE GENOME TO Y hard masked!!!
rule alignment_SCC_females_females:
    input:
        fq1 = os.path.join(config["proj_path"], "trimmed_fastqs/AFR/females/{sample}_NEAT_trimmed_read1.fastq.gz"),
        fq2 = os.path.join(config["proj_path"], "trimmed_fastqs/AFR/females/{sample}_NEAT_trimmed_read2.fastq.gz")
    output:
        os.path.join(config["proj_path"], "bams/AFR/females/{sample}_GRCh38_YHardMasked_sorted.bam")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_HardMasked"],
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
rule index_bam_females_females:
    input:
        os.path.join(config["proj_path"], "bams/AFR/females/{sample}_GRCh38_YHardMasked_sorted.bam")
    output:
        os.path.join(config["proj_path"], "bams/AFR/females/{sample}_GRCh38_YHardMasked_sorted.bam.bai")
    shell:
        "samtools index {input}"


#------------------------------------------------------------------------------#
# Step: Mark duplicates
#------------------------------------------------------------------------------#
rule MarkDups_females_females:
    input:
        bam = os.path.join(config["proj_path"], "bams/AFR/females/{sample}_GRCh38_YHardMasked_sorted.bam"),
        bai = os.path.join(config["proj_path"], "bams/AFR/females/{sample}_GRCh38_YHardMasked_sorted.bam.bai")
    output:
        bam = os.path.join(config["proj_path"], "bams/AFR/females/{sample}_GRCh38_YHardMasked_sorted_mkdups.bam"),
        metrics = os.path.join(config["proj_path"], "stats/AFR/females/{sample}_GRCh38_YHardMasked_sorted_mkdups_metrics.txt")
    shell:
        """
        picard -Xmx14g MarkDuplicates I={input.bam} O={output.bam} M={output.metrics} VALIDATION_STRINGENCY=LENIENT
        """


#------------------------------------------------------------------------------#
# Step: Index bams from last step
#------------------------------------------------------------------------------#
rule index_MarkDups_bam_females_females:
    input:
        os.path.join(config["proj_path"], "bams/AFR/females/{sample}_GRCh38_YHardMasked_sorted_mkdups.bam")
    output:
        os.path.join(config["proj_path"], "bams/AFR/females/{sample}_GRCh38_YHardMasked_sorted_mkdups.bam.bai")
    shell:
        "samtools index {input}"


#------------------------------------------------------------------------------#
# Step: Call variants - diploid mode
#------------------------------------------------------------------------------#
# TO DO: Once scratch is back up:
# - run this step with EMIT_ALL_ACTIVE_SITES so that I can get ref sites that are
#   called. This will help me calculate true negatives

rule make_gvcfs_females_diploid_females:
    input:
        bam = os.path.join(config["proj_path"], "bams/AFR/females/{sample}_GRCh38_YHardMasked_sorted_mkdups.bam"),
        bai = os.path.join(config["proj_path"], "bams/AFR/females/{sample}_GRCh38_YHardMasked_sorted_mkdups.bam.bai")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_HardMasked"],
        ploidy = 2,
        chrms = "{chrms}"
    output:
        gvcf = os.path.join(config["proj_path"], "gvcfs/AFR/females/{sample}_{chrms}_GRCh38_YHardMasked_diploid.g.vcf.gz")
        #gvcf = os.path.join(config["proj_path"], "gvcfs/AFR/females/{sample}_GRCh38_YHardMasked_diploid.g.vcf.gz")
    shell:
        """
        gatk --java-options '-Xmx4g' HaplotypeCaller -R {params.ref} -I {input.bam} -L {params.chrms} -ploidy {params.ploidy} -O {output.gvcf} -ERC GVCF
        """



#------------------------------------------------------------------------------#
# Step: Combine gvcfs and joint genotype - all together, diploid
#------------------------------------------------------------------------------#
rule CombineGVCFs_females:
    input:
        gvcf =  expand(os.path.join(config["proj_path"], "gvcfs/AFR/females/{sample}_{chrms}_GRCh38_YHardMasked_diploid.g.vcf.gz"), sample=config["females"], chrms=config["chromosomes"])
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_HardMasked"],
        gvcfs = expand(("-V " + config["proj_path"] + "/gvcfs/AFR/females/{sample}_{chrms}_GRCh38_YHardMasked_diploid.g.vcf.gz"), sample=config["females"], chrms=config["chromosomes"]),
        chrms = "{chrms}"
    output:
        os.path.join(config["proj_path"], "combine_g_vcfs/AFR/females/{chrms}_GRCh38_YHardMasked_diploid_combined.g.vcf.gz")
    shell:
        """
        gatk CombineGVCFs -R {params.ref} {params.gvcfs} -L {params.chrms} -O {output}
        """


# this sites: https://pixy.readthedocs.io/en/latest/generating_invar/generating_invar.html
# suggests that you can use a flag called -all-sites
# https://gatk.broadinstitute.org/hc/en-us/articles/4404607598875-GenotypeGVCFs
# --include-non-variant-sites / -all-sites. Include loci found to be non-variant
# after genotyping. boolean  false
rule GenotypeGVCFs_females:
    input:
        os.path.join(config["proj_path"], "combine_g_vcfs/AFR/females/{chrms}_GRCh38_YHardMasked_diploid_combined.g.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_HardMasked"],
        chrms = "{chrms}"
    output:
        os.path.join(config["proj_path"], "joint_called_vcfs/AFR/females/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_raw.vcf.gz")
    shell:
        """
        gatk GenotypeGVCFs -R {params.ref} -V {input} -L {params.chrms} -O {output}
        """

#------------------------------------------------------------------------------#
# Step: Extract SNPs
#------------------------------------------------------------------------------#
rule select_SNPs_females:
    input:
        os.path.join(config["proj_path"], "joint_called_vcfs/AFR/females/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_raw.vcf.gz")
    output:
        os.path.join(config["proj_path"], "joint_called_vcfs/AFR/females/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_HardMasked"],
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
# TO DO - Do each hard filter separately.
rule hard_filter_variant_filtration_all_females:
    input:
        vcf = os.path.join(config["proj_path"], "joint_called_vcfs/AFR/females/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    params:
        chrms = "{chrms}",
        ref = config["genome_paths"]["Ref_GRCh38_Y_HardMasked"],
        qd = lambda wildcards: config[wildcards.filtering_options]["QD"],
        qual = lambda wildcards: config[wildcards.filtering_options]["QUAL"],
        sor = lambda wildcards: config[wildcards.filtering_options]["SOR"],
        fs = lambda wildcards: config[wildcards.filtering_options]["FS"],
        mq = lambda wildcards: config[wildcards.filtering_options]["MQ"],
        mqranksum = lambda wildcards: config[wildcards.filtering_options]["MQRankSum"],
        readposranksum = lambda wildcards: config[wildcards.filtering_options]["ReadPosRankSum"],
        intermediatevcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/all/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs_{filtering_options}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/all/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "QD < {params.qd}" --filter-name "QD{params.qd}" -filter "QUAL < {params.qual}" --filter-name "QUAL{params.qual}" -filter "SOR > {params.sor}" --filter-name "SOR{params.sor}" -filter "FS > {params.fs}" --filter-name "FS{params.fs}" -filter "MQ < {params.mq}" --filter-name "MQ{params.mq}" -filter "MQRankSum < {params.mqranksum}" --filter-name "MQRankSum{params.mqranksum}" -filter "ReadPosRankSum < {params.readposranksum}" --filter-name "ReadPosRankSum{params.readposranksum}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """


#------------------------------------------------------------------------------#
# Step: For X chromosome, extract non PARS (exclude PARs)
#------------------------------------------------------------------------------#
rule extract_X_nonPARs_females:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/all/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered.vcf.gz")
    output:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/all/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered_nonPARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_HardMasked"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.chrms} -XL {params.ilist} """
        """-O {output} """


#------------------------------------------------------------------------------#
# TODO: Step: For X chromosome, extract PARS (we will also analyze PARs)
#------------------------------------------------------------------------------#
rule extract_X_PARs_females:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/all/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered.vcf.gz")
    output:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/all/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered_PARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_HardMasked"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.ilist} """
        """-O {output} """


#------------------------------------------------------------------------------#
# Step: For Y and X chromosomes (non PARS) get per sample stats
#------------------------------------------------------------------------------#

rule bcftools_PSCs_Y_females:
    input:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/all/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered.vcf.gz")
    output:
        os.path.join(config["proj_path"], "bcftools_stats/AFR/females/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_{filtering_options}_filtered.bcftools.stats.PSC.txt")
    shell:
        """
        bcftools stats -s - {input} | grep PSC > {output}
        """

'''
#------------------------------------------------------------------------------#
# Step: Each filter separately
#------------------------------------------------------------------------------#
##############
# Filter: QD #
##############
rule female_hard_filter_variant_filtration_QD:
    input:
        vcf = os.path.join(config["proj_path"], "joint_called_vcfs/AFR/females/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    params:
        chrms = "{chrms}",
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        qd = "{filter}",
        intermediatevcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/QD/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs_QD_{filter}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/QD/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "QD < {params.qd}" --filter-name "QD{params.qd}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """

rule female_extract_X_nonPARs_QD:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/QD/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/QD/chrX_nonPARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}_nonPARs.vcf.gz")
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


rule female_extract_X_PARs_QD:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/QD/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/QD/chrX_PARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}_PARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.ilist} """
        """-O {output} """


################
# Filter: QUAL #
################
rule female_hard_filter_variant_filtration_QUAL:
    input:
        vcf = os.path.join(config["proj_path"], "joint_called_vcfs/AFR/females/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    params:
        chrms = "{chrms}",
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        qual = "{filter}",
        intermediatevcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/QUAL/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs_QUAL_{filter}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/QUAL/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "QUAL < {params.qual}" --filter-name "QUAL{params.qual}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """

rule female_extract_X_nonPARs_QUAL:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/QUAL/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/QUAL/chrX_nonPARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}_nonPARs.vcf.gz")
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

rule female_extract_X_PARs_QUAL:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/QUAL/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/QUAL/chrX_PARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}_PARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.ilist} """
        """-O {output} """

###############
# Filter: SOR #
###############
rule female_hard_filter_variant_filtration_SOR:
    input:
        vcf = os.path.join(config["proj_path"], "joint_called_vcfs/AFR/females/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    params:
        chrms = "{chrms}",
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        sor = "{filter}",
        intermediatevcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/SOR/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs_SOR_{filter}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/SOR/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "SOR > {params.sor}" --filter-name "SOR{params.sor}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """

rule female_extract_X_nonPARs_SOR:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/SOR/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/SOR/chrX_nonPARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}_nonPARs.vcf.gz")
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

rule female_extract_X_PARs_SOR:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/SOR/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/SOR/chrX_PARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}_PARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.ilist} """
        """-O {output} """

##############
# Filter: FS #
##############
rule female_hard_filter_variant_filtration_FS:
    input:
        vcf = os.path.join(config["proj_path"], "joint_called_vcfs/AFR/females/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    params:
        chrms = "{chrms}",
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        fs = "{filter}",
        intermediatevcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/FS/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs_FS_{filter}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/FS/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "FS > {params.fs}" --filter-name "FS{params.fs}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """

rule female_extract_X_nonPARs_FS:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/FS/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/FS/chrX_nonPARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}_nonPARs.vcf.gz")
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

rule female_extract_X_PARs_FS:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/FS/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/FS/chrX_PARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}_PARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.ilist} """
        """-O {output} """

##############
# Filter: MQ #
##############
rule female_hard_filter_variant_filtration_MQ:
    input:
        vcf = os.path.join(config["proj_path"], "joint_called_vcfs/AFR/females/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    params:
        chrms = "{chrms}",
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        mq = "{filter}",
        intermediatevcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/MQ/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs_MQ_{filter}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/MQ/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "MQ < {params.mq}" --filter-name "MQ{params.mq}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """

rule female_extract_X_nonPARs_MQ:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/MQ/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/MQ/chrX_nonPARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}_nonPARs.vcf.gz")
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

rule female_extract_X_PARs_MQ:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/MQ/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/MQ/chrX_PARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}_PARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.ilist} """
        """-O {output} """

#####################
# Filter: MQRankSum #
#####################
rule female_hard_filter_variant_filtration_MQRankSum:
    input:
        vcf = os.path.join(config["proj_path"], "joint_called_vcfs/AFR/females/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    params:
        chrms = "{chrms}",
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        mqranksum = "{filter}",
        intermediatevcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/MQRankSum/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs_MQRankSum_{filter}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/MQRankSum/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "MQRankSum < {params.mqranksum}" --filter-name "MQRankSum{params.mqranksum}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """

rule female_extract_X_nonPARs_MQRankSum:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/MQRankSum/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/MQRankSum/chrX_nonPARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}_nonPARs.vcf.gz")
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

rule female_extract_X_PARs_MQRankSum:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/MQRankSum/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/MQRankSum/chrX_PARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}_PARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.ilist} """
        """-O {output} """

##########################
# Filter: ReadPosRankSum #
##########################
rule female_hard_filter_variant_filtration_ReadPosRankSum:
    input:
        vcf = os.path.join(config["proj_path"], "joint_called_vcfs/AFR/females/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    params:
        chrms = "{chrms}",
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        readposranksum = "{filter}",
        intermediatevcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/ReadPosRankSum/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs_ReadPosRankSum_{filter}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/ReadPosRankSum/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "ReadPosRankSum < {params.readposranksum}" --filter-name "ReadPosRankSum{params.readposranksum}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """

rule female_extract_X_nonPARs_ReadPosRankSum:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/ReadPosRankSum/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/ReadPosRankSum/chrX_nonPARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}_nonPARs.vcf.gz")
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

rule female_extract_X_PARs_ReadPosRankSum:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/ReadPosRankSum/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}.vcf.gz")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/ReadPosRankSum/chrX_PARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}_PARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.ilist} """
        """-O {output} """

##########################
# Filter: AN #
##########################
rule female_hard_filter_variant_filtration_AN:
    input:
        vcf = os.path.join(config["proj_path"], "joint_called_vcfs/AFR/females/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    params:
        chrms = "{chrms}",
        AN = "{filter}"
    output:
        ovcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/AN/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}.vcf.gz"),
        ovcfidx = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/AN/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}.vcf.gz.tbi")
    shell:
        """
        bcftools filter -e'INFO/AN<{params.AN}' {input.vcf} -Oz -o {output.ovcf};
        tabix -p vcf {output.ovcf}
        """

rule female_extract_X_nonPARs_AN:
    input:
         ovcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/AN/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}.vcf.gz"),
         ovcfidx = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/AN/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}.vcf.gz.tbi")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/AN/chrX_nonPARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}_nonPARs.vcf.gz")
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

rule female_extract_X_PARs_AN:
    input:
         ovcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/AN/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}.vcf.gz"),
         ovcfidx = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/AN/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}.vcf.gz.tbi")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/AN/chrX_PARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}_PARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input.ovcf} """
        """-L {params.ilist} """
        """-O {output} """

##########################
# Filter: DP #
##########################
rule female_hard_filter_variant_filtration_DP:
    input:
        vcf = os.path.join(config["proj_path"], "joint_called_vcfs/AFR/females/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    params:
        chrms = "{chrms}",
        DP = "{filter}"
    output:
        ovcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/DP/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}.vcf.gz"),
        ovcfidx = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/DP/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}.vcf.gz.tbi")
    shell:
        """
        bcftools filter -e'FORMAT/DP<{params.DP}' -S . {input.vcf} -Oz -o {output.ovcf};
        tabix -p vcf {output.ovcf}
        """

rule female_extract_X_nonPARs_DP:
    input:
         ovcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/DP/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}.vcf.gz"),
         ovcfidx = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/DP/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}.vcf.gz.tbi")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/DP/chrX_nonPARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}_nonPARs.vcf.gz")
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

rule female_extract_X_PARs_DP:
    input:
         ovcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/DP/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}.vcf.gz"),
         ovcfidx = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/DP/{chrms}/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}.vcf.gz.tbi")
    output:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/DP/chrX_PARs/{chrms}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}_PARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input.ovcf} """
        """-L {params.ilist} """
        """-O {output} """

'''



#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Also align FEMALE data to default ref to test between scc and default on X
#                                  DEFAULT                                     #
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Step: Alignment (DEFAULT)
#------------------------------------------------------------------------------#
rule default_alignment_females:
    input:
        fq1 = os.path.join(config["proj_path"], "trimmed_fastqs/AFR/females/{sample}_NEAT_trimmed_read1.fastq.gz"),
        fq2 = os.path.join(config["proj_path"], "trimmed_fastqs/AFR/females/{sample}_NEAT_trimmed_read2.fastq.gz")
    output:
        os.path.join(config["proj_path"], "bams/AFR/females/default/{sample}_GRCh38_default_sorted.bam")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
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
rule index_bam_females_default:
    input:
        os.path.join(config["proj_path"], "bams/AFR/females/default/{sample}_GRCh38_default_sorted.bam")
    output:
        os.path.join(config["proj_path"], "bams/AFR/females/default/{sample}_GRCh38_default_sorted.bam.bai")
    shell:
        "samtools index {input}"


#------------------------------------------------------------------------------#
# Step: Mark duplicates
#------------------------------------------------------------------------------#
rule MarkDups_females_default:
    input:
        bam = os.path.join(config["proj_path"], "bams/AFR/females/default/{sample}_GRCh38_default_sorted.bam"),
        bai = os.path.join(config["proj_path"], "bams/AFR/females/default/{sample}_GRCh38_default_sorted.bam.bai")
    output:
        bam = os.path.join(config["proj_path"], "bams/AFR/females/default/{sample}_GRCh38_default_sorted_mkdups.bam"),
        metrics = os.path.join(config["proj_path"], "stats/AFR/females/default/{sample}_GRCh38_default_sorted_mkdups_metrics.txt")
    shell:
        """
        picard -Xmx14g MarkDuplicates I={input.bam} O={output.bam} M={output.metrics} VALIDATION_STRINGENCY=LENIENT
        """


#------------------------------------------------------------------------------#
# Step: Index bams from last step
#------------------------------------------------------------------------------#
rule index_MarkDups_bam_females_default:
    input:
        os.path.join(config["proj_path"], "bams/AFR/females/default/{sample}_GRCh38_default_sorted_mkdups.bam")
    output:
        os.path.join(config["proj_path"], "bams/AFR/females/default/{sample}_GRCh38_default_sorted_mkdups.bam.bai")
    shell:
        "samtools index {input}"


#------------------------------------------------------------------------------#
# Step: Call variants - diploid mode
#------------------------------------------------------------------------------#
# TO DO: Once scratch is back up:
# - run this step with EMIT_ALL_ACTIVE_SITES so that I can get ref sites that are
#   called. This will help me calculate true negatives

rule make_gvcfs_females_diploid_default:
    input:
        bam = os.path.join(config["proj_path"], "bams/AFR/females/default/{sample}_GRCh38_default_sorted_mkdups.bam"),
        bai = os.path.join(config["proj_path"], "bams/AFR/females/default/{sample}_GRCh38_default_sorted_mkdups.bam.bai")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        ploidy = 2,
        chrms = "{chrms}"
    output:
        gvcf = os.path.join(config["proj_path"], "gvcfs/AFR/females/default/{sample}_{chrms}_GRCh38_default_diploid.g.vcf.gz")
        #gvcf = os.path.join(config["proj_path"], "gvcfs/AFR/females/default/{sample}_GRCh38_default_diploid.g.vcf.gz")
    shell:
        """
        gatk --java-options '-Xmx4g' HaplotypeCaller -R {params.ref} -I {input.bam} -L {params.chrms} -ploidy {params.ploidy} -O {output.gvcf} -ERC GVCF
        """



#------------------------------------------------------------------------------#
# Step: Combine gvcfs and joint genotype - all together, diploid
#------------------------------------------------------------------------------#
rule CombineGVCFs_females_default:
    input:
        gvcf =  expand(os.path.join(config["proj_path"], "gvcfs/AFR/females/default/{sample}_{chrms}_GRCh38_default_diploid.g.vcf.gz"), sample=config["females"], chrms=config["chromosomes"])
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        gvcfs = expand(("-V " + config["proj_path"] + "/gvcfs/AFR/females/default/{sample}_{chrms}_GRCh38_default_diploid.g.vcf.gz"), sample=config["females"], chrms=config["chromosomes"]),
        chrms = "{chrms}"
    output:
        os.path.join(config["proj_path"], "combine_g_vcfs/AFR/females/default/{chrms}_GRCh38_default_diploid_combined.g.vcf.gz")
    shell:
        """
        gatk CombineGVCFs -R {params.ref} {params.gvcfs} -L {params.chrms} -O {output}
        """


# this sites: https://pixy.readthedocs.io/en/latest/generating_invar/generating_invar.html
# suggests that you can use a flag called -all-sites
# https://gatk.broadinstitute.org/hc/en-us/articles/4404607598875-GenotypeGVCFs
# --include-non-variant-sites / -all-sites. Include loci found to be non-variant
# after genotyping. boolean  false
rule GenotypeGVCFs_females_default:
    input:
        os.path.join(config["proj_path"], "combine_g_vcfs/AFR/females/default/{chrms}_GRCh38_default_diploid_combined.g.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        chrms = "{chrms}"
    output:
        os.path.join(config["proj_path"], "joint_called_vcfs/AFR/females/default/{chrms}_GRCh38_default_gatk_diploid_called_raw.vcf.gz")
    shell:
        """
        gatk GenotypeGVCFs -R {params.ref} -V {input} -L {params.chrms} -O {output}
        """

#------------------------------------------------------------------------------#
# Step: Extract SNPs
#------------------------------------------------------------------------------#
rule select_SNPs_females_default:
    input:
        os.path.join(config["proj_path"], "joint_called_vcfs/AFR/females/default/{chrms}_GRCh38_default_gatk_diploid_called_raw.vcf.gz")
    output:
        os.path.join(config["proj_path"], "joint_called_vcfs/AFR/females/default/{chrms}_GRCh38_default_gatk_diploid_called_raw_SNPs.vcf.gz")
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
# TO DO - Do each hard filter separately.
rule hard_filter_variant_filtration_all_females_default:
    input:
        vcf = os.path.join(config["proj_path"], "joint_called_vcfs/AFR/females/default/{chrms}_GRCh38_default_gatk_diploid_called_raw_SNPs.vcf.gz")
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
        intermediatevcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/default/all/{chrms}_GRCh38_default_gatk_diploid_called_raw_SNPs_{filtering_options}.vcf.gz")
    output:
        ovcf = os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/default/all/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "QD < {params.qd}" --filter-name "QD{params.qd}" -filter "QUAL < {params.qual}" --filter-name "QUAL{params.qual}" -filter "SOR > {params.sor}" --filter-name "SOR{params.sor}" -filter "FS > {params.fs}" --filter-name "FS{params.fs}" -filter "MQ < {params.mq}" --filter-name "MQ{params.mq}" -filter "MQRankSum < {params.mqranksum}" --filter-name "MQRankSum{params.mqranksum}" -filter "ReadPosRankSum < {params.readposranksum}" --filter-name "ReadPosRankSum{params.readposranksum}" -O {params.intermediatevcf};
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf} -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """


#------------------------------------------------------------------------------#
# Step: For X chromosome, extract non PARS (exclude PARs)
#------------------------------------------------------------------------------#
rule extract_X_nonPARs_females_default:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/default/all/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered.vcf.gz")
    output:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/default/all/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered_nonPARs.vcf.gz")
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


#------------------------------------------------------------------------------#
# TODO: Step: For X chromosome, extract PARS (we will also analyze PARs)
#------------------------------------------------------------------------------#
rule extract_X_PARs_females_default:
    input:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/default/all/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered.vcf.gz")
    output:
         os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/default/all/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered_PARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Default"],
        chrms = "{chrms}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.ilist} """
        """-O {output} """


#------------------------------------------------------------------------------#
# Step: For Y and X chromosomes (non PARS) get per sample stats
#------------------------------------------------------------------------------#
rule bcftools_PSCs_Y_females_default:
    input:
        os.path.join(config["proj_path"], "hard_filtered_vcfs/AFR/females/default/all/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered.vcf.gz")
    output:
        os.path.join(config["proj_path"], "bcftools_stats/AFR/females/default/{chrms}_GRCh38_default_gatk_diploid_called_SNPs_{filtering_options}_filtered.bcftools.stats.PSC.txt")
    shell:
        """
        bcftools stats -s - {input} | grep PSC > {output}
        """
