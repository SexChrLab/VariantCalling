import os

# Environment: variant_calling_simulations_project

configfile: "ASN_xy_var_calling_pilot_04_haploid.config.json"

wildcard_constraints:
    chrmsXnonPARs= '|'.join([re.escape(x) for x in config["chrX"]]),
    chrmsY= '|'.join([re.escape(x) for x in config["chrY"]]),
    chrmsM= '|'.join([re.escape(x) for x in config["chrM"]])


rule all:
    input:
        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf"), sample=config["samples"], chrmsXnonPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrmsY}.recode.vcf"), sample=config["samples"], chrmsY=config["chrY"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrmsM}.recode.vcf"), sample=config["samples"], chrmsM=config["chrM"]),

        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrms=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/chrY/{sample}_{chrms}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrms=config["chrY"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/chrM/{sample}_{chrmsM}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrmsM=config["chrM"]),

        # MALES - DEFAULT #
        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf"), sample=config["samples"], chrmsXnonPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrmsY}.recode.vcf"), sample=config["samples"], chrmsY=config["chrY"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrmsM}.recode.vcf"), sample=config["samples"], chrmsM=config["chrM"]),

        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/default/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrms=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/default/chrY/{sample}_{chrms}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrms=config["chrY"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/default/chrM/{sample}_{chrmsM}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrmsM=config["chrM"]),

'''
        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/QD/chrX_nonPARs/by_sample/{sample}_{chrms}_QD_{filter}_nonPARs.recode.vcf"), sample=config["samples"], chrms=config["chrX"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/QD/chrY/by_sample/{sample}_{chrms}_QD_{filter}.recode.vcf"), sample=config["samples"], chrms=config["chrY"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/QD/chrX_nonPARs/{sample}_{chrms}_nonPARs_QD_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrms=config["chrX"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/QD/chrY/{sample}_{chrms}_QD_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrms=config["chrY"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/QD/chrM/by_sample/{sample}_{chrmsM}_QD_{filter}.recode.vcf"), sample=config["samples"], chrmsM=config["chrM"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/QD/chrM/{sample}_{chrmsM}_QD_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrmsM=config["chrM"], filter=config["QD_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/QUAL/chrX_nonPARs/by_sample/{sample}_{chrms}_QUAL_{filter}_nonPARs.recode.vcf"), sample=config["samples"], chrms=config["chrX"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/QUAL/chrY/by_sample/{sample}_{chrms}_QUAL_{filter}.recode.vcf"), sample=config["samples"], chrms=config["chrY"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/QUAL/chrX_nonPARs/{sample}_{chrms}_nonPARs_QUAL_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrms=config["chrX"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/QUAL/chrY/{sample}_{chrms}_QUAL_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrms=config["chrY"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/QUAL/chrM/by_sample/{sample}_{chrmsM}_QUAL_{filter}.recode.vcf"), sample=config["samples"], chrmsM=config["chrM"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/QUAL/chrM/{sample}_{chrmsM}_QUAL_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrmsM=config["chrM"], filter=config["QUAL_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/SOR/chrX_nonPARs/by_sample/{sample}_{chrms}_SOR_{filter}_nonPARs.recode.vcf"), sample=config["samples"], chrms=config["chrX"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/SOR/chrY/by_sample/{sample}_{chrms}_SOR_{filter}.recode.vcf"), sample=config["samples"], chrms=config["chrY"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/SOR/chrX_nonPARs/{sample}_{chrms}_nonPARs_SOR_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrms=config["chrX"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/SOR/chrY/{sample}_{chrms}_SOR_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrms=config["chrY"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/SOR/chrM/by_sample/{sample}_{chrmsM}_SOR_{filter}.recode.vcf"), sample=config["samples"], chrmsM=config["chrM"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/SOR/chrM/{sample}_{chrmsM}_SOR_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrmsM=config["chrM"], filter=config["SOR_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/FS/chrX_nonPARs/by_sample/{sample}_{chrms}_FS_{filter}_nonPARs.recode.vcf"), sample=config["samples"], chrms=config["chrX"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/FS/chrY/by_sample/{sample}_{chrms}_FS_{filter}.recode.vcf"), sample=config["samples"], chrms=config["chrY"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/FS/chrX_nonPARs/{sample}_{chrms}_nonPARs_FS_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrms=config["chrX"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/FS/chrY/{sample}_{chrms}_FS_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrms=config["chrY"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/FS/chrM/by_sample/{sample}_{chrmsM}_FS_{filter}.recode.vcf"), sample=config["samples"], chrmsM=config["chrM"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/FS/chrM/{sample}_{chrmsM}_FS_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrmsM=config["chrM"], filter=config["FS_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/MQ/chrX_nonPARs/by_sample/{sample}_{chrms}_MQ_{filter}_nonPARs.recode.vcf"), sample=config["samples"], chrms=config["chrX"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/MQ/chrY/by_sample/{sample}_{chrms}_MQ_{filter}.recode.vcf"), sample=config["samples"], chrms=config["chrY"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/MQ/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQ_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrms=config["chrX"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/MQ/chrY/{sample}_{chrms}_MQ_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrms=config["chrY"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/MQ/chrM/by_sample/{sample}_{chrmsM}_MQ_{filter}.recode.vcf"), sample=config["samples"], chrmsM=config["chrM"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/MQ/chrM/{sample}_{chrmsM}_MQ_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrmsM=config["chrM"], filter=config["MQ_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/MQRankSum/chrX_nonPARs/by_sample/{sample}_{chrms}_MQRankSum_{filter}_nonPARs.recode.vcf"), sample=config["samples"], chrms=config["chrX"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/MQRankSum/chrY/by_sample/{sample}_{chrms}_MQRankSum_{filter}.recode.vcf"), sample=config["samples"], chrms=config["chrY"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/MQRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrms=config["chrX"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/MQRankSum/chrY/{sample}_{chrms}_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrms=config["chrY"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/MQRankSum/chrM/by_sample/{sample}_{chrmsM}_MQRankSum_{filter}.recode.vcf"), sample=config["samples"], chrmsM=config["chrM"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/MQRankSum/chrM/{sample}_{chrmsM}_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrmsM=config["chrM"], filter=config["MQRankSum_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrX_nonPARs/by_sample/{sample}_{chrms}_ReadPosRankSum_{filter}_nonPARs.recode.vcf"), sample=config["samples"], chrms=config["chrX"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrY/by_sample/{sample}_{chrms}_ReadPosRankSum_{filter}.recode.vcf"), sample=config["samples"], chrms=config["chrY"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/ReadPosRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrms=config["chrX"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/ReadPosRankSum/chrY/{sample}_{chrms}_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrms=config["chrY"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrM/by_sample/{sample}_{chrmsM}_ReadPosRankSum_{filter}.recode.vcf"), sample=config["samples"], chrmsM=config["chrM"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/ReadPosRankSum/chrM/{sample}_{chrmsM}_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrmsM=config["chrM"], filter=config["ReadPosRankSum_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/DP/chrX_nonPARs/by_sample/{sample}_{chrms}_DP_{filter}_nonPARs.recode.vcf"), sample=config["samples"], chrms=config["chrX"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/DP/chrY/by_sample/{sample}_{chrms}_DP_{filter}.recode.vcf"), sample=config["samples"], chrms=config["chrY"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/DP/chrX_nonPARs/{sample}_{chrms}_nonPARs_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrms=config["chrX"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/DP/chrY/{sample}_{chrms}_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrms=config["chrY"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/DP/chrM/by_sample/{sample}_{chrmsM}_DP_{filter}.recode.vcf"), sample=config["samples"], chrmsM=config["chrM"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/DP/chrM/{sample}_{chrmsM}_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrmsM=config["chrM"], filter=config["DP_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/AN/chrX_nonPARs/by_sample/{sample}_{chrms}_AN_{filter}_nonPARs.recode.vcf"), sample=config["samples"], chrms=config["chrX"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/AN/chrY/by_sample/{sample}_{chrms}_AN_{filter}.recode.vcf"), sample=config["samples"], chrms=config["chrY"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/AN/chrX_nonPARs/{sample}_{chrms}_nonPARs_AN_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrms=config["chrX"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/AN/chrY/{sample}_{chrms}_AN_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrms=config["chrY"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/AN/chrM/by_sample/{sample}_{chrmsM}_AN_{filter}.recode.vcf"), sample=config["samples"], chrmsM=config["chrM"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/AN/chrM/{sample}_{chrmsM}_AN_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrmsM=config["chrM"], filter=config["AN_filters"]),

        # MALES - DEFAULT #
        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf"), sample=config["samples"], chrmsXnonPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrmsY}.recode.vcf"), sample=config["samples"], chrmsY=config["chrY"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrmsM}.recode.vcf"), sample=config["samples"], chrmsM=config["chrM"]),

        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/default/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrms=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/default/chrY/{sample}_{chrms}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrms=config["chrY"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/default/chrM/{sample}_{chrmsM}_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrmsM=config["chrM"]),

'''



'''
# No filter, not sure if doing this...so keep here for now
        expand(os.path.join(config["scratch_proj_path"], "haploid/joint_called_vcfs/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs_nonPARs.vcf.gz"), chrmsXnonPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/joint_called_vcfs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_no_filter.recode.vcf"), sample=config["samples"], chrmsXnonPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/joint_called_vcfs/by_sample/{sample}_{chrmsY}_no_filter.recode.vcf"), sample=config["samples"], chrmsY=config["chrY"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/no_filter/chrX_nonPARs/{sample}_{chrms}_nonPARs_no_filter_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrms=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/no_filter/chrY/{sample}_{chrms}_no_filter_golden_vs_called_performance_metrics.txt"), sample=config["samples"], chrms=config["chrY"])
'''

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - ALL
#-------------------------------------------------------------------------------#
rule extract_samples_call_VCF_XnonPARs:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/all/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_gatkHardFilter_filtered_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_samples_call_VCF_Y:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/all/{chrmsY}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_gatkHardFilter_filtered.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrmsY}")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrmsY}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule get_performance_metrics_XnonPARs:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrms}_nonPARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrms}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/chrY/{sample}_{chrms}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/chrY/{sample}_{chrms}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule extract_samples_call_VCF_M:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/all/{chrmsM}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_gatkHardFilter_filtered.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrmsM}")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrmsM}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_M:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrmsM}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrmsM}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/chrM/{sample}_{chrmsM}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/chrM/{sample}_{chrmsM}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

'''
#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - QD
#-------------------------------------------------------------------------------#
rule extract_samples_call_VCF_XnonPARs_QD:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/QD/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_QD_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/QD/chrX_nonPARs/by_sample/{sample}_{chrms}_QD_{filter}_nonPARs")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/QD/chrX_nonPARs/by_sample/{sample}_{chrms}_QD_{filter}_nonPARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule extract_samples_call_VCF_Y_QD:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/QD/chrY/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_QD_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/QD/chrY/by_sample/{sample}_{chrms}_QD_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/QD/chrY/by_sample/{sample}_{chrms}_QD_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule get_performance_metrics_XnonPARs_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/QD/chrX_nonPARs/by_sample/{sample}_{chrms}_QD_{filter}_nonPARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/QD/chrX_nonPARs/{sample}_{chrms}_nonPARs_QD_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/QD/chrX_nonPARs/{sample}_{chrms}_nonPARs_QD_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/QD/{chrms}/by_sample/{sample}_{chrms}_QD_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/QD/chrY/{sample}_{chrms}_QD_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/QD/chrY/{sample}_{chrms}_QD_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule extract_samples_call_VCF_M_QD:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/QD/chrM/{chrmsM}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_QD_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/QD/chrM/by_sample/{sample}_{chrmsM}_QD_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/QD/chrM/by_sample/{sample}_{chrmsM}_QD_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_M_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrmsM}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/QD/{chrmsM}/by_sample/{sample}_{chrmsM}_QD_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/QD/chrM/{sample}_{chrmsM}_QD_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/QD/chrM/{sample}_{chrmsM}_QD_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - QUAL
#-------------------------------------------------------------------------------#
rule extract_samples_call_VCF_XnonPARs_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/QUAL/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_QUAL_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/QUAL/chrX_nonPARs/by_sample/{sample}_{chrms}_QUAL_{filter}_nonPARs")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/QUAL/chrX_nonPARs/by_sample/{sample}_{chrms}_QUAL_{filter}_nonPARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule extract_samples_call_VCF_Y_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/QUAL/chrY/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_QUAL_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/QUAL/chrY/by_sample/{sample}_{chrms}_QUAL_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/QUAL/chrY/by_sample/{sample}_{chrms}_QUAL_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule get_performance_metrics_XnonPARs_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/QUAL/chrX_nonPARs/by_sample/{sample}_{chrms}_QUAL_{filter}_nonPARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/QUAL/chrX_nonPARs/{sample}_{chrms}_nonPARs_QUAL_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/QUAL/chrX_nonPARs/{sample}_{chrms}_nonPARs_QUAL_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/QUAL/chrY/by_sample/{sample}_{chrms}_QUAL_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/QUAL/chrY/{sample}_{chrms}_QUAL_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/QUAL/chrY/{sample}_{chrms}_QUAL_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule extract_samples_call_VCF_M_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/QUAL/chrM/{chrmsM}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_QUAL_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/QUAL/chrM/by_sample/{sample}_{chrmsM}_QUAL_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/QUAL/chrM/by_sample/{sample}_{chrmsM}_QUAL_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_M_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrmsM}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/QUAL/{chrmsM}/by_sample/{sample}_{chrmsM}_QUAL_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/QUAL/chrM/{sample}_{chrmsM}_QUAL_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/QUAL/chrM/{sample}_{chrmsM}_QUAL_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - SOR
#-------------------------------------------------------------------------------#
rule extract_samples_call_VCF_XnonPARs_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/SOR/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_SOR_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/SOR/chrX_nonPARs/by_sample/{sample}_{chrms}_SOR_{filter}_nonPARs")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/SOR/chrX_nonPARs/by_sample/{sample}_{chrms}_SOR_{filter}_nonPARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule extract_samples_call_VCF_Y_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/SOR/chrY/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_SOR_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/SOR/chrY/by_sample/{sample}_{chrms}_SOR_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/SOR/chrY/by_sample/{sample}_{chrms}_SOR_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule get_performance_metrics_XnonPARs_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/SOR/chrX_nonPARs/by_sample/{sample}_{chrms}_SOR_{filter}_nonPARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/SOR/chrX_nonPARs/{sample}_{chrms}_nonPARs_SOR_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/SOR/chrX_nonPARs/{sample}_{chrms}_nonPARs_SOR_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/SOR/chrY/by_sample/{sample}_{chrms}_SOR_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/SOR/chrY/{sample}_{chrms}_SOR_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/SOR/chrY/{sample}_{chrms}_SOR_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule extract_samples_call_VCF_M_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/SOR/chrM/{chrmsM}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_SOR_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/SOR/chrM/by_sample/{sample}_{chrmsM}_SOR_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/SOR/chrM/by_sample/{sample}_{chrmsM}_SOR_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_M_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrmsM}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/SOR/{chrmsM}/by_sample/{sample}_{chrmsM}_SOR_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/SOR/chrM/{sample}_{chrmsM}_SOR_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/SOR/chrM/{sample}_{chrmsM}_SOR_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - FS
#-------------------------------------------------------------------------------#
rule extract_samples_call_VCF_XnonPARs_FS:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/FS/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_FS_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/FS/chrX_nonPARs/by_sample/{sample}_{chrms}_FS_{filter}_nonPARs")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/FS/chrX_nonPARs/by_sample/{sample}_{chrms}_FS_{filter}_nonPARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule extract_samples_call_VCF_Y_FS:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/FS/chrY/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_FS_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/FS/chrY/by_sample/{sample}_{chrms}_FS_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/FS/chrY/by_sample/{sample}_{chrms}_FS_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule get_performance_metrics_XnonPARs_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/FS/chrX_nonPARs/by_sample/{sample}_{chrms}_FS_{filter}_nonPARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/FS/chrX_nonPARs/{sample}_{chrms}_nonPARs_FS_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/FS/chrX_nonPARs/{sample}_{chrms}_nonPARs_FS_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/FS/chrY/by_sample/{sample}_{chrms}_FS_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/FS/chrY/{sample}_{chrms}_FS_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/FS/chrY/{sample}_{chrms}_FS_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule extract_samples_call_VCF_M_FS:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/FS/chrM/{chrmsM}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_FS_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/FS/chrM/by_sample/{sample}_{chrmsM}_FS_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/FS/chrM/by_sample/{sample}_{chrmsM}_FS_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_M_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrmsM}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/FS/{chrmsM}/by_sample/{sample}_{chrmsM}_FS_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/FS/chrM/{sample}_{chrmsM}_FS_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/FS/chrM/{sample}_{chrmsM}_FS_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - MQ
#-------------------------------------------------------------------------------#
rule extract_samples_call_VCF_XnonPARs_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/MQ/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_MQ_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/MQ/chrX_nonPARs/by_sample/{sample}_{chrms}_MQ_{filter}_nonPARs")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/MQ/chrX_nonPARs/by_sample/{sample}_{chrms}_MQ_{filter}_nonPARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule extract_samples_call_VCF_Y_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/MQ/chrY/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_MQ_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/MQ/chrY/by_sample/{sample}_{chrms}_MQ_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/MQ/chrY/by_sample/{sample}_{chrms}_MQ_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule get_performance_metrics_XnonPARs_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/MQ/chrX_nonPARs/by_sample/{sample}_{chrms}_MQ_{filter}_nonPARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/MQ/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQ_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/MQ/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQ_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/MQ/chrY/by_sample/{sample}_{chrms}_MQ_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/MQ/chrY/{sample}_{chrms}_MQ_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/MQ/chrY/{sample}_{chrms}_MQ_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule extract_samples_call_VCF_M_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/MQ/chrM/{chrmsM}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_MQ_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/MQ/chrM/by_sample/{sample}_{chrmsM}_MQ_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/MQ/chrM/by_sample/{sample}_{chrmsM}_MQ_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_M_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrmsM}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/MQ/{chrmsM}/by_sample/{sample}_{chrmsM}_MQ_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/MQ/chrM/{sample}_{chrmsM}_MQ_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/MQ/chrM/{sample}_{chrmsM}_MQ_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - MQRankSum
#-------------------------------------------------------------------------------#
rule extract_samples_call_VCF_XnonPARs_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/MQRankSum/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_MQRankSum_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/MQRankSum/chrX_nonPARs/by_sample/{sample}_{chrms}_MQRankSum_{filter}_nonPARs")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/MQRankSum/chrX_nonPARs/by_sample/{sample}_{chrms}_MQRankSum_{filter}_nonPARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule extract_samples_call_VCF_Y_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/MQRankSum/chrY/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_MQRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/MQRankSum/chrY/by_sample/{sample}_{chrms}_MQRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/MQRankSum/chrY/by_sample/{sample}_{chrms}_MQRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule get_performance_metrics_XnonPARs_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/MQRankSum/chrX_nonPARs/by_sample/{sample}_{chrms}_MQRankSum_{filter}_nonPARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/MQRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/MQRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/MQRankSum/chrY/by_sample/{sample}_{chrms}_MQRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/MQRankSum/chrY/{sample}_{chrms}_MQRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/MQRankSum/chrY/{sample}_{chrms}_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule extract_samples_call_VCF_M_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/MQRankSum/chrM/{chrmsM}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_MQRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/MQRankSum/chrM/by_sample/{sample}_{chrmsM}_MQRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/MQRankSum/chrM/by_sample/{sample}_{chrmsM}_MQRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_M_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrmsM}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/MQRankSum/{chrmsM}/by_sample/{sample}_{chrmsM}_MQRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/MQRankSum/chrM/{sample}_{chrmsM}_MQRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/MQRankSum/chrM/{sample}_{chrmsM}_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - ReadPosRankSum
#-------------------------------------------------------------------------------#
rule extract_samples_call_VCF_XnonPARs_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_ReadPosRankSum_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrX_nonPARs/by_sample/{sample}_{chrms}_ReadPosRankSum_{filter}_nonPARs")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrX_nonPARs/by_sample/{sample}_{chrms}_ReadPosRankSum_{filter}_nonPARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule extract_samples_call_VCF_Y_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrY/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_ReadPosRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrY/by_sample/{sample}_{chrms}_ReadPosRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrY/by_sample/{sample}_{chrms}_ReadPosRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule get_performance_metrics_XnonPARs_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrX_nonPARs/by_sample/{sample}_{chrms}_ReadPosRankSum_{filter}_nonPARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/ReadPosRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_ReadPosRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/ReadPosRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrY/by_sample/{sample}_{chrms}_ReadPosRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/ReadPosRankSum/chrY/{sample}_{chrms}_ReadPosRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/ReadPosRankSum/chrY/{sample}_{chrms}_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule extract_samples_call_VCF_M_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrM/{chrmsM}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_ReadPosRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrM/by_sample/{sample}_{chrmsM}_ReadPosRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrM/by_sample/{sample}_{chrmsM}_ReadPosRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_M_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrmsM}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/ReadPosRankSum/{chrmsM}/by_sample/{sample}_{chrmsM}_ReadPosRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/ReadPosRankSum/chrM/{sample}_{chrmsM}_ReadPosRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/ReadPosRankSum/chrM/{sample}_{chrmsM}_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - DP
#-------------------------------------------------------------------------------#
rule extract_samples_call_VCF_XnonPARs_DP:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/DP/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_DP_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/DP/chrX_nonPARs/by_sample/{sample}_{chrms}_DP_{filter}_nonPARs")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/DP/chrX_nonPARs/by_sample/{sample}_{chrms}_DP_{filter}_nonPARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule extract_samples_call_VCF_Y_DP:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/DP/chrY/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_DP_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/DP/chrY/by_sample/{sample}_{chrms}_DP_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/DP/chrY/by_sample/{sample}_{chrms}_DP_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule get_performance_metrics_XnonPARs_DP:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/DP/chrX_nonPARs/by_sample/{sample}_{chrms}_DP_{filter}_nonPARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/DP/chrX_nonPARs/{sample}_{chrms}_nonPARs_DP_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/DP/chrX_nonPARs/{sample}_{chrms}_nonPARs_DP_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y_DP:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/DP/chrY/by_sample/{sample}_{chrms}_DP_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/DP/chrY/{sample}_{chrms}_DP_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/DP/chrY/{sample}_{chrms}_DP_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule extract_samples_call_VCF_M_DP:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/DP/chrM/{chrmsM}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_DP_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/DP/chrM/by_sample/{sample}_{chrmsM}_DP_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/DP/chrM/by_sample/{sample}_{chrmsM}_DP_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_M_DP:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrmsM}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/DP/{chrmsM}/by_sample/{sample}_{chrmsM}_DP_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/DP/chrM/{sample}_{chrmsM}_DP_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/DP/chrM/{sample}_{chrmsM}_DP_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - AN
#-------------------------------------------------------------------------------#
rule extract_samples_call_VCF_XnonPARs_AN:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/AN/chrX_nonPARs/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_AN_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/AN/chrX_nonPARs/by_sample/{sample}_{chrms}_AN_{filter}_nonPARs")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/AN/chrX_nonPARs/by_sample/{sample}_{chrms}_AN_{filter}_nonPARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule extract_samples_call_VCF_Y_AN:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/AN/chrY/{chrms}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_AN_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/AN/chrY/by_sample/{sample}_{chrms}_AN_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/AN/chrY/by_sample/{sample}_{chrms}_AN_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule get_performance_metrics_XnonPARs_AN:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/AN/chrX_nonPARs/by_sample/{sample}_{chrms}_AN_{filter}_nonPARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/AN/chrX_nonPARs/{sample}_{chrms}_nonPARs_AN_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/AN/chrX_nonPARs/{sample}_{chrms}_nonPARs_AN_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y_AN:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/AN/chrY/by_sample/{sample}_{chrms}_AN_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/AN/chrY/{sample}_{chrms}_AN_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/AN/chrY/{sample}_{chrms}_AN_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule extract_samples_call_VCF_M_AN:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/AN/chrM/{chrmsM}_GRCh38_YPARsMasked_gatk_haploid_called_SNPs_filtered_AN_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/AN/chrM/by_sample/{sample}_{chrmsM}_AN_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/AN/chrM/by_sample/{sample}_{chrmsM}_AN_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_M_AN:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrmsM}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/AN/{chrmsM}/by_sample/{sample}_{chrmsM}_AN_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/AN/chrM/{sample}_{chrmsM}_AN_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/AN/chrM/{sample}_{chrmsM}_AN_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

'''

'''

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - No filters
#-------------------------------------------------------------------------------#
rule extract_nonPARs_nofilter:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/joint_called_vcfs/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/joint_called_vcfs/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs_nonPARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrmsXnonPARs}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.chrms} -XL {params.ilist} """
        """-O {output} """

rule extract_samples_call_VCF_XnonPARs_no_filter:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/joint_called_vcfs/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_haploid_called_raw_SNPs_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/joint_called_vcfs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_no_filter")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/joint_called_vcfs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_no_filter.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_samples_call_VCF_Y_no_filter:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/joint_called_vcfs/{chrmsY}_GRCh38_YPARsMasked_gatk_haploid_called_raw.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/joint_called_vcfs/by_sample/{sample}_{chrmsY}_no_filter")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/joint_called_vcfs/by_sample/{sample}_{chrmsY}_no_filter.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_XnonPARs_no_filter:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/joint_called_vcfs/by_sample/{sample}_{chrms}_nonPARs_no_filter.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/no_filter/chrX_nonPARs/{sample}_{chrms}_nonPARs_no_filter_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/no_filter/chrX_nonPARs/{sample}_{chrms}_nonPARs_no_filter_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y_no_filter:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/joint_called_vcfs/by_sample/{sample}_{chrms}_no_filter.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/no_filter/chrY/{sample}_{chrms}_no_filter_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/no_filter/chrY/{sample}_{chrms}_no_filter_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """
'''

#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#
# MALES - DEFAULT
#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - ALL
#-------------------------------------------------------------------------------#
rule default_extract_samples_call_VCF_XnonPARs:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/default/all/{chrmsXnonPARs}_GRCh38_default_gatk_haploid_called_SNPs_gatkHardFilter_filtered_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_samples_call_VCF_Y:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/default/all/{chrmsY}_GRCh38_default_gatk_haploid_called_SNPs_gatkHardFilter_filtered.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrmsY}")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrmsY}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule default_get_performance_metrics_XnonPARs:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrms}_nonPARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/default/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/default/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_Y:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrms}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/default/chrY/{sample}_{chrms}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/default/chrY/{sample}_{chrms}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule default_extract_samples_call_VCF_M:
    input:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/default/all/{chrmsM}_GRCh38_default_gatk_haploid_called_SNPs_gatkHardFilter_filtered.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrmsM}")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrmsM}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_get_performance_metrics_M:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrmsM}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "haploid/hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrmsM}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/default/chrM/{sample}_{chrmsM}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "haploid/compare_VCFs/ASN/males/default/chrM/{sample}_{chrmsM}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """
