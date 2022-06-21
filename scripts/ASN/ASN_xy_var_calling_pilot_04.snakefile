import os

# Environment: variant_calling_simulations_project

configfile: "ASN_xy_var_calling_pilot_04.config.json"


wildcard_constraints:
    chrmsXnonPARs= '|'.join([re.escape(x) for x in config["chrX"]]),
    chrmsXPARs= '|'.join([re.escape(x) for x in config["chrX"]]),
    chrmsY= '|'.join([re.escape(x) for x in config["chrY"]]),
    chrmsA= '|'.join([re.escape(x) for x in config["autosomes"]]),
    chrmsM= '|'.join([re.escape(x) for x in config["chrM"]])



rule all:
    input:
        # MALES #
        # SCC #
        expand(os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsA}_golden.vcf"), sample=config["males"], chrmsA=config["autosomes"]),
        expand(os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrmsXnonPARs}_nonPARs_golden.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsXPARs}_PARs_golden.vcf"), sample=config["males"], chrmsXPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrmsY}_golden.vcf"), sample=config["males"], chrmsY=config["chrY"]),
        expand(os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrmsM}_golden.vcf"), sample=config["males"], chrmsM=config["chrM"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrmsA}_autos.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrmsXPARs}_PARs.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrmsY}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrmsM}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"]),

        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/autos/{sample}_{chrms}_autos_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/chrX_PARs/{sample}_{chrms}_PARs_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/chrY/{sample}_{chrms}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/chrM/{sample}_{chrms}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"]),

        # FEMALES #
        # SCC #
        expand(os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsA}_golden.vcf"), sample=config["females"], chrmsA=config["autosomes"]),
        expand(os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsXnonPARs}_nonPARs_golden.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsXPARs}_PARs_golden.vcf"), sample=config["females"], chrmsXPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrmsM}_golden.vcf"), sample=config["females"], chrmsM=config["chrM"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/all/by_sample/{sample}_{chrmsA}_autos.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/all/by_sample/{sample}_{chrmsXPARs}_PARs.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/all/by_sample/{sample}_{chrmsM}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"]),

        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/autos/{sample}_{chrms}_autos_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/chrX_PARs/{sample}_{chrms}_PARs_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/chrM/{sample}_{chrms}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"]),

        # FEMALES #
        # DEFAULT #
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/all/by_sample/{sample}_{chrmsA}_autos.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/all/by_sample/{sample}_{chrmsXPARs}_PARs.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/all/by_sample/{sample}_{chrmsM}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"]),

        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/autos/{sample}_{chrms}_autos_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/chrX_PARs/{sample}_{chrms}_PARs_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/chrM/{sample}_{chrms}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"]),

        # MALES DEFAULT #
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrmsA}_autos.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrmsXPARs}_PARs.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrmsY}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrmsM}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"]),

        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/default/autos/{sample}_{chrms}_autos_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/default/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/default/chrX_PARs/{sample}_{chrms}_PARs_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/default/chrY/{sample}_{chrms}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/default/chrM/{sample}_{chrms}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"]),


'''
        # hard filters
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QD/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QD_{filter}.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QD/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QD_{filter}.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QD/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QD_{filter}.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QD/chrY/by_sample/{sample}_{chrmsY}_QD_{filter}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/QD/autos/{sample}_{chrms}_autos_QD_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/QD/chrX_nonPARs/{sample}_{chrms}_nonPARs_QD_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/QD/chrX_PARs/{sample}_{chrms}_PARs_QD_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/QD/chrY/{sample}_{chrms}_QD_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QD/chrM/by_sample/{sample}_{chrmsM}_QD_{filter}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/QD/chrM/{sample}_{chrms}_QD_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"], filter=config["QD_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QUAL/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QUAL_{filter}.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QUAL/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QUAL_{filter}.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QUAL/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QUAL_{filter}.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QUAL/chrY/by_sample/{sample}_{chrmsY}_QUAL_{filter}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/QUAL/autos/{sample}_{chrms}_autos_QUAL_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/QUAL/chrX_nonPARs/{sample}_{chrms}_nonPARs_QUAL_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/QUAL/chrX_PARs/{sample}_{chrms}_PARs_QUAL_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/QUAL/chrY/{sample}_{chrms}_QUAL_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QUAL/chrM/by_sample/{sample}_{chrmsM}_QUAL_{filter}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/QUAL/chrM/{sample}_{chrms}_QUAL_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"], filter=config["QUAL_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/SOR/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_SOR_{filter}.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/SOR/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_SOR_{filter}.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/SOR/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_SOR_{filter}.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/SOR/chrY/by_sample/{sample}_{chrmsY}_SOR_{filter}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/SOR/autos/{sample}_{chrms}_autos_SOR_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/SOR/chrX_nonPARs/{sample}_{chrms}_nonPARs_SOR_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/SOR/chrX_PARs/{sample}_{chrms}_PARs_SOR_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/SOR/chrY/{sample}_{chrms}_SOR_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/SOR/chrM/by_sample/{sample}_{chrmsM}_SOR_{filter}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/SOR/chrM/{sample}_{chrms}_SOR_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"], filter=config["SOR_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/FS/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_FS_{filter}.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/FS/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_FS_{filter}.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/FS/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_FS_{filter}.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/FS/chrY/by_sample/{sample}_{chrmsY}_FS_{filter}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/FS/autos/{sample}_{chrms}_autos_FS_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/FS/chrX_nonPARs/{sample}_{chrms}_nonPARs_FS_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/FS/chrX_PARs/{sample}_{chrms}_PARs_FS_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/FS/chrY/{sample}_{chrms}_FS_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/FS/chrM/by_sample/{sample}_{chrmsM}_FS_{filter}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/FS/chrM/{sample}_{chrms}_FS_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"], filter=config["FS_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQ/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQ_{filter}.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQ/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQ_{filter}.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQ/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQ_{filter}.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQ/chrY/by_sample/{sample}_{chrmsY}_MQ_{filter}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/MQ/autos/{sample}_{chrms}_autos_MQ_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/MQ/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQ_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/MQ/chrX_PARs/{sample}_{chrms}_PARs_MQ_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/MQ/chrY/{sample}_{chrms}_MQ_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQ/chrM/by_sample/{sample}_{chrmsM}_MQ_{filter}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/MQ/chrM/{sample}_{chrms}_MQ_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"], filter=config["MQ_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQRankSum_{filter}.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQRankSum_{filter}.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQRankSum_{filter}.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQRankSum/chrY/by_sample/{sample}_{chrmsY}_MQRankSum_{filter}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/MQRankSum/autos/{sample}_{chrms}_autos_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/MQRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/MQRankSum/chrX_PARs/{sample}_{chrms}_PARs_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/MQRankSum/chrY/{sample}_{chrms}_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQRankSum/chrM/by_sample/{sample}_{chrmsM}_MQRankSum_{filter}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/MQRankSum/chrM/{sample}_{chrms}_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"], filter=config["MQRankSum_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/ReadPosRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_ReadPosRankSum_{filter}.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_ReadPosRankSum_{filter}.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_ReadPosRankSum_{filter}.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrY/by_sample/{sample}_{chrmsY}_ReadPosRankSum_{filter}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/ReadPosRankSum/autos/{sample}_{chrms}_autos_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/ReadPosRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/ReadPosRankSum/chrX_PARs/{sample}_{chrms}_PARs_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/ReadPosRankSum/chrY/{sample}_{chrms}_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrM/by_sample/{sample}_{chrmsM}_ReadPosRankSum_{filter}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/ReadPosRankSum/chrM/{sample}_{chrms}_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"], filter=config["ReadPosRankSum_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/DP/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_DP_{filter}.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/DP/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_DP_{filter}.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/DP/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_DP_{filter}.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/DP/chrY/by_sample/{sample}_{chrmsY}_DP_{filter}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/DP/autos/{sample}_{chrms}_autos_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/DP/chrX_nonPARs/{sample}_{chrms}_nonPARs_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/DP/chrX_PARs/{sample}_{chrms}_PARs_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/DP/chrY/{sample}_{chrms}_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/DP/chrM/by_sample/{sample}_{chrmsM}_DP_{filter}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/DP/chrM/{sample}_{chrms}_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"], filter=config["DP_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/AN/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_AN_{filter}.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/AN/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_AN_{filter}.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/AN/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_AN_{filter}.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/AN/chrY/by_sample/{sample}_{chrmsY}_AN_{filter}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/AN/autos/{sample}_{chrms}_autos_AN_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/AN/chrX_nonPARs/{sample}_{chrms}_nonPARs_AN_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/AN/chrX_PARs/{sample}_{chrms}_PARs_AN_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/AN/chrY/{sample}_{chrms}_AN_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/AN/chrM/by_sample/{sample}_{chrmsM}_AN_{filter}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/AN/chrM/{sample}_{chrms}_AN_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"], filter=config["AN_filters"]),


        # FEMALES #
        # SCC #
        expand(os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsA}_golden.vcf"), sample=config["females"], chrmsA=config["autosomes"]),
        expand(os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsXnonPARs}_nonPARs_golden.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsXPARs}_PARs_golden.vcf"), sample=config["females"], chrmsXPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrmsM}_golden.vcf"), sample=config["females"], chrmsM=config["chrM"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/all/by_sample/{sample}_{chrmsA}_autos.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/all/by_sample/{sample}_{chrmsXPARs}_PARs.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/all/by_sample/{sample}_{chrmsM}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"]),

        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/autos/{sample}_{chrms}_autos_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/chrX_PARs/{sample}_{chrms}_PARs_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/chrM/{sample}_{chrms}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"]),

        # hard filters
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QD/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QD_{filter}.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QD/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QD_{filter}.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QD/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QD_{filter}.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/QD/autos/{sample}_{chrms}_autos_QD_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/QD/chrX_nonPARs/{sample}_{chrms}_nonPARs_QD_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/QD/chrX_PARs/{sample}_{chrms}_PARs_QD_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QD/chrM/by_sample/{sample}_{chrmsM}_QD_{filter}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/QD/chrM/{sample}_{chrms}_QD_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"], filter=config["QD_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QUAL/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QUAL_{filter}.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QUAL/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QUAL_{filter}.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QUAL/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QUAL_{filter}.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/QUAL/autos/{sample}_{chrms}_autos_QUAL_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/QUAL/chrX_nonPARs/{sample}_{chrms}_nonPARs_QUAL_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/QUAL/chrX_PARs/{sample}_{chrms}_PARs_QUAL_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QUAL/chrM/by_sample/{sample}_{chrmsM}_QUAL_{filter}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/QUAL/chrM/{sample}_{chrms}_QUAL_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"], filter=config["QUAL_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/SOR/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_SOR_{filter}.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/SOR/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_SOR_{filter}.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/SOR/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_SOR_{filter}.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/SOR/autos/{sample}_{chrms}_autos_SOR_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/SOR/chrX_nonPARs/{sample}_{chrms}_nonPARs_SOR_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/SOR/chrX_PARs/{sample}_{chrms}_PARs_SOR_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/SOR/chrM/by_sample/{sample}_{chrmsM}_SOR_{filter}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/SOR/chrM/{sample}_{chrms}_SOR_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"], filter=config["SOR_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/FS/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_FS_{filter}.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/FS/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_FS_{filter}.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/FS/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_FS_{filter}.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/FS/autos/{sample}_{chrms}_autos_FS_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/FS/chrX_nonPARs/{sample}_{chrms}_nonPARs_FS_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/FS/chrX_PARs/{sample}_{chrms}_PARs_FS_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/FS/chrM/by_sample/{sample}_{chrmsM}_FS_{filter}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/FS/chrM/{sample}_{chrms}_FS_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"], filter=config["FS_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQ/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQ_{filter}.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQ/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQ_{filter}.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQ/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQ_{filter}.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/MQ/autos/{sample}_{chrms}_autos_MQ_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/MQ/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQ_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/MQ/chrX_PARs/{sample}_{chrms}_PARs_MQ_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQ/chrM/by_sample/{sample}_{chrmsM}_MQ_{filter}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/MQ/chrM/{sample}_{chrms}_MQ_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"], filter=config["MQ_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQRankSum_{filter}.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQRankSum_{filter}.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQRankSum_{filter}.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/MQRankSum/autos/{sample}_{chrms}_autos_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/MQRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/MQRankSum/chrX_PARs/{sample}_{chrms}_PARs_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQRankSum/chrM/by_sample/{sample}_{chrmsM}_MQRankSum_{filter}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/MQRankSum/chrM/{sample}_{chrms}_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"], filter=config["MQRankSum_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/ReadPosRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_ReadPosRankSum_{filter}.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/ReadPosRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_ReadPosRankSum_{filter}.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/ReadPosRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_ReadPosRankSum_{filter}.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/ReadPosRankSum/autos/{sample}_{chrms}_autos_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/ReadPosRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/ReadPosRankSum/chrX_PARs/{sample}_{chrms}_PARs_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/ReadPosRankSum/chrM/by_sample/{sample}_{chrmsM}_ReadPosRankSum_{filter}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/ReadPosRankSum/chrM/{sample}_{chrms}_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"], filter=config["ReadPosRankSum_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/DP/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_DP_{filter}.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/DP/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_DP_{filter}.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/DP/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_DP_{filter}.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/DP/autos/{sample}_{chrms}_autos_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/DP/chrX_nonPARs/{sample}_{chrms}_nonPARs_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/DP/chrX_PARs/{sample}_{chrms}_PARs_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/DP/chrM/by_sample/{sample}_{chrmsM}_DP_{filter}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/DP/chrM/{sample}_{chrms}_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"], filter=config["DP_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/AN/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_AN_{filter}.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/AN/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_AN_{filter}.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/AN/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_AN_{filter}.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/AN/autos/{sample}_{chrms}_autos_AN_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/AN/chrX_nonPARs/{sample}_{chrms}_nonPARs_AN_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/AN/chrX_PARs/{sample}_{chrms}_PARs_AN_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/AN/chrM/by_sample/{sample}_{chrmsM}_AN_{filter}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/AN/chrM/{sample}_{chrms}_AN_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"], filter=config["AN_filters"]),

'''

#-------------------------------------------------------------------------------#
# Step: Prep golden (simulated) VCFs
#-------------------------------------------------------------------------------#
rule unzip_golden_VCFs_autos:
    input:
        os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsA}_golden.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsA}_golden.vcf")
    shell:
        """
        gunzip {input}
        """

rule unzip_golden_VCFs_X:
    input:
        os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrmsXnonPARs}_nonPARs_golden.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrmsXnonPARs}_nonPARs_golden.vcf")
    shell:
        """
        gunzip {input}
        """

rule unzip_golden_VCFs_X_PARs:
    input:
        os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsXPARs}_PARs_golden.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsXPARs}_PARs_golden.vcf")
    shell:
        """
        gunzip {input}
        """

rule unzip_golden_VCFs_Y:
    input:
        os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrmsY}_golden.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrmsY}_golden.vcf")
    shell:
        """
        gunzip {input}
        """

rule unzip_golden_VCFs_M:
    input:
        os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrmsM}_golden.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrmsM}_golden.vcf")
    shell:
        """
        gunzip {input}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - ALL
#-------------------------------------------------------------------------------#
rule extract_males_call_VCF_autos:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/all/{chrmsA}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrmsA}_autos")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrmsA}_autos.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XnonPARs:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/all/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XPARs:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/all/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrmsXPARs}_PARs")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrmsXPARs}_PARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_Y:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/all/{chrmsY}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrmsY}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrmsY}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_M:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/all/{chrmsM}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrmsM}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrmsM}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_autos:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrms}_autos.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/autos/{sample}_{chrms}_autos_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/autos/{sample}_{chrms}_autos_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XnonPARs:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrms}_nonPARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XPARs:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrms}_PARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/chrX_PARs/{sample}_{chrms}_PARs_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/chrX_PARs/{sample}_{chrms}_PARs_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrms}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/chrY/{sample}_{chrms}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/chrY/{sample}_{chrms}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_M:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/all/by_sample/{sample}_{chrms}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/chrM/{sample}_{chrms}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/chrM/{sample}_{chrms}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

'''
#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - QD
#-------------------------------------------------------------------------------#
rule extract_males_call_VCF_autos_QD:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QD/{chrmsA}/{chrmsA}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QD/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QD_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QD/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QD_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XnonPARs_QD:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QD/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QD/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QD_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QD/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QD_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XPARs_QD:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QD/chrX_PARs/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QD/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QD_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QD/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QD_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_Y_QD:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QD/chrY/{chrmsY}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QD/chrY/by_sample/{sample}_{chrmsY}_QD_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QD/chrY/by_sample/{sample}_{chrmsY}_QD_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule get_performance_metrics_autos_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QD/{chrms}/by_sample/{sample}_{chrms}_autos_QD_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/QD/autos/{sample}_{chrms}_autos_QD_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/QD/autos/{sample}_{chrms}_autos_QD_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XnonPARs_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QD/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_QD_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/QD/chrX_nonPARs/{sample}_{chrms}_nonPARs_QD_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/QD/chrX_nonPARs/{sample}_{chrms}_nonPARs_QD_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XPARs_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QD/chrX_PARs/by_sample/{sample}_{chrms}_PARs_QD_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/QD/chrX_PARs/{sample}_{chrms}_PARs_QD_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/QD/chrX_PARs/{sample}_{chrms}_PARs_QD_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QD/chrY/by_sample/{sample}_{chrms}_QD_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/QD/chrY/{sample}_{chrms}_QD_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/QD/chrY/{sample}_{chrms}_QD_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule extract_males_call_VCF_M_QD:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QD/chrM/{chrmsM}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QD/chrM/by_sample/{sample}_{chrmsM}_QD_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QD/chrM/by_sample/{sample}_{chrmsM}_QD_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_M_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QD/chrM/by_sample/{sample}_{chrms}_QD_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/QD/chrM/{sample}_{chrms}_QD_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/QD/chrM/{sample}_{chrms}_QD_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - QUAL
#-------------------------------------------------------------------------------#
rule extract_males_call_VCF_autos_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QUAL/{chrmsA}/{chrmsA}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QUAL/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QUAL_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QUAL/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QUAL_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XnonPARs_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QUAL/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QUAL/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QUAL_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QUAL/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QUAL_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XPARs_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QUAL/chrX_PARs/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QUAL/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QUAL_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QUAL/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QUAL_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_Y_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QUAL/chrY/{chrmsY}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QUAL/chrY/by_sample/{sample}_{chrmsY}_QUAL_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QUAL/chrY/by_sample/{sample}_{chrmsY}_QUAL_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule get_performance_metrics_autos_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QUAL/{chrms}/by_sample/{sample}_{chrms}_autos_QUAL_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/QUAL/autos/{sample}_{chrms}_autos_QUAL_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/QUAL/autos/{sample}_{chrms}_autos_QUAL_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XnonPARs_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QUAL/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_QUAL_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/QUAL/chrX_nonPARs/{sample}_{chrms}_nonPARs_QUAL_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/QUAL/chrX_nonPARs/{sample}_{chrms}_nonPARs_QUAL_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XPARs_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QUAL/chrX_PARs/by_sample/{sample}_{chrms}_PARs_QUAL_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/QUAL/chrX_PARs/{sample}_{chrms}_PARs_QUAL_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/QUAL/chrX_PARs/{sample}_{chrms}_PARs_QUAL_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QUAL/chrY/by_sample/{sample}_{chrms}_QUAL_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/QUAL/chrY/{sample}_{chrms}_QUAL_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/QUAL/chrY/{sample}_{chrms}_QUAL_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule extract_males_call_VCF_M_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QUAL/chrM/{chrmsM}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QUAL/chrM/by_sample/{sample}_{chrmsM}_QUAL_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QUAL/chrM/by_sample/{sample}_{chrmsM}_QUAL_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_M_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/QUAL/chrM/by_sample/{sample}_{chrms}_QUAL_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/QUAL/chrM/{sample}_{chrms}_QUAL_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/QUAL/chrM/{sample}_{chrms}_QUAL_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - SOR
#-------------------------------------------------------------------------------#
rule extract_males_call_VCF_autos_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/SOR/{chrmsA}/{chrmsA}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/SOR/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_SOR_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/SOR/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_SOR_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XnonPARs_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/SOR/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/SOR/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_SOR_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/SOR/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_SOR_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XPARs_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/SOR/chrX_PARs/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/SOR/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_SOR_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/SOR/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_SOR_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_Y_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/SOR/chrY/{chrmsY}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/SOR/chrY/by_sample/{sample}_{chrmsY}_SOR_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/SOR/chrY/by_sample/{sample}_{chrmsY}_SOR_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule get_performance_metrics_autos_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/SOR/{chrms}/by_sample/{sample}_{chrms}_autos_SOR_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/SOR/autos/{sample}_{chrms}_autos_SOR_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/SOR/autos/{sample}_{chrms}_autos_SOR_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XnonPARs_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/SOR/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_SOR_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/SOR/chrX_nonPARs/{sample}_{chrms}_nonPARs_SOR_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/SOR/chrX_nonPARs/{sample}_{chrms}_nonPARs_SOR_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XPARs_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/SOR/chrX_PARs/by_sample/{sample}_{chrms}_PARs_SOR_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/SOR/chrX_PARs/{sample}_{chrms}_PARs_SOR_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/SOR/chrX_PARs/{sample}_{chrms}_PARs_SOR_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/SOR/chrY/by_sample/{sample}_{chrms}_SOR_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/SOR/chrY/{sample}_{chrms}_SOR_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/SOR/chrY/{sample}_{chrms}_SOR_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule extract_males_call_VCF_M_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/SOR/chrM/{chrmsM}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/SOR/chrM/by_sample/{sample}_{chrmsM}_SOR_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/SOR/chrM/by_sample/{sample}_{chrmsM}_SOR_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_M_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/SOR/chrM/by_sample/{sample}_{chrms}_SOR_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/SOR/chrM/{sample}_{chrms}_SOR_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/SOR/chrM/{sample}_{chrms}_SOR_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - FS
#-------------------------------------------------------------------------------#
rule extract_males_call_VCF_autos_FS:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/FS/{chrmsA}/{chrmsA}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/FS/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_FS_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/FS/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_FS_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XnonPARs_FS:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/FS/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/FS/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_FS_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/FS/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_FS_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XPARs_FS:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/FS/chrX_PARs/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/FS/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_FS_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/FS/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_FS_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_Y_FS:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/FS/chrY/{chrmsY}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/FS/chrY/by_sample/{sample}_{chrmsY}_FS_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/FS/chrY/by_sample/{sample}_{chrmsY}_FS_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule get_performance_metrics_autos_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/FS/{chrms}/by_sample/{sample}_{chrms}_autos_FS_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/FS/autos/{sample}_{chrms}_autos_FS_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/FS/autos/{sample}_{chrms}_autos_FS_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XnonPARs_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/FS/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_FS_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/FS/chrX_nonPARs/{sample}_{chrms}_nonPARs_FS_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/FS/chrX_nonPARs/{sample}_{chrms}_nonPARs_FS_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XPARs_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/FS/chrX_PARs/by_sample/{sample}_{chrms}_PARs_FS_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/FS/chrX_PARs/{sample}_{chrms}_PARs_FS_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/FS/chrX_PARs/{sample}_{chrms}_PARs_FS_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/FS/chrY/by_sample/{sample}_{chrms}_FS_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/FS/chrY/{sample}_{chrms}_FS_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/FS/chrY/{sample}_{chrms}_FS_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule extract_males_call_VCF_M_FS:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/FS/chrM/{chrmsM}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/FS/chrM/by_sample/{sample}_{chrmsM}_FS_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/FS/chrM/by_sample/{sample}_{chrmsM}_FS_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_M_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/FS/chrM/by_sample/{sample}_{chrms}_FS_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/FS/chrM/{sample}_{chrms}_FS_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/FS/chrM/{sample}_{chrms}_FS_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - MQ
#-------------------------------------------------------------------------------#
rule extract_males_call_VCF_autos_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQ/{chrmsA}/{chrmsA}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQ/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQ_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQ/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQ_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XnonPARs_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQ/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQ/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQ_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQ/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQ_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XPARs_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQ/chrX_PARs/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQ/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQ_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQ/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQ_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_Y_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQ/chrY/{chrmsY}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQ/chrY/by_sample/{sample}_{chrmsY}_MQ_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQ/chrY/by_sample/{sample}_{chrmsY}_MQ_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule get_performance_metrics_autos_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQ/{chrms}/by_sample/{sample}_{chrms}_autos_MQ_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/MQ/autos/{sample}_{chrms}_autos_MQ_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/MQ/autos/{sample}_{chrms}_autos_MQ_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XnonPARs_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQ/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_MQ_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/MQ/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQ_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/MQ/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQ_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XPARs_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQ/chrX_PARs/by_sample/{sample}_{chrms}_PARs_MQ_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/MQ/chrX_PARs/{sample}_{chrms}_PARs_MQ_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/MQ/chrX_PARs/{sample}_{chrms}_PARs_MQ_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQ/chrY/by_sample/{sample}_{chrms}_MQ_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/MQ/chrY/{sample}_{chrms}_MQ_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/MQ/chrY/{sample}_{chrms}_MQ_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule extract_males_call_VCF_M_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQ/chrM/{chrmsM}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQ/chrM/by_sample/{sample}_{chrmsM}_MQ_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQ/chrM/by_sample/{sample}_{chrmsM}_MQ_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_M_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQ/chrM/by_sample/{sample}_{chrms}_MQ_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/MQ/chrM/{sample}_{chrms}_MQ_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/MQ/chrM/{sample}_{chrms}_MQ_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - MQRankSum
#-------------------------------------------------------------------------------#
rule extract_males_call_VCF_autos_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQRankSum/{chrmsA}/{chrmsA}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XnonPARs_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQRankSum/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XPARs_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQRankSum/chrX_PARs/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_Y_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQRankSum/chrY/{chrmsY}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQRankSum/chrY/by_sample/{sample}_{chrmsY}_MQRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQRankSum/chrY/by_sample/{sample}_{chrmsY}_MQRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule get_performance_metrics_autos_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQRankSum/{chrms}/by_sample/{sample}_{chrms}_autos_MQRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/MQRankSum/autos/{sample}_{chrms}_autos_MQRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/MQRankSum/autos/{sample}_{chrms}_autos_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XnonPARs_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQRankSum/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_MQRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/MQRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/MQRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XPARs_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQRankSum/chrX_PARs/by_sample/{sample}_{chrms}_PARs_MQRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/MQRankSum/chrX_PARs/{sample}_{chrms}_PARs_MQRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/MQRankSum/chrX_PARs/{sample}_{chrms}_PARs_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQRankSum/chrY/by_sample/{sample}_{chrms}_MQRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/MQRankSum/chrY/{sample}_{chrms}_MQRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/MQRankSum/chrY/{sample}_{chrms}_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule extract_males_call_VCF_M_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQRankSum/chrM/{chrmsM}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQRankSum/chrM/by_sample/{sample}_{chrmsM}_MQRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQRankSum/chrM/by_sample/{sample}_{chrmsM}_MQRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_M_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/MQRankSum/chrM/by_sample/{sample}_{chrms}_MQRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/MQRankSum/chrM/{sample}_{chrms}_MQRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/MQRankSum/chrM/{sample}_{chrms}_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - ReadPosRankSum
#-------------------------------------------------------------------------------#
rule extract_males_call_VCF_autos_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/ReadPosRankSum/{chrmsA}/{chrmsA}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/ReadPosRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_ReadPosRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/ReadPosRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_ReadPosRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XnonPARs_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_ReadPosRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_ReadPosRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XPARs_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrX_PARs/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_ReadPosRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_ReadPosRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_Y_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrY/{chrmsY}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrY/by_sample/{sample}_{chrmsY}_ReadPosRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrY/by_sample/{sample}_{chrmsY}_ReadPosRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule get_performance_metrics_autos_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/ReadPosRankSum/{chrms}/by_sample/{sample}_{chrms}_autos_ReadPosRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/ReadPosRankSum/autos/{sample}_{chrms}_autos_ReadPosRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/ReadPosRankSum/autos/{sample}_{chrms}_autos_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XnonPARs_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_ReadPosRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/ReadPosRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_ReadPosRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/ReadPosRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XPARs_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrX_PARs/by_sample/{sample}_{chrms}_PARs_ReadPosRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/ReadPosRankSum/chrX_PARs/{sample}_{chrms}_PARs_ReadPosRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/ReadPosRankSum/chrX_PARs/{sample}_{chrms}_PARs_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrY/by_sample/{sample}_{chrms}_ReadPosRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/ReadPosRankSum/chrY/{sample}_{chrms}_ReadPosRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/ReadPosRankSum/chrY/{sample}_{chrms}_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule extract_males_call_VCF_M_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrM/{chrmsM}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrM/by_sample/{sample}_{chrmsM}_ReadPosRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrM/by_sample/{sample}_{chrmsM}_ReadPosRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_M_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/ReadPosRankSum/chrM/by_sample/{sample}_{chrms}_ReadPosRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/ReadPosRankSum/chrM/{sample}_{chrms}_ReadPosRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/ReadPosRankSum/chrM/{sample}_{chrms}_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - DP
#-------------------------------------------------------------------------------#
rule extract_males_call_VCF_autos_DP:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/DP/{chrmsA}/{chrmsA}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/DP/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_DP_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/DP/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_DP_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XnonPARs_DP:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/DP/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/DP/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_DP_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/DP/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_DP_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XPARs_DP:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/DP/chrX_PARs/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/DP/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_DP_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/DP/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_DP_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_Y_DP:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/DP/chrY/{chrmsY}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/DP/chrY/by_sample/{sample}_{chrmsY}_DP_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/DP/chrY/by_sample/{sample}_{chrmsY}_DP_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule get_performance_metrics_autos_DP:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/DP/{chrms}/by_sample/{sample}_{chrms}_autos_DP_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/DP/autos/{sample}_{chrms}_autos_DP_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/DP/autos/{sample}_{chrms}_autos_DP_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XnonPARs_DP:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/DP/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_DP_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/DP/chrX_nonPARs/{sample}_{chrms}_nonPARs_DP_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/DP/chrX_nonPARs/{sample}_{chrms}_nonPARs_DP_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XPARs_DP:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/DP/chrX_PARs/by_sample/{sample}_{chrms}_PARs_DP_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/DP/chrX_PARs/{sample}_{chrms}_PARs_DP_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/DP/chrX_PARs/{sample}_{chrms}_PARs_DP_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y_DP:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/DP/chrY/by_sample/{sample}_{chrms}_DP_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/DP/chrY/{sample}_{chrms}_DP_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/DP/chrY/{sample}_{chrms}_DP_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule extract_males_call_VCF_M_DP:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/DP/chrM/{chrmsM}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/DP/chrM/by_sample/{sample}_{chrmsM}_DP_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/DP/chrM/by_sample/{sample}_{chrmsM}_DP_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_M_DP:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/DP/chrM/by_sample/{sample}_{chrms}_DP_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/DP/chrM/{sample}_{chrms}_DP_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/DP/chrM/{sample}_{chrms}_DP_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - AN
#-------------------------------------------------------------------------------#
rule extract_males_call_VCF_autos_AN:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/AN/{chrmsA}/{chrmsA}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/AN/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_AN_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/AN/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_AN_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XnonPARs_AN:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/AN/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/AN/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_AN_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/AN/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_AN_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XPARs_AN:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/AN/chrX_PARs/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/AN/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_AN_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/AN/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_AN_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_Y_AN:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/AN/chrY/{chrmsY}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/AN/chrY/by_sample/{sample}_{chrmsY}_AN_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/AN/chrY/by_sample/{sample}_{chrmsY}_AN_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule get_performance_metrics_autos_AN:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/AN/{chrms}/by_sample/{sample}_{chrms}_autos_AN_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/AN/autos/{sample}_{chrms}_autos_AN_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/AN/autos/{sample}_{chrms}_autos_AN_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XnonPARs_AN:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/AN/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_AN_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/AN/chrX_nonPARs/{sample}_{chrms}_nonPARs_AN_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/AN/chrX_nonPARs/{sample}_{chrms}_nonPARs_AN_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XPARs_AN:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/AN/chrX_PARs/by_sample/{sample}_{chrms}_PARs_AN_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/AN/chrX_PARs/{sample}_{chrms}_PARs_AN_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/AN/chrX_PARs/{sample}_{chrms}_PARs_AN_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y_AN:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/AN/chrY/by_sample/{sample}_{chrms}_AN_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/AN/chrY/{sample}_{chrms}_AN_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/AN/chrY/{sample}_{chrms}_AN_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule extract_males_call_VCF_M_AN:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/AN/chrM/{chrmsM}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/AN/chrM/by_sample/{sample}_{chrmsM}_AN_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/AN/chrM/by_sample/{sample}_{chrmsM}_AN_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_M_AN:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/AN/chrM/by_sample/{sample}_{chrms}_AN_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/AN/chrM/{sample}_{chrms}_AN_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/AN/chrM/{sample}_{chrms}_AN_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """
'''


'''
#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - No filters
#-------------------------------------------------------------------------------#
rule extract_males_call_VCF_autos_no_filter:
    input:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/males/{chrmsA}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/males/by_sample/{sample}_{chrmsA}_autos_no_filter")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/males/by_sample/{sample}_{chrmsA}_autos_no_filter.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_nonPARs_nofilter:
    input:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/males/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/males/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_nonPARs.vcf.gz")
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

rule extract_PARs_nofilter:
    input:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/males/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/males/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_PARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrmsXPARs}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.ilist} """
        """-O {output} """
#"""-L {params.chrms} -XL {params.ilist} """

rule extract_males_call_VCF_XnonPARs_no_filter:
    input:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/males/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/males/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_no_filter")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/males/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_no_filter.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XPARs_no_filter:
    input:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/males/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/males/by_sample/{sample}_{chrmsXPARs}_PARs_no_filter")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/males/by_sample/{sample}_{chrmsXPARs}_PARs_no_filter.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_Y_no_filter:
    input:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/males/{chrmsY}_GRCh38_YPARsMasked_gatk_diploid_called_raw.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/males/by_sample/{sample}_{chrmsY}_no_filter")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/males/by_sample/{sample}_{chrmsY}_no_filter.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_autos_no_filter:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/males/by_sample/{sample}_{chrms}_autos_no_filter.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/no_filter/autos/{sample}_{chrms}_autos_no_filter_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/no_filter/autos/{sample}_{chrms}_autos_no_filter_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XnonPARs_no_filter:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/males/by_sample/{sample}_{chrms}_nonPARs_no_filter.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/no_filter/chrX_nonPARs/{sample}_{chrms}_nonPARs_no_filter_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/no_filter/chrX_nonPARs/{sample}_{chrms}_nonPARs_no_filter_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XPARs_no_filter:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/males/by_sample/{sample}_{chrms}_PARs_no_filter.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/no_filter/chrX_PARs/{sample}_{chrms}_PARs_no_filter_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/no_filter/chrX_PARs/{sample}_{chrms}_PARs_no_filter_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y_no_filter:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/males/by_sample/{sample}_{chrms}_no_filter.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/no_filter/chrY/{sample}_{chrms}_no_filter_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/no_filter/chrY/{sample}_{chrms}_no_filter_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """
'''



#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#
# MALES DEFAULT
#-------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - ALL
#-------------------------------------------------------------------------------#
rule default_extract_males_call_VCF_autos:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/default/all/{chrmsA}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrmsA}_autos")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrmsA}_autos.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_XnonPARs:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/default/all/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_XPARs:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/default/all/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrmsXPARs}_PARs")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrmsXPARs}_PARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_Y:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/default/all/{chrmsY}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrmsY}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrmsY}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_M:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/default/all/{chrmsM}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrmsM}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrmsM}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_get_performance_metrics_autos:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrms}_autos.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/default/autos/{sample}_{chrms}_autos_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/default/autos/{sample}_{chrms}_autos_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_XnonPARs:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrms}_nonPARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/default/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/default/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_XPARs:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrms}_PARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/default/chrX_PARs/{sample}_{chrms}_PARs_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/default/chrX_PARs/{sample}_{chrms}_PARs_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_Y:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrms}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/default/chrY/{sample}_{chrms}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/default/chrY/{sample}_{chrms}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_M:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/males/default/all/by_sample/{sample}_{chrms}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/default/chrM/{sample}_{chrms}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/males/default/chrM/{sample}_{chrms}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """



#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#                                   FEMALES                                    #
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Step: Prep golden (simulated) VCFs
#------------------------------------------------------------------------------#
rule females_unzip_golden_VCFs_autos:
    input:
        os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsA}_golden.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsA}_golden.vcf")
    shell:
        """
        gunzip {input}
        """

rule females_unzip_golden_VCFs_X:
    input:
        os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsXnonPARs}_nonPARs_golden.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsXnonPARs}_nonPARs_golden.vcf")
    shell:
        """
        gunzip {input}
        """

rule females_unzip_golden_VCFs_X_PARs:
    input:
        os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsXPARs}_PARs_golden.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsXPARs}_PARs_golden.vcf")
    shell:
        """
        gunzip {input}
        """

rule females_unzip_golden_VCFs_M:
    input:
        os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrmsM}_golden.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrmsM}_golden.vcf")
    shell:
        """
        gunzip {input}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - ALL
#-------------------------------------------------------------------------------#
rule females_extract_females_call_VCF_autos:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/all/{chrmsA}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/all/by_sample/{sample}_{chrmsA}_autos")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/all/by_sample/{sample}_{chrmsA}_autos.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_extract_females_call_VCF_XnonPARs:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/all/{chrmsXnonPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_extract_females_call_VCF_XPARs:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/all/{chrmsXPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/all/by_sample/{sample}_{chrmsXPARs}_PARs")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/all/by_sample/{sample}_{chrmsXPARs}_PARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_extract_females_call_VCF_M:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/all/{chrmsM}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/all/by_sample/{sample}_{chrmsM}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/all/by_sample/{sample}_{chrmsM}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule females_get_performance_metrics_autos:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/all/by_sample/{sample}_{chrms}_autos.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/autos/{sample}_{chrms}_autos_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/autos/{sample}_{chrms}_autos_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_get_performance_metrics_XnonPARs:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/all/by_sample/{sample}_{chrms}_nonPARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_get_performance_metrics_XPARs:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/all/by_sample/{sample}_{chrms}_PARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/chrX_PARs/{sample}_{chrms}_PARs_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/chrX_PARs/{sample}_{chrms}_PARs_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_get_performance_metrics_M:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/all/by_sample/{sample}_{chrms}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/chrM/{sample}_{chrms}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/chrM/{sample}_{chrms}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

'''
# FILTERS #
#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - QD
#-------------------------------------------------------------------------------#
rule female_extract_call_VCF_autos_QD:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QD/{chrmsA}/{chrmsA}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QD/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QD_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QD/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QD_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XnonPARs_QD:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QD/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QD/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QD_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QD/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QD_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XPARs_QD:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QD/chrX_PARs/{chrmsXPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QD/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QD_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QD/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QD_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """



rule female_get_performance_metrics_autos_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QD/{chrms}/by_sample/{sample}_{chrms}_autos_QD_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/QD/autos/{sample}_{chrms}_autos_QD_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/QD/autos/{sample}_{chrms}_autos_QD_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XnonPARs_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QD/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_QD_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/QD/chrX_nonPARs/{sample}_{chrms}_nonPARs_QD_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/QD/chrX_nonPARs/{sample}_{chrms}_nonPARs_QD_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XPARs_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QD/chrX_PARs/by_sample/{sample}_{chrms}_PARs_QD_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/QD/chrX_PARs/{sample}_{chrms}_PARs_QD_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/QD/chrX_PARs/{sample}_{chrms}_PARs_QD_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule female_extract_call_VCF_M_QD:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QD/chrM/{chrmsM}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QD/chrM/by_sample/{sample}_{chrmsM}_QD_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QD/chrM/by_sample/{sample}_{chrmsM}_QD_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_get_performance_metrics_M_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QD/chrM/by_sample/{sample}_{chrms}_QD_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/QD/chrM/{sample}_{chrms}_QD_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/QD/chrM/{sample}_{chrms}_QD_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - QUAL
#-------------------------------------------------------------------------------#
rule female_extract_call_VCF_autos_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QUAL/{chrmsA}/{chrmsA}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QUAL/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QUAL_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QUAL/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QUAL_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XnonPARs_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QUAL/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QUAL/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QUAL_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QUAL/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QUAL_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XPARs_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QUAL/chrX_PARs/{chrmsXPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QUAL/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QUAL_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QUAL/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QUAL_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule female_get_performance_metrics_autos_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QUAL/{chrms}/by_sample/{sample}_{chrms}_autos_QUAL_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/QUAL/autos/{sample}_{chrms}_autos_QUAL_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/QUAL/autos/{sample}_{chrms}_autos_QUAL_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XnonPARs_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QUAL/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_QUAL_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/QUAL/chrX_nonPARs/{sample}_{chrms}_nonPARs_QUAL_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/QUAL/chrX_nonPARs/{sample}_{chrms}_nonPARs_QUAL_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XPARs_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QUAL/chrX_PARs/by_sample/{sample}_{chrms}_PARs_QUAL_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/QUAL/chrX_PARs/{sample}_{chrms}_PARs_QUAL_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/QUAL/chrX_PARs/{sample}_{chrms}_PARs_QUAL_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule female_extract_call_VCF_M_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QUAL/chrM/{chrmsM}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QUAL/chrM/by_sample/{sample}_{chrmsM}_QUAL_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QUAL/chrM/by_sample/{sample}_{chrmsM}_QUAL_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_get_performance_metrics_M_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/QUAL/chrM/by_sample/{sample}_{chrms}_QUAL_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/QUAL/chrM/{sample}_{chrms}_QUAL_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/QUAL/chrM/{sample}_{chrms}_QUAL_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - SOR
#-------------------------------------------------------------------------------#
rule female_extract_call_VCF_autos_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/SOR/{chrmsA}/{chrmsA}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/SOR/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_SOR_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/SOR/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_SOR_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XnonPARs_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/SOR/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/SOR/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_SOR_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/SOR/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_SOR_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XPARs_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/SOR/chrX_PARs/{chrmsXPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/SOR/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_SOR_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/SOR/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_SOR_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """



rule female_get_performance_metrics_autos_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/SOR/{chrms}/by_sample/{sample}_{chrms}_autos_SOR_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/SOR/autos/{sample}_{chrms}_autos_SOR_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/SOR/autos/{sample}_{chrms}_autos_SOR_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XnonPARs_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/SOR/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_SOR_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/SOR/chrX_nonPARs/{sample}_{chrms}_nonPARs_SOR_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/SOR/chrX_nonPARs/{sample}_{chrms}_nonPARs_SOR_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XPARs_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/SOR/chrX_PARs/by_sample/{sample}_{chrms}_PARs_SOR_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/SOR/chrX_PARs/{sample}_{chrms}_PARs_SOR_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/SOR/chrX_PARs/{sample}_{chrms}_PARs_SOR_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule female_extract_call_VCF_M_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/SOR/chrM/{chrmsM}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/SOR/chrM/by_sample/{sample}_{chrmsM}_SOR_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/SOR/chrM/by_sample/{sample}_{chrmsM}_SOR_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_get_performance_metrics_M_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/SOR/chrM/by_sample/{sample}_{chrms}_SOR_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/SOR/chrM/{sample}_{chrms}_SOR_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/SOR/chrM/{sample}_{chrms}_SOR_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - FS
#-------------------------------------------------------------------------------#
rule female_extract_call_VCF_autos_FS:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/FS/{chrmsA}/{chrmsA}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/FS/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_FS_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/FS/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_FS_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XnonPARs_FS:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/FS/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/FS/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_FS_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/FS/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_FS_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XPARs_FS:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/FS/chrX_PARs/{chrmsXPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/FS/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_FS_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/FS/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_FS_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule female_get_performance_metrics_autos_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/FS/{chrms}/by_sample/{sample}_{chrms}_autos_FS_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/FS/autos/{sample}_{chrms}_autos_FS_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/FS/autos/{sample}_{chrms}_autos_FS_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XnonPARs_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/FS/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_FS_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/FS/chrX_nonPARs/{sample}_{chrms}_nonPARs_FS_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/FS/chrX_nonPARs/{sample}_{chrms}_nonPARs_FS_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XPARs_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/FS/chrX_PARs/by_sample/{sample}_{chrms}_PARs_FS_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/FS/chrX_PARs/{sample}_{chrms}_PARs_FS_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/FS/chrX_PARs/{sample}_{chrms}_PARs_FS_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule female_extract_call_VCF_M_FS:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/FS/chrM/{chrmsM}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/FS/chrM/by_sample/{sample}_{chrmsM}_FS_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/FS/chrM/by_sample/{sample}_{chrmsM}_FS_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_get_performance_metrics_M_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/FS/chrM/by_sample/{sample}_{chrms}_FS_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/FS/chrM/{sample}_{chrms}_FS_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/FS/chrM/{sample}_{chrms}_FS_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - MQ
#-------------------------------------------------------------------------------#
rule female_extract_call_VCF_autos_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQ/{chrmsA}/{chrmsA}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQ/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQ_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQ/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQ_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XnonPARs_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQ/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQ/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQ_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQ/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQ_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XPARs_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQ/chrX_PARs/{chrmsXPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQ/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQ_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQ/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQ_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """



rule female_get_performance_metrics_autos_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQ/{chrms}/by_sample/{sample}_{chrms}_autos_MQ_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/MQ/autos/{sample}_{chrms}_autos_MQ_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/MQ/autos/{sample}_{chrms}_autos_MQ_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XnonPARs_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQ/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_MQ_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/MQ/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQ_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/MQ/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQ_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XPARs_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQ/chrX_PARs/by_sample/{sample}_{chrms}_PARs_MQ_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/MQ/chrX_PARs/{sample}_{chrms}_PARs_MQ_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/MQ/chrX_PARs/{sample}_{chrms}_PARs_MQ_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule female_extract_call_VCF_M_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQ/chrM/{chrmsM}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQ/chrM/by_sample/{sample}_{chrmsM}_MQ_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQ/chrM/by_sample/{sample}_{chrmsM}_MQ_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_get_performance_metrics_M_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQ/chrM/by_sample/{sample}_{chrms}_MQ_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/MQ/chrM/{sample}_{chrms}_MQ_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/MQ/chrM/{sample}_{chrms}_MQ_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - MQRankSum
#-------------------------------------------------------------------------------#
rule female_extract_call_VCF_autos_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQRankSum/{chrmsA}/{chrmsA}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XnonPARs_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQRankSum/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XPARs_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQRankSum/chrX_PARs/{chrmsXPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule female_get_performance_metrics_autos_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQRankSum/{chrms}/by_sample/{sample}_{chrms}_autos_MQRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/MQRankSum/autos/{sample}_{chrms}_autos_MQRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/MQRankSum/autos/{sample}_{chrms}_autos_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XnonPARs_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQRankSum/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_MQRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/MQRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/MQRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XPARs_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQRankSum/chrX_PARs/by_sample/{sample}_{chrms}_PARs_MQRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/MQRankSum/chrX_PARs/{sample}_{chrms}_PARs_MQRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/MQRankSum/chrX_PARs/{sample}_{chrms}_PARs_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule female_extract_call_VCF_M_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQRankSum/chrM/{chrmsM}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQRankSum/chrM/by_sample/{sample}_{chrmsM}_MQRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQRankSum/chrM/by_sample/{sample}_{chrmsM}_MQRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_get_performance_metrics_M_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/MQRankSum/chrM/by_sample/{sample}_{chrms}_MQRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/MQRankSum/chrM/{sample}_{chrms}_MQRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/MQRankSum/chrM/{sample}_{chrms}_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - ReadPosRankSum
#-------------------------------------------------------------------------------#
rule female_extract_call_VCF_autos_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/ReadPosRankSum/{chrmsA}/{chrmsA}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/ReadPosRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_ReadPosRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/ReadPosRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_ReadPosRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XnonPARs_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/ReadPosRankSum/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/ReadPosRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_ReadPosRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/ReadPosRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_ReadPosRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XPARs_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/ReadPosRankSum/chrX_PARs/{chrmsXPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/ReadPosRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_ReadPosRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/ReadPosRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_ReadPosRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule female_get_performance_metrics_autos_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/ReadPosRankSum/{chrms}/by_sample/{sample}_{chrms}_autos_ReadPosRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/ReadPosRankSum/autos/{sample}_{chrms}_autos_ReadPosRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/ReadPosRankSum/autos/{sample}_{chrms}_autos_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XnonPARs_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/ReadPosRankSum/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_ReadPosRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/ReadPosRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_ReadPosRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/ReadPosRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XPARs_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/ReadPosRankSum/chrX_PARs/by_sample/{sample}_{chrms}_PARs_ReadPosRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/ReadPosRankSum/chrX_PARs/{sample}_{chrms}_PARs_ReadPosRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/ReadPosRankSum/chrX_PARs/{sample}_{chrms}_PARs_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule female_extract_call_VCF_M_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/ReadPosRankSum/chrM/{chrmsM}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/ReadPosRankSum/chrM/by_sample/{sample}_{chrmsM}_ReadPosRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/ReadPosRankSum/chrM/by_sample/{sample}_{chrmsM}_ReadPosRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_get_performance_metrics_M_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/ReadPosRankSum/chrM/by_sample/{sample}_{chrms}_ReadPosRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/ReadPosRankSum/chrM/{sample}_{chrms}_ReadPosRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/ReadPosRankSum/chrM/{sample}_{chrms}_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - DP
#-------------------------------------------------------------------------------#
rule female_extract_call_VCF_autos_DP:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/DP/{chrmsA}/{chrmsA}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/DP/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_DP_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/DP/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_DP_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XnonPARs_DP:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/DP/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/DP/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_DP_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/DP/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_DP_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XPARs_DP:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/DP/chrX_PARs/{chrmsXPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/DP/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_DP_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/DP/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_DP_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """



rule female_get_performance_metrics_autos_DP:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/DP/{chrms}/by_sample/{sample}_{chrms}_autos_DP_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/DP/autos/{sample}_{chrms}_autos_DP_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/DP/autos/{sample}_{chrms}_autos_DP_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XnonPARs_DP:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/DP/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_DP_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/DP/chrX_nonPARs/{sample}_{chrms}_nonPARs_DP_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/DP/chrX_nonPARs/{sample}_{chrms}_nonPARs_DP_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XPARs_DP:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/DP/chrX_PARs/by_sample/{sample}_{chrms}_PARs_DP_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/DP/chrX_PARs/{sample}_{chrms}_PARs_DP_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/DP/chrX_PARs/{sample}_{chrms}_PARs_DP_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule female_extract_call_VCF_M_DP:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/DP/chrM/{chrmsM}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/DP/chrM/by_sample/{sample}_{chrmsM}_DP_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/DP/chrM/by_sample/{sample}_{chrmsM}_DP_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_get_performance_metrics_M_DP:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/DP/chrM/by_sample/{sample}_{chrms}_DP_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/DP/chrM/{sample}_{chrms}_DP_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/DP/chrM/{sample}_{chrms}_DP_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - AN
#-------------------------------------------------------------------------------#
rule female_extract_call_VCF_autos_AN:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/AN/{chrmsA}/{chrmsA}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/AN/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_AN_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/AN/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_AN_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XnonPARs_AN:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/AN/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/AN/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_AN_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/AN/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_AN_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XPARs_AN:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/AN/chrX_PARs/{chrmsXPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/AN/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_AN_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/AN/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_AN_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule female_get_performance_metrics_autos_AN:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/AN/{chrms}/by_sample/{sample}_{chrms}_autos_AN_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/AN/autos/{sample}_{chrms}_autos_AN_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/AN/autos/{sample}_{chrms}_autos_AN_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XnonPARs_AN:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/AN/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_AN_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/AN/chrX_nonPARs/{sample}_{chrms}_nonPARs_AN_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/AN/chrX_nonPARs/{sample}_{chrms}_nonPARs_AN_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XPARs_AN:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/AN/chrX_PARs/by_sample/{sample}_{chrms}_PARs_AN_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/AN/chrX_PARs/{sample}_{chrms}_PARs_AN_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/AN/chrX_PARs/{sample}_{chrms}_PARs_AN_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule female_extract_call_VCF_M_AN:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/AN/chrM/{chrmsM}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/AN/chrM/by_sample/{sample}_{chrmsM}_AN_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/AN/chrM/by_sample/{sample}_{chrmsM}_AN_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_get_performance_metrics_M_AN:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/AN/chrM/by_sample/{sample}_{chrms}_AN_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/AN/chrM/{sample}_{chrms}_AN_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/AN/chrM/{sample}_{chrms}_AN_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """
'''


'''

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - No filters
#-------------------------------------------------------------------------------#
rule females_extract_females_call_VCF_autos_no_filter:
    input:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/{chrmsA}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/by_sample/{sample}_{chrmsA}_autos_no_filter")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/by_sample/{sample}_{chrmsA}_autos_no_filter.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_extract_nonPARs_nofilter:
    input:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_nonPARs.vcf.gz")
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

rule females_extract_PARs_nofilter:
    input:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_PARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrmsXPARs}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.ilist} """
        """-O {output} """
#"""-L {params.chrms} -XL {params.ilist} """

rule females_extract_females_call_VCF_XnonPARs_no_filter:
    input:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_no_filter")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_no_filter.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_extract_females_call_VCF_XPARs_no_filter:
    input:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/by_sample/{sample}_{chrmsXPARs}_PARs_no_filter")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/by_sample/{sample}_{chrmsXPARs}_PARs_no_filter.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_get_performance_metrics_autos_no_filter:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/by_sample/{sample}_{chrms}_autos_no_filter.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/no_filter/autos/{sample}_{chrms}_autos_no_filter_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/no_filter/autos/{sample}_{chrms}_autos_no_filter_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_get_performance_metrics_XnonPARs_no_filter:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/by_sample/{sample}_{chrms}_nonPARs_no_filter.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/no_filter/chrX_nonPARs/{sample}_{chrms}_nonPARs_no_filter_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/no_filter/chrX_nonPARs/{sample}_{chrms}_nonPARs_no_filter_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_get_performance_metrics_XPARs_no_filter:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/by_sample/{sample}_{chrms}_PARs_no_filter.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/no_filter/chrX_PARs/{sample}_{chrms}_PARs_no_filter_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/no_filter/chrX_PARs/{sample}_{chrms}_PARs_no_filter_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """
'''

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#                                   FEMALES                                    #
#                                   DEFAULT                                    #
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#


#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - ALL
#-------------------------------------------------------------------------------#
rule females_default_extract_females_call_VCF_autos:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/all/{chrmsA}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/all/by_sample/{sample}_{chrmsA}_autos")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/all/by_sample/{sample}_{chrmsA}_autos.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_default_extract_females_call_VCF_XnonPARs:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/all/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_default_extract_females_call_VCF_XPARs:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/all/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/all/by_sample/{sample}_{chrmsXPARs}_PARs")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/all/by_sample/{sample}_{chrmsXPARs}_PARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_default_extract_females_call_VCF_M:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/all/{chrmsM}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/all/by_sample/{sample}_{chrmsM}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/all/by_sample/{sample}_{chrmsM}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_default_get_performance_metrics_autos:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/all/by_sample/{sample}_{chrms}_autos.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/autos/{sample}_{chrms}_autos_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/autos/{sample}_{chrms}_autos_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_default_get_performance_metrics_XnonPARs:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/all/by_sample/{sample}_{chrms}_nonPARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_default_get_performance_metrics_XPARs:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/all/by_sample/{sample}_{chrms}_PARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/chrX_PARs/{sample}_{chrms}_PARs_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/chrX_PARs/{sample}_{chrms}_PARs_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_default_get_performance_metrics_M:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/all/by_sample/{sample}_{chrms}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/chrM/{sample}_{chrms}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/chrM/{sample}_{chrms}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """



# For females, I dont think we are testing different filters. We just want to
# compare SCC to default

'''
#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - QD
#-------------------------------------------------------------------------------#
rule females_default_extract_females_call_VCF_autos_QD:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/QD/{chrmsA}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_QD.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/QD/by_sample/{sample}_{chrmsA}_autos_QD")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/QD/by_sample/{sample}_{chrmsA}_autos_QD.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_default_extract_females_call_VCF_XnonPARs_QD:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/QD/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_QD_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/QD/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QD")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/QD/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QD.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_default_extract_females_call_VCF_XPARs_QD:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/QD/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_QD_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/QD/by_sample/{sample}_{chrmsXPARs}_PARs_QD")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/QD/by_sample/{sample}_{chrmsXPARs}_PARs_QD.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule females_default_get_performance_metrics_autos_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/QD/by_sample/{sample}_{chrms}_autos_QD.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/QD/autos/{sample}_{chrms}_autos_QD_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/QD/autos/{sample}_{chrms}_autos_QD_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_default_get_performance_metrics_XnonPARs_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/QD/by_sample/{sample}_{chrms}_nonPARs_QD.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/QD/chrX_nonPARs/{sample}_{chrms}_nonPARs_QD_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/QD/chrX_nonPARs/{sample}_{chrms}_nonPARs_QD_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_default_get_performance_metrics_XPARs_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/QD/by_sample/{sample}_{chrms}_PARs_QD.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/QD/chrX_PARs/{sample}_{chrms}_PARs_QD_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/QD/chrX_PARs/{sample}_{chrms}_PARs_QD_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - QUAL
#-------------------------------------------------------------------------------#
rule females_default_extract_females_call_VCF_autos_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/QUAL/{chrmsA}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_QUAL.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/QUAL/by_sample/{sample}_{chrmsA}_autos_QUAL")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/QUAL/by_sample/{sample}_{chrmsA}_autos_QUAL.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_default_extract_females_call_VCF_XnonPARs_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/QUAL/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_QUAL_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/QUAL/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QUAL")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/QUAL/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QUAL.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_default_extract_females_call_VCF_XPARs_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/QUAL/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_QUAL_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/QUAL/by_sample/{sample}_{chrmsXPARs}_PARs_QUAL")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/QUAL/by_sample/{sample}_{chrmsXPARs}_PARs_QUAL.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule females_default_get_performance_metrics_autos_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/QUAL/by_sample/{sample}_{chrms}_autos_QUAL.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/QUAL/autos/{sample}_{chrms}_autos_QUAL_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/QUAL/autos/{sample}_{chrms}_autos_QUAL_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_default_get_performance_metrics_XnonPARs_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/QUAL/by_sample/{sample}_{chrms}_nonPARs_QUAL.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/QUAL/chrX_nonPARs/{sample}_{chrms}_nonPARs_QUAL_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/QUAL/chrX_nonPARs/{sample}_{chrms}_nonPARs_QUAL_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_default_get_performance_metrics_XPARs_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/QUAL/by_sample/{sample}_{chrms}_PARs_QUAL.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/QUAL/chrX_PARs/{sample}_{chrms}_PARs_QUAL_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/QUAL/chrX_PARs/{sample}_{chrms}_PARs_QUAL_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - SOR
#-------------------------------------------------------------------------------#
rule females_default_extract_females_call_VCF_autos_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/SOR/{chrmsA}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_SOR.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/SOR/by_sample/{sample}_{chrmsA}_autos_SOR")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/SOR/by_sample/{sample}_{chrmsA}_autos_SOR.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_default_extract_females_call_VCF_XnonPARs_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/SOR/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_SOR_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/SOR/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_SOR")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/SOR/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_SOR.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_default_extract_females_call_VCF_XPARs_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/SOR/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_SOR_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/SOR/by_sample/{sample}_{chrmsXPARs}_PARs_SOR")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/SOR/by_sample/{sample}_{chrmsXPARs}_PARs_SOR.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule females_default_get_performance_metrics_autos_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/SOR/by_sample/{sample}_{chrms}_autos_SOR.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/SOR/autos/{sample}_{chrms}_autos_SOR_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/SOR/autos/{sample}_{chrms}_autos_SOR_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_default_get_performance_metrics_XnonPARs_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/SOR/by_sample/{sample}_{chrms}_nonPARs_SOR.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/SOR/chrX_nonPARs/{sample}_{chrms}_nonPARs_SOR_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/SOR/chrX_nonPARs/{sample}_{chrms}_nonPARs_SOR_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_default_get_performance_metrics_XPARs_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/SOR/by_sample/{sample}_{chrms}_PARs_SOR.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/SOR/chrX_PARs/{sample}_{chrms}_PARs_SOR_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/SOR/chrX_PARs/{sample}_{chrms}_PARs_SOR_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - FS
#-------------------------------------------------------------------------------#
rule females_default_extract_females_call_VCF_autos_FS:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/FS/{chrmsA}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_FS.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/FS/by_sample/{sample}_{chrmsA}_autos_FS")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/FS/by_sample/{sample}_{chrmsA}_autos_FS.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_default_extract_females_call_VCF_XnonPARs_FS:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/FS/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_FS_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/FS/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_FS")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/FS/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_FS.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_default_extract_females_call_VCF_XPARs_FS:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/FS/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_FS_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/FS/by_sample/{sample}_{chrmsXPARs}_PARs_FS")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/FS/by_sample/{sample}_{chrmsXPARs}_PARs_FS.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule females_default_get_performance_metrics_autos_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/FS/by_sample/{sample}_{chrms}_autos_FS.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/FS/autos/{sample}_{chrms}_autos_FS_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/FS/autos/{sample}_{chrms}_autos_FS_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_default_get_performance_metrics_XnonPARs_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/FS/by_sample/{sample}_{chrms}_nonPARs_FS.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/FS/chrX_nonPARs/{sample}_{chrms}_nonPARs_FS_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/FS/chrX_nonPARs/{sample}_{chrms}_nonPARs_FS_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_default_get_performance_metrics_XPARs_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/FS/by_sample/{sample}_{chrms}_PARs_FS.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/FS/chrX_PARs/{sample}_{chrms}_PARs_FS_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/FS/chrX_PARs/{sample}_{chrms}_PARs_FS_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - MQ
#-------------------------------------------------------------------------------#
rule females_default_extract_females_call_VCF_autos_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/MQ/{chrmsA}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_MQ.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/MQ/by_sample/{sample}_{chrmsA}_autos_MQ")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/MQ/by_sample/{sample}_{chrmsA}_autos_MQ.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_default_extract_females_call_VCF_XnonPARs_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/MQ/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_MQ_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/MQ/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQ")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/MQ/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQ.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_default_extract_females_call_VCF_XnonPARs_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/MQ/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_MQ_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/MQ/by_sample/{sample}_{chrmsXPARs}_PARs_MQ")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/MQ/by_sample/{sample}_{chrmsXPARs}_PARs_MQ.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule females_default_get_performance_metrics_autos_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/MQ/by_sample/{sample}_{chrms}_autos_MQ.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/MQ/autos/{sample}_{chrms}_autos_MQ_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/MQ/autos/{sample}_{chrms}_autos_MQ_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_default_get_performance_metrics_XnonPARs_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/MQ/by_sample/{sample}_{chrms}_nonPARs_MQ.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/MQ/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQ_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/MQ/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQ_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_default_get_performance_metrics_XPARs_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/MQ/by_sample/{sample}_{chrms}_PARs_MQ.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/MQ/chrX_PARs/{sample}_{chrms}_PARs_MQ_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/MQ/chrX_PARs/{sample}_{chrms}_PARs_MQ_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - MQRankSum
#-------------------------------------------------------------------------------#
rule females_default_extract_females_call_VCF_autos_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/MQRankSum/{chrmsA}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_MQRankSum.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/MQRankSum/by_sample/{sample}_{chrmsA}_autos_MQRankSum")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/MQRankSum/by_sample/{sample}_{chrmsA}_autos_MQRankSum.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_default_extract_females_call_VCF_XnonPARs_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/MQRankSum/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_MQRankSum_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/MQRankSum/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQRankSum")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/MQRankSum/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQRankSum.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_default_extract_females_call_VCF_XPARs_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/MQRankSum/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_MQRankSum_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/MQRankSum/by_sample/{sample}_{chrmsXPARs}_PARs_MQRankSum")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/MQRankSum/by_sample/{sample}_{chrmsXPARs}_PARs_MQRankSum.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule females_default_get_performance_metrics_autos_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/MQRankSum/by_sample/{sample}_{chrms}_autos_MQRankSum.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/MQRankSum/autos/{sample}_{chrms}_autos_MQRankSum_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/MQRankSum/autos/{sample}_{chrms}_autos_MQRankSum_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_default_get_performance_metrics_XnonPARs_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/MQRankSum/by_sample/{sample}_{chrms}_nonPARs_MQRankSum.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/MQRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQRankSum_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/MQRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQRankSum_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_default_get_performance_metrics_XPARs_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/MQRankSum/by_sample/{sample}_{chrms}_PARs_MQRankSum.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/MQRankSum/chrX_PARs/{sample}_{chrms}_PARs_MQRankSum_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/MQRankSum/chrX_PARs/{sample}_{chrms}_PARs_MQRankSum_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - ReadPosRankSum
#-------------------------------------------------------------------------------#
rule females_default_extract_females_call_VCF_autos_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/ReadPosRankSum/{chrmsA}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_ReadPosRankSum.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/ReadPosRankSum/by_sample/{sample}_{chrmsA}_autos_ReadPosRankSum")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/ReadPosRankSum/by_sample/{sample}_{chrmsA}_autos_ReadPosRankSum.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_default_extract_females_call_VCF_XnonPARs_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/ReadPosRankSum/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_ReadPosRankSum_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/ReadPosRankSum/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_ReadPosRankSum")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/ReadPosRankSum/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_ReadPosRankSum.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_default_extract_females_call_VCF_XPARs_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/ReadPosRankSum/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_ReadPosRankSum_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/ReadPosRankSum/by_sample/{sample}_{chrmsXPARs}_PARs_ReadPosRankSum")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/ReadPosRankSum/by_sample/{sample}_{chrmsXPARs}_PARs_ReadPosRankSum.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule females_default_get_performance_metrics_autos_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/ReadPosRankSum/by_sample/{sample}_{chrms}_autos_ReadPosRankSum.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/ReadPosRankSum/autos/{sample}_{chrms}_autos_ReadPosRankSum_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/ReadPosRankSum/autos/{sample}_{chrms}_autos_ReadPosRankSum_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_default_get_performance_metrics_XnonPARs_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/ReadPosRankSum/by_sample/{sample}_{chrms}_nonPARs_ReadPosRankSum.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/ReadPosRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_ReadPosRankSum_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/ReadPosRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_ReadPosRankSum_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_default_get_performance_metrics_XPARs_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/ASN/females/default/ReadPosRankSum/by_sample/{sample}_{chrms}_PARs_ReadPosRankSum.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/ReadPosRankSum/chrX_PARs/{sample}_{chrms}_PARs_ReadPosRankSum_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/ReadPosRankSum/chrX_PARs/{sample}_{chrms}_PARs_ReadPosRankSum_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - No filters
#-------------------------------------------------------------------------------#
rule females_default_extract_females_call_VCF_autos_no_filter:
    input:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/default/{chrmsA}_GRCh38_default_gatk_diploid_called_raw_SNPs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/default/by_sample/{sample}_{chrmsA}_autos_no_filter")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/default/by_sample/{sample}_{chrmsA}_autos_no_filter.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_default_extract_nonPARs_nofilter:
    input:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/default/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_raw_SNPs.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/default/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_raw_SNPs_nonPARs.vcf.gz")
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

rule females_default_extract_PARs_nofilter:
    input:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/default/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_raw_SNPs.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/default/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_raw_SNPs_PARs.vcf.gz")
    params:
        ref = config["genome_paths"]["Ref_GRCh38_Y_PARsMasked"],
        chrms = "{chrmsXPARs}",
        ilist = config["par_intervals"]
    shell:
        """gatk SelectVariants """
        """-R {params.ref} """
        """-V {input} """
        """-L {params.ilist} """
        """-O {output} """
#"""-L {params.chrms} -XL {params.ilist} """

rule females_default_extract_females_call_VCF_XnonPARs_no_filter:
    input:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/default/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_raw_SNPs_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/default/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_no_filter")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/default/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_no_filter.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_default_extract_females_call_VCF_XPARs_no_filter:
    input:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/default/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_raw_SNPs_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/default/by_sample/{sample}_{chrmsXPARs}_PARs_no_filter")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/default/by_sample/{sample}_{chrmsXPARs}_PARs_no_filter.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_default_get_performance_metrics_autos_no_filter:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/default/by_sample/{sample}_{chrms}_autos_no_filter.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/no_filter/autos/{sample}_{chrms}_autos_no_filter_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/no_filter/autos/{sample}_{chrms}_autos_no_filter_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_default_get_performance_metrics_XnonPARs_no_filter:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/default/by_sample/{sample}_{chrms}_nonPARs_no_filter.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/no_filter/chrX_nonPARs/{sample}_{chrms}_nonPARs_no_filter_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/no_filter/chrX_nonPARs/{sample}_{chrms}_nonPARs_no_filter_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_default_get_performance_metrics_XPARs_no_filter:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/ASN/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/ASN/females/default/by_sample/{sample}_{chrms}_PARs_no_filter.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/no_filter/chrX_PARs/{sample}_{chrms}_PARs_no_filter_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/ASN/females/default/no_filter/chrX_PARs/{sample}_{chrms}_PARs_no_filter_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """
'''
