import os

# Environment: variant_calling_simulations_project

configfile: "xy_var_calling_pilot_04.config.json"


wildcard_constraints:
    chrmsXnonPARs= '|'.join([re.escape(x) for x in config["chrX"]]),
    chrmsXPARs= '|'.join([re.escape(x) for x in config["chrX"]]),
    chrmsY= '|'.join([re.escape(x) for x in config["chrY"]]),
    #chrmsA= '|'.join([re.escape(x) for x in config["autosomes_new"]]),
    chrmsA= '|'.join([re.escape(x) for x in config["autosomes"]]),
    chrmsM= '|'.join([re.escape(x) for x in config["chrM"]])



rule all:
    input:
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_DP_{filter}.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"], filter=config["DP_filters_new"]), # males - scc - hard filters
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_DP_{filter}.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"], filter=config["DP_filters_new"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_DP_{filter}.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"], filter=config["DP_filters_new"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP/chrY/by_sample/{sample}_{chrmsY}_DP_{filter}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"], filter=config["DP_filters_new"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP/autos/{sample}_{chrms}_autos_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"], filter=config["DP_filters_new"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP/chrX_nonPARs/{sample}_{chrms}_nonPARs_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["DP_filters_new"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP/chrX_PARs/{sample}_{chrms}_PARs_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["DP_filters_new"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP/chrY/{sample}_{chrms}_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"], filter=config["DP_filters_new"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP/chrM/by_sample/{sample}_{chrmsM}_DP_{filter}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"], filter=config["DP_filters_new"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP/chrM/{sample}_{chrms}_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"], filter=config["DP_filters_new"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP_INFO/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_DP_INFO_{filter}.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"], filter=config["DP_INFO_filters"]), # males - scc - hard filters
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP_INFO/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_DP_INFO_{filter}.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"], filter=config["DP_INFO_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP_INFO/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_DP_INFO_{filter}.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"], filter=config["DP_INFO_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP_INFO/chrY/by_sample/{sample}_{chrmsY}_DP_INFO_{filter}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"], filter=config["DP_INFO_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP_INFO/autos/{sample}_{chrms}_autos_DP_INFO_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"], filter=config["DP_INFO_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP_INFO/chrX_nonPARs/{sample}_{chrms}_nonPARs_DP_INFO_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["DP_INFO_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP_INFO/chrX_PARs/{sample}_{chrms}_PARs_DP_INFO_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["DP_INFO_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP_INFO/chrY/{sample}_{chrms}_DP_INFO_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"], filter=config["DP_INFO_filters"]),
        #expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP_INFO/chrM/by_sample/{sample}_{chrmsM}_DP_INFO_{filter}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"], filter=config["DP_filters_new"]),
        #expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP_INFO/chrM/{sample}_{chrms}_DP_INFO_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"], filter=config["DP_filters_new"]),

        #females -scc -  hard filters
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_DP_{filter}.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"], filter=config["DP_filters_new"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_DP_{filter}.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"], filter=config["DP_filters_new"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_DP_{filter}.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"], filter=config["DP_filters_new"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/DP/autos/{sample}_{chrms}_autos_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"], filter=config["DP_filters_new"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/DP/chrX_nonPARs/{sample}_{chrms}_nonPARs_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["DP_filters_new"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/DP/chrX_PARs/{sample}_{chrms}_PARs_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["DP_filters_new"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP/chrM/by_sample/{sample}_{chrmsM}_DP_{filter}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"], filter=config["DP_filters_new"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/DP/chrM/{sample}_{chrms}_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"], filter=config["DP_filters_new"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP_INFO/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_DP_INFO_{filter}.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"], filter=config["DP_INFO_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP_INFO/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_DP_INFO_{filter}.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"], filter=config["DP_INFO_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP_INFO/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_DP_INFO_{filter}.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"], filter=config["DP_INFO_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/DP_INFO/autos/{sample}_{chrms}_autos_DP_INFO_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"], filter=config["DP_INFO_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/DP_INFO/chrX_nonPARs/{sample}_{chrms}_nonPARs_DP_INFO_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["DP_INFO_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/DP_INFO/chrX_PARs/{sample}_{chrms}_PARs_DP_INFO_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["DP_INFO_filters"]),
        #expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP_INFO/chrM/by_sample/{sample}_{chrmsM}_DP_INFO_{filter}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"], filter=config["DP_INFO_filters"]),
        #expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/DP_INFO/chrM/{sample}_{chrms}_DP_INFO_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"], filter=config["DP_INFO_filters"]),

'''
        # FEMALES #
        # DEFAULT #
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/all/by_sample/{sample}_{chrmsA}_autos.recode.vcf"), sample=config["females"], chrmsA=config["autosomes_new"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/all/by_sample/{sample}_{chrmsXPARs}_PARs.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/all/by_sample/{sample}_{chrmsM}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"]),

        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/autos/{sample}_{chrms}_autos_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes_new"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/chrX_PARs/{sample}_{chrms}_PARs_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/chrM/{sample}_{chrms}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"]),

        # MALES #
        # DEFAULT #
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/all/by_sample/{sample}_{chrmsA}_autos.recode.vcf"), sample=config["males"], chrmsA=config["autosomes_new"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/all/by_sample/{sample}_{chrmsXPARs}_PARs.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/all/by_sample/{sample}_{chrmsY}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/all/by_sample/{sample}_{chrmsM}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"]),

        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/autos/{sample}_{chrms}_autos_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes_new"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/chrX_PARs/{sample}_{chrms}_PARs_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/chrY/{sample}_{chrms}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/chrM/{sample}_{chrms}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"]),


        # MALES #
        # SCC #
        expand(os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsA}_golden.vcf"), sample=config["males"], chrmsA=config["autosomes_new"]),
        expand(os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrmsXnonPARs}_nonPARs_golden.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsXPARs}_PARs_golden.vcf"), sample=config["males"], chrmsXPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrmsY}_golden.vcf"), sample=config["males"], chrmsY=config["chrY"]),
        expand(os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrmsM}_golden.vcf"), sample=config["males"], chrmsM=config["chrM"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrmsA}_autos.recode.vcf"), sample=config["males"], chrmsA=config["autosomes_new"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrmsXPARs}_PARs.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrmsY}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrmsM}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"]),

        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/autos/{sample}_{chrms}_autos_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes_new"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/chrX_PARs/{sample}_{chrms}_PARs_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/chrY/{sample}_{chrms}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/chrM/{sample}_{chrms}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"]),


        # FEMALES #
        # SCC #
        expand(os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsA}_golden.vcf"), sample=config["females"], chrmsA=config["autosomes_new"]),
        expand(os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsXnonPARs}_nonPARs_golden.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsXPARs}_PARs_golden.vcf"), sample=config["females"], chrmsXPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrmsM}_golden.vcf"), sample=config["females"], chrmsM=config["chrM"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/all/by_sample/{sample}_{chrmsA}_autos.recode.vcf"), sample=config["females"], chrmsA=config["autosomes_new"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/all/by_sample/{sample}_{chrmsXPARs}_PARs.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/all/by_sample/{sample}_{chrmsM}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"]),

        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/autos/{sample}_{chrms}_autos_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes_new"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/chrX_PARs/{sample}_{chrms}_PARs_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/chrM/{sample}_{chrms}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"]),



        # males - scc - hard filters
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QD/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QD_{filter}.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QD/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QD_{filter}.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QD/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QD_{filter}.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QD/chrY/by_sample/{sample}_{chrmsY}_QD_{filter}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/QD/autos/{sample}_{chrms}_autos_QD_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/QD/chrX_nonPARs/{sample}_{chrms}_nonPARs_QD_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/QD/chrX_PARs/{sample}_{chrms}_PARs_QD_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/QD/chrY/{sample}_{chrms}_QD_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QD/chrM/by_sample/{sample}_{chrmsM}_QD_{filter}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/QD/chrM/{sample}_{chrms}_QD_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"], filter=config["QD_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QUAL/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QUAL_{filter}.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QUAL/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QUAL_{filter}.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QUAL/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QUAL_{filter}.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QUAL/chrY/by_sample/{sample}_{chrmsY}_QUAL_{filter}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/QUAL/autos/{sample}_{chrms}_autos_QUAL_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/QUAL/chrX_nonPARs/{sample}_{chrms}_nonPARs_QUAL_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/QUAL/chrX_PARs/{sample}_{chrms}_PARs_QUAL_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/QUAL/chrY/{sample}_{chrms}_QUAL_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QUAL/chrM/by_sample/{sample}_{chrmsM}_QUAL_{filter}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/QUAL/chrM/{sample}_{chrms}_QUAL_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"], filter=config["QUAL_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/SOR/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_SOR_{filter}.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/SOR/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_SOR_{filter}.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/SOR/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_SOR_{filter}.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/SOR/chrY/by_sample/{sample}_{chrmsY}_SOR_{filter}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/SOR/autos/{sample}_{chrms}_autos_SOR_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/SOR/chrX_nonPARs/{sample}_{chrms}_nonPARs_SOR_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/SOR/chrX_PARs/{sample}_{chrms}_PARs_SOR_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/SOR/chrY/{sample}_{chrms}_SOR_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/SOR/chrM/by_sample/{sample}_{chrmsM}_SOR_{filter}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/SOR/chrM/{sample}_{chrms}_SOR_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"], filter=config["SOR_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/FS/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_FS_{filter}.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/FS/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_FS_{filter}.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/FS/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_FS_{filter}.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/FS/chrY/by_sample/{sample}_{chrmsY}_FS_{filter}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/FS/autos/{sample}_{chrms}_autos_FS_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/FS/chrX_nonPARs/{sample}_{chrms}_nonPARs_FS_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/FS/chrX_PARs/{sample}_{chrms}_PARs_FS_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/FS/chrY/{sample}_{chrms}_FS_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/FS/chrM/by_sample/{sample}_{chrmsM}_FS_{filter}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/FS/chrM/{sample}_{chrms}_FS_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"], filter=config["FS_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQ/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQ_{filter}.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQ/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQ_{filter}.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQ/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQ_{filter}.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQ/chrY/by_sample/{sample}_{chrmsY}_MQ_{filter}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/MQ/autos/{sample}_{chrms}_autos_MQ_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/MQ/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQ_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/MQ/chrX_PARs/{sample}_{chrms}_PARs_MQ_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/MQ/chrY/{sample}_{chrms}_MQ_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQ/chrM/by_sample/{sample}_{chrmsM}_MQ_{filter}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/MQ/chrM/{sample}_{chrms}_MQ_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"], filter=config["MQ_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQRankSum_{filter}.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQRankSum_{filter}.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQRankSum_{filter}.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQRankSum/chrY/by_sample/{sample}_{chrmsY}_MQRankSum_{filter}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/MQRankSum/autos/{sample}_{chrms}_autos_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/MQRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/MQRankSum/chrX_PARs/{sample}_{chrms}_PARs_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/MQRankSum/chrY/{sample}_{chrms}_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQRankSum/chrM/by_sample/{sample}_{chrmsM}_MQRankSum_{filter}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/MQRankSum/chrM/{sample}_{chrms}_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"], filter=config["MQRankSum_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/ReadPosRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_ReadPosRankSum_{filter}.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/ReadPosRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_ReadPosRankSum_{filter}.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/ReadPosRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_ReadPosRankSum_{filter}.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/ReadPosRankSum/chrY/by_sample/{sample}_{chrmsY}_ReadPosRankSum_{filter}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/ReadPosRankSum/autos/{sample}_{chrms}_autos_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/ReadPosRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/ReadPosRankSum/chrX_PARs/{sample}_{chrms}_PARs_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/ReadPosRankSum/chrY/{sample}_{chrms}_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/ReadPosRankSum/chrM/by_sample/{sample}_{chrmsM}_ReadPosRankSum_{filter}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/ReadPosRankSum/chrM/{sample}_{chrms}_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"], filter=config["ReadPosRankSum_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_DP_{filter}.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_DP_{filter}.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_DP_{filter}.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP/chrY/by_sample/{sample}_{chrmsY}_DP_{filter}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP/autos/{sample}_{chrms}_autos_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP/chrX_nonPARs/{sample}_{chrms}_nonPARs_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP/chrX_PARs/{sample}_{chrms}_PARs_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP/chrY/{sample}_{chrms}_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP/chrM/by_sample/{sample}_{chrmsM}_DP_{filter}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP/chrM/{sample}_{chrms}_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"], filter=config["DP_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/AN/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_AN_{filter}.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/AN/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_AN_{filter}.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/AN/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_AN_{filter}.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/AN/chrY/by_sample/{sample}_{chrmsY}_AN_{filter}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/AN/autos/{sample}_{chrms}_autos_AN_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/AN/chrX_nonPARs/{sample}_{chrms}_nonPARs_AN_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/AN/chrX_PARs/{sample}_{chrms}_PARs_AN_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/AN/chrY/{sample}_{chrms}_AN_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/AN/chrM/by_sample/{sample}_{chrmsM}_AN_{filter}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/AN/chrM/{sample}_{chrms}_AN_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"], filter=config["AN_filters"]),


        #females -scc -  hard filters
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QD/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QD_{filter}.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QD/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QD_{filter}.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QD/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QD_{filter}.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/QD/autos/{sample}_{chrms}_autos_QD_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/QD/chrX_nonPARs/{sample}_{chrms}_nonPARs_QD_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/QD/chrX_PARs/{sample}_{chrms}_PARs_QD_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QD/chrM/by_sample/{sample}_{chrmsM}_QD_{filter}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/QD/chrM/{sample}_{chrms}_QD_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"], filter=config["QD_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QUAL/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QUAL_{filter}.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QUAL/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QUAL_{filter}.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QUAL/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QUAL_{filter}.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/QUAL/autos/{sample}_{chrms}_autos_QUAL_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/QUAL/chrX_nonPARs/{sample}_{chrms}_nonPARs_QUAL_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/QUAL/chrX_PARs/{sample}_{chrms}_PARs_QUAL_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QUAL/chrM/by_sample/{sample}_{chrmsM}_QUAL_{filter}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/QUAL/chrM/{sample}_{chrms}_QUAL_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"], filter=config["QUAL_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/SOR/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_SOR_{filter}.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/SOR/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_SOR_{filter}.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/SOR/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_SOR_{filter}.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/SOR/autos/{sample}_{chrms}_autos_SOR_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/SOR/chrX_nonPARs/{sample}_{chrms}_nonPARs_SOR_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/SOR/chrX_PARs/{sample}_{chrms}_PARs_SOR_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/SOR/chrM/by_sample/{sample}_{chrmsM}_SOR_{filter}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/SOR/chrM/{sample}_{chrms}_SOR_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"], filter=config["SOR_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/FS/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_FS_{filter}.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/FS/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_FS_{filter}.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/FS/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_FS_{filter}.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/FS/autos/{sample}_{chrms}_autos_FS_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/FS/chrX_nonPARs/{sample}_{chrms}_nonPARs_FS_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/FS/chrX_PARs/{sample}_{chrms}_PARs_FS_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/FS/chrM/by_sample/{sample}_{chrmsM}_FS_{filter}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/FS/chrM/{sample}_{chrms}_FS_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"], filter=config["FS_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQ/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQ_{filter}.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQ/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQ_{filter}.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQ/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQ_{filter}.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/MQ/autos/{sample}_{chrms}_autos_MQ_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/MQ/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQ_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/MQ/chrX_PARs/{sample}_{chrms}_PARs_MQ_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQ/chrM/by_sample/{sample}_{chrmsM}_MQ_{filter}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/MQ/chrM/{sample}_{chrms}_MQ_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"], filter=config["MQ_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQRankSum_{filter}.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQRankSum_{filter}.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQRankSum_{filter}.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/MQRankSum/autos/{sample}_{chrms}_autos_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/MQRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/MQRankSum/chrX_PARs/{sample}_{chrms}_PARs_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQRankSum/chrM/by_sample/{sample}_{chrmsM}_MQRankSum_{filter}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/MQRankSum/chrM/{sample}_{chrms}_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"], filter=config["MQRankSum_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/ReadPosRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_ReadPosRankSum_{filter}.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/ReadPosRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_ReadPosRankSum_{filter}.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/ReadPosRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_ReadPosRankSum_{filter}.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/ReadPosRankSum/autos/{sample}_{chrms}_autos_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/ReadPosRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/ReadPosRankSum/chrX_PARs/{sample}_{chrms}_PARs_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/ReadPosRankSum/chrM/by_sample/{sample}_{chrmsM}_ReadPosRankSum_{filter}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/ReadPosRankSum/chrM/{sample}_{chrms}_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"], filter=config["ReadPosRankSum_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_DP_{filter}.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_DP_{filter}.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_DP_{filter}.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/DP/autos/{sample}_{chrms}_autos_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/DP/chrX_nonPARs/{sample}_{chrms}_nonPARs_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/DP/chrX_PARs/{sample}_{chrms}_PARs_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP/chrM/by_sample/{sample}_{chrmsM}_DP_{filter}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/DP/chrM/{sample}_{chrms}_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"], filter=config["DP_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/AN/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_AN_{filter}.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/AN/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_AN_{filter}.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/AN/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_AN_{filter}.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/AN/autos/{sample}_{chrms}_autos_AN_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/AN/chrX_nonPARs/{sample}_{chrms}_nonPARs_AN_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/AN/chrX_PARs/{sample}_{chrms}_PARs_AN_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/AN/chrM/by_sample/{sample}_{chrmsM}_AN_{filter}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/AN/chrM/{sample}_{chrms}_AN_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"], filter=config["AN_filters"]),

        # FEMALES #
        # DEFAULT #
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/all/by_sample/{sample}_{chrmsA}_autos.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/all/by_sample/{sample}_{chrmsXPARs}_PARs.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/all/by_sample/{sample}_{chrmsM}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"]),

        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/autos/{sample}_{chrms}_autos_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/chrX_PARs/{sample}_{chrms}_PARs_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/chrM/{sample}_{chrms}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QD/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QD_{filter}.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QD/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QD_{filter}.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QD/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QD_{filter}.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/QD/autos/{sample}_{chrms}_autos_QD_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/QD/chrX_nonPARs/{sample}_{chrms}_nonPARs_QD_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/QD/chrX_PARs/{sample}_{chrms}_PARs_QD_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QD/chrM/by_sample/{sample}_{chrmsM}_QD_{filter}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/QD/chrM/{sample}_{chrms}_QD_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"], filter=config["QD_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QUAL/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QUAL_{filter}.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QUAL/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QUAL_{filter}.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QUAL/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QUAL_{filter}.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/QUAL/autos/{sample}_{chrms}_autos_QUAL_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/QUAL/chrX_nonPARs/{sample}_{chrms}_nonPARs_QUAL_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/QUAL/chrX_PARs/{sample}_{chrms}_PARs_QUAL_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QUAL/chrM/by_sample/{sample}_{chrmsM}_QUAL_{filter}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/QUAL/chrM/{sample}_{chrms}_QUAL_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"], filter=config["QUAL_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/SOR/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_SOR_{filter}.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/SOR/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_SOR_{filter}.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/SOR/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_SOR_{filter}.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/SOR/autos/{sample}_{chrms}_autos_SOR_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/SOR/chrX_nonPARs/{sample}_{chrms}_nonPARs_SOR_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/SOR/chrX_PARs/{sample}_{chrms}_PARs_SOR_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/SOR/chrM/by_sample/{sample}_{chrmsM}_SOR_{filter}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/SOR/chrM/{sample}_{chrms}_SOR_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"], filter=config["SOR_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/FS/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_FS_{filter}.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/FS/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_FS_{filter}.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/FS/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_FS_{filter}.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/FS/autos/{sample}_{chrms}_autos_FS_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/FS/chrX_nonPARs/{sample}_{chrms}_nonPARs_FS_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/FS/chrX_PARs/{sample}_{chrms}_PARs_FS_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/FS/chrM/by_sample/{sample}_{chrmsM}_FS_{filter}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/FS/chrM/{sample}_{chrms}_FS_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"], filter=config["FS_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQ/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQ_{filter}.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQ/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQ_{filter}.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQ/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQ_{filter}.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/MQ/autos/{sample}_{chrms}_autos_MQ_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/MQ/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQ_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/MQ/chrX_PARs/{sample}_{chrms}_PARs_MQ_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQ/chrM/by_sample/{sample}_{chrmsM}_MQ_{filter}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/MQ/chrM/{sample}_{chrms}_MQ_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"], filter=config["MQ_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQRankSum_{filter}.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQRankSum_{filter}.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQRankSum_{filter}.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/MQRankSum/autos/{sample}_{chrms}_autos_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/MQRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/MQRankSum/chrX_PARs/{sample}_{chrms}_PARs_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQRankSum/chrM/by_sample/{sample}_{chrmsM}_MQRankSum_{filter}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/MQRankSum/chrM/{sample}_{chrms}_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"], filter=config["MQRankSum_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/ReadPosRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_ReadPosRankSum_{filter}.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/ReadPosRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_ReadPosRankSum_{filter}.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/ReadPosRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_ReadPosRankSum_{filter}.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/ReadPosRankSum/autos/{sample}_{chrms}_autos_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/ReadPosRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/ReadPosRankSum/chrX_PARs/{sample}_{chrms}_PARs_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/ReadPosRankSum/chrM/by_sample/{sample}_{chrmsM}_ReadPosRankSum_{filter}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/ReadPosRankSum/chrM/{sample}_{chrms}_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"], filter=config["ReadPosRankSum_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/DP/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_DP_{filter}.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/DP/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_DP_{filter}.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/DP/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_DP_{filter}.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/DP/autos/{sample}_{chrms}_autos_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/DP/chrX_nonPARs/{sample}_{chrms}_nonPARs_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/DP/chrX_PARs/{sample}_{chrms}_PARs_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/DP/chrM/by_sample/{sample}_{chrmsM}_DP_{filter}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/DP/chrM/{sample}_{chrms}_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"], filter=config["DP_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/AN/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_AN_{filter}.recode.vcf"), sample=config["females"], chrmsA=config["autosomes"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/AN/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_AN_{filter}.recode.vcf"), sample=config["females"], chrmsXnonPARs=config["chrX"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/AN/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_AN_{filter}.recode.vcf"), sample=config["females"], chrmsXPARs=config["chrX"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/AN/autos/{sample}_{chrms}_autos_AN_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["autosomes"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/AN/chrX_nonPARs/{sample}_{chrms}_nonPARs_AN_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/AN/chrX_PARs/{sample}_{chrms}_PARs_AN_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrX"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/AN/chrM/by_sample/{sample}_{chrmsM}_AN_{filter}.recode.vcf"), sample=config["females"], chrmsM=config["chrM"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/AN/chrM/{sample}_{chrms}_AN_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["females"], chrms=config["chrM"], filter=config["AN_filters"]),


        # MALES DEFAULT #
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/all/by_sample/{sample}_{chrmsA}_autos.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/all/by_sample/{sample}_{chrmsXPARs}_PARs.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/all/by_sample/{sample}_{chrmsY}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/all/by_sample/{sample}_{chrmsM}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"]),

        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/autos/{sample}_{chrms}_autos_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/chrX_PARs/{sample}_{chrms}_PARs_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/chrY/{sample}_{chrms}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/chrM/{sample}_{chrms}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QD/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QD_{filter}.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QD/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QD_{filter}.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QD/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QD_{filter}.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QD/chrY/by_sample/{sample}_{chrmsY}_QD_{filter}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/QD/autos/{sample}_{chrms}_autos_QD_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/QD/chrX_nonPARs/{sample}_{chrms}_nonPARs_QD_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/QD/chrX_PARs/{sample}_{chrms}_PARs_QD_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/QD/chrY/{sample}_{chrms}_QD_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QD/chrM/by_sample/{sample}_{chrmsM}_QD_{filter}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"], filter=config["QD_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/QD/chrM/{sample}_{chrms}_QD_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"], filter=config["QD_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QUAL/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QUAL_{filter}.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QUAL/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QUAL_{filter}.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QUAL/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QUAL_{filter}.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QUAL/chrY/by_sample/{sample}_{chrmsY}_QUAL_{filter}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/QUAL/autos/{sample}_{chrms}_autos_QUAL_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/QUAL/chrX_nonPARs/{sample}_{chrms}_nonPARs_QUAL_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/QUAL/chrX_PARs/{sample}_{chrms}_PARs_QUAL_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/QUAL/chrY/{sample}_{chrms}_QUAL_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QUAL/chrM/by_sample/{sample}_{chrmsM}_QUAL_{filter}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"], filter=config["QUAL_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/QUAL/chrM/{sample}_{chrms}_QUAL_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"], filter=config["QUAL_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/SOR/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_SOR_{filter}.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/SOR/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_SOR_{filter}.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/SOR/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_SOR_{filter}.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/SOR/chrY/by_sample/{sample}_{chrmsY}_SOR_{filter}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/SOR/autos/{sample}_{chrms}_autos_SOR_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/SOR/chrX_nonPARs/{sample}_{chrms}_nonPARs_SOR_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/SOR/chrX_PARs/{sample}_{chrms}_PARs_SOR_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/SOR/chrY/{sample}_{chrms}_SOR_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/SOR/chrM/by_sample/{sample}_{chrmsM}_SOR_{filter}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"], filter=config["SOR_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/SOR/chrM/{sample}_{chrms}_SOR_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"], filter=config["SOR_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/FS/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_FS_{filter}.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/FS/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_FS_{filter}.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/FS/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_FS_{filter}.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/FS/chrY/by_sample/{sample}_{chrmsY}_FS_{filter}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/FS/autos/{sample}_{chrms}_autos_FS_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/FS/chrX_nonPARs/{sample}_{chrms}_nonPARs_FS_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/FS/chrX_PARs/{sample}_{chrms}_PARs_FS_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/FS/chrY/{sample}_{chrms}_FS_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/FS/chrM/by_sample/{sample}_{chrmsM}_FS_{filter}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"], filter=config["FS_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/FS/chrM/{sample}_{chrms}_FS_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"], filter=config["FS_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQ/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQ_{filter}.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQ/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQ_{filter}.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQ/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQ_{filter}.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQ/chrY/by_sample/{sample}_{chrmsY}_MQ_{filter}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/MQ/autos/{sample}_{chrms}_autos_MQ_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/MQ/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQ_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/MQ/chrX_PARs/{sample}_{chrms}_PARs_MQ_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/MQ/chrY/{sample}_{chrms}_MQ_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQ/chrM/by_sample/{sample}_{chrmsM}_MQ_{filter}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"], filter=config["MQ_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/MQ/chrM/{sample}_{chrms}_MQ_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"], filter=config["MQ_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQRankSum_{filter}.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQRankSum_{filter}.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQRankSum_{filter}.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQRankSum/chrY/by_sample/{sample}_{chrmsY}_MQRankSum_{filter}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/MQRankSum/autos/{sample}_{chrms}_autos_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/MQRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/MQRankSum/chrX_PARs/{sample}_{chrms}_PARs_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/MQRankSum/chrY/{sample}_{chrms}_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQRankSum/chrM/by_sample/{sample}_{chrmsM}_MQRankSum_{filter}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"], filter=config["MQRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/MQRankSum/chrM/{sample}_{chrms}_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"], filter=config["MQRankSum_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/ReadPosRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_ReadPosRankSum_{filter}.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/ReadPosRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_ReadPosRankSum_{filter}.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/ReadPosRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_ReadPosRankSum_{filter}.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/ReadPosRankSum/chrY/by_sample/{sample}_{chrmsY}_ReadPosRankSum_{filter}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/ReadPosRankSum/autos/{sample}_{chrms}_autos_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/ReadPosRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/ReadPosRankSum/chrX_PARs/{sample}_{chrms}_PARs_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/ReadPosRankSum/chrY/{sample}_{chrms}_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/ReadPosRankSum/chrM/by_sample/{sample}_{chrmsM}_ReadPosRankSum_{filter}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"], filter=config["ReadPosRankSum_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/ReadPosRankSum/chrM/{sample}_{chrms}_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"], filter=config["ReadPosRankSum_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/DP/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_DP_{filter}.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/DP/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_DP_{filter}.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/DP/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_DP_{filter}.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/DP/chrY/by_sample/{sample}_{chrmsY}_DP_{filter}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/DP/autos/{sample}_{chrms}_autos_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/DP/chrX_nonPARs/{sample}_{chrms}_nonPARs_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/DP/chrX_PARs/{sample}_{chrms}_PARs_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/DP/chrY/{sample}_{chrms}_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/DP/chrM/by_sample/{sample}_{chrmsM}_DP_{filter}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"], filter=config["DP_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/DP/chrM/{sample}_{chrms}_DP_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"], filter=config["DP_filters"]),

        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/AN/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_AN_{filter}.recode.vcf"), sample=config["males"], chrmsA=config["autosomes"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/AN/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_AN_{filter}.recode.vcf"), sample=config["males"], chrmsXnonPARs=config["chrX"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/AN/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_AN_{filter}.recode.vcf"), sample=config["males"], chrmsXPARs=config["chrX"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/AN/chrY/by_sample/{sample}_{chrmsY}_AN_{filter}.recode.vcf"), sample=config["males"], chrmsY=config["chrY"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/AN/autos/{sample}_{chrms}_autos_AN_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["autosomes"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/AN/chrX_nonPARs/{sample}_{chrms}_nonPARs_AN_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/AN/chrX_PARs/{sample}_{chrms}_PARs_AN_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrX"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/AN/chrY/{sample}_{chrms}_AN_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrY"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/AN/chrM/by_sample/{sample}_{chrmsM}_AN_{filter}.recode.vcf"), sample=config["males"], chrmsM=config["chrM"], filter=config["AN_filters"]),
        expand(os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/AN/chrM/{sample}_{chrms}_AN_{filter}_golden_vs_called_performance_metrics.txt"), sample=config["males"], chrms=config["chrM"], filter=config["AN_filters"]),

'''


#-------------------------------------------------------------------------------#
# Step: Prep golden (simulated) VCFs
#-------------------------------------------------------------------------------#
rule unzip_golden_VCFs_autos:
    input:
        os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsA}_golden.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsA}_golden.vcf")
    shell:
        """
        gunzip {input}
        """

rule unzip_golden_VCFs_X:
    input:
        os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrmsXnonPARs}_nonPARs_golden.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrmsXnonPARs}_nonPARs_golden.vcf")
    shell:
        """
        gunzip {input}
        """

rule unzip_golden_VCFs_X_PARs:
    input:
        os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsXPARs}_PARs_golden.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsXPARs}_PARs_golden.vcf")
    shell:
        """
        gunzip {input}
        """

rule unzip_golden_VCFs_Y:
    input:
        os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrmsY}_golden.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrmsY}_golden.vcf")
    shell:
        """
        gunzip {input}
        """

rule unzip_golden_VCFs_M:
    input:
        os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrmsM}_golden.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrmsM}_golden.vcf")
    shell:
        """
        gunzip {input}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - ALL
#-------------------------------------------------------------------------------#
rule extract_males_call_VCF_autos:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/{chrmsA}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrmsA}_autos")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrmsA}_autos.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XnonPARs:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XPARs:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrmsXPARs}_PARs")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrmsXPARs}_PARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_Y:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/{chrmsY}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrmsY}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrmsY}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_M:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/{chrmsM}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrmsM}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrmsM}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_autos:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrms}_autos.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/autos/{sample}_{chrms}_autos_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/autos/{sample}_{chrms}_autos_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XnonPARs:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrms}_nonPARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XPARs:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrms}_PARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/chrX_PARs/{sample}_{chrms}_PARs_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/chrX_PARs/{sample}_{chrms}_PARs_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrms}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/chrY/{sample}_{chrms}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/chrY/{sample}_{chrms}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_M:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/all/by_sample/{sample}_{chrms}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/chrM/{sample}_{chrms}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/chrM/{sample}_{chrms}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - QD
#-------------------------------------------------------------------------------#
rule extract_males_call_VCF_autos_QD:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QD/{chrmsA}/{chrmsA}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QD/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QD_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QD/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QD_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XnonPARs_QD:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QD/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QD/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QD_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QD/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QD_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XPARs_QD:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QD/chrX_PARs/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QD/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QD_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QD/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QD_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_Y_QD:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QD/chrY/{chrmsY}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QD/chrY/by_sample/{sample}_{chrmsY}_QD_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QD/chrY/by_sample/{sample}_{chrmsY}_QD_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule get_performance_metrics_autos_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QD/{chrms}/by_sample/{sample}_{chrms}_autos_QD_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/QD/autos/{sample}_{chrms}_autos_QD_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/QD/autos/{sample}_{chrms}_autos_QD_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XnonPARs_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QD/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_QD_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/QD/chrX_nonPARs/{sample}_{chrms}_nonPARs_QD_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/QD/chrX_nonPARs/{sample}_{chrms}_nonPARs_QD_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XPARs_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QD/chrX_PARs/by_sample/{sample}_{chrms}_PARs_QD_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/QD/chrX_PARs/{sample}_{chrms}_PARs_QD_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/QD/chrX_PARs/{sample}_{chrms}_PARs_QD_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QD/chrY/by_sample/{sample}_{chrms}_QD_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/QD/chrY/{sample}_{chrms}_QD_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/QD/chrY/{sample}_{chrms}_QD_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule extract_males_call_VCF_M_QD:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QD/chrM/{chrmsM}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QD/chrM/by_sample/{sample}_{chrmsM}_QD_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QD/chrM/by_sample/{sample}_{chrmsM}_QD_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_M_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QD/chrM/by_sample/{sample}_{chrms}_QD_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/QD/chrM/{sample}_{chrms}_QD_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/QD/chrM/{sample}_{chrms}_QD_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - QUAL
#-------------------------------------------------------------------------------#
rule extract_males_call_VCF_autos_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QUAL/{chrmsA}/{chrmsA}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QUAL/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QUAL_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QUAL/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QUAL_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XnonPARs_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QUAL/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QUAL/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QUAL_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QUAL/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QUAL_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XPARs_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QUAL/chrX_PARs/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QUAL/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QUAL_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QUAL/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QUAL_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_Y_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QUAL/chrY/{chrmsY}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QUAL/chrY/by_sample/{sample}_{chrmsY}_QUAL_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QUAL/chrY/by_sample/{sample}_{chrmsY}_QUAL_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule get_performance_metrics_autos_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QUAL/{chrms}/by_sample/{sample}_{chrms}_autos_QUAL_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/QUAL/autos/{sample}_{chrms}_autos_QUAL_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/QUAL/autos/{sample}_{chrms}_autos_QUAL_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XnonPARs_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QUAL/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_QUAL_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/QUAL/chrX_nonPARs/{sample}_{chrms}_nonPARs_QUAL_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/QUAL/chrX_nonPARs/{sample}_{chrms}_nonPARs_QUAL_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XPARs_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QUAL/chrX_PARs/by_sample/{sample}_{chrms}_PARs_QUAL_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/QUAL/chrX_PARs/{sample}_{chrms}_PARs_QUAL_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/QUAL/chrX_PARs/{sample}_{chrms}_PARs_QUAL_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QUAL/chrY/by_sample/{sample}_{chrms}_QUAL_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/QUAL/chrY/{sample}_{chrms}_QUAL_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/QUAL/chrY/{sample}_{chrms}_QUAL_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule extract_males_call_VCF_M_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QUAL/chrM/{chrmsM}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QUAL/chrM/by_sample/{sample}_{chrmsM}_QUAL_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QUAL/chrM/by_sample/{sample}_{chrmsM}_QUAL_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_M_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/QUAL/chrM/by_sample/{sample}_{chrms}_QUAL_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/QUAL/chrM/{sample}_{chrms}_QUAL_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/QUAL/chrM/{sample}_{chrms}_QUAL_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - SOR
#-------------------------------------------------------------------------------#
rule extract_males_call_VCF_autos_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/SOR/{chrmsA}/{chrmsA}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/SOR/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_SOR_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/SOR/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_SOR_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XnonPARs_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/SOR/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/SOR/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_SOR_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/SOR/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_SOR_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XPARs_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/SOR/chrX_PARs/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/SOR/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_SOR_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/SOR/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_SOR_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_Y_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/SOR/chrY/{chrmsY}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/SOR/chrY/by_sample/{sample}_{chrmsY}_SOR_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/SOR/chrY/by_sample/{sample}_{chrmsY}_SOR_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule get_performance_metrics_autos_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/SOR/{chrms}/by_sample/{sample}_{chrms}_autos_SOR_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/SOR/autos/{sample}_{chrms}_autos_SOR_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/SOR/autos/{sample}_{chrms}_autos_SOR_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XnonPARs_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/SOR/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_SOR_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/SOR/chrX_nonPARs/{sample}_{chrms}_nonPARs_SOR_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/SOR/chrX_nonPARs/{sample}_{chrms}_nonPARs_SOR_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XPARs_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/SOR/chrX_PARs/by_sample/{sample}_{chrms}_PARs_SOR_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/SOR/chrX_PARs/{sample}_{chrms}_PARs_SOR_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/SOR/chrX_PARs/{sample}_{chrms}_PARs_SOR_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/SOR/chrY/by_sample/{sample}_{chrms}_SOR_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/SOR/chrY/{sample}_{chrms}_SOR_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/SOR/chrY/{sample}_{chrms}_SOR_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule extract_males_call_VCF_M_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/SOR/chrM/{chrmsM}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/SOR/chrM/by_sample/{sample}_{chrmsM}_SOR_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/SOR/chrM/by_sample/{sample}_{chrmsM}_SOR_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_M_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/SOR/chrM/by_sample/{sample}_{chrms}_SOR_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/SOR/chrM/{sample}_{chrms}_SOR_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/SOR/chrM/{sample}_{chrms}_SOR_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - FS
#-------------------------------------------------------------------------------#
rule extract_males_call_VCF_autos_FS:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/FS/{chrmsA}/{chrmsA}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/FS/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_FS_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/FS/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_FS_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XnonPARs_FS:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/FS/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/FS/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_FS_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/FS/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_FS_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XPARs_FS:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/FS/chrX_PARs/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/FS/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_FS_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/FS/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_FS_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_Y_FS:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/FS/chrY/{chrmsY}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/FS/chrY/by_sample/{sample}_{chrmsY}_FS_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/FS/chrY/by_sample/{sample}_{chrmsY}_FS_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule get_performance_metrics_autos_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/FS/{chrms}/by_sample/{sample}_{chrms}_autos_FS_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/FS/autos/{sample}_{chrms}_autos_FS_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/FS/autos/{sample}_{chrms}_autos_FS_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XnonPARs_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/FS/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_FS_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/FS/chrX_nonPARs/{sample}_{chrms}_nonPARs_FS_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/FS/chrX_nonPARs/{sample}_{chrms}_nonPARs_FS_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XPARs_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/FS/chrX_PARs/by_sample/{sample}_{chrms}_PARs_FS_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/FS/chrX_PARs/{sample}_{chrms}_PARs_FS_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/FS/chrX_PARs/{sample}_{chrms}_PARs_FS_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/FS/chrY/by_sample/{sample}_{chrms}_FS_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/FS/chrY/{sample}_{chrms}_FS_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/FS/chrY/{sample}_{chrms}_FS_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule extract_males_call_VCF_M_FS:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/FS/chrM/{chrmsM}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/FS/chrM/by_sample/{sample}_{chrmsM}_FS_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/FS/chrM/by_sample/{sample}_{chrmsM}_FS_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_M_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/FS/chrM/by_sample/{sample}_{chrms}_FS_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/FS/chrM/{sample}_{chrms}_FS_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/FS/chrM/{sample}_{chrms}_FS_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - MQ
#-------------------------------------------------------------------------------#
rule extract_males_call_VCF_autos_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQ/{chrmsA}/{chrmsA}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQ/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQ_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQ/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQ_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XnonPARs_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQ/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQ/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQ_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQ/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQ_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XPARs_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQ/chrX_PARs/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQ/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQ_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQ/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQ_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_Y_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQ/chrY/{chrmsY}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQ/chrY/by_sample/{sample}_{chrmsY}_MQ_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQ/chrY/by_sample/{sample}_{chrmsY}_MQ_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule get_performance_metrics_autos_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQ/{chrms}/by_sample/{sample}_{chrms}_autos_MQ_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/MQ/autos/{sample}_{chrms}_autos_MQ_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/MQ/autos/{sample}_{chrms}_autos_MQ_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XnonPARs_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQ/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_MQ_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/MQ/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQ_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/MQ/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQ_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XPARs_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQ/chrX_PARs/by_sample/{sample}_{chrms}_PARs_MQ_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/MQ/chrX_PARs/{sample}_{chrms}_PARs_MQ_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/MQ/chrX_PARs/{sample}_{chrms}_PARs_MQ_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQ/chrY/by_sample/{sample}_{chrms}_MQ_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/MQ/chrY/{sample}_{chrms}_MQ_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/MQ/chrY/{sample}_{chrms}_MQ_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule extract_males_call_VCF_M_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQ/chrM/{chrmsM}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQ/chrM/by_sample/{sample}_{chrmsM}_MQ_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQ/chrM/by_sample/{sample}_{chrmsM}_MQ_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_M_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQ/chrM/by_sample/{sample}_{chrms}_MQ_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/MQ/chrM/{sample}_{chrms}_MQ_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/MQ/chrM/{sample}_{chrms}_MQ_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - MQRankSum
#-------------------------------------------------------------------------------#
rule extract_males_call_VCF_autos_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQRankSum/{chrmsA}/{chrmsA}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XnonPARs_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQRankSum/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XPARs_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQRankSum/chrX_PARs/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_Y_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQRankSum/chrY/{chrmsY}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQRankSum/chrY/by_sample/{sample}_{chrmsY}_MQRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQRankSum/chrY/by_sample/{sample}_{chrmsY}_MQRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule get_performance_metrics_autos_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQRankSum/{chrms}/by_sample/{sample}_{chrms}_autos_MQRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/MQRankSum/autos/{sample}_{chrms}_autos_MQRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/MQRankSum/autos/{sample}_{chrms}_autos_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XnonPARs_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQRankSum/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_MQRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/MQRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/MQRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XPARs_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQRankSum/chrX_PARs/by_sample/{sample}_{chrms}_PARs_MQRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/MQRankSum/chrX_PARs/{sample}_{chrms}_PARs_MQRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/MQRankSum/chrX_PARs/{sample}_{chrms}_PARs_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQRankSum/chrY/by_sample/{sample}_{chrms}_MQRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/MQRankSum/chrY/{sample}_{chrms}_MQRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/MQRankSum/chrY/{sample}_{chrms}_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule extract_males_call_VCF_M_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQRankSum/chrM/{chrmsM}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQRankSum/chrM/by_sample/{sample}_{chrmsM}_MQRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQRankSum/chrM/by_sample/{sample}_{chrmsM}_MQRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_M_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/MQRankSum/chrM/by_sample/{sample}_{chrms}_MQRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/MQRankSum/chrM/{sample}_{chrms}_MQRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/MQRankSum/chrM/{sample}_{chrms}_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - ReadPosRankSum
#-------------------------------------------------------------------------------#
rule extract_males_call_VCF_autos_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/ReadPosRankSum/{chrmsA}/{chrmsA}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/ReadPosRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_ReadPosRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/ReadPosRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_ReadPosRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XnonPARs_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/ReadPosRankSum/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/ReadPosRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_ReadPosRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/ReadPosRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_ReadPosRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XPARs_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/ReadPosRankSum/chrX_PARs/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/ReadPosRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_ReadPosRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/ReadPosRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_ReadPosRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_Y_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/ReadPosRankSum/chrY/{chrmsY}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/ReadPosRankSum/chrY/by_sample/{sample}_{chrmsY}_ReadPosRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/ReadPosRankSum/chrY/by_sample/{sample}_{chrmsY}_ReadPosRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule get_performance_metrics_autos_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/ReadPosRankSum/{chrms}/by_sample/{sample}_{chrms}_autos_ReadPosRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/ReadPosRankSum/autos/{sample}_{chrms}_autos_ReadPosRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/ReadPosRankSum/autos/{sample}_{chrms}_autos_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XnonPARs_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/ReadPosRankSum/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_ReadPosRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/ReadPosRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_ReadPosRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/ReadPosRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XPARs_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/ReadPosRankSum/chrX_PARs/by_sample/{sample}_{chrms}_PARs_ReadPosRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/ReadPosRankSum/chrX_PARs/{sample}_{chrms}_PARs_ReadPosRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/ReadPosRankSum/chrX_PARs/{sample}_{chrms}_PARs_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/ReadPosRankSum/chrY/by_sample/{sample}_{chrms}_ReadPosRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/ReadPosRankSum/chrY/{sample}_{chrms}_ReadPosRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/ReadPosRankSum/chrY/{sample}_{chrms}_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule extract_males_call_VCF_M_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/ReadPosRankSum/chrM/{chrmsM}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/ReadPosRankSum/chrM/by_sample/{sample}_{chrmsM}_ReadPosRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/ReadPosRankSum/chrM/by_sample/{sample}_{chrmsM}_ReadPosRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_M_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/ReadPosRankSum/chrM/by_sample/{sample}_{chrms}_ReadPosRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/ReadPosRankSum/chrM/{sample}_{chrms}_ReadPosRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/ReadPosRankSum/chrM/{sample}_{chrms}_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - DP
#-------------------------------------------------------------------------------#
rule extract_males_call_VCF_autos_DP:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP/{chrmsA}/{chrmsA}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_DP_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_DP_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XnonPARs_DP:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_DP_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_DP_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XPARs_DP:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP/chrX_PARs/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_DP_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_DP_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_Y_DP:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP/chrY/{chrmsY}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP/chrY/by_sample/{sample}_{chrmsY}_DP_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP/chrY/by_sample/{sample}_{chrmsY}_DP_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule get_performance_metrics_autos_DP:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP/{chrms}/by_sample/{sample}_{chrms}_autos_DP_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP/autos/{sample}_{chrms}_autos_DP_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP/autos/{sample}_{chrms}_autos_DP_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XnonPARs_DP:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_DP_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP/chrX_nonPARs/{sample}_{chrms}_nonPARs_DP_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP/chrX_nonPARs/{sample}_{chrms}_nonPARs_DP_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XPARs_DP:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP/chrX_PARs/by_sample/{sample}_{chrms}_PARs_DP_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP/chrX_PARs/{sample}_{chrms}_PARs_DP_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP/chrX_PARs/{sample}_{chrms}_PARs_DP_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y_DP:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP/chrY/by_sample/{sample}_{chrms}_DP_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP/chrY/{sample}_{chrms}_DP_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP/chrY/{sample}_{chrms}_DP_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule extract_males_call_VCF_M_DP:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP/chrM/{chrmsM}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP/chrM/by_sample/{sample}_{chrmsM}_DP_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP/chrM/by_sample/{sample}_{chrmsM}_DP_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_M_DP:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP/chrM/by_sample/{sample}_{chrms}_DP_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP/chrM/{sample}_{chrms}_DP_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP/chrM/{sample}_{chrms}_DP_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """



rule extract_males_call_VCF_autos_DP_INFO:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP_INFO/{chrmsA}/{chrmsA}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_DP_INFO_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP_INFO/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_DP_INFO_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP_INFO/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_DP_INFO_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XnonPARs_DP_INFO:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP_INFO/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_DP_INFO_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP_INFO/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_DP_INFO_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP_INFO/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_DP_INFO_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XPARs_DP_INFO:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP_INFO/chrX_PARs/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_DP_INFO_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP_INFO/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_DP_INFO_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP_INFO/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_DP_INFO_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_Y_DP_INFO:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP_INFO/chrY/{chrmsY}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_DP_INFO_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP_INFO/chrY/by_sample/{sample}_{chrmsY}_DP_INFO_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP_INFO/chrY/by_sample/{sample}_{chrmsY}_DP_INFO_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule get_performance_metrics_autos_DP_INFO:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP_INFO/{chrms}/by_sample/{sample}_{chrms}_autos_DP_INFO_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP_INFO/autos/{sample}_{chrms}_autos_DP_INFO_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP_INFO/autos/{sample}_{chrms}_autos_DP_INFO_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XnonPARs_DP_INFO:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP_INFO/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_DP_INFO_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP_INFO/chrX_nonPARs/{sample}_{chrms}_nonPARs_DP_INFO_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP_INFO/chrX_nonPARs/{sample}_{chrms}_nonPARs_DP_INFO_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XPARs_DP_INFO:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP_INFO/chrX_PARs/by_sample/{sample}_{chrms}_PARs_DP_INFO_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP_INFO/chrX_PARs/{sample}_{chrms}_PARs_DP_INFO_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP_INFO/chrX_PARs/{sample}_{chrms}_PARs_DP_INFO_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y_DP_INFO:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP_INFO/chrY/by_sample/{sample}_{chrms}_DP_INFO_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP_INFO/chrY/{sample}_{chrms}_DP_INFO_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP_INFO/chrY/{sample}_{chrms}_DP_INFO_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule extract_males_call_VCF_M_DP_INFO:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP_INFO/chrM/{chrmsM}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_DP_INFO_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP_INFO/chrM/by_sample/{sample}_{chrmsM}_DP_INFO_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP_INFO/chrM/by_sample/{sample}_{chrmsM}_DP_INFO_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_M_DP_INFO:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/DP_INFO/chrM/by_sample/{sample}_{chrms}_DP_INFO_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP_INFO/chrM/{sample}_{chrms}_DP_INFO_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/DP_INFO/chrM/{sample}_{chrms}_DP_INFO_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - AN
#-------------------------------------------------------------------------------#
rule extract_males_call_VCF_autos_AN:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/AN/{chrmsA}/{chrmsA}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/AN/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_AN_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/AN/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_AN_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XnonPARs_AN:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/AN/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/AN/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_AN_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/AN/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_AN_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XPARs_AN:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/AN/chrX_PARs/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/AN/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_AN_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/AN/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_AN_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_Y_AN:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/AN/chrY/{chrmsY}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/AN/chrY/by_sample/{sample}_{chrmsY}_AN_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/AN/chrY/by_sample/{sample}_{chrmsY}_AN_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule get_performance_metrics_autos_AN:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/AN/{chrms}/by_sample/{sample}_{chrms}_autos_AN_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/AN/autos/{sample}_{chrms}_autos_AN_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/AN/autos/{sample}_{chrms}_autos_AN_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XnonPARs_AN:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/AN/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_AN_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/AN/chrX_nonPARs/{sample}_{chrms}_nonPARs_AN_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/AN/chrX_nonPARs/{sample}_{chrms}_nonPARs_AN_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XPARs_AN:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/AN/chrX_PARs/by_sample/{sample}_{chrms}_PARs_AN_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/AN/chrX_PARs/{sample}_{chrms}_PARs_AN_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/AN/chrX_PARs/{sample}_{chrms}_PARs_AN_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y_AN:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/AN/chrY/by_sample/{sample}_{chrms}_AN_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/AN/chrY/{sample}_{chrms}_AN_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/AN/chrY/{sample}_{chrms}_AN_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule extract_males_call_VCF_M_AN:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/AN/chrM/{chrmsM}_GRCh38_YPARsMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/AN/chrM/by_sample/{sample}_{chrmsM}_AN_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/AN/chrM/by_sample/{sample}_{chrmsM}_AN_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_M_AN:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/AN/chrM/by_sample/{sample}_{chrms}_AN_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/AN/chrM/{sample}_{chrms}_AN_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/AN/chrM/{sample}_{chrms}_AN_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """




'''
#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - No filters
#-------------------------------------------------------------------------------#
rule extract_males_call_VCF_autos_no_filter:
    input:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/males/{chrmsA}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/males/by_sample/{sample}_{chrmsA}_autos_no_filter")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/males/by_sample/{sample}_{chrmsA}_autos_no_filter.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_nonPARs_nofilter:
    input:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/males/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/males/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_nonPARs.vcf.gz")
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
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/males/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/males/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_PARs.vcf.gz")
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
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/males/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/males/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_no_filter")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/males/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_no_filter.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_XPARs_no_filter:
    input:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/males/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/males/by_sample/{sample}_{chrmsXPARs}_PARs_no_filter")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/males/by_sample/{sample}_{chrmsXPARs}_PARs_no_filter.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule extract_males_call_VCF_Y_no_filter:
    input:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/males/{chrmsY}_GRCh38_YPARsMasked_gatk_diploid_called_raw.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/males/by_sample/{sample}_{chrmsY}_no_filter")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/males/by_sample/{sample}_{chrmsY}_no_filter.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule get_performance_metrics_autos_no_filter:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/males/by_sample/{sample}_{chrms}_autos_no_filter.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/no_filter/autos/{sample}_{chrms}_autos_no_filter_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/no_filter/autos/{sample}_{chrms}_autos_no_filter_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XnonPARs_no_filter:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/males/by_sample/{sample}_{chrms}_nonPARs_no_filter.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/no_filter/chrX_nonPARs/{sample}_{chrms}_nonPARs_no_filter_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/no_filter/chrX_nonPARs/{sample}_{chrms}_nonPARs_no_filter_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_XPARs_no_filter:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/males/by_sample/{sample}_{chrms}_PARs_no_filter.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/no_filter/chrX_PARs/{sample}_{chrms}_PARs_no_filter_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/no_filter/chrX_PARs/{sample}_{chrms}_PARs_no_filter_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule get_performance_metrics_Y_no_filter:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/males/by_sample/{sample}_{chrms}_no_filter.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/no_filter/chrY/{sample}_{chrms}_no_filter_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/no_filter/chrY/{sample}_{chrms}_no_filter_golden_vs_called_performance_metrics.txt")
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
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/all/{chrmsA}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/all/by_sample/{sample}_{chrmsA}_autos")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/all/by_sample/{sample}_{chrmsA}_autos.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_XnonPARs:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/all/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_XPARs:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/all/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/all/by_sample/{sample}_{chrmsXPARs}_PARs")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/all/by_sample/{sample}_{chrmsXPARs}_PARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_Y:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/all/{chrmsY}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/all/by_sample/{sample}_{chrmsY}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/all/by_sample/{sample}_{chrmsY}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_M:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/all/{chrmsM}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/all/by_sample/{sample}_{chrmsM}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/all/by_sample/{sample}_{chrmsM}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_get_performance_metrics_autos:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/all/by_sample/{sample}_{chrms}_autos.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/autos/{sample}_{chrms}_autos_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/autos/{sample}_{chrms}_autos_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_XnonPARs:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/all/by_sample/{sample}_{chrms}_nonPARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_XPARs:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/all/by_sample/{sample}_{chrms}_PARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/chrX_PARs/{sample}_{chrms}_PARs_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/chrX_PARs/{sample}_{chrms}_PARs_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_Y:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/all/by_sample/{sample}_{chrms}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/chrY/{sample}_{chrms}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/chrY/{sample}_{chrms}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_M:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/all/by_sample/{sample}_{chrms}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/chrM/{sample}_{chrms}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/chrM/{sample}_{chrms}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

###########
# FILTERS #
###########

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - QD
#-------------------------------------------------------------------------------#
rule default_extract_males_call_VCF_autos_QD:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QD/{chrmsA}/{chrmsA}_GRCh38_default_gatk_diploid_called_SNPs_filtered_QD_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QD/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QD_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QD/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QD_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_XnonPARs_QD:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QD/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_QD_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QD/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QD_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QD/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QD_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_XPARs_QD:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QD/chrX_PARs/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_QD_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QD/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QD_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QD/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QD_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_Y_QD:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QD/chrY/{chrmsY}_GRCh38_default_gatk_diploid_called_SNPs_filtered_QD_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QD/chrY/by_sample/{sample}_{chrmsY}_QD_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QD/chrY/by_sample/{sample}_{chrmsY}_QD_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule default_get_performance_metrics_autos_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QD/{chrms}/by_sample/{sample}_{chrms}_autos_QD_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/QD/autos/{sample}_{chrms}_autos_QD_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/QD/autos/{sample}_{chrms}_autos_QD_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_XnonPARs_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QD/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_QD_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/QD/chrX_nonPARs/{sample}_{chrms}_nonPARs_QD_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/QD/chrX_nonPARs/{sample}_{chrms}_nonPARs_QD_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_XPARs_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QD/chrX_PARs/by_sample/{sample}_{chrms}_PARs_QD_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/QD/chrX_PARs/{sample}_{chrms}_PARs_QD_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/QD/chrX_PARs/{sample}_{chrms}_PARs_QD_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_Y_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QD/chrY/by_sample/{sample}_{chrms}_QD_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/QD/chrY/{sample}_{chrms}_QD_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/QD/chrY/{sample}_{chrms}_QD_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_extract_males_call_VCF_M_QD:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QD/chrM/{chrmsM}_GRCh38_default_gatk_diploid_called_SNPs_filtered_QD_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QD/chrM/by_sample/{sample}_{chrmsM}_QD_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QD/chrM/by_sample/{sample}_{chrmsM}_QD_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_get_performance_metrics_M_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QD/chrM/by_sample/{sample}_{chrms}_QD_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/QD/chrM/{sample}_{chrms}_QD_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/QD/chrM/{sample}_{chrms}_QD_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - QUAL
#-------------------------------------------------------------------------------#
rule default_extract_males_call_VCF_autos_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QUAL/{chrmsA}/{chrmsA}_GRCh38_default_gatk_diploid_called_SNPs_filtered_QUAL_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QUAL/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QUAL_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QUAL/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QUAL_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_XnonPARs_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QUAL/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_QUAL_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QUAL/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QUAL_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QUAL/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QUAL_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_XPARs_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QUAL/chrX_PARs/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_QUAL_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QUAL/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QUAL_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QUAL/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QUAL_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_Y_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QUAL/chrY/{chrmsY}_GRCh38_default_gatk_diploid_called_SNPs_filtered_QUAL_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QUAL/chrY/by_sample/{sample}_{chrmsY}_QUAL_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QUAL/chrY/by_sample/{sample}_{chrmsY}_QUAL_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule default_get_performance_metrics_autos_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QUAL/{chrms}/by_sample/{sample}_{chrms}_autos_QUAL_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/QUAL/autos/{sample}_{chrms}_autos_QUAL_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/QUAL/autos/{sample}_{chrms}_autos_QUAL_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_XnonPARs_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QUAL/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_QUAL_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/QUAL/chrX_nonPARs/{sample}_{chrms}_nonPARs_QUAL_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/QUAL/chrX_nonPARs/{sample}_{chrms}_nonPARs_QUAL_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_XPARs_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QUAL/chrX_PARs/by_sample/{sample}_{chrms}_PARs_QUAL_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/QUAL/chrX_PARs/{sample}_{chrms}_PARs_QUAL_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/QUAL/chrX_PARs/{sample}_{chrms}_PARs_QUAL_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_Y_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QUAL/chrY/by_sample/{sample}_{chrms}_QUAL_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/QUAL/chrY/{sample}_{chrms}_QUAL_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/QUAL/chrY/{sample}_{chrms}_QUAL_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_extract_males_call_VCF_M_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QUAL/chrM/{chrmsM}_GRCh38_default_gatk_diploid_called_SNPs_filtered_QUAL_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QUAL/chrM/by_sample/{sample}_{chrmsM}_QUAL_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QUAL/chrM/by_sample/{sample}_{chrmsM}_QUAL_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_get_performance_metrics_M_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/QUAL/chrM/by_sample/{sample}_{chrms}_QUAL_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/QUAL/chrM/{sample}_{chrms}_QUAL_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/QUAL/chrM/{sample}_{chrms}_QUAL_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - SOR
#-------------------------------------------------------------------------------#
rule default_extract_males_call_VCF_autos_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/SOR/{chrmsA}/{chrmsA}_GRCh38_default_gatk_diploid_called_SNPs_filtered_SOR_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/SOR/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_SOR_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/SOR/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_SOR_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_XnonPARs_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/SOR/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_SOR_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/SOR/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_SOR_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/SOR/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_SOR_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_XPARs_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/SOR/chrX_PARs/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_SOR_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/SOR/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_SOR_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/SOR/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_SOR_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_Y_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/SOR/chrY/{chrmsY}_GRCh38_default_gatk_diploid_called_SNPs_filtered_SOR_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/SOR/chrY/by_sample/{sample}_{chrmsY}_SOR_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/SOR/chrY/by_sample/{sample}_{chrmsY}_SOR_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule default_get_performance_metrics_autos_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/SOR/{chrms}/by_sample/{sample}_{chrms}_autos_SOR_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/SOR/autos/{sample}_{chrms}_autos_SOR_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/SOR/autos/{sample}_{chrms}_autos_SOR_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_XnonPARs_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/SOR/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_SOR_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/SOR/chrX_nonPARs/{sample}_{chrms}_nonPARs_SOR_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/SOR/chrX_nonPARs/{sample}_{chrms}_nonPARs_SOR_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_XPARs_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/SOR/chrX_PARs/by_sample/{sample}_{chrms}_PARs_SOR_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/SOR/chrX_PARs/{sample}_{chrms}_PARs_SOR_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/SOR/chrX_PARs/{sample}_{chrms}_PARs_SOR_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_Y_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/SOR/chrY/by_sample/{sample}_{chrms}_SOR_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/SOR/chrY/{sample}_{chrms}_SOR_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/SOR/chrY/{sample}_{chrms}_SOR_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_extract_males_call_VCF_M_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/SOR/chrM/{chrmsM}_GRCh38_default_gatk_diploid_called_SNPs_filtered_SOR_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/SOR/chrM/by_sample/{sample}_{chrmsM}_SOR_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/SOR/chrM/by_sample/{sample}_{chrmsM}_SOR_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_get_performance_metrics_M_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/SOR/chrM/by_sample/{sample}_{chrms}_SOR_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/SOR/chrM/{sample}_{chrms}_SOR_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/SOR/chrM/{sample}_{chrms}_SOR_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - FS
#-------------------------------------------------------------------------------#
rule default_extract_males_call_VCF_autos_FS:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/FS/{chrmsA}/{chrmsA}_GRCh38_default_gatk_diploid_called_SNPs_filtered_FS_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/FS/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_FS_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/FS/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_FS_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_XnonPARs_FS:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/FS/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_FS_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/FS/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_FS_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/FS/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_FS_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_XPARs_FS:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/FS/chrX_PARs/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_FS_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/FS/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_FS_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/FS/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_FS_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_Y_FS:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/FS/chrY/{chrmsY}_GRCh38_default_gatk_diploid_called_SNPs_filtered_FS_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/FS/chrY/by_sample/{sample}_{chrmsY}_FS_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/FS/chrY/by_sample/{sample}_{chrmsY}_FS_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule default_get_performance_metrics_autos_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/FS/{chrms}/by_sample/{sample}_{chrms}_autos_FS_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/FS/autos/{sample}_{chrms}_autos_FS_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/FS/autos/{sample}_{chrms}_autos_FS_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_XnonPARs_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/FS/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_FS_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/FS/chrX_nonPARs/{sample}_{chrms}_nonPARs_FS_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/FS/chrX_nonPARs/{sample}_{chrms}_nonPARs_FS_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_XPARs_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/FS/chrX_PARs/by_sample/{sample}_{chrms}_PARs_FS_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/FS/chrX_PARs/{sample}_{chrms}_PARs_FS_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/FS/chrX_PARs/{sample}_{chrms}_PARs_FS_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_Y_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/FS/chrY/by_sample/{sample}_{chrms}_FS_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/FS/chrY/{sample}_{chrms}_FS_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/FS/chrY/{sample}_{chrms}_FS_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_extract_males_call_VCF_M_FS:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/FS/chrM/{chrmsM}_GRCh38_default_gatk_diploid_called_SNPs_filtered_FS_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/FS/chrM/by_sample/{sample}_{chrmsM}_FS_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/FS/chrM/by_sample/{sample}_{chrmsM}_FS_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_get_performance_metrics_M_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/FS/chrM/by_sample/{sample}_{chrms}_FS_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/FS/chrM/{sample}_{chrms}_FS_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/FS/chrM/{sample}_{chrms}_FS_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - MQ
#-------------------------------------------------------------------------------#
rule default_extract_males_call_VCF_autos_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQ/{chrmsA}/{chrmsA}_GRCh38_default_gatk_diploid_called_SNPs_filtered_MQ_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQ/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQ_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQ/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQ_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_XnonPARs_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQ/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_MQ_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQ/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQ_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQ/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQ_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_XPARs_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQ/chrX_PARs/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_MQ_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQ/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQ_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQ/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQ_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_Y_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQ/chrY/{chrmsY}_GRCh38_default_gatk_diploid_called_SNPs_filtered_MQ_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQ/chrY/by_sample/{sample}_{chrmsY}_MQ_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQ/chrY/by_sample/{sample}_{chrmsY}_MQ_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule default_get_performance_metrics_autos_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQ/{chrms}/by_sample/{sample}_{chrms}_autos_MQ_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/MQ/autos/{sample}_{chrms}_autos_MQ_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/MQ/autos/{sample}_{chrms}_autos_MQ_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_XnonPARs_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQ/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_MQ_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/MQ/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQ_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/MQ/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQ_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_XPARs_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQ/chrX_PARs/by_sample/{sample}_{chrms}_PARs_MQ_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/MQ/chrX_PARs/{sample}_{chrms}_PARs_MQ_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/MQ/chrX_PARs/{sample}_{chrms}_PARs_MQ_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_Y_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQ/chrY/by_sample/{sample}_{chrms}_MQ_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/MQ/chrY/{sample}_{chrms}_MQ_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/MQ/chrY/{sample}_{chrms}_MQ_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_extract_males_call_VCF_M_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQ/chrM/{chrmsM}_GRCh38_default_gatk_diploid_called_SNPs_filtered_MQ_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQ/chrM/by_sample/{sample}_{chrmsM}_MQ_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQ/chrM/by_sample/{sample}_{chrmsM}_MQ_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_get_performance_metrics_M_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQ/chrM/by_sample/{sample}_{chrms}_MQ_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/MQ/chrM/{sample}_{chrms}_MQ_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/MQ/chrM/{sample}_{chrms}_MQ_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - MQRankSum
#-------------------------------------------------------------------------------#
rule default_extract_males_call_VCF_autos_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQRankSum/{chrmsA}/{chrmsA}_GRCh38_default_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_XnonPARs_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQRankSum/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_XPARs_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQRankSum/chrX_PARs/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_Y_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQRankSum/chrY/{chrmsY}_GRCh38_default_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQRankSum/chrY/by_sample/{sample}_{chrmsY}_MQRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQRankSum/chrY/by_sample/{sample}_{chrmsY}_MQRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule default_get_performance_metrics_autos_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQRankSum/{chrms}/by_sample/{sample}_{chrms}_autos_MQRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/MQRankSum/autos/{sample}_{chrms}_autos_MQRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/MQRankSum/autos/{sample}_{chrms}_autos_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_XnonPARs_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQRankSum/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_MQRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/MQRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/MQRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_XPARs_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQRankSum/chrX_PARs/by_sample/{sample}_{chrms}_PARs_MQRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/MQRankSum/chrX_PARs/{sample}_{chrms}_PARs_MQRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/MQRankSum/chrX_PARs/{sample}_{chrms}_PARs_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_Y_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQRankSum/chrY/by_sample/{sample}_{chrms}_MQRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/MQRankSum/chrY/{sample}_{chrms}_MQRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/MQRankSum/chrY/{sample}_{chrms}_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_extract_males_call_VCF_M_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQRankSum/chrM/{chrmsM}_GRCh38_default_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQRankSum/chrM/by_sample/{sample}_{chrmsM}_MQRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQRankSum/chrM/by_sample/{sample}_{chrmsM}_MQRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_get_performance_metrics_M_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/MQRankSum/chrM/by_sample/{sample}_{chrms}_MQRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/MQRankSum/chrM/{sample}_{chrms}_MQRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/MQRankSum/chrM/{sample}_{chrms}_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - ReadPosRankSum
#-------------------------------------------------------------------------------#
rule default_extract_males_call_VCF_autos_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/ReadPosRankSum/{chrmsA}/{chrmsA}_GRCh38_default_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/ReadPosRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_ReadPosRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/ReadPosRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_ReadPosRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_XnonPARs_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/ReadPosRankSum/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/ReadPosRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_ReadPosRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/ReadPosRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_ReadPosRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_XPARs_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/ReadPosRankSum/chrX_PARs/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/ReadPosRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_ReadPosRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/ReadPosRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_ReadPosRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_Y_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/ReadPosRankSum/chrY/{chrmsY}_GRCh38_default_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/ReadPosRankSum/chrY/by_sample/{sample}_{chrmsY}_ReadPosRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/ReadPosRankSum/chrY/by_sample/{sample}_{chrmsY}_ReadPosRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule default_get_performance_metrics_autos_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/ReadPosRankSum/{chrms}/by_sample/{sample}_{chrms}_autos_ReadPosRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/ReadPosRankSum/autos/{sample}_{chrms}_autos_ReadPosRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/ReadPosRankSum/autos/{sample}_{chrms}_autos_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_XnonPARs_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/ReadPosRankSum/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_ReadPosRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/ReadPosRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_ReadPosRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/ReadPosRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_XPARs_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/ReadPosRankSum/chrX_PARs/by_sample/{sample}_{chrms}_PARs_ReadPosRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/ReadPosRankSum/chrX_PARs/{sample}_{chrms}_PARs_ReadPosRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/ReadPosRankSum/chrX_PARs/{sample}_{chrms}_PARs_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_Y_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/ReadPosRankSum/chrY/by_sample/{sample}_{chrms}_ReadPosRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/ReadPosRankSum/chrY/{sample}_{chrms}_ReadPosRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/ReadPosRankSum/chrY/{sample}_{chrms}_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_extract_males_call_VCF_M_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/ReadPosRankSum/chrM/{chrmsM}_GRCh38_default_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/ReadPosRankSum/chrM/by_sample/{sample}_{chrmsM}_ReadPosRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/ReadPosRankSum/chrM/by_sample/{sample}_{chrmsM}_ReadPosRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_get_performance_metrics_M_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/ReadPosRankSum/chrM/by_sample/{sample}_{chrms}_ReadPosRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/ReadPosRankSum/chrM/{sample}_{chrms}_ReadPosRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/ReadPosRankSum/chrM/{sample}_{chrms}_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - DP
#-------------------------------------------------------------------------------#
rule default_extract_males_call_VCF_autos_DP:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/DP/{chrmsA}/{chrmsA}_GRCh38_default_gatk_diploid_called_SNPs_filtered_DP_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/DP/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_DP_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/DP/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_DP_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_XnonPARs_DP:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/DP/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_DP_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/DP/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_DP_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/DP/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_DP_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_XPARs_DP:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/DP/chrX_PARs/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_DP_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/DP/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_DP_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/DP/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_DP_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_Y_DP:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/DP/chrY/{chrmsY}_GRCh38_default_gatk_diploid_called_SNPs_filtered_DP_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/DP/chrY/by_sample/{sample}_{chrmsY}_DP_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/DP/chrY/by_sample/{sample}_{chrmsY}_DP_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule default_get_performance_metrics_autos_DP:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/DP/{chrms}/by_sample/{sample}_{chrms}_autos_DP_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/DP/autos/{sample}_{chrms}_autos_DP_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/DP/autos/{sample}_{chrms}_autos_DP_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_XnonPARs_DP:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/DP/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_DP_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/DP/chrX_nonPARs/{sample}_{chrms}_nonPARs_DP_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/DP/chrX_nonPARs/{sample}_{chrms}_nonPARs_DP_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_XPARs_DP:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/DP/chrX_PARs/by_sample/{sample}_{chrms}_PARs_DP_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/DP/chrX_PARs/{sample}_{chrms}_PARs_DP_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/DP/chrX_PARs/{sample}_{chrms}_PARs_DP_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_Y_DP:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/DP/chrY/by_sample/{sample}_{chrms}_DP_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/DP/chrY/{sample}_{chrms}_DP_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/DP/chrY/{sample}_{chrms}_DP_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_extract_males_call_VCF_M_DP:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/DP/chrM/{chrmsM}_GRCh38_default_gatk_diploid_called_SNPs_filtered_DP_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/DP/chrM/by_sample/{sample}_{chrmsM}_DP_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/DP/chrM/by_sample/{sample}_{chrmsM}_DP_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_get_performance_metrics_M_DP:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/DP/chrM/by_sample/{sample}_{chrms}_DP_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/DP/chrM/{sample}_{chrms}_DP_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/DP/chrM/{sample}_{chrms}_DP_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - AN
#-------------------------------------------------------------------------------#
rule default_extract_males_call_VCF_autos_AN:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/AN/{chrmsA}/{chrmsA}_GRCh38_default_gatk_diploid_called_SNPs_filtered_AN_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/AN/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_AN_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/AN/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_AN_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_XnonPARs_AN:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/AN/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_AN_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/AN/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_AN_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/AN/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_AN_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_XPARs_AN:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/AN/chrX_PARs/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_AN_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/AN/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_AN_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/AN/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_AN_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_extract_males_call_VCF_Y_AN:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/AN/chrY/{chrmsY}_GRCh38_default_gatk_diploid_called_SNPs_filtered_AN_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/AN/chrY/by_sample/{sample}_{chrmsY}_AN_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/AN/chrY/by_sample/{sample}_{chrmsY}_AN_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule default_get_performance_metrics_autos_AN:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/AN/{chrms}/by_sample/{sample}_{chrms}_autos_AN_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/AN/autos/{sample}_{chrms}_autos_AN_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/AN/autos/{sample}_{chrms}_autos_AN_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_XnonPARs_AN:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/AN/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_AN_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/AN/chrX_nonPARs/{sample}_{chrms}_nonPARs_AN_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/AN/chrX_nonPARs/{sample}_{chrms}_nonPARs_AN_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_XPARs_AN:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/AN/chrX_PARs/by_sample/{sample}_{chrms}_PARs_AN_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/AN/chrX_PARs/{sample}_{chrms}_PARs_AN_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/AN/chrX_PARs/{sample}_{chrms}_PARs_AN_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_get_performance_metrics_Y_AN:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/AN/chrY/by_sample/{sample}_{chrms}_AN_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/AN/chrY/{sample}_{chrms}_AN_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/AN/chrY/{sample}_{chrms}_AN_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_extract_males_call_VCF_M_AN:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/AN/chrM/{chrmsM}_GRCh38_default_gatk_diploid_called_SNPs_filtered_AN_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/AN/chrM/by_sample/{sample}_{chrmsM}_AN_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/AN/chrM/by_sample/{sample}_{chrmsM}_AN_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_get_performance_metrics_M_AN:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/males/default/AN/chrM/by_sample/{sample}_{chrms}_AN_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/AN/chrM/{sample}_{chrms}_AN_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/males/default/AN/chrM/{sample}_{chrms}_AN_{filter}_golden_vs_called_performance_metrics.txt")
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
        os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsA}_golden.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsA}_golden.vcf")
    shell:
        """
        gunzip {input}
        """

rule females_unzip_golden_VCFs_X:
    input:
        os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsXnonPARs}_nonPARs_golden.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsXnonPARs}_nonPARs_golden.vcf")
    shell:
        """
        gunzip {input}
        """

rule females_unzip_golden_VCFs_X_PARs:
    input:
        os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsXPARs}_PARs_golden.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrmsXPARs}_PARs_golden.vcf")
    shell:
        """
        gunzip {input}
        """

rule females_unzip_golden_VCFs_M:
    input:
        os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrmsM}_golden.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrmsM}_golden.vcf")
    shell:
        """
        gunzip {input}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - ALL
#-------------------------------------------------------------------------------#
rule females_extract_females_call_VCF_autos:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/all/{chrmsA}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/all/by_sample/{sample}_{chrmsA}_autos")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/all/by_sample/{sample}_{chrmsA}_autos.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_extract_females_call_VCF_XnonPARs:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/all/{chrmsXnonPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_extract_females_call_VCF_XPARs:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/all/{chrmsXPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/all/by_sample/{sample}_{chrmsXPARs}_PARs")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/all/by_sample/{sample}_{chrmsXPARs}_PARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_extract_females_call_VCF_M:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/all/{chrmsM}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/all/by_sample/{sample}_{chrmsM}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/all/by_sample/{sample}_{chrmsM}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule females_get_performance_metrics_autos:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/all/by_sample/{sample}_{chrms}_autos.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/autos/{sample}_{chrms}_autos_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/autos/{sample}_{chrms}_autos_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_get_performance_metrics_XnonPARs:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/all/by_sample/{sample}_{chrms}_nonPARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_get_performance_metrics_XPARs:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/all/by_sample/{sample}_{chrms}_PARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/chrX_PARs/{sample}_{chrms}_PARs_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/chrX_PARs/{sample}_{chrms}_PARs_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_get_performance_metrics_M:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/all/by_sample/{sample}_{chrms}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/chrM/{sample}_{chrms}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/chrM/{sample}_{chrms}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

# FILTERS #
#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - QD
#-------------------------------------------------------------------------------#
rule female_extract_call_VCF_autos_QD:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QD/{chrmsA}/{chrmsA}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QD/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QD_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QD/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QD_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XnonPARs_QD:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QD/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QD/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QD_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QD/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QD_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XPARs_QD:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QD/chrX_PARs/{chrmsXPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QD/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QD_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QD/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QD_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """



rule female_get_performance_metrics_autos_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QD/{chrms}/by_sample/{sample}_{chrms}_autos_QD_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/QD/autos/{sample}_{chrms}_autos_QD_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/QD/autos/{sample}_{chrms}_autos_QD_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XnonPARs_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QD/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_QD_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/QD/chrX_nonPARs/{sample}_{chrms}_nonPARs_QD_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/QD/chrX_nonPARs/{sample}_{chrms}_nonPARs_QD_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XPARs_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QD/chrX_PARs/by_sample/{sample}_{chrms}_PARs_QD_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/QD/chrX_PARs/{sample}_{chrms}_PARs_QD_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/QD/chrX_PARs/{sample}_{chrms}_PARs_QD_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule female_extract_call_VCF_M_QD:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QD/chrM/{chrmsM}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QD_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QD/chrM/by_sample/{sample}_{chrmsM}_QD_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QD/chrM/by_sample/{sample}_{chrmsM}_QD_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_get_performance_metrics_M_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QD/chrM/by_sample/{sample}_{chrms}_QD_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/QD/chrM/{sample}_{chrms}_QD_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/QD/chrM/{sample}_{chrms}_QD_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - QUAL
#-------------------------------------------------------------------------------#
rule female_extract_call_VCF_autos_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QUAL/{chrmsA}/{chrmsA}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QUAL/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QUAL_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QUAL/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QUAL_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XnonPARs_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QUAL/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QUAL/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QUAL_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QUAL/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QUAL_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XPARs_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QUAL/chrX_PARs/{chrmsXPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QUAL/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QUAL_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QUAL/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QUAL_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule female_get_performance_metrics_autos_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QUAL/{chrms}/by_sample/{sample}_{chrms}_autos_QUAL_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/QUAL/autos/{sample}_{chrms}_autos_QUAL_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/QUAL/autos/{sample}_{chrms}_autos_QUAL_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XnonPARs_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QUAL/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_QUAL_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/QUAL/chrX_nonPARs/{sample}_{chrms}_nonPARs_QUAL_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/QUAL/chrX_nonPARs/{sample}_{chrms}_nonPARs_QUAL_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XPARs_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QUAL/chrX_PARs/by_sample/{sample}_{chrms}_PARs_QUAL_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/QUAL/chrX_PARs/{sample}_{chrms}_PARs_QUAL_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/QUAL/chrX_PARs/{sample}_{chrms}_PARs_QUAL_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule female_extract_call_VCF_M_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QUAL/chrM/{chrmsM}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_QUAL_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QUAL/chrM/by_sample/{sample}_{chrmsM}_QUAL_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QUAL/chrM/by_sample/{sample}_{chrmsM}_QUAL_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_get_performance_metrics_M_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/QUAL/chrM/by_sample/{sample}_{chrms}_QUAL_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/QUAL/chrM/{sample}_{chrms}_QUAL_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/QUAL/chrM/{sample}_{chrms}_QUAL_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - SOR
#-------------------------------------------------------------------------------#
rule female_extract_call_VCF_autos_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/SOR/{chrmsA}/{chrmsA}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/SOR/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_SOR_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/SOR/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_SOR_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XnonPARs_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/SOR/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/SOR/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_SOR_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/SOR/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_SOR_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XPARs_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/SOR/chrX_PARs/{chrmsXPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/SOR/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_SOR_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/SOR/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_SOR_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """



rule female_get_performance_metrics_autos_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/SOR/{chrms}/by_sample/{sample}_{chrms}_autos_SOR_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/SOR/autos/{sample}_{chrms}_autos_SOR_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/SOR/autos/{sample}_{chrms}_autos_SOR_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XnonPARs_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/SOR/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_SOR_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/SOR/chrX_nonPARs/{sample}_{chrms}_nonPARs_SOR_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/SOR/chrX_nonPARs/{sample}_{chrms}_nonPARs_SOR_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XPARs_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/SOR/chrX_PARs/by_sample/{sample}_{chrms}_PARs_SOR_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/SOR/chrX_PARs/{sample}_{chrms}_PARs_SOR_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/SOR/chrX_PARs/{sample}_{chrms}_PARs_SOR_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule female_extract_call_VCF_M_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/SOR/chrM/{chrmsM}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_SOR_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/SOR/chrM/by_sample/{sample}_{chrmsM}_SOR_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/SOR/chrM/by_sample/{sample}_{chrmsM}_SOR_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_get_performance_metrics_M_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/SOR/chrM/by_sample/{sample}_{chrms}_SOR_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/SOR/chrM/{sample}_{chrms}_SOR_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/SOR/chrM/{sample}_{chrms}_SOR_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - FS
#-------------------------------------------------------------------------------#
rule female_extract_call_VCF_autos_FS:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/FS/{chrmsA}/{chrmsA}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/FS/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_FS_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/FS/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_FS_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XnonPARs_FS:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/FS/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/FS/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_FS_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/FS/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_FS_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XPARs_FS:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/FS/chrX_PARs/{chrmsXPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/FS/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_FS_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/FS/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_FS_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule female_get_performance_metrics_autos_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/FS/{chrms}/by_sample/{sample}_{chrms}_autos_FS_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/FS/autos/{sample}_{chrms}_autos_FS_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/FS/autos/{sample}_{chrms}_autos_FS_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XnonPARs_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/FS/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_FS_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/FS/chrX_nonPARs/{sample}_{chrms}_nonPARs_FS_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/FS/chrX_nonPARs/{sample}_{chrms}_nonPARs_FS_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XPARs_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/FS/chrX_PARs/by_sample/{sample}_{chrms}_PARs_FS_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/FS/chrX_PARs/{sample}_{chrms}_PARs_FS_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/FS/chrX_PARs/{sample}_{chrms}_PARs_FS_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule female_extract_call_VCF_M_FS:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/FS/chrM/{chrmsM}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_FS_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/FS/chrM/by_sample/{sample}_{chrmsM}_FS_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/FS/chrM/by_sample/{sample}_{chrmsM}_FS_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_get_performance_metrics_M_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/FS/chrM/by_sample/{sample}_{chrms}_FS_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/FS/chrM/{sample}_{chrms}_FS_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/FS/chrM/{sample}_{chrms}_FS_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - MQ
#-------------------------------------------------------------------------------#
rule female_extract_call_VCF_autos_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQ/{chrmsA}/{chrmsA}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQ/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQ_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQ/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQ_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XnonPARs_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQ/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQ/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQ_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQ/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQ_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XPARs_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQ/chrX_PARs/{chrmsXPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQ/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQ_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQ/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQ_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """



rule female_get_performance_metrics_autos_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQ/{chrms}/by_sample/{sample}_{chrms}_autos_MQ_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/MQ/autos/{sample}_{chrms}_autos_MQ_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/MQ/autos/{sample}_{chrms}_autos_MQ_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XnonPARs_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQ/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_MQ_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/MQ/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQ_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/MQ/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQ_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XPARs_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQ/chrX_PARs/by_sample/{sample}_{chrms}_PARs_MQ_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/MQ/chrX_PARs/{sample}_{chrms}_PARs_MQ_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/MQ/chrX_PARs/{sample}_{chrms}_PARs_MQ_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule female_extract_call_VCF_M_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQ/chrM/{chrmsM}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQ_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQ/chrM/by_sample/{sample}_{chrmsM}_MQ_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQ/chrM/by_sample/{sample}_{chrmsM}_MQ_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_get_performance_metrics_M_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQ/chrM/by_sample/{sample}_{chrms}_MQ_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/MQ/chrM/{sample}_{chrms}_MQ_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/MQ/chrM/{sample}_{chrms}_MQ_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - MQRankSum
#-------------------------------------------------------------------------------#
rule female_extract_call_VCF_autos_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQRankSum/{chrmsA}/{chrmsA}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XnonPARs_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQRankSum/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XPARs_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQRankSum/chrX_PARs/{chrmsXPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule female_get_performance_metrics_autos_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQRankSum/{chrms}/by_sample/{sample}_{chrms}_autos_MQRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/MQRankSum/autos/{sample}_{chrms}_autos_MQRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/MQRankSum/autos/{sample}_{chrms}_autos_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XnonPARs_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQRankSum/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_MQRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/MQRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/MQRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XPARs_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQRankSum/chrX_PARs/by_sample/{sample}_{chrms}_PARs_MQRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/MQRankSum/chrX_PARs/{sample}_{chrms}_PARs_MQRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/MQRankSum/chrX_PARs/{sample}_{chrms}_PARs_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule female_extract_call_VCF_M_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQRankSum/chrM/{chrmsM}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQRankSum/chrM/by_sample/{sample}_{chrmsM}_MQRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQRankSum/chrM/by_sample/{sample}_{chrmsM}_MQRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_get_performance_metrics_M_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/MQRankSum/chrM/by_sample/{sample}_{chrms}_MQRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/MQRankSum/chrM/{sample}_{chrms}_MQRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/MQRankSum/chrM/{sample}_{chrms}_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - ReadPosRankSum
#-------------------------------------------------------------------------------#
rule female_extract_call_VCF_autos_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/ReadPosRankSum/{chrmsA}/{chrmsA}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/ReadPosRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_ReadPosRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/ReadPosRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_ReadPosRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XnonPARs_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/ReadPosRankSum/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/ReadPosRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_ReadPosRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/ReadPosRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_ReadPosRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XPARs_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/ReadPosRankSum/chrX_PARs/{chrmsXPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/ReadPosRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_ReadPosRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/ReadPosRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_ReadPosRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule female_get_performance_metrics_autos_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/ReadPosRankSum/{chrms}/by_sample/{sample}_{chrms}_autos_ReadPosRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/ReadPosRankSum/autos/{sample}_{chrms}_autos_ReadPosRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/ReadPosRankSum/autos/{sample}_{chrms}_autos_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XnonPARs_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/ReadPosRankSum/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_ReadPosRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/ReadPosRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_ReadPosRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/ReadPosRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XPARs_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/ReadPosRankSum/chrX_PARs/by_sample/{sample}_{chrms}_PARs_ReadPosRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/ReadPosRankSum/chrX_PARs/{sample}_{chrms}_PARs_ReadPosRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/ReadPosRankSum/chrX_PARs/{sample}_{chrms}_PARs_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule female_extract_call_VCF_M_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/ReadPosRankSum/chrM/{chrmsM}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/ReadPosRankSum/chrM/by_sample/{sample}_{chrmsM}_ReadPosRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/ReadPosRankSum/chrM/by_sample/{sample}_{chrmsM}_ReadPosRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_get_performance_metrics_M_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/ReadPosRankSum/chrM/by_sample/{sample}_{chrms}_ReadPosRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/ReadPosRankSum/chrM/{sample}_{chrms}_ReadPosRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/ReadPosRankSum/chrM/{sample}_{chrms}_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - DP
#-------------------------------------------------------------------------------#
rule female_extract_call_VCF_autos_DP:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP/{chrmsA}/{chrmsA}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_DP_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_DP_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XnonPARs_DP:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_DP_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_DP_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XPARs_DP:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP/chrX_PARs/{chrmsXPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_DP_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_DP_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """



rule female_get_performance_metrics_autos_DP:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP/{chrms}/by_sample/{sample}_{chrms}_autos_DP_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/DP/autos/{sample}_{chrms}_autos_DP_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/DP/autos/{sample}_{chrms}_autos_DP_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XnonPARs_DP:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_DP_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/DP/chrX_nonPARs/{sample}_{chrms}_nonPARs_DP_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/DP/chrX_nonPARs/{sample}_{chrms}_nonPARs_DP_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XPARs_DP:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP/chrX_PARs/by_sample/{sample}_{chrms}_PARs_DP_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/DP/chrX_PARs/{sample}_{chrms}_PARs_DP_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/DP/chrX_PARs/{sample}_{chrms}_PARs_DP_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule female_extract_call_VCF_M_DP:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP/chrM/{chrmsM}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_DP_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP/chrM/by_sample/{sample}_{chrmsM}_DP_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP/chrM/by_sample/{sample}_{chrmsM}_DP_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_get_performance_metrics_M_DP:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP/chrM/by_sample/{sample}_{chrms}_DP_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/DP/chrM/{sample}_{chrms}_DP_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/DP/chrM/{sample}_{chrms}_DP_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - DP
#-------------------------------------------------------------------------------#
rule female_extract_call_VCF_autos_DP_INFO:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP_INFO/{chrmsA}/{chrmsA}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_DP_INFO_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP_INFO/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_DP_INFO_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP_INFO/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_DP_INFO_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XnonPARs_DP_INFO:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP_INFO/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_DP_INFO_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP_INFO/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_DP_INFO_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP_INFO/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_DP_INFO_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XPARs_DP_INFO:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP_INFO/chrX_PARs/{chrmsXPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_DP_INFO_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP_INFO/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_DP_INFO_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP_INFO/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_DP_INFO_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """



rule female_get_performance_metrics_autos_DP_INFO:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP_INFO/{chrms}/by_sample/{sample}_{chrms}_autos_DP_INFO_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/DP_INFO/autos/{sample}_{chrms}_autos_DP_INFO_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/DP_INFO/autos/{sample}_{chrms}_autos_DP_INFO_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XnonPARs_DP_INFO:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP_INFO/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_DP_INFO_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/DP_INFO/chrX_nonPARs/{sample}_{chrms}_nonPARs_DP_INFO_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/DP_INFO/chrX_nonPARs/{sample}_{chrms}_nonPARs_DP_INFO_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XPARs_DP_INFO:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP_INFO/chrX_PARs/by_sample/{sample}_{chrms}_PARs_DP_INFO_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/DP_INFO/chrX_PARs/{sample}_{chrms}_PARs_DP_INFO_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/DP_INFO/chrX_PARs/{sample}_{chrms}_PARs_DP_INFO_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule female_extract_call_VCF_M_DP_INFO:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP_INFO/chrM/{chrmsM}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_DP_INFO_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP_INFO/chrM/by_sample/{sample}_{chrmsM}_DP_INFO_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP_INFO/chrM/by_sample/{sample}_{chrmsM}_DP_INFO_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_get_performance_metrics_M_DP_INFO:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/DP_INFO/chrM/by_sample/{sample}_{chrms}_DP_INFO_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/DP_INFO/chrM/{sample}_{chrms}_DP_INFO_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/DP_INFO/chrM/{sample}_{chrms}_DP_INFO_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """



#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - AN
#-------------------------------------------------------------------------------#
rule female_extract_call_VCF_autos_AN:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/AN/{chrmsA}/{chrmsA}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/AN/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_AN_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/AN/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_AN_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XnonPARs_AN:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/AN/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/AN/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_AN_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/AN/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_AN_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_extract_call_VCF_XPARs_AN:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/AN/chrX_PARs/{chrmsXPARs}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/AN/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_AN_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/AN/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_AN_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule female_get_performance_metrics_autos_AN:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/AN/{chrms}/by_sample/{sample}_{chrms}_autos_AN_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/AN/autos/{sample}_{chrms}_autos_AN_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/AN/autos/{sample}_{chrms}_autos_AN_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XnonPARs_AN:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/AN/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_AN_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/AN/chrX_nonPARs/{sample}_{chrms}_nonPARs_AN_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/AN/chrX_nonPARs/{sample}_{chrms}_nonPARs_AN_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule female_get_performance_metrics_XPARs_AN:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/AN/chrX_PARs/by_sample/{sample}_{chrms}_PARs_AN_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/AN/chrX_PARs/{sample}_{chrms}_PARs_AN_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/AN/chrX_PARs/{sample}_{chrms}_PARs_AN_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule female_extract_call_VCF_M_AN:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/AN/chrM/{chrmsM}_GRCh38_YHardMasked_gatk_diploid_called_SNPs_filtered_AN_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/AN/chrM/by_sample/{sample}_{chrmsM}_AN_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/AN/chrM/by_sample/{sample}_{chrmsM}_AN_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule female_get_performance_metrics_M_AN:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/AN/chrM/by_sample/{sample}_{chrms}_AN_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/AN/chrM/{sample}_{chrms}_AN_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/AN/chrM/{sample}_{chrms}_AN_{filter}_golden_vs_called_performance_metrics.txt")
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
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/{chrmsA}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/by_sample/{sample}_{chrmsA}_autos_no_filter")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/by_sample/{sample}_{chrmsA}_autos_no_filter.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_extract_nonPARs_nofilter:
    input:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_nonPARs.vcf.gz")
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
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_PARs.vcf.gz")
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
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/{chrmsXnonPARs}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_no_filter")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_no_filter.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_extract_females_call_VCF_XPARs_no_filter:
    input:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/{chrmsXPARs}_GRCh38_YPARsMasked_gatk_diploid_called_raw_SNPs_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/by_sample/{sample}_{chrmsXPARs}_PARs_no_filter")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/by_sample/{sample}_{chrmsXPARs}_PARs_no_filter.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_get_performance_metrics_autos_no_filter:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/by_sample/{sample}_{chrms}_autos_no_filter.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/no_filter/autos/{sample}_{chrms}_autos_no_filter_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/no_filter/autos/{sample}_{chrms}_autos_no_filter_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_get_performance_metrics_XnonPARs_no_filter:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/by_sample/{sample}_{chrms}_nonPARs_no_filter.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/no_filter/chrX_nonPARs/{sample}_{chrms}_nonPARs_no_filter_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/no_filter/chrX_nonPARs/{sample}_{chrms}_nonPARs_no_filter_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_get_performance_metrics_XPARs_no_filter:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/by_sample/{sample}_{chrms}_PARs_no_filter.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/no_filter/chrX_PARs/{sample}_{chrms}_PARs_no_filter_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/no_filter/chrX_PARs/{sample}_{chrms}_PARs_no_filter_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """




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
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/all/{chrmsA}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/all/by_sample/{sample}_{chrmsA}_autos")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/all/by_sample/{sample}_{chrmsA}_autos.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_default_extract_females_call_VCF_XnonPARs:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/all/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/all/by_sample/{sample}_{chrmsXnonPARs}_nonPARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_default_extract_females_call_VCF_XPARs:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/all/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/all/by_sample/{sample}_{chrmsXPARs}_PARs")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/all/by_sample/{sample}_{chrmsXPARs}_PARs.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_default_extract_females_call_VCF_M:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/all/{chrmsM}_GRCh38_default_gatk_diploid_called_SNPs_gatkHardFilter_filtered.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/all/by_sample/{sample}_{chrmsM}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/all/by_sample/{sample}_{chrmsM}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_default_get_performance_metrics_autos:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/all/by_sample/{sample}_{chrms}_autos.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/autos/{sample}_{chrms}_autos_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/autos/{sample}_{chrms}_autos_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_default_get_performance_metrics_XnonPARs:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/all/by_sample/{sample}_{chrms}_nonPARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/chrX_nonPARs/{sample}_{chrms}_nonPARs_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_default_get_performance_metrics_XPARs:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/all/by_sample/{sample}_{chrms}_PARs.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/chrX_PARs/{sample}_{chrms}_PARs_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/chrX_PARs/{sample}_{chrms}_PARs_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_default_get_performance_metrics_M:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/all/by_sample/{sample}_{chrms}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/chrM/{sample}_{chrms}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/chrM/{sample}_{chrms}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

# For females, I dont think we are testing different filters. We just want to
# compare SCC to default
############
# FILTERS #
############
#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - QD
#-------------------------------------------------------------------------------#
rule default_female_extract_call_VCF_autos_QD:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QD/{chrmsA}/{chrmsA}_GRCh38_default_gatk_diploid_called_SNPs_filtered_QD_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QD/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QD_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QD/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QD_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_female_extract_call_VCF_XnonPARs_QD:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QD/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_QD_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QD/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QD_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QD/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QD_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_female_extract_call_VCF_XPARs_QD:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QD/chrX_PARs/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_QD_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QD/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QD_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QD/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QD_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """



rule default_female_get_performance_metrics_autos_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QD/{chrms}/by_sample/{sample}_{chrms}_autos_QD_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/QD/autos/{sample}_{chrms}_autos_QD_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/QD/autos/{sample}_{chrms}_autos_QD_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_female_get_performance_metrics_XnonPARs_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QD/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_QD_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/QD/chrX_nonPARs/{sample}_{chrms}_nonPARs_QD_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/QD/chrX_nonPARs/{sample}_{chrms}_nonPARs_QD_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_female_get_performance_metrics_XPARs_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QD/chrX_PARs/by_sample/{sample}_{chrms}_PARs_QD_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/QD/chrX_PARs/{sample}_{chrms}_PARs_QD_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/QD/chrX_PARs/{sample}_{chrms}_PARs_QD_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule default_female_extract_call_VCF_M_QD:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QD/chrM/{chrmsM}_GRCh38_default_gatk_diploid_called_SNPs_filtered_QD_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QD/chrM/by_sample/{sample}_{chrmsM}_QD_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QD/chrM/by_sample/{sample}_{chrmsM}_QD_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_female_get_performance_metrics_M_QD:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QD/chrM/by_sample/{sample}_{chrms}_QD_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/QD/chrM/{sample}_{chrms}_QD_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/QD/chrM/{sample}_{chrms}_QD_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - QUAL
#-------------------------------------------------------------------------------#
rule default_female_extract_call_VCF_autos_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QUAL/{chrmsA}/{chrmsA}_GRCh38_default_gatk_diploid_called_SNPs_filtered_QUAL_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QUAL/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QUAL_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QUAL/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_QUAL_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_female_extract_call_VCF_XnonPARs_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QUAL/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_QUAL_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QUAL/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QUAL_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QUAL/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_QUAL_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_female_extract_call_VCF_XPARs_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QUAL/chrX_PARs/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_QUAL_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QUAL/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QUAL_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QUAL/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_QUAL_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule default_female_get_performance_metrics_autos_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QUAL/{chrms}/by_sample/{sample}_{chrms}_autos_QUAL_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/QUAL/autos/{sample}_{chrms}_autos_QUAL_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/QUAL/autos/{sample}_{chrms}_autos_QUAL_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_female_get_performance_metrics_XnonPARs_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QUAL/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_QUAL_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/QUAL/chrX_nonPARs/{sample}_{chrms}_nonPARs_QUAL_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/QUAL/chrX_nonPARs/{sample}_{chrms}_nonPARs_QUAL_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_female_get_performance_metrics_XPARs_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QUAL/chrX_PARs/by_sample/{sample}_{chrms}_PARs_QUAL_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/QUAL/chrX_PARs/{sample}_{chrms}_PARs_QUAL_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/QUAL/chrX_PARs/{sample}_{chrms}_PARs_QUAL_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule default_female_extract_call_VCF_M_QUAL:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QUAL/chrM/{chrmsM}_GRCh38_default_gatk_diploid_called_SNPs_filtered_QUAL_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QUAL/chrM/by_sample/{sample}_{chrmsM}_QUAL_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QUAL/chrM/by_sample/{sample}_{chrmsM}_QUAL_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_female_get_performance_metrics_M_QUAL:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/QUAL/chrM/by_sample/{sample}_{chrms}_QUAL_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/QUAL/chrM/{sample}_{chrms}_QUAL_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/QUAL/chrM/{sample}_{chrms}_QUAL_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - SOR
#-------------------------------------------------------------------------------#
rule default_female_extract_call_VCF_autos_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/SOR/{chrmsA}/{chrmsA}_GRCh38_default_gatk_diploid_called_SNPs_filtered_SOR_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/SOR/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_SOR_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/SOR/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_SOR_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_female_extract_call_VCF_XnonPARs_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/SOR/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_SOR_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/SOR/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_SOR_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/SOR/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_SOR_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_female_extract_call_VCF_XPARs_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/SOR/chrX_PARs/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_SOR_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/SOR/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_SOR_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/SOR/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_SOR_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """



rule default_female_get_performance_metrics_autos_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/SOR/{chrms}/by_sample/{sample}_{chrms}_autos_SOR_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/SOR/autos/{sample}_{chrms}_autos_SOR_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/SOR/autos/{sample}_{chrms}_autos_SOR_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_female_get_performance_metrics_XnonPARs_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/SOR/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_SOR_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/SOR/chrX_nonPARs/{sample}_{chrms}_nonPARs_SOR_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/SOR/chrX_nonPARs/{sample}_{chrms}_nonPARs_SOR_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_female_get_performance_metrics_XPARs_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/SOR/chrX_PARs/by_sample/{sample}_{chrms}_PARs_SOR_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/SOR/chrX_PARs/{sample}_{chrms}_PARs_SOR_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/SOR/chrX_PARs/{sample}_{chrms}_PARs_SOR_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule default_female_extract_call_VCF_M_SOR:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/SOR/chrM/{chrmsM}_GRCh38_default_gatk_diploid_called_SNPs_filtered_SOR_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/SOR/chrM/by_sample/{sample}_{chrmsM}_SOR_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/SOR/chrM/by_sample/{sample}_{chrmsM}_SOR_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_female_get_performance_metrics_M_SOR:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/SOR/chrM/by_sample/{sample}_{chrms}_SOR_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/SOR/chrM/{sample}_{chrms}_SOR_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/SOR/chrM/{sample}_{chrms}_SOR_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - FS
#-------------------------------------------------------------------------------#
rule default_female_extract_call_VCF_autos_FS:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/FS/{chrmsA}/{chrmsA}_GRCh38_default_gatk_diploid_called_SNPs_filtered_FS_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/FS/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_FS_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/FS/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_FS_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_female_extract_call_VCF_XnonPARs_FS:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/FS/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_FS_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/FS/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_FS_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/FS/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_FS_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_female_extract_call_VCF_XPARs_FS:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/FS/chrX_PARs/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_FS_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/FS/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_FS_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/FS/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_FS_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule default_female_get_performance_metrics_autos_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/FS/{chrms}/by_sample/{sample}_{chrms}_autos_FS_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/FS/autos/{sample}_{chrms}_autos_FS_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/FS/autos/{sample}_{chrms}_autos_FS_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_female_get_performance_metrics_XnonPARs_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/FS/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_FS_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/FS/chrX_nonPARs/{sample}_{chrms}_nonPARs_FS_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/FS/chrX_nonPARs/{sample}_{chrms}_nonPARs_FS_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_female_get_performance_metrics_XPARs_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/FS/chrX_PARs/by_sample/{sample}_{chrms}_PARs_FS_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/FS/chrX_PARs/{sample}_{chrms}_PARs_FS_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/FS/chrX_PARs/{sample}_{chrms}_PARs_FS_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule default_female_extract_call_VCF_M_FS:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/FS/chrM/{chrmsM}_GRCh38_default_gatk_diploid_called_SNPs_filtered_FS_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/FS/chrM/by_sample/{sample}_{chrmsM}_FS_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/FS/chrM/by_sample/{sample}_{chrmsM}_FS_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_female_get_performance_metrics_M_FS:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/FS/chrM/by_sample/{sample}_{chrms}_FS_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/FS/chrM/{sample}_{chrms}_FS_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/FS/chrM/{sample}_{chrms}_FS_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - MQ
#-------------------------------------------------------------------------------#
rule default_female_extract_call_VCF_autos_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQ/{chrmsA}/{chrmsA}_GRCh38_default_gatk_diploid_called_SNPs_filtered_MQ_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQ/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQ_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQ/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQ_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_female_extract_call_VCF_XnonPARs_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQ/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_MQ_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQ/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQ_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQ/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQ_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_female_extract_call_VCF_XPARs_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQ/chrX_PARs/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_MQ_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQ/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQ_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQ/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQ_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """



rule default_female_get_performance_metrics_autos_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQ/{chrms}/by_sample/{sample}_{chrms}_autos_MQ_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/MQ/autos/{sample}_{chrms}_autos_MQ_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/MQ/autos/{sample}_{chrms}_autos_MQ_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_female_get_performance_metrics_XnonPARs_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQ/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_MQ_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/MQ/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQ_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/MQ/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQ_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_female_get_performance_metrics_XPARs_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQ/chrX_PARs/by_sample/{sample}_{chrms}_PARs_MQ_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/MQ/chrX_PARs/{sample}_{chrms}_PARs_MQ_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/MQ/chrX_PARs/{sample}_{chrms}_PARs_MQ_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule default_female_extract_call_VCF_M_MQ:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQ/chrM/{chrmsM}_GRCh38_default_gatk_diploid_called_SNPs_filtered_MQ_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQ/chrM/by_sample/{sample}_{chrmsM}_MQ_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQ/chrM/by_sample/{sample}_{chrmsM}_MQ_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_female_get_performance_metrics_M_MQ:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQ/chrM/by_sample/{sample}_{chrms}_MQ_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/MQ/chrM/{sample}_{chrms}_MQ_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/MQ/chrM/{sample}_{chrms}_MQ_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - MQRankSum
#-------------------------------------------------------------------------------#
rule default_female_extract_call_VCF_autos_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQRankSum/{chrmsA}/{chrmsA}_GRCh38_default_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_MQRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_female_extract_call_VCF_XnonPARs_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQRankSum/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_MQRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_female_extract_call_VCF_XPARs_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQRankSum/chrX_PARs/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_MQRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule default_female_get_performance_metrics_autos_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQRankSum/{chrms}/by_sample/{sample}_{chrms}_autos_MQRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/MQRankSum/autos/{sample}_{chrms}_autos_MQRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/MQRankSum/autos/{sample}_{chrms}_autos_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_female_get_performance_metrics_XnonPARs_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQRankSum/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_MQRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/MQRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/MQRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_female_get_performance_metrics_XPARs_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQRankSum/chrX_PARs/by_sample/{sample}_{chrms}_PARs_MQRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/MQRankSum/chrX_PARs/{sample}_{chrms}_PARs_MQRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/MQRankSum/chrX_PARs/{sample}_{chrms}_PARs_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule default_female_extract_call_VCF_M_MQRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQRankSum/chrM/{chrmsM}_GRCh38_default_gatk_diploid_called_SNPs_filtered_MQRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQRankSum/chrM/by_sample/{sample}_{chrmsM}_MQRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQRankSum/chrM/by_sample/{sample}_{chrmsM}_MQRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_female_get_performance_metrics_M_MQRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/MQRankSum/chrM/by_sample/{sample}_{chrms}_MQRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/MQRankSum/chrM/{sample}_{chrms}_MQRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/MQRankSum/chrM/{sample}_{chrms}_MQRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - ReadPosRankSum
#-------------------------------------------------------------------------------#
rule default_female_extract_call_VCF_autos_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/ReadPosRankSum/{chrmsA}/{chrmsA}_GRCh38_default_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/ReadPosRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_ReadPosRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/ReadPosRankSum/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_ReadPosRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_female_extract_call_VCF_XnonPARs_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/ReadPosRankSum/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/ReadPosRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_ReadPosRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/ReadPosRankSum/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_ReadPosRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_female_extract_call_VCF_XPARs_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/ReadPosRankSum/chrX_PARs/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/ReadPosRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_ReadPosRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/ReadPosRankSum/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_ReadPosRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule default_female_get_performance_metrics_autos_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/ReadPosRankSum/{chrms}/by_sample/{sample}_{chrms}_autos_ReadPosRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/ReadPosRankSum/autos/{sample}_{chrms}_autos_ReadPosRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/ReadPosRankSum/autos/{sample}_{chrms}_autos_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_female_get_performance_metrics_XnonPARs_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/ReadPosRankSum/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_ReadPosRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/ReadPosRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_ReadPosRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/ReadPosRankSum/chrX_nonPARs/{sample}_{chrms}_nonPARs_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_female_get_performance_metrics_XPARs_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/ReadPosRankSum/chrX_PARs/by_sample/{sample}_{chrms}_PARs_ReadPosRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/ReadPosRankSum/chrX_PARs/{sample}_{chrms}_PARs_ReadPosRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/ReadPosRankSum/chrX_PARs/{sample}_{chrms}_PARs_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule default_female_extract_call_VCF_M_ReadPosRankSum:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/ReadPosRankSum/chrM/{chrmsM}_GRCh38_default_gatk_diploid_called_SNPs_filtered_ReadPosRankSum_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/ReadPosRankSum/chrM/by_sample/{sample}_{chrmsM}_ReadPosRankSum_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/ReadPosRankSum/chrM/by_sample/{sample}_{chrmsM}_ReadPosRankSum_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_female_get_performance_metrics_M_ReadPosRankSum:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/ReadPosRankSum/chrM/by_sample/{sample}_{chrms}_ReadPosRankSum_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/ReadPosRankSum/chrM/{sample}_{chrms}_ReadPosRankSum_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/ReadPosRankSum/chrM/{sample}_{chrms}_ReadPosRankSum_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - DP
#-------------------------------------------------------------------------------#
rule default_female_extract_call_VCF_autos_DP:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/DP/{chrmsA}/{chrmsA}_GRCh38_default_gatk_diploid_called_SNPs_filtered_DP_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/DP/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_DP_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/DP/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_DP_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_female_extract_call_VCF_XnonPARs_DP:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/DP/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_DP_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/DP/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_DP_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/DP/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_DP_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_female_extract_call_VCF_XPARs_DP:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/DP/chrX_PARs/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_DP_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/DP/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_DP_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/DP/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_DP_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """



rule default_female_get_performance_metrics_autos_DP:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/DP/{chrms}/by_sample/{sample}_{chrms}_autos_DP_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/DP/autos/{sample}_{chrms}_autos_DP_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/DP/autos/{sample}_{chrms}_autos_DP_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_female_get_performance_metrics_XnonPARs_DP:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/DP/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_DP_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/DP/chrX_nonPARs/{sample}_{chrms}_nonPARs_DP_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/DP/chrX_nonPARs/{sample}_{chrms}_nonPARs_DP_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_female_get_performance_metrics_XPARs_DP:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/DP/chrX_PARs/by_sample/{sample}_{chrms}_PARs_DP_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/DP/chrX_PARs/{sample}_{chrms}_PARs_DP_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/DP/chrX_PARs/{sample}_{chrms}_PARs_DP_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule default_female_extract_call_VCF_M_DP:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/DP/chrM/{chrmsM}_GRCh38_default_gatk_diploid_called_SNPs_filtered_DP_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/DP/chrM/by_sample/{sample}_{chrmsM}_DP_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/DP/chrM/by_sample/{sample}_{chrmsM}_DP_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_female_get_performance_metrics_M_DP:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/DP/chrM/by_sample/{sample}_{chrms}_DP_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/DP/chrM/{sample}_{chrms}_DP_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/DP/chrM/{sample}_{chrms}_DP_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - AN
#-------------------------------------------------------------------------------#
rule default_female_extract_call_VCF_autos_AN:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/AN/{chrmsA}/{chrmsA}_GRCh38_default_gatk_diploid_called_SNPs_filtered_AN_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/AN/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_AN_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/AN/{chrmsA}/by_sample/{sample}_{chrmsA}_autos_AN_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_female_extract_call_VCF_XnonPARs_AN:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/AN/chrX_nonPARs/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_AN_{filter}_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/AN/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_AN_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/AN/chrX_nonPARs/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_AN_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_female_extract_call_VCF_XPARs_AN:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/AN/chrX_PARs/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_SNPs_filtered_AN_{filter}_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/AN/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_AN_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/AN/chrX_PARs/by_sample/{sample}_{chrmsXPARs}_PARs_AN_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """


rule default_female_get_performance_metrics_autos_AN:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/AN/{chrms}/by_sample/{sample}_{chrms}_autos_AN_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/AN/autos/{sample}_{chrms}_autos_AN_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/AN/autos/{sample}_{chrms}_autos_AN_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_female_get_performance_metrics_XnonPARs_AN:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/AN/chrX_nonPARs/by_sample/{sample}_{chrms}_nonPARs_AN_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/AN/chrX_nonPARs/{sample}_{chrms}_nonPARs_AN_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/AN/chrX_nonPARs/{sample}_{chrms}_nonPARs_AN_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule default_female_get_performance_metrics_XPARs_AN:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/AN/chrX_PARs/by_sample/{sample}_{chrms}_PARs_AN_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/AN/chrX_PARs/{sample}_{chrms}_PARs_AN_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/AN/chrX_PARs/{sample}_{chrms}_PARs_AN_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """


rule default_female_extract_call_VCF_M_AN:
    input:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/AN/chrM/{chrmsM}_GRCh38_default_gatk_diploid_called_SNPs_filtered_AN_{filter}.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/AN/chrM/by_sample/{sample}_{chrmsM}_AN_{filter}")
    output:
        os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/AN/chrM/by_sample/{sample}_{chrmsM}_AN_{filter}.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule default_female_get_performance_metrics_M_AN:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "hard_filtered_vcfs/EUR/females/default/AN/chrM/by_sample/{sample}_{chrms}_AN_{filter}.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "haploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/AN/chrM/{sample}_{chrms}_AN_{filter}_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/AN/chrM/{sample}_{chrms}_AN_{filter}_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """



'''
#-------------------------------------------------------------------------------#
# Step: Get performance metrics on filtered VCFs - No filters
#-------------------------------------------------------------------------------#
rule females_default_extract_females_call_VCF_autos_no_filter:
    input:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/default/{chrmsA}_GRCh38_default_gatk_diploid_called_raw_SNPs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/default/by_sample/{sample}_{chrmsA}_autos_no_filter")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/default/by_sample/{sample}_{chrmsA}_autos_no_filter.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_default_extract_nonPARs_nofilter:
    input:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/default/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_raw_SNPs.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/default/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_raw_SNPs_nonPARs.vcf.gz")
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
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/default/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_raw_SNPs.vcf.gz")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/default/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_raw_SNPs_PARs.vcf.gz")
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
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/default/{chrmsXnonPARs}_GRCh38_default_gatk_diploid_called_raw_SNPs_nonPARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/default/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_no_filter")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/default/by_sample/{sample}_{chrmsXnonPARs}_nonPARs_no_filter.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_default_extract_females_call_VCF_XPARs_no_filter:
    input:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/default/{chrmsXPARs}_GRCh38_default_gatk_diploid_called_raw_SNPs_PARs.vcf.gz")
    params:
        smpl = "{sample}",
        vcf = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/default/by_sample/{sample}_{chrmsXPARs}_PARs_no_filter")
    output:
        os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/default/by_sample/{sample}_{chrmsXPARs}_PARs_no_filter.recode.vcf")
    shell:
        """
        vcftools --gzvcf {input} --indv {params.smpl} --recode --out {params.vcf}
        """

rule females_default_get_performance_metrics_autos_no_filter:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/default/by_sample/{sample}_{chrms}_autos_no_filter.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/no_filter/autos/{sample}_{chrms}_autos_no_filter_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/no_filter/autos/{sample}_{chrms}_autos_no_filter_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_default_get_performance_metrics_XnonPARs_no_filter:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/default/by_sample/{sample}_{chrms}_nonPARs_no_filter.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/no_filter/chrX_nonPARs/{sample}_{chrms}_nonPARs_no_filter_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/no_filter/chrX_nonPARs/{sample}_{chrms}_nonPARs_no_filter_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """

rule females_default_get_performance_metrics_XPARs_no_filter:
    input:
        golden = os.path.join(config["scratch_proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf"),
        called = os.path.join(config["scratch_proj_path"], "joint_called_vcfs/EUR/females/default/by_sample/{sample}_{chrms}_PARs_no_filter.recode.vcf")
    params:
        scriptpth = config["compare_VCFs_path"],
        gploidy = "diploid",
        outstm = os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/no_filter/chrX_PARs/{sample}_{chrms}_PARs_no_filter_golden_vs_called")
    output:
        os.path.join(config["scratch_proj_path"], "compare_VCFs/EUR/females/default/no_filter/chrX_PARs/{sample}_{chrms}_PARs_no_filter_golden_vs_called_performance_metrics.txt")
    shell:
        """
        python3 {params.scriptpth}compare_VCFs_fix.py --golden {input.golden} --called {input.called} --gploidy {params.gploidy} --out {params.outstm}
        """
'''
