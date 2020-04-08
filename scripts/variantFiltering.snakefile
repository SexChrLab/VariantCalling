# Environment: varCalling

configfile: "variantFiltering.config.json"


rule all:
    input:
        #expand(os.path.join(config["out_dir"], "stats_pre_filter/autosomes", "{chrs}.gatk.called.raw.{vcf_options_auto}", "{chrs}.gatk.called.raw.{vcf_options_auto}_annotations.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        #expand(os.path.join(config["out_dir"], "stats_pre_filter/sex_chromosomes", "{chrs}.gatk.called.raw.{vcf_options_x}", "{chrs}.gatk.called.raw.{vcf_options_x}_annotations.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "stats_pre_filter/autosomes", "{chrs}.gatk.called.raw.{vcf_options_auto}", "{chrs}.gatk.called.raw.{vcf_options_auto}_AN_plot.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "stats_pre_filter/autosomes", "{chrs}.gatk.called.raw.{vcf_options_auto}", "{chrs}.gatk.called.raw.{vcf_options_auto}_QD_plot.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "stats_pre_filter/autosomes", "{chrs}.gatk.called.raw.{vcf_options_auto}", "{chrs}.gatk.called.raw.{vcf_options_auto}_DP_plot.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "stats_pre_filter/autosomes", "{chrs}.gatk.called.raw.{vcf_options_auto}", "{chrs}.gatk.called.raw.{vcf_options_auto}_MQ_plot.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "stats_pre_filter/autosomes", "{chrs}.gatk.called.raw.{vcf_options_auto}", "{chrs}.gatk.called.raw.{vcf_options_auto}_DP_statistics.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "stats_pre_filter/sex_chromosomes", "{chrs}.gatk.called.raw.{vcf_options_x}", "{chrs}.gatk.called.raw.{vcf_options_x}_AN_plot.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "stats_pre_filter/sex_chromosomes", "{chrs}.gatk.called.raw.{vcf_options_x}", "{chrs}.gatk.called.raw.{vcf_options_x}_QD_plot.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "stats_pre_filter/sex_chromosomes", "{chrs}.gatk.called.raw.{vcf_options_x}", "{chrs}.gatk.called.raw.{vcf_options_x}_DP_plot.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "stats_pre_filter/sex_chromosomes", "{chrs}.gatk.called.raw.{vcf_options_x}", "{chrs}.gatk.called.raw.{vcf_options_x}_MQ_plot.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "stats_pre_filter/sex_chromosomes", "{chrs}.gatk.called.raw.{vcf_options_x}", "{chrs}.gatk.called.raw.{vcf_options_x}_DP_statistics.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "stats_pre_filter/autosomes", "{chrs}.gatk.called.raw.{vcf_options_auto}", "{chrs}.gatk.called.raw.{vcf_options_auto}_num_sites.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "stats_pre_filter/sex_chromosomes", "{chrs}.gatk.called.raw.{vcf_options_x}", "{chrs}.gatk.called.raw.{vcf_options_x}_num_sites.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        #expand(os.path.join(config["out_dir"], "vqsr/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}_output.recal"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        #expand(os.path.join(config["out_dir"], "vqsr/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}_output.tranches"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        #expand(os.path.join(config["out_dir"], "vqsr/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}_output.plots.R"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        #expand(os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}_output.recal"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        #expand(os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}_output.tranches"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        #expand(os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}_output.plots.R"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        #expand(os.path.join(config["out_dir"], "vqsr/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.vcf.gz"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        #expand(os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.vcf.gz"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "vqsr/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.vcf.gz"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.vcf.gz"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        #expand(os.path.join(config["out_dir"], "stats_post_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv_annotations.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        #expand(os.path.join(config["out_dir"], "stats_post_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv_annotations.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "stats_post_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv_AN_plot.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "stats_post_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv_QD_plot.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "stats_post_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv_DP_plot.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "stats_post_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv_MQ_plot.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "stats_post_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv_DP_statistics.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "stats_post_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv_AN_plot.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "stats_post_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv_QD_plot.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "stats_post_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv_DP_plot.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "stats_post_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv_MQ_plot.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "stats_post_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv_DP_statistics.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "stats_post_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv_num_sites.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "stats_post_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv_num_sites.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        #expand(os.path.join(config["out_dir"], "hard_filter/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.snps.vcf.gz"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        #expand(os.path.join(config["out_dir"], "hard_filter/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.snps.vcf.gz"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        #expand(os.path.join(config["out_dir"], "hard_filter/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.vcf.gz"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        #expand(os.path.join(config["out_dir"], "hard_filter/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.vcf.gz"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "hard_filter/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.vcf.gz"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "hard_filter/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.vcf.gz"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "stats_post_hard_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv_annotations.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "stats_post_hard_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv_annotations.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "stats_post_hard_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv_AN_plot.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "stats_post_hard_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv_QD_plot.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "stats_post_hard_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv_DP_plot.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "stats_post_hard_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv_MQ_plot.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "stats_post_hard_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv_DP_statistics.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "stats_post_hard_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv_AN_plot.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "stats_post_hard_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv_QD_plot.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "stats_post_hard_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv_DP_plot.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "stats_post_hard_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv_MQ_plot.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "stats_post_hard_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv_DP_statistics.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "stats_post_hard_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv_num_sites.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "stats_post_hard_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv_num_sites.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "raw/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.array.sites.vcf.gz"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "raw/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites.vcf.gz"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "vqsr/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.vcf.gz"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.vcf.gz"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "hard_filter/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.vcf.gz"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "hard_filter/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.vcf.gz"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"])


#-------------------------------------------------------------------------------
# Get statistics of raw VCFs #
# Get statistics of raw VCFs: Rule 1
rule getStatsPreFilterAutos:
    input:
        vcf = os.path.join(config["in_vcf_dir"], "{chrs}.gatk.called.raw.{vcf_options_auto}.vcf.gz")
    params:
        script = os.path.join(config["scripts_dir"], "extract_stats_from_vcf.py"),
        AN = "AN",
        QD = "QD",
        MQ = "MQ",
        DP = "DP"
    output:
        stats = os.path.join(config["out_dir"], "stats_pre_filter/autosomes", "{chrs}.gatk.called.raw.{vcf_options_auto}", "{chrs}.gatk.called.raw.{vcf_options_auto}_annotations.txt")
    shell:
        """
        python {params.script} {params.AN} {params.QD} {params.MQ} {params.DP} --vcf {input.vcf} --outfile {output.stats}
        """

# Get statistics of raw VCFs: Rule 2
rule getStatsPreFilterXchr:
    input:
        vcf = os.path.join(config["in_vcf_dir"],"{chrs}.gatk.called.raw.{vcf_options_x}.vcf.gz")
    params:
        script = os.path.join(config["scripts_dir"], "extract_stats_from_vcf.py"),
        AN = "AN",
        QD = "QD",
        MQ = "MQ",
        DP = "DP"
    output:
        stats = os.path.join(config["out_dir"], "stats_pre_filter/sex_chromosomes", "{chrs}.gatk.called.raw.{vcf_options_x}", "{chrs}.gatk.called.raw.{vcf_options_x}_annotations.txt")
    shell:
        """
        python {params.script} {params.AN} {params.QD} {params.MQ} {params.DP} --vcf {input.vcf} --outfile {output.stats}
        """

# Get statistics of raw VCFs: Rule 3
rule plotDensityPreFilterAutos:
    input:
        os.path.join(config["out_dir"], "stats_pre_filter/autosomes", "{chrs}.gatk.called.raw.{vcf_options_auto}", "{chrs}.gatk.called.raw.{vcf_options_auto}_annotations.txt")
    params:
        script = os.path.join(config["scripts_dir"], "plot_stats.R")
    output:
        AN_plot = os.path.join(config["out_dir"], "stats_pre_filter/autosomes", "{chrs}.gatk.called.raw.{vcf_options_auto}", "{chrs}.gatk.called.raw.{vcf_options_auto}_AN_plot.png"),
        QD_plot = os.path.join(config["out_dir"], "stats_pre_filter/autosomes", "{chrs}.gatk.called.raw.{vcf_options_auto}", "{chrs}.gatk.called.raw.{vcf_options_auto}_QD_plot.png"),
        DP_plot = os.path.join(config["out_dir"], "stats_pre_filter/autosomes", "{chrs}.gatk.called.raw.{vcf_options_auto}", "{chrs}.gatk.called.raw.{vcf_options_auto}_DP_plot.png"),
        MQ_plot = os.path.join(config["out_dir"], "stats_pre_filter/autosomes", "{chrs}.gatk.called.raw.{vcf_options_auto}", "{chrs}.gatk.called.raw.{vcf_options_auto}_MQ_plot.png"),
        DP_statistics = os.path.join(config["out_dir"], "stats_pre_filter/autosomes", "{chrs}.gatk.called.raw.{vcf_options_auto}", "{chrs}.gatk.called.raw.{vcf_options_auto}_DP_statistics.txt")
    shell:
        """
        Rscript {params.script} {input} {output.AN_plot} {output.QD_plot} {output.DP_plot} {output.MQ_plot} {output.DP_statistics}
        """

# Get statistics of raw VCFs: Rule 4
rule plotDensityPreFilterXchr:
    input:
        os.path.join(config["out_dir"], "stats_pre_filter/sex_chromosomes", "{chrs}.gatk.called.raw.{vcf_options_x}", "{chrs}.gatk.called.raw.{vcf_options_x}_annotations.txt")
    params:
        script = os.path.join(config["scripts_dir"], "plot_stats.R")
    output:
        AN_plot = os.path.join(config["out_dir"], "stats_pre_filter/sex_chromosomes", "{chrs}.gatk.called.raw.{vcf_options_x}", "{chrs}.gatk.called.raw.{vcf_options_x}_AN_plot.png"),
        QD_plot = os.path.join(config["out_dir"], "stats_pre_filter/sex_chromosomes", "{chrs}.gatk.called.raw.{vcf_options_x}", "{chrs}.gatk.called.raw.{vcf_options_x}_QD_plot.png"),
        DP_plot = os.path.join(config["out_dir"], "stats_pre_filter/sex_chromosomes", "{chrs}.gatk.called.raw.{vcf_options_x}", "{chrs}.gatk.called.raw.{vcf_options_x}_DP_plot.png"),
        MQ_plot = os.path.join(config["out_dir"], "stats_pre_filter/sex_chromosomes", "{chrs}.gatk.called.raw.{vcf_options_x}", "{chrs}.gatk.called.raw.{vcf_options_x}_MQ_plot.png"),
        DP_statistics = os.path.join(config["out_dir"], "stats_pre_filter/sex_chromosomes", "{chrs}.gatk.called.raw.{vcf_options_x}", "{chrs}.gatk.called.raw.{vcf_options_x}_DP_statistics.txt")
    shell:
        """
        Rscript {params.script} {input} {output.AN_plot} {output.QD_plot} {output.DP_plot} {output.MQ_plot} {output.DP_statistics}
        """

# Get statistics of raw VCFs: Rule 5
rule countNumSitesPreFilterAutos:
    input:
        vcf = os.path.join(config["in_vcf_dir"], "{chrs}.gatk.called.raw.{vcf_options_auto}.vcf.gz")
    output:
        siteinfo = os.path.join(config["out_dir"], "stats_pre_filter/autosomes", "{chrs}.gatk.called.raw.{vcf_options_auto}", "{chrs}.gatk.called.raw.{vcf_options_auto}_num_sites.txt")
    shell:
        """
        zcat {input.vcf} | grep -v '#' | wc -l > {output.siteinfo}
        """

# Get statistics of raw VCFs: Rule 6
rule countNumSitesPreFilterXchr:
    input:
        vcf = os.path.join(config["in_vcf_dir"], "{chrs}.gatk.called.raw.{vcf_options_x}.vcf.gz")
    output:
        siteinfo = os.path.join(config["out_dir"], "stats_pre_filter/sex_chromosomes", "{chrs}.gatk.called.raw.{vcf_options_x}", "{chrs}.gatk.called.raw.{vcf_options_x}_num_sites.txt")
    shell:
        """
        zcat {input.vcf} | grep -v '#' | wc -l > {output.siteinfo}
        """
# Finished getting statistics of raw VCFs.

#-------------------------------------------------------------------------------
# VQSR Steps #
# These rules with perform VQSR filtering on the autosomes and x chromosome VCFs
# VQSR Step - Variant Recalibrator. Autosomes
rule VariantRecalibratorAutos:
    input:
        vcf = os.path.join(config["in_vcf_dir"], "{chrs}.gatk.called.raw.{vcf_options_auto}.vcf.gz"),
        hapmap = config["hapmap"],
        omni = config["omni"],
        thousandG = config["thousandG"],
        dbsnp = config["dbsnp"]
    params:
        chrms = "{chrs}",
        ref = lambda wildcards: config[wildcards.vcf_options_auto]["refgenome"]
    output:
        recal = os.path.join(config["out_dir"], "vqsr/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}_output.recal"),
        tranches = os.path.join(config["out_dir"], "vqsr/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}_output.tranches"),
        rplots = os.path.join(config["out_dir"], "vqsr/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}_output.plots.R")
    shell:
        """gatk --java-options "-Xmx16g" VariantRecalibrator """
        """-R {params.ref} -V {input.vcf} -L {params.chrms} """
        """--resource:hapmap,known=false,training=true,truth=true,prior=15.0 {input.hapmap} """
        """--resource:omni,known=false,training=true,truth=false,prior=12.0 {input.omni} """
        """--resource:1000G,known=false,training=true,truth=false,prior=10.0 {input.thousandG} """
        """--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp} """
        """-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR """
        """-mode SNP """
        """-O {output.recal} """
        """--tranches-file {output.tranches} """
        """--rscript-file {output.rplots} """
#"""-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff """
'''
rule VariantRecalibratorXchr:
    input:
        vcf = os.path.join(config["in_vcf_dir"], "{chrs}.gatk.called.raw.{vcf_options_x}.vcf.gz"),
        hapmap = config["hapmap"],
        omni = config["omni"],
        thousandG = config["thousandG"],
        dbsnp = config["dbsnp"]
    params:
        chrms = "{chrs}",
        ref = lambda wildcards: config[wildcards.vcf_options_x]["refgenome"]
    output:
        recal = os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}_output.recal"),
        tranches = os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}_output.tranches"),
        rplots = os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}_output.plots.R")
    shell:
        """gatk --java-options "-Xmx16g" VariantRecalibrator """
        """-R {params.ref} -V {input.vcf}  -L {params.chrms} """
        """--resource:hapmap,known=false,training=true,truth=true,prior=15.0 {input.hapmap} """
        """--resource:omni,known=false,training=true,truth=false,prior=12.0 {input.omni} """
        """--resource:1000G,known=false,training=true,truth=false,prior=10.0 {input.thousandG} """
        """--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp} """
        """-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR """
        """-mode SNP """
        """-O {output.recal} """
        """--tranches-file {output.tranches} """
        """--rscript-file {output.rplots} """
'''

# VQSR Step - Variant Recalibrator. X chromosome
rule VariantRecalibratorXchr:
    input:
        vcf = os.path.join(config["in_vcf_dir"], "{chrs}.gatk.called.raw.{vcf_options_x}.vcf.gz"),
        hapmap = config["hapmap"],
        omni = config["omni"],
        thousandG = config["thousandG"],
        dbsnp = config["dbsnp"]
    params:
        chrms = "{chrs}",
        ref = lambda wildcards: config[wildcards.vcf_options_x]["refgenome"],
    output:
        recal = os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}_output.recal"),
        tranches = os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}_output.tranches"),
        rplots = os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}_output.plots.R")
    run:
        if wildcards.vcf_options_x == "haploid":
            shell("gatk VariantRecalibrator -R {params.ref} -V {input.vcf}  -L {params.chrms} --resource:hapmap,known=false,training=true,truth=true,prior=15.0 {input.hapmap} --resource:omni,known=false,training=true,truth=false,prior=12.0 {input.omni} --resource:1000G,known=false,training=true,truth=false,prior=10.0 {input.thousandG} --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp} -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP -O {output.recal} --max-gaussians 4 --tranches-file {output.tranches} --rscript-file {output.rplots}")
            print("haploid used")
        else:
            shell("gatk VariantRecalibrator -R {params.ref} -V {input.vcf}  -L {params.chrms} --resource:hapmap,known=false,training=true,truth=true,prior=15.0 {input.hapmap} --resource:omni,known=false,training=true,truth=false,prior=12.0 {input.omni} --resource:1000G,known=false,training=true,truth=false,prior=10.0 {input.thousandG} --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp} -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP -O {output.recal} --tranches-file {output.tranches} --rscript-file {output.rplots}")

# VQSR Step - Apply VQSR. Autosomes
rule ApplyVQSRAutos:
    input:
        vcf = os.path.join(config["in_vcf_dir"], "{chrs}.gatk.called.raw.{vcf_options_auto}.vcf.gz"),
        recal = os.path.join(config["out_dir"], "vqsr/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}_output.recal"),
        tranches = os.path.join(config["out_dir"], "vqsr/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}_output.tranches")
    params:
        chrms = "{chrs}",
        ref = lambda wildcards: config[wildcards.vcf_options_auto]["refgenome"]
    output:
        outvcf = os.path.join(config["out_dir"], "vqsr/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.vcf.gz")
    shell:
        """gatk --java-options "-Xmx16g" ApplyVQSR """
        """-R {params.ref} """
        """-V {input.vcf} """
        """-O {output.outvcf} """
        """-L {params.chrms} """
        """--truth-sensitivity-filter-level 99.0 """
        """--tranches-file {input.tranches} """
        """--recal-file {input.recal} """
        """-mode SNP """

# VQSR Step - Apply VQSR. X chromosome
rule ApplyVQSRXchr:
    input:
        vcf = os.path.join(config["in_vcf_dir"], "{chrs}.gatk.called.raw.{vcf_options_x}.vcf.gz"),
        recal = os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}_output.recal"),
        tranches = os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}_output.tranches")
    params:
        chrms = "{chrs}",
        ref = lambda wildcards: config[wildcards.vcf_options_x]["refgenome"]
    output:
        outvcf = os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.vcf.gz")
    shell:
        """gatk --java-options "-Xmx16g" ApplyVQSR """
        """-R {params.ref} """
        """-V {input.vcf} """
        """-O {output.outvcf} """
        """-L {params.chrms} """
        """--truth-sensitivity-filter-level 99.0 """
        """--tranches-file {input.tranches} """
        """--recal-file {input.recal} """
        """-mode SNP """

# VQSR Step - Select Variants. Autosomes
rule SelectVariantsAutos:
    input:
        vcf = os.path.join(config["out_dir"], "vqsr/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.vcf.gz")
    params:
        chrms = "{chrs}",
        ref = lambda wildcards: config[wildcards.vcf_options_auto]["refgenome"]
    output:
        outvcf = os.path.join(config["out_dir"], "vqsr/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.vcf.gz")
    shell:
        """gatk --java-options "-Xmx16g" SelectVariants """
        """-R {params.ref} """
        """-V {input.vcf} """
        """ -L {params.chrms} """
        """--exclude-filtered """
        """-O {output.outvcf} """

# VQSR Step - Select Variants. X chromosome
rule SelectVariantsXchr:
    input:
        vcf = os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.vcf.gz")
    params:
        chrms = "{chrs}",
        ref = lambda wildcards: config[wildcards.vcf_options_x]["refgenome"]
    output:
        outvcf = os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.vcf.gz")
    shell:
        """gatk --java-options "-Xmx16g" SelectVariants """
        """-R {params.ref} """
        """-V {input.vcf} """
        """ -L {params.chrms} """
        """--exclude-filtered """
        """-O {output.outvcf} """

#-------------------------------------------------------------------------------
# Get statistics for VCFs post VQSR #
# Get statistics for VCFs post VQSR: Rule 1
rule getStatsPostFilterAutos:
    input:
        vcfauto = os.path.join(config["out_dir"], "vqsr/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.vcf.gz")
    params:
        script = os.path.join(config["scripts_dir"], "extract_stats_from_vcf.py"),
        AN = "AN",
        QD = "QD",
        MQ = "MQ",
        DP = "DP"
    output:
        statsauto = os.path.join(config["out_dir"], "stats_post_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv_annotations.txt")
    shell:
        """
        python {params.script} {params.AN} {params.QD} {params.MQ} {params.DP} --vcf {input.vcfauto} --outfile {output.statsauto}
        """

# Get statistics for VCFs post VQSR: Rule 2
rule getStatsPostFilterXchr:
    input:
        vcfx = os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.vcf.gz")
    params:
        script = os.path.join(config["scripts_dir"], "extract_stats_from_vcf.py"),
        AN = "AN",
        QD = "QD",
        MQ = "MQ",
        DP = "DP"
    output:
        statsx = os.path.join(config["out_dir"], "stats_post_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv_annotations.txt")
    shell:
        """
        python {params.script} {params.AN} {params.QD} {params.MQ} {params.DP} --vcf {input.vcfx} --outfile {output.statsx}
        """

# Get statistics for VCFs post VQSR: Rule 3
rule plotDensityPostFilterAutos:
    input:
        os.path.join(config["out_dir"], "stats_post_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv_annotations.txt"),
    params:
        script = os.path.join(config["scripts_dir"], "plot_stats.R")
    output:
        AN_plot = os.path.join(config["out_dir"], "stats_post_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv_AN_plot.png"),
        QD_plot = os.path.join(config["out_dir"], "stats_post_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv_QD_plot.png"),
        DP_plot = os.path.join(config["out_dir"], "stats_post_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv_DP_plot.png"),
        MQ_plot = os.path.join(config["out_dir"], "stats_post_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv_MQ_plot.png"),
        DP_statistics = os.path.join(config["out_dir"], "stats_post_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv_DP_statistics.txt")
    shell:
        """
        Rscript {params.script} {input} {output.AN_plot} {output.QD_plot} {output.DP_plot} {output.MQ_plot} {output.DP_statistics}
        """

# Get statistics for VCFs post VQSR: Rule 4
rule plotDensityPostFilterXchr:
    input:
        os.path.join(config["out_dir"], "stats_post_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv_annotations.txt"),
    params:
        script = os.path.join(config["scripts_dir"], "plot_stats.R")
    output:
        AN_plot = os.path.join(config["out_dir"], "stats_post_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv_AN_plot.png"),
        QD_plot = os.path.join(config["out_dir"], "stats_post_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv_QD_plot.png"),
        DP_plot = os.path.join(config["out_dir"], "stats_post_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv_DP_plot.png"),
        MQ_plot = os.path.join(config["out_dir"], "stats_post_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv_MQ_plot.png"),
        DP_statistics = os.path.join(config["out_dir"], "stats_post_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv_DP_statistics.txt")
    shell:
        """
        Rscript {params.script} {input} {output.AN_plot} {output.QD_plot} {output.DP_plot} {output.MQ_plot} {output.DP_statistics}
        """

# Get statistics for VCFs post VQSR: Rule 5
rule countNumSitesPostFilterAutos:
    input:
        vcf = os.path.join(config["out_dir"], "vqsr/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.vcf.gz")
    output:
        siteinfo = os.path.join(config["out_dir"], "stats_post_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv_num_sites.txt")
    shell:
        """
        zcat {input.vcf} | grep -v '#' | wc -l > {output.siteinfo}
        """

# Get statistics for VCFs post VQSR: Rule 6
rule countNumSitesPostFilterXchr:
    input:
        vcf = os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.vcf.gz")
    output:
        siteinfo = os.path.join(config["out_dir"], "stats_post_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv_num_sites.txt")
    shell:
        """
        zcat {input.vcf} | grep -v '#' | wc -l > {output.siteinfo}
        """

#-------------------------------------------------------------------------------
# Hard Filtering Steps #
# These rules with perform VQSR filtering on the autosomes and x chromosome VCFs
# Hard Filtering Step - Select Variants. Autosomes
rule SelectVariantsHardFilterAutos:
    input:
        vcf = os.path.join(config["in_vcf_dir"], "{chrs}.gatk.called.raw.{vcf_options_auto}.vcf.gz"),
    params:
        chrms = "{chrs}",
        ref = lambda wildcards: config[wildcards.vcf_options_auto]["refgenome"],
    output:
        ovcf = os.path.join(config["out_dir"], "hard_filter/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.snps.vcf.gz")
    shell:
        """
        gatk SelectVariants -R {params.ref} -V {input.vcf} -L {params.chrms} -select-type SNP -O {output.ovcf}
        """

# Hard Filtering Step - Select Variants. X chromosome
rule SelectVariantsHardFilterXchr:
    input:
        vcf = os.path.join(config["in_vcf_dir"], "{chrs}.gatk.called.raw.{vcf_options_x}.vcf.gz"),
    params:
        chrms = "{chrs}",
        ref = lambda wildcards: config[wildcards.vcf_options_x]["refgenome"]
    output:
        ovcf = os.path.join(config["out_dir"], "hard_filter/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.snps.vcf.gz")
    shell:
        """
        gatk SelectVariants -R {params.ref} -V {input.vcf} -L {params.chrms} -select-type SNP -O {output.ovcf}
        """

# Hard Filtering Step - Variant Filtration. Autosomes
rule VariantFiltrationHardFilterAutos:
    input:
        vcf = os.path.join(config["out_dir"], "hard_filter/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.snps.vcf.gz")
    params:
        chrms = "{chrs}",
        ref = lambda wildcards: config[wildcards.vcf_options_auto]["refgenome"],
        qd = lambda wildcards: config[wildcards.filtering_options]["QD"],
        qual = lambda wildcards: config[wildcards.filtering_options]["QUAL"],
        sor = lambda wildcards: config[wildcards.filtering_options]["SOR"],
        fs = lambda wildcards: config[wildcards.filtering_options]["FS"],
        mq = lambda wildcards: config[wildcards.filtering_options]["MQ"],
        mqranksum = lambda wildcards: config[wildcards.filtering_options]["MQRankSum"],
        readposranksum = lambda wildcards: config[wildcards.filtering_options]["ReadPosRankSum"],
        intermediatevcf = os.path.join(config["out_dir"], "hard_filter/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}")
    output:
        ovcf = os.path.join(config["out_dir"], "hard_filter/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.vcf.gz")
    shell:
        """
        gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "QD < {params.qd}" --filter-name "QD{params.qd}" -filter "QUAL < {params.qual}" --filter-name "QUAL{params.qual}" -filter "SOR > {params.sor}" --filter-name "SOR{params.sor}" -filter "FS > {params.fs}" --filter-name "FS{params.fs}" -filter "MQ < {params.mq}" --filter-name "MQ{params.mq}" -filter "MQRankSum < {params.mqranksum}" --filter-name "MQRankSum{params.mqranksum}" -filter "ReadPosRankSum < {params.readposranksum}" --filter-name "ReadPosRankSum{params.readposranksum}" -O {params.intermediatevcf}.vcf.gz;
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf}.vcf.gz -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """

# Hard Filtering Step - Variant Filtration. X chromosome
rule VariantFiltrationHardFilterXchr:
    input:
        vcf = os.path.join(config["out_dir"], "hard_filter/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.snps.vcf.gz")
    params:
        chrms = "{chrs}",
        ref = lambda wildcards: config[wildcards.vcf_options_x]["refgenome"],
        qd = lambda wildcards: config[wildcards.filtering_options]["QD"],
        qual = lambda wildcards: config[wildcards.filtering_options]["QUAL"],
        sor = lambda wildcards: config[wildcards.filtering_options]["SOR"],
        fs = lambda wildcards: config[wildcards.filtering_options]["FS"],
        mq = lambda wildcards: config[wildcards.filtering_options]["MQ"],
        mqranksum = lambda wildcards: config[wildcards.filtering_options]["MQRankSum"],
        readposranksum = lambda wildcards: config[wildcards.filtering_options]["ReadPosRankSum"],
        intermediatevcf = os.path.join(config["out_dir"], "hard_filter/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}")
    output:
        ovcf = os.path.join(config["out_dir"], "hard_filter/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.vcf.gz")
    shell:
        """gatk VariantFiltration -R {params.ref} -V {input.vcf} -L {params.chrms} -filter "QD < {params.qd}" --filter-name "QD{params.qd}" -filter "QUAL < {params.qual}" --filter-name "QUAL{params.qual}" -filter "SOR > {params.sor}" --filter-name "SOR{params.sor}" -filter "FS > {params.fs}" --filter-name "FS{params.fs}" -filter "MQ < {params.mq}" --filter-name "MQ{params.mq}" -filter "MQRankSum < {params.mqranksum}" --filter-name "MQRankSum{params.mqranksum}" -filter "ReadPosRankSum < {params.readposranksum}" --filter-name "ReadPosRankSum{params.readposranksum}" -O {params.intermediatevcf}.vcf.gz;
        gatk --java-options "-Xmx16g" SelectVariants -R {params.ref} -V {params.intermediatevcf}.vcf.gz -L {params.chrms} --exclude-filtered -O {output.ovcf}
        """

#-------------------------------------------------------------------------------
# Get statistics for VCFs post Hard Filtering #
# Get statistics for VCFs post Hard Filtering: Rule 1
rule getStatsPostHardFilterAutos:
    input:
        vcf = os.path.join(config["out_dir"], "hard_filter/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.vcf.gz")
    params:
        script = os.path.join(config["scripts_dir"], "extract_stats_from_vcf.py"),
        AN = "AN",
        QD = "QD",
        MQ = "MQ",
        DP = "DP"
    output:
        stats = os.path.join(config["out_dir"], "stats_post_hard_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv_annotations.txt")
    shell:
        """
        python {params.script} {params.AN} {params.QD} {params.MQ} {params.DP} --vcf {input.vcf} --outfile {output.stats}
        """

# Get statistics for VCFs post Hard Filtering: Rule 2
rule getStatsPostHardFilterXchr:
    input:
        vcf = os.path.join(config["out_dir"], "hard_filter/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.vcf.gz")
    params:
        script = os.path.join(config["scripts_dir"], "extract_stats_from_vcf.py"),
        AN = "AN",
        QD = "QD",
        MQ = "MQ",
        DP = "DP"
    output:
        stats = os.path.join(config["out_dir"], "stats_post_hard_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv_annotations.txt")
    shell:
        """
        python {params.script} {params.AN} {params.QD} {params.MQ} {params.DP} --vcf {input.vcf} --outfile {output.stats}
        """

# Get statistics for VCFs post Hard Filtering: Rule 3
rule plotDensityPostHardFilterAutos:
    input:
        os.path.join(config["out_dir"], "stats_post_hard_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv_annotations.txt")
    params:
        script = os.path.join(config["scripts_dir"], "plot_stats.R")
    output:
        AN_plot = os.path.join(config["out_dir"], "stats_post_hard_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv_AN_plot.png"),
        QD_plot = os.path.join(config["out_dir"], "stats_post_hard_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv_QD_plot.png"),
        DP_plot = os.path.join(config["out_dir"], "stats_post_hard_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv_DP_plot.png"),
        MQ_plot = os.path.join(config["out_dir"], "stats_post_hard_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv_MQ_plot.png"),
        DP_statistics = os.path.join(config["out_dir"], "stats_post_hard_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv_DP_statistics.txt")
    shell:
        """
        Rscript {params.script} {input} {output.AN_plot} {output.QD_plot} {output.DP_plot} {output.MQ_plot} {output.DP_statistics}
        """

# Get statistics for VCFs post Hard Filtering: Rule 4
rule plotDensityPostHardFilterXchr:
    input:
        os.path.join(config["out_dir"], "stats_post_hard_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv_annotations.txt")
    params:
        script = os.path.join(config["scripts_dir"], "plot_stats.R")
    output:
        AN_plot = os.path.join(config["out_dir"], "stats_post_hard_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv_AN_plot.png"),
        QD_plot = os.path.join(config["out_dir"], "stats_post_hard_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv_QD_plot.png"),
        DP_plot = os.path.join(config["out_dir"], "stats_post_hard_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv_DP_plot.png"),
        MQ_plot = os.path.join(config["out_dir"], "stats_post_hard_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv_MQ_plot.png"),
        DP_statistics = os.path.join(config["out_dir"], "stats_post_hard_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv_DP_statistics.txt")
    shell:
        """
        Rscript {params.script} {input} {output.AN_plot} {output.QD_plot} {output.DP_plot} {output.MQ_plot} {output.DP_statistics}
        """

# Get statistics for VCFs post Hard Filtering: Rule 5
rule countNumSitesPostHardFilterAutos:
    input:
        vcf = os.path.join(config["out_dir"], "hard_filter/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.vcf.gz")
    output:
        siteinfo = os.path.join(config["out_dir"], "stats_post_hard_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv_num_sites.txt")
    shell:
        """
        zcat {input.vcf} | grep -v '#' | wc -l > {output.siteinfo}
        """

# Get statistics for VCFs post Hard Filtering: Rule 6
rule countNumSitesPostHardFilterXchr:
    input:
        vcf = os.path.join(config["out_dir"], "hard_filter/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.vcf.gz")
    output:
        siteinfo = os.path.join(config["out_dir"], "stats_post_hard_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv_num_sites.txt")
    shell:
        """
        zcat {input.vcf} | grep -v '#' | wc -l > {output.siteinfo}
        """

#-------------------------------------------------------------------------------
# Extract array site from VCFs #
# For the raw VCFs, extract the sites from the array. We will eventually calculate
# SFS and pi for unfiltered and filtered VCFs
# Raw VCFs #
# Autosomes
rule extractArraySitesPreFilterAutos:
    input:
        vcf = os.path.join(config["in_vcf_dir"], "{chrs}.gatk.called.raw.{vcf_options_auto}.vcf.gz")
    params:
        posfn = config["array_pos"]
    output:
        os.path.join(config["out_dir"], "raw/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.array.sites.vcf.gz")
    shell:
        "vcftools --gzvcf {input.vcf} --positions {params.posfn} --recode --stdout | bgzip -c > {output}"

# Sex chromosomes
rule extractArraySitesPreFilterXchr:
    input:
        vcf = os.path.join(config["in_vcf_dir"], "{chrs}.gatk.called.raw.{vcf_options_x}.vcf.gz")
    params:
        posfn = config["array_pos"]
    output:
        os.path.join(config["out_dir"], "raw/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites.vcf.gz")
    shell:
        "vcftools --gzvcf {input.vcf} --positions {params.posfn} --recode --stdout | bgzip -c > {output}"

# VQSR VCFs #
# Autosomes
rule extractArraySitesVQSRAutos:
    input:
        vcf = os.path.join(config["out_dir"], "vqsr/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.vcf.gz")
    params:
        posfn = config["array_pos"]
    output:
        os.path.join(config["out_dir"], "vqsr/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.vcf.gz")
    shell:
        "vcftools --gzvcf {input.vcf} --positions {params.posfn} --recode --stdout | bgzip -c > {output}"

# X chromosome
rule extractArraySitesVQSRXchr:
    input:
        vcf = os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.vcf.gz")
    params:
        posfn = config["array_pos"]
    output:
        os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.vcf.gz")
    shell:
        "vcftools --gzvcf {input.vcf} --positions {params.posfn} --recode --stdout | bgzip -c > {output}"

# Hard Filtered VCFs#
# Autosomes
rule extractArraySitesHardFilterAutos:
    input:
        vcf = os.path.join(config["out_dir"], "hard_filter/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.vcf.gz")
    params:
        posfn = config["array_pos"]
    output:
        os.path.join(config["out_dir"], "hard_filter/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.vcf.gz")
    shell:
        "vcftools --gzvcf {input.vcf} --positions {params.posfn} --recode --stdout | bgzip -c > {output}"

# X chromosome
rule extractArraySitesHardFilterXchr:
    input:
        vcf = os.path.join(config["out_dir"], "hard_filter/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.vcf.gz")
    params:
        posfn = config["array_pos"]
    output:
        os.path.join(config["out_dir"], "hard_filter/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.vcf.gz")
    shell:
        "vcftools --gzvcf {input.vcf} --positions {params.posfn} --recode --stdout | bgzip -c > {output}"
