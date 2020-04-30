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
        expand(os.path.join(config["out_dir"], "hard_filter/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.vcf.gz"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/pre_filter/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.array.sites.SFS.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "results/pre_filter/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites.SFS.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "results/vqsr/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.SFS.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "results/vqsr/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.SFS.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "results/hard_filter/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.SFS.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/hard_filter/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.SFS.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/merged/array_sites/sfs/autosomes", "{chrs}", "{chrs}.{vcf_options_auto}.raw.vqsr.{filtering_options}.SFS.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/merged/array_sites/sfs/sex_chromosomes", "{chrs}", "{chrs}.{vcf_options_x}.raw.vqsr.{filtering_options}.SFS.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "stats_pre_filter_array_sites/autosomes", "{chrs}.gatk.called.raw.{vcf_options_auto}.array.sites", "{chrs}.gatk.called.raw.{vcf_options_auto}.array.sites_num_sites.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "stats_pre_filter_array_sites/sex_chromosomes", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites_num_sites.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "stats_post_filter_array_sites/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites_num_sites.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "stats_post_filter_array_sites/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites_num_sites.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "stats_post_hard_filter_array_sites/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites_num_sites.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "stats_post_hard_filter_array_sites/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites_num_sites.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/pre_filter/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.SFS.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "results/pre_filter/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.SFS.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "results/vqsr/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.SFS.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "results/vqsr/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.SFS.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "results/hard_filter/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.SFS.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/hard_filter/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.SFS.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/merged/all_variable_sites/sfs/autosomes", "{chrs}", "{chrs}.{vcf_options_auto}.raw.vqsr.{filtering_options}.SFS.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/merged/all_variable_sites/sfs/sex_chromosomes", "{chrs}", "{chrs}.{vcf_options_x}.raw.vqsr.{filtering_options}.SFS.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/pre_filter/pi/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.pi.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "results/pre_filter/pi/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.pi.per.site.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "raw/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.noPARs.noXTR.vcf.gz"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "results/pre_filter/pi/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.noPARs.noXTR.pi.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "results/pre_filter/pi/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.noPARs.noXTR.pi.per.site.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "results/vqsr/pi/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.pi.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "results/vqsr/pi/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.pi.per.site.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.noPARs.noXTR.vcf.gz"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "results/vqsr/pi/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.noPARs.noXTR.pi.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "results/vqsr/pi/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.noPARs.noXTR.pi.per.site.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "results/hard_filter/pi/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.pi.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/hard_filter/pi/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.pi.per.site.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "hard_filter/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.noPARs.noXTR.vcf.gz"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/hard_filter/pi/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.noPARs.noXTR.pi.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/hard_filter/pi/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.noPARs.noXTR.pi.per.site.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/pre_filter/pi/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.array.sites.pi.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "results/pre_filter/pi/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.array.sites.pi.per.site.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "raw/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites.noPARs.noXTR.vcf.gz"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "results/pre_filter/pi/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites.noPARs.noXTR.pi.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "results/pre_filter/pi/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites.noPARs.noXTR.pi.per.site.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "results/vqsr/pi/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.pi.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "results/vqsr/pi/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.pi.per.site.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.noPARs.noXTR.vcf.gz"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "results/vqsr/pi/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.noPARs.noXTR.pi.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "results/vqsr/pi/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.noPARs.noXTR.pi.per.site.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "results/hard_filter/pi/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.pi.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/hard_filter/pi/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.pi.per.site.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "hard_filter/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.noPARs.noXTR.vcf.gz"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/hard_filter/pi/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.noPARs.noXTR.pi.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/hard_filter/pi/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.noPARs.noXTR.pi.per.site.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.chr.pos.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.chr.pos.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites/autosomes", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.venn.diagram.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.chr.pos.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.chr.pos.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites/sex_chromosomes", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.venn.diagram.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites_array_sites/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.chr.pos.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites_array_sites/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.chr.pos.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites_array_sites/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.chr.pos.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites_array_sites/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.chr.pos.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites_array_sites/autosomes", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.array.sites.venn.diagram.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites_array_sites/sex_chromosomes", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.array.sites.venn.diagram.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "raw/autosomes/filters", "{chrs}.gatk.called.raw.{vcf_options_auto}.high.qual.neutral.sites.vcf.gz"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "raw/sex_chromosomes/filters", "{chrs}.gatk.called.raw.{vcf_options_x}.high.qual.neutral.sites.noPARsXTR.vcf.gz"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "vqsr/autosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.high.qual.neutral.sites.vcf.gz"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "vqsr/sex_chromosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.high.qual.neutral.sites.noPARsXTR.vcf.gz"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "hard_filter/autosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.high.qual.neutral.sites.vcf.gz"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "hard_filter/sex_chromosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.high.qual.neutral.sites.noPARsXTR.vcf.gz"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "raw/autosomes/filters", "{chrs}.gatk.called.raw.{vcf_options_auto}.array.sites.high.qual.neutral.sites.vcf.gz"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "raw/sex_chromosomes/filters", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites.high.qual.neutral.sites.noPARsXTR.vcf.gz"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "vqsr/autosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.high.qual.neutral.sites.vcf.gz"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "vqsr/sex_chromosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.high.qual.neutral.sites.noPARsXTR.vcf.gz"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "hard_filter/autosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.vcf.gz"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "hard_filter/sex_chromosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.noPARsXTR.vcf.gz"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "stats_pre_filter/autosomes", "{chrs}.gatk.called.raw.{vcf_options_auto}.high.qual.neutral.sites", "{chrs}.gatk.called.raw.{vcf_options_auto}.high.qual.neutral.sites_num_sites.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "stats_pre_filter/sex_chromosomes", "{chrs}.gatk.called.raw.{vcf_options_x}.high.qual.neutral.sites.noPARsXTR", "{chrs}.gatk.called.raw.{vcf_options_x}.high.qual.neutral.noPARsXTR.sites_num_sites.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "stats_post_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.high.qual.neutral.sites", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.high.qual.neutral.sites_num_sites.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "stats_post_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.high.qual.neutral.sites.noPARsXTR", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.high.qual.neutral.sites.noPARsXTR_num_sites.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "stats_post_hard_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.high.qual.neutral.sites", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.high.qual.neutral.sites_num_sites.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "stats_post_hard_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.high.qual.neutral.sites.noPARsXTR", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.high.qual.neutral.sites.noPARsXTR_num_sites.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "stats_pre_filter_array_sites/autosomes", "{chrs}.gatk.called.raw.{vcf_options_auto}.array.sites.high.qual.neutral.sites", "{chrs}.gatk.called.raw.{vcf_options_auto}.array.sites.high.qual.neutral.sites_num_sites.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "stats_pre_filter_array_sites/sex_chromosomes", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites.high.qual.neutral.sites.noPARsXTR", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites.high.qual.neutral.noPARsXTR.sites_num_sites.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "stats_post_filter_array_sites/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.high.qual.neutral.sites", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.high.qual.neutral.sites_num_sites.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "stats_post_filter_array_sites/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.high.qual.neutral.sites.noPARsXTR", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.high.qual.neutral.sites.noPARsXTR_num_sites.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "stats_post_hard_filter_array_sites/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites_num_sites.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "stats_post_hard_filter_array_sites/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.noPARsXTR", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.noPARsXTR_num_sites.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites_filters/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.high.qual.neutral.sites.chr.pos.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites_filters/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.high.qual.neutral.sites.chr.pos.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites_filters/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.high.qual.neutral.sites.noPARsXTR.chr.pos.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites_filters/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.high.qual.neutral.sites.noPARsXTR.chr.pos.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites_filters/autosomes", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.high.qual.neutral.sites.venn.diagram.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites_filters/sex_chromosomes", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.high.qual.neutral.sites.noPARsXTR.venn.diagram.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.high.qual.neutral.sites.chr.pos.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.chr.pos.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.high.qual.neutral.sites.noPARsXTR.chr.pos.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.noPARsXTR.chr.pos.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/autosomes", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.array.sites.high.qual.neutral.sites.venn.diagram.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/sex_chromosomes", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.array.sites.high.qual.neutral.sites.noPARsXTR.venn.diagram.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/pre_filter_region_filters/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.high.qual.neutral.sites.SFS.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "results/pre_filter_region_filters/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.high.qual.neutral.sites.noPARsXTR.SFS.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "results/vqsr_region_filters/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.high.qual.neutral.sites.SFS.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "results/vqsr_region_filters/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.high.qual.neutral.sites.noPARsXTR.SFS.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "results/hard_filter_region_filters/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.high.qual.neutral.sites.SFS.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/hard_filter_region_filters/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.high.qual.neutral.sites.noPARsXTR.SFS.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/merged/all_variable_sites_region_filters/sfs/autosomes", "{chrs}", "{chrs}.{vcf_options_auto}.raw.vqsr.{filtering_options}.high.qual.neutral.sites.SFS.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/merged/all_variable_sites_region_filters/sfs/sex_chromosomes", "{chrs}", "{chrs}.{vcf_options_x}.raw.vqsr.{filtering_options}.high.qual.neutral.sites.noPARs.noXTR.SFS.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/pre_filter_array_sites_region_filters/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.array.sites.high.qual.neutral.sites.SFS.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "results/pre_filter_array_sites_region_filters/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites.high.qual.neutral.sites.noPARsXTR.SFS.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "results/vqsr_array_sites_region_filters/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.high.qual.neutral.sites.SFS.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"]),
        expand(os.path.join(config["out_dir"], "results/vqsr_array_sites_region_filters/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.high.qual.neutral.sites.noPARsXTR.SFS.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"]),
        expand(os.path.join(config["out_dir"], "results/hard_filter_array_sites_region_filters/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.SFS.txt"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/hard_filter_array_sites_region_filters/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.noPARsXTR.SFS.txt"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/merged/array_sites_region_filters/sfs/autosomes", "{chrs}", "{chrs}.{vcf_options_auto}.raw.vqsr.{filtering_options}.array.sites.high.qual.neutral.sites.SFS.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/merged/array_sites_region_filters/sfs/sex_chromosomes", "{chrs}", "{chrs}.{vcf_options_x}.raw.vqsr.{filtering_options}.array.sites.high.qual.neutral.sites.noPARsXTR.SFS.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites/autosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.DP.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites/autosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.QD.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites/sex_chromosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.DP.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites/sex_chromosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.QD.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites_array_sites/autosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.array.sites.DP.png"),  chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites_array_sites/autosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.array.sites.QD.png"),  chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites_array_sites/sex_chromosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.array.sites.DP.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites_array_sites/sex_chromosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.array.sites.QD.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites_filters/autosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.high.qual.neutral.sites.DP.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites_filters/autosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.high.qual.neutral.sites.QD.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites_filters/sex_chromosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.high.qual.neutral.sites.noPARsXTR.DP.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites_filters/sex_chromosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.high.qual.neutral.sites.noPARsXTR.QD.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/autosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.array.sites.high.qual.neutral.sites.DP.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/autosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.array.sites.high.qual.neutral.sites.QD.png"), chrs=config["autos"], vcf_options_auto=config["vcf_options_auto"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/sex_chromosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.array.sites.high.qual.neutral.sites.noPARsXTR.DP.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"]),
        expand(os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/sex_chromosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.array.sites.high.qual.neutral.sites.noPARsXTR.QD.png"), chrs=config["x"], vcf_options_x=config["vcf_options_x"], filtering_options=config["filtering_options"])


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
# Venn diagram of overlapping sites - all called variants VQSR and hard filter #
# 1. Make a list per vcf with chr_pos
# Autosomes
rule chrPosListsVQSRAutos:
    input:
        os.path.join(config["out_dir"], "vqsr/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.vcf.gz")
    output:
        os.path.join(config["out_dir"], "results/overlap_sites/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.chr.pos.txt")
    shell:
        """
        zcat {input} | grep -v '#' | awk '{{ print $1'_'$2 }}' >  {output}
        """

rule chrPosListsHardFilterAutos:
    input:
        os.path.join(config["out_dir"], "hard_filter/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.vcf.gz")
    output:
        os.path.join(config["out_dir"], "results/overlap_sites/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.chr.pos.txt")
    shell:
        """
        zcat {input} | grep -v '#' | awk '{{ print $1'_'$2 }}' >  {output}
        """

# X chromosome
rule chrPosListsVQSRXchr:
    input:
        os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.vcf.gz")
    output:
        os.path.join(config["out_dir"], "results/overlap_sites/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.chr.pos.txt")
    shell:
        """
        zcat {input} | grep -v '#' | awk '{{ print $1'_'$2 }}' >  {output}
        """

rule chrPosListsHardFilterXchr:
    input:
        os.path.join(config["out_dir"], "hard_filter/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.vcf.gz")
    output:
        os.path.join(config["out_dir"], "results/overlap_sites/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.chr.pos.txt")
    shell:
        """
        zcat {input} | grep -v '#' | awk '{{ print $1'_'$2 }}' >  {output}
        """

# 2. Plot results
# Autosomes
rule plotOverlapAutos:
    input:
        vqsr = os.path.join(config["out_dir"], "results/overlap_sites/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.chr.pos.txt"),
        hardf = os.path.join(config["out_dir"], "results/overlap_sites/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.chr.pos.txt")
    params:
        srpt = config["scripts_dir"]
    output:
        os.path.join(config["out_dir"], "results/overlap_sites/autosomes", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.venn.diagram.png")
    shell:
        "Rscript {params.srpt}overlappingSites.R {input.vqsr} {input.hardf} {output}"

# X chromosome
rule plotOverlapXchr:
    input:
        vqsr = os.path.join(config["out_dir"], "results/overlap_sites/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.chr.pos.txt"),
        hardf = os.path.join(config["out_dir"], "results/overlap_sites/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.chr.pos.txt")
    params:
        srpt = config["scripts_dir"]
    output:
        os.path.join(config["out_dir"], "results/overlap_sites/sex_chromosomes", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.venn.diagram.png")
    shell:
        "Rscript {params.srpt}overlappingSites.R {input.vqsr} {input.hardf} {output}"

#-------------------------------------------------------------------------------
# Annotation plots for overlapping sites and sites uniquely retained post VQSR
#                                                           and hard filtering #
#-------------------------------------------------------------------------------
# All sites - no additional filters
# Autosomes #
rule densityPlotsAutos:
    input:
        vqsr = os.path.join(config["out_dir"], "results/overlap_sites/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.chr.pos.txt"),
        hardf = os.path.join(config["out_dir"], "results/overlap_sites/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.chr.pos.txt"),
        vqsrvcf = os.path.join(config["out_dir"], "vqsr/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.vcf.gz"),
        hardfvcf = os.path.join(config["out_dir"], "hard_filter/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.vcf.gz")
    params:
        srpt = config["scripts_dir"],
        chrm = "{chrs}",
        out1 = os.path.join(config["out_dir"], "results/overlap_sites/autosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_auto}.vqsr.unique.txt"),
        out2 = os.path.join(config["out_dir"], "results/overlap_sites/autosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_auto}.{filtering_options}.unique.txt"),
        out3 = os.path.join(config["out_dir"], "results/overlap_sites/autosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.overlap.txt"),
        vcf1 = os.path.join(config["out_dir"], "results/overlap_sites/autosomes", "{chrs}", "density_plots", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.unique.vcf.gz"),
        vcf2 = os.path.join(config["out_dir"], "results/overlap_sites/autosomes", "{chrs}", "density_plots", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.unique.vcf.gz"),
        vcf3 = os.path.join(config["out_dir"], "results/overlap_sites/autosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.sv.overlap.vcf.gz"),
        annotation1 = os.path.join(config["out_dir"], "results/overlap_sites/autosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.unique_annotations.txt"),
        annotation2 = os.path.join(config["out_dir"], "results/overlap_sites/autosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_auto}.{filtering_options}.unique_annotations.txt"),
        annotation3 = os.path.join(config["out_dir"], "results/overlap_sites/autosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.overlap_annotations.txt")
    output:
        outDP = os.path.join(config["out_dir"], "results/overlap_sites/autosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.DP.png"),
        outQD = os.path.join(config["out_dir"], "results/overlap_sites/autosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.QD.png")
    shell:
        """
        awk 'FNR==NR {{a[$0]++; next}} !a[$0]' {input.hardf} {input.vqsr} | sed -e 's/{params.chrm}/{params.chrm}\t/g' > {params.out1};
        awk 'FNR==NR {{a[$0]++; next}} !a[$0]' {input.vqsr} {input.hardf} | sed -e 's/{params.chrm}/{params.chrm}\t/g' > {params.out2};
        awk 'NR==FNR {{ lines[$0]=1; next }} $0 in lines' {input.hardf} {input.vqsr} | sed -e 's/{params.chrm}/{params.chrm}\t/g' > {params.out3};
        bcftools view -O z  {input.vqsrvcf} -T {params.out1} -o {params.vcf1};
        bcftools view -O z  {input.hardfvcf} -T {params.out2} -o {params.vcf2};
        bcftools view -O z  {input.vqsrvcf} -T {params.out3} -o {params.vcf3};
        python {params.srpt}extract_stats_from_vcf.py AN QD MQ DP --vcf {params.vcf1} --outfile {params.annotation1};
        python {params.srpt}extract_stats_from_vcf.py AN QD MQ DP --vcf {params.vcf2} --outfile {params.annotation2};
        python {params.srpt}extract_stats_from_vcf.py AN QD MQ DP --vcf {params.vcf3} --outfile {params.annotation3};
        Rscript {params.srpt}annotation_density_plots.R {params.annotation1} {params.annotation2} {params.annotation3} {output.outDP} {output.outQD}
        """

# X chromosome #
rule densityPlotsXchr:
    input:
        vqsr = os.path.join(config["out_dir"], "results/overlap_sites/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.chr.pos.txt"),
        hardf = os.path.join(config["out_dir"], "results/overlap_sites/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.chr.pos.txt"),
        vqsrvcf = os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.vcf.gz"),
        hardfvcf = os.path.join(config["out_dir"], "hard_filter/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.vcf.gz")
    params:
        srpt = config["scripts_dir"],
        chrm = "{chrs}",
        out1 = os.path.join(config["out_dir"], "results/overlap_sites/sex_chromosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_x}.vqsr.unique.txt"),
        out2 = os.path.join(config["out_dir"], "results/overlap_sites/sex_chromosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_x}.{filtering_options}.unique.txt"),
        out3 = os.path.join(config["out_dir"], "results/overlap_sites/sex_chromosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.overlap.txt"),
        vcf1 = os.path.join(config["out_dir"], "results/overlap_sites/sex_chromosomes", "{chrs}", "density_plots", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.unique.vcf.gz"),
        vcf2 = os.path.join(config["out_dir"], "results/overlap_sites/sex_chromosomes", "{chrs}", "density_plots", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.unique.vcf.gz"),
        vcf3 = os.path.join(config["out_dir"], "results/overlap_sites/sex_chromosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.sv.overlap.vcf.gz"),
        annotation1 = os.path.join(config["out_dir"], "results/overlap_sites/sex_chromosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.unique_annotations.txt"),
        annotation2 = os.path.join(config["out_dir"], "results/overlap_sites/sex_chromosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_x}.{filtering_options}.unique_annotations.txt"),
        annotation3 = os.path.join(config["out_dir"], "results/overlap_sites/sex_chromosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.overlap_annotations.txt")
    output:
        outDP = os.path.join(config["out_dir"], "results/overlap_sites/sex_chromosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.DP.png"),
        outQD = os.path.join(config["out_dir"], "results/overlap_sites/sex_chromosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.QD.png")
    shell:
        """
        awk 'FNR==NR {{a[$0]++; next}} !a[$0]' {input.hardf} {input.vqsr} | sed -e 's/{params.chrm}/{params.chrm}\t/g' > {params.out1};
        awk 'FNR==NR {{a[$0]++; next}} !a[$0]' {input.vqsr} {input.hardf} | sed -e 's/{params.chrm}/{params.chrm}\t/g' > {params.out2};
        awk 'NR==FNR {{ lines[$0]=1; next }} $0 in lines' {input.hardf} {input.vqsr} | sed -e 's/{params.chrm}/{params.chrm}\t/g' > {params.out3};
        bcftools view -O z  {input.vqsrvcf} -T {params.out1} -o {params.vcf1};
        bcftools view -O z  {input.hardfvcf} -T {params.out2} -o {params.vcf2};
        bcftools view -O z  {input.vqsrvcf} -T {params.out3} -o {params.vcf3};
        python {params.srpt}extract_stats_from_vcf.py AN QD MQ DP --vcf {params.vcf1} --outfile {params.annotation1};
        python {params.srpt}extract_stats_from_vcf.py AN QD MQ DP --vcf {params.vcf2} --outfile {params.annotation2};
        python {params.srpt}extract_stats_from_vcf.py AN QD MQ DP --vcf {params.vcf3} --outfile {params.annotation3};
        Rscript {params.srpt}annotation_density_plots.R {params.annotation1} {params.annotation2} {params.annotation3} {output.outDP} {output.outQD}
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


#-------------------------------------------------------------------------------
# Get number of sites post array site extraction from VQSR and hard filtered
# VCFs.
# Raw vcfs array sites autosomes
rule countNumSitesRawArraySitesAutos:
    input:
        vcf = os.path.join(config["out_dir"], "raw/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.array.sites.vcf.gz")
    output:
        siteinfo = os.path.join(config["out_dir"], "stats_pre_filter_array_sites/autosomes", "{chrs}.gatk.called.raw.{vcf_options_auto}.array.sites", "{chrs}.gatk.called.raw.{vcf_options_auto}.array.sites_num_sites.txt")
    shell:
        """
        zcat {input.vcf} | grep -v '#' | wc -l > {output.siteinfo}
        """

# Raw vcfs array sites X chromosome
rule countNumSitesRawArraySitesXchr:
    input:
        vcf = os.path.join(config["out_dir"], "raw/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites.vcf.gz")
    output:
        siteinfo = os.path.join(config["out_dir"], "stats_pre_filter_array_sites/sex_chromosomes", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites_num_sites.txt")
    shell:
        """
        zcat {input.vcf} | grep -v '#' | wc -l > {output.siteinfo}
        """

# VQSR array sites autosomes
rule countNumSitesVQSRArraySitesAutos:
    input:
        vcf = os.path.join(config["out_dir"], "vqsr/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.vcf.gz")
    output:
        siteinfo = os.path.join(config["out_dir"], "stats_post_filter_array_sites/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites_num_sites.txt")
    shell:
        """
        zcat {input.vcf} | grep -v '#' | wc -l > {output.siteinfo}
        """

# VQSR array sites X chromosome
rule countNumSitesVQSRArraySitesXchr:
    input:
        vcf = os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.vcf.gz")
    output:
        siteinfo = os.path.join(config["out_dir"], "stats_post_filter_array_sites/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites_num_sites.txt")
    shell:
        """
        zcat {input.vcf} | grep -v '#' | wc -l > {output.siteinfo}
        """

# Hard filter array sites autosomes
rule countNumSitesHardFilterArraySitesAutos:
    input:
        vcf = os.path.join(config["out_dir"], "hard_filter/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.vcf.gz")
    output:
        siteinfo = os.path.join(config["out_dir"], "stats_post_hard_filter_array_sites/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites_num_sites.txt")
    shell:
        """
        zcat {input.vcf} | grep -v '#' | wc -l > {output.siteinfo}
        """

# Hard filter array sites X chromosome
rule countNumSitesHardFilterArraySitesXchr:
    input:
        vcf = os.path.join(config["out_dir"], "hard_filter/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.vcf.gz")
    output:
        siteinfo = os.path.join(config["out_dir"], "stats_post_hard_filter_array_sites/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites_num_sites.txt")
    shell:
        """
        zcat {input.vcf} | grep -v '#' | wc -l > {output.siteinfo}
        """


#-------------------------------------------------------------------------------
# Venn diagram of overlapping sites - array sites VQSR and hard filter #
# 1. Make a list per vcf with chr_pos
# Autosomes
rule chrPosListsVQSRArraySitesAutos:
    input:
        os.path.join(config["out_dir"], "vqsr/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.vcf.gz")
    output:
        os.path.join(config["out_dir"], "results/overlap_sites_array_sites/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.chr.pos.txt")
    shell:
        """
        zcat {input} | grep -v '#' | awk '{{ print $1'_'$2 }}' >  {output}
        """

rule chrPosListsHardFilterArraySitesAutos:
    input:
        os.path.join(config["out_dir"], "hard_filter/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.vcf.gz")
    output:
        os.path.join(config["out_dir"], "results/overlap_sites_array_sites/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.chr.pos.txt")
    shell:
        """
        zcat {input} | grep -v '#' | awk '{{ print $1'_'$2 }}' >  {output}
        """

# X chromosome
rule chrPosListsVQSRArraySitesXchr:
    input:
        os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.vcf.gz")
    output:
        os.path.join(config["out_dir"], "results/overlap_sites_array_sites/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.chr.pos.txt")
    shell:
        """
        zcat {input} | grep -v '#' | awk '{{ print $1'_'$2 }}' >  {output}
        """

rule chrPosListsHardFilterArraySitesXchr:
    input:
        os.path.join(config["out_dir"], "hard_filter/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.vcf.gz")
    output:
        os.path.join(config["out_dir"], "results/overlap_sites_array_sites/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.chr.pos.txt")
    shell:
        """
        zcat {input} | grep -v '#' | awk '{{ print $1'_'$2 }}' >  {output}
        """

# 2. Plot results
# Autosomes
rule plotOverlapArraySitesAutos:
    input:
        vqsr = os.path.join(config["out_dir"], "results/overlap_sites_array_sites/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.chr.pos.txt"),
        hardf = os.path.join(config["out_dir"], "results/overlap_sites_array_sites/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.chr.pos.txt")
    params:
        srpt = config["scripts_dir"]
    output:
        os.path.join(config["out_dir"], "results/overlap_sites_array_sites/autosomes", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.array.sites.venn.diagram.png")
    shell:
        "Rscript {params.srpt}overlappingSites.R {input.vqsr} {input.hardf} {output}"

# X chromosome
rule plotOverlapArraySitesXchr:
    input:
        vqsr = os.path.join(config["out_dir"], "results/overlap_sites_array_sites/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.chr.pos.txt"),
        hardf = os.path.join(config["out_dir"], "results/overlap_sites_array_sites/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.chr.pos.txt")
    params:
        srpt = config["scripts_dir"]
    output:
        os.path.join(config["out_dir"], "results/overlap_sites_array_sites/sex_chromosomes", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.array.sites.venn.diagram.png")
    shell:
        "Rscript {params.srpt}overlappingSites.R {input.vqsr} {input.hardf} {output}"


#-------------------------------------------------------------------------------
# Annotation plots for overlapping sites and sites uniquely retained post VQSR
#                                                           and hard filtering #
#-------------------------------------------------------------------------------
# Array sites - no additional filters
# Autosomes #
rule densityPlotsArraySitesAutos:
    input:
        vqsr = os.path.join(config["out_dir"], "results/overlap_sites_array_sites/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.chr.pos.txt"),
        hardf = os.path.join(config["out_dir"], "results/overlap_sites_array_sites/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.chr.pos.txt"),
        vqsrvcf = os.path.join(config["out_dir"], "vqsr/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.vcf.gz"),
        hardfvcf = os.path.join(config["out_dir"], "hard_filter/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.vcf.gz"),
        orgvqsrvcf = os.path.join(config["out_dir"], "vqsr/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.vcf.gz"),
        orghardfvcf = os.path.join(config["out_dir"], "hard_filter/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.vcf.gz"),
    params:
        srpt = config["scripts_dir"],
        chrm = "{chrs}",
        out1 = os.path.join(config["out_dir"], "results/overlap_sites_array_sites/autosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_auto}.vqsr.array.sites.unique.txt"),
        out2 = os.path.join(config["out_dir"], "results/overlap_sites_array_sites/autosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_auto}.{filtering_options}.array.sites.unique.txt"),
        out3 = os.path.join(config["out_dir"], "results/overlap_sites_array_sites/autosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.array.sites.overlap.txt"),
        vcf1 = os.path.join(config["out_dir"], "results/overlap_sites_array_sites/autosomes", "{chrs}", "density_plots", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.unique.vcf.gz"),
        vcf2 = os.path.join(config["out_dir"], "results/overlap_sites_array_sites/autosomes", "{chrs}", "density_plots", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.unique.vcf.gz"),
        vcf3 = os.path.join(config["out_dir"], "results/overlap_sites_array_sites/autosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.sv.array.sites.overlap.vcf.gz"),
        annotation1 = os.path.join(config["out_dir"], "results/overlap_sites_array_sites/autosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.array.sites.unique_annotations.txt"),
        annotation2 = os.path.join(config["out_dir"], "results/overlap_sites_array_sites/autosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_auto}.{filtering_options}.array.sites.unique_annotations.txt"),
        annotation3 = os.path.join(config["out_dir"], "results/overlap_sites_array_sites/autosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.array.sites.overlap_annotations.txt")
    output:
        outDP = os.path.join(config["out_dir"], "results/overlap_sites_array_sites/autosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.array.sites.DP.png"),
        outQD = os.path.join(config["out_dir"], "results/overlap_sites_array_sites/autosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.array.sites.QD.png")
    shell:
        """
        awk 'FNR==NR {{a[$0]++; next}} !a[$0]' {input.hardf} {input.vqsr} | sed -e 's/{params.chrm}/{params.chrm}\t/g' > {params.out1};
        awk 'FNR==NR {{a[$0]++; next}} !a[$0]' {input.vqsr} {input.hardf} | sed -e 's/{params.chrm}/{params.chrm}\t/g' > {params.out2};
        awk 'NR==FNR {{ lines[$0]=1; next }} $0 in lines' {input.hardf} {input.vqsr} | sed -e 's/{params.chrm}/{params.chrm}\t/g' > {params.out3};
        bcftools view -O z  {input.orgvqsrvcf} -T {params.out1} -o {params.vcf1};
        bcftools view -O z  {input.orghardfvcf} -T {params.out2} -o {params.vcf2};
        bcftools view -O z  {input.orgvqsrvcf} -T {params.out3} -o {params.vcf3};
        python {params.srpt}extract_stats_from_vcf.py AN QD MQ DP --vcf {params.vcf1} --outfile {params.annotation1};
        python {params.srpt}extract_stats_from_vcf.py AN QD MQ DP --vcf {params.vcf2} --outfile {params.annotation2};
        python {params.srpt}extract_stats_from_vcf.py AN QD MQ DP --vcf {params.vcf3} --outfile {params.annotation3};
        Rscript {params.srpt}annotation_density_plots.R {params.annotation1} {params.annotation2} {params.annotation3} {output.outDP} {output.outQD}
        """

# X chromosome #
rule densityPlotsArraySitesXchr:
    input:
        vqsr = os.path.join(config["out_dir"], "results/overlap_sites_array_sites/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.chr.pos.txt"),
        hardf = os.path.join(config["out_dir"], "results/overlap_sites_array_sites/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.chr.pos.txt"),
        vqsrvcf = os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.vcf.gz"),
        hardfvcf = os.path.join(config["out_dir"], "hard_filter/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.vcf.gz"),
        orgvqsrvcf = os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.vcf.gz"),
        orghardfvcf = os.path.join(config["out_dir"], "hard_filter/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.vcf.gz")
    params:
        srpt = config["scripts_dir"],
        chrm = "{chrs}",
        out1 = os.path.join(config["out_dir"], "results/overlap_sites_array_sites/sex_chromosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_x}.vqsr.array.sites.unique.txt"),
        out2 = os.path.join(config["out_dir"], "results/overlap_sites_array_sites/sex_chromosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_x}.{filtering_options}.array.sites.unique.txt"),
        out3 = os.path.join(config["out_dir"], "results/overlap_sites_array_sites/sex_chromosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.array.sites.overlap.txt"),
        vcf1 = os.path.join(config["out_dir"], "results/overlap_sites_array_sites/sex_chromosomes", "{chrs}", "density_plots", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.unique.vcf.gz"),
        vcf2 = os.path.join(config["out_dir"], "results/overlap_sites_array_sites/sex_chromosomes", "{chrs}", "density_plots", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.unique.vcf.gz"),
        vcf3 = os.path.join(config["out_dir"], "results/overlap_sites_array_sites/sex_chromosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.sv.array.sites.overlap.vcf.gz"),
        annotation1 = os.path.join(config["out_dir"], "results/overlap_sites_array_sites/sex_chromosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.array.sites.unique_annotations.txt"),
        annotation2 = os.path.join(config["out_dir"], "results/overlap_sites_array_sites/sex_chromosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_x}.{filtering_options}.array.sites.unique_annotations.txt"),
        annotation3 = os.path.join(config["out_dir"], "results/overlap_sites_array_sites/sex_chromosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.array.sites.overlap_annotations.txt")
    output:
        outDP = os.path.join(config["out_dir"], "results/overlap_sites_array_sites/sex_chromosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.array.sites.DP.png"),
        outQD = os.path.join(config["out_dir"], "results/overlap_sites_array_sites/sex_chromosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.array.sites.QD.png")
    shell:
        """
        awk 'FNR==NR {{a[$0]++; next}} !a[$0]' {input.hardf} {input.vqsr} | sed -e 's/{params.chrm}/{params.chrm}\t/g' > {params.out1};
        awk 'FNR==NR {{a[$0]++; next}} !a[$0]' {input.vqsr} {input.hardf} | sed -e 's/{params.chrm}/{params.chrm}\t/g' > {params.out2};
        awk 'NR==FNR {{ lines[$0]=1; next }} $0 in lines' {input.hardf} {input.vqsr} | sed -e 's/{params.chrm}/{params.chrm}\t/g' > {params.out3};
        bcftools view -O z  {input.orgvqsrvcf} -T {params.out1} -o {params.vcf1};
        bcftools view -O z  {input.orghardfvcf} -T {params.out2} -o {params.vcf2};
        bcftools view -O z  {input.orgvqsrvcf} -T {params.out3} -o {params.vcf3};
        python {params.srpt}extract_stats_from_vcf.py AN QD MQ DP --vcf {params.vcf1} --outfile {params.annotation1};
        python {params.srpt}extract_stats_from_vcf.py AN QD MQ DP --vcf {params.vcf2} --outfile {params.annotation2};
        python {params.srpt}extract_stats_from_vcf.py AN QD MQ DP --vcf {params.vcf3} --outfile {params.annotation3};
        Rscript {params.srpt}annotation_density_plots.R {params.annotation1} {params.annotation2} {params.annotation3} {output.outDP} {output.outQD}
        """


#-------------------------------------------------------------------------------
# SITES RESTRICTED TO THE ARRAY #
# Generate SFS and plot results #
# Pre filter #
# Autosomes
rule SFSPreFilterArraySitesAutos:
    input:
        vcf = os.path.join(config["out_dir"], "raw/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.array.sites.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"],
        plotscriptdir = config["popgen_plotting_scripts_dir"],
        sfstxt = os.path.join(config["out_dir"], "results/pre_filter/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.array.sites.SFS.txt")
    output:
        sfstxt = os.path.join(config["out_dir"], "results/pre_filter/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.array.sites.SFS.txt"),
        png = os.path.join(config["out_dir"], "results/pre_filter/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.array.sites.SFS.png")
    shell:
        """
        python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {params.sfstxt} --ploidy diploid;
        Rscript {params.plotscriptdir}plot_sfs.R {params.sfstxt} {output.png}
        """

# X chromosome
rule SFSPreFilterArraySitesXchr:
    input:
        vcf = os.path.join(config["out_dir"], "raw/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"],
        plotscriptdir = config["popgen_plotting_scripts_dir"],
        sfstxt = os.path.join(config["out_dir"], "results/pre_filter/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites.SFS.txt")
    output:
        sfstxt = os.path.join(config["out_dir"], "results/pre_filter/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites.SFS.txt"),
        png = os.path.join(config["out_dir"], "results/pre_filter/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites.SFS.png")
    run:
        if wildcards.vcf_options_x == "haploid":
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {params.sfstxt} --ploidy haploid")
            shell("Rscript {params.plotscriptdir}plot_sfs.R {params.sfstxt} {output.png}")
        else:
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {params.sfstxt} --ploidy diploid")
            shell("Rscript {params.plotscriptdir}plot_sfs.R {params.sfstxt} {output.png}")

# Post VQSR #
# Autosomes
rule SFSVQSRArraySitesAutos:
    input:
        vcf = os.path.join(config["out_dir"], "vqsr/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"],
        plotscriptdir = config["popgen_plotting_scripts_dir"],
        sfstxt = os.path.join(config["out_dir"], "results/vqsr/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.SFS.txt")
    output:
        sfstxt = os.path.join(config["out_dir"], "results/vqsr/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.SFS.txt"),
        png = os.path.join(config["out_dir"], "results/vqsr/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.SFS.png")
    shell:
        """
        python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {params.sfstxt} --ploidy diploid;
        Rscript {params.plotscriptdir}plot_sfs.R {params.sfstxt} {output.png}
        """

# X chromosome
rule SFSVQSRArraySitesXchr:
    input:
        vcf = os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"],
        plotscriptdir = config["popgen_plotting_scripts_dir"],
        sfstxt = os.path.join(config["out_dir"], "results/vqsr/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.SFS.txt")
    output:
        sfstxt = os.path.join(config["out_dir"], "results/vqsr/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.SFS.txt"),
        png = os.path.join(config["out_dir"], "results/vqsr/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.SFS.png")
    run:
        if wildcards.vcf_options_x == "haploid":
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {params.sfstxt} --ploidy haploid")
            shell("Rscript {params.plotscriptdir}plot_sfs.R {params.sfstxt} {output.png}")
        else:
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {params.sfstxt} --ploidy diploid")
            shell("Rscript {params.plotscriptdir}plot_sfs.R {params.sfstxt} {output.png}")

# Post Hard Filter #
# Autosomes
rule SFSHardFilterArraySitesAutos:
    input:
        vcf = os.path.join(config["out_dir"], "hard_filter/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"],
        plotscriptdir = config["popgen_plotting_scripts_dir"],
        sfstxt = os.path.join(config["out_dir"], "results/hard_filter/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.SFS.txt")
    output:
        sfstxt = os.path.join(config["out_dir"], "results/hard_filter/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.SFS.txt"),
        png = os.path.join(config["out_dir"], "results/hard_filter/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.SFS.png")
    shell:
        """
        python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {params.sfstxt} --ploidy diploid;
        Rscript {params.plotscriptdir}plot_sfs.R {params.sfstxt} {output.png}
        """

# X chromosome
rule SFSHardFilterArraySitesXchr:
    input:
        vcf = os.path.join(config["out_dir"], "hard_filter/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"],
        plotscriptdir = config["popgen_plotting_scripts_dir"],
        sfstxt = os.path.join(config["out_dir"], "results/hard_filter/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.SFS.txt")
    output:
        sfstxt = os.path.join(config["out_dir"], "results/hard_filter/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.SFS.txt"),
        png = os.path.join(config["out_dir"], "results/hard_filter/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.SFS.png")
    run:
        if wildcards.vcf_options_x == "haploid":
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {params.sfstxt} --ploidy haploid")
            shell("Rscript {params.plotscriptdir}plot_sfs.R {params.sfstxt} {output}")
        else:
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {params.sfstxt} --ploidy diploid")
            shell("Rscript {params.plotscriptdir}plot_sfs.R {params.sfstxt} {output.png}")


#-------------------------------------------------------------------------------
# SITES RESTRICTED TO THE ARRAY #
# Plot SFSs in groups: raw, vqsr, and hard filter SFSs on 1 plot
# Autosomes
# NOTE: For R script order of input is important - Raw first, than VQSR, than hard filter
rule SFSAllArraySitesAutos:
    input:
        rawfn = os.path.join(config["out_dir"], "results/pre_filter/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.array.sites.SFS.txt"),
        vqsrfn = os.path.join(config["out_dir"], "results/vqsr/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.SFS.txt"),
        hfilterfn = os.path.join(config["out_dir"], "results/hard_filter/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.SFS.txt")
    params:
        scrpt = config["scripts_dir"]
    output:
        os.path.join(config["out_dir"], "results/merged/array_sites/sfs/autosomes", "{chrs}", "{chrs}.{vcf_options_auto}.raw.vqsr.{filtering_options}.SFS.png")
    shell:
        "Rscript {params.scrpt}sfs_plotting.R {input.rawfn} {input.vqsrfn} {input.hfilterfn} {output}"

# X chromosome
rule SFSAllArraySitesXchr:
    input:
        rawfn = os.path.join(config["out_dir"], "results/pre_filter/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites.SFS.txt"),
        vqsrfn = os.path.join(config["out_dir"], "results/vqsr/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.SFS.txt"),
        hfilterfn = os.path.join(config["out_dir"], "results/hard_filter/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.SFS.txt")
    params:
        scrpt = config["scripts_dir"]
    output:
        os.path.join(config["out_dir"], "results/merged/array_sites/sfs/sex_chromosomes", "{chrs}", "{chrs}.{vcf_options_x}.raw.vqsr.{filtering_options}.SFS.png")
    shell:
        "Rscript {params.scrpt}sfs_plotting.R {input.rawfn} {input.vqsrfn} {input.hfilterfn} {output}"


#-------------------------------------------------------------------------------
# ALL VARIANT SITES #
# Generate SFS and plot results #
# Pre filter #
# Autosomes
rule SFSPreFilterAutos:
    input:
        vcf = os.path.join(config["in_vcf_dir"], "{chrs}.gatk.called.raw.{vcf_options_auto}.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"],
        plotscriptdir = config["popgen_plotting_scripts_dir"],
        sfstxt = os.path.join(config["out_dir"], "results/pre_filter/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.SFS.txt")
    output:
        sfstxt = os.path.join(config["out_dir"], "results/pre_filter/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.SFS.txt"),
        png = os.path.join(config["out_dir"], "results/pre_filter/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.SFS.png")
    shell:
        """
        python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {params.sfstxt} --ploidy diploid;
        Rscript {params.plotscriptdir}plot_sfs.R {params.sfstxt} {output.png}
        """

# X chromosome
rule SFSPreFilterXchr:
    input:
        vcf = os.path.join(config["in_vcf_dir"], "{chrs}.gatk.called.raw.{vcf_options_x}.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"],
        plotscriptdir = config["popgen_plotting_scripts_dir"],
        sfstxt = os.path.join(config["out_dir"], "results/pre_filter/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.SFS.txt")
    output:
        sfstxt = os.path.join(config["out_dir"], "results/pre_filter/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.SFS.txt"),
        png = os.path.join(config["out_dir"], "results/pre_filter/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.SFS.png")
    run:
        if wildcards.vcf_options_x == "haploid":
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {params.sfstxt} --ploidy haploid")
            shell("Rscript {params.plotscriptdir}plot_sfs.R {params.sfstxt} {output.png}")
        else:
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {params.sfstxt} --ploidy diploid")
            shell("Rscript {params.plotscriptdir}plot_sfs.R {params.sfstxt} {output.png}")

# Post VQSR #
# Autosomes
rule SFSVQSRAutos:
    input:
        vcf = os.path.join(config["out_dir"], "vqsr/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"],
        plotscriptdir = config["popgen_plotting_scripts_dir"],
        sfstxt = os.path.join(config["out_dir"], "results/vqsr/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.SFS.txt")
    output:
        sfstxt = os.path.join(config["out_dir"], "results/vqsr/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.SFS.txt"),
        png = os.path.join(config["out_dir"], "results/vqsr/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.SFS.png")
    shell:
        """
        python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {params.sfstxt} --ploidy diploid;
        Rscript {params.plotscriptdir}plot_sfs.R {params.sfstxt} {output.png}
        """

# X chromosome
rule SFSVQSRXchr:
    input:
        vcf = os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"],
        plotscriptdir = config["popgen_plotting_scripts_dir"],
        sfstxt = os.path.join(config["out_dir"], "results/vqsr/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.SFS.txt")
    output:
        sfstxt = os.path.join(config["out_dir"], "results/vqsr/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.SFS.txt"),
        png = os.path.join(config["out_dir"], "results/vqsr/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.SFS.png")
    run:
        if wildcards.vcf_options_x == "haploid":
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {params.sfstxt} --ploidy haploid")
            shell("Rscript {params.plotscriptdir}plot_sfs.R {params.sfstxt} {output.png}")
        else:
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {params.sfstxt} --ploidy diploid")
            shell("Rscript {params.plotscriptdir}plot_sfs.R {params.sfstxt} {output.png}")

# Post Hard Filter #
# Autosomes
rule SFSHardFilterAutos:
    input:
        vcf = os.path.join(config["out_dir"], "hard_filter/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"],
        plotscriptdir = config["popgen_plotting_scripts_dir"],
        sfstxt = os.path.join(config["out_dir"], "results/hard_filter/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.SFS.txt")
    output:
        sfstxt = os.path.join(config["out_dir"], "results/hard_filter/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.SFS.txt"),
        png = os.path.join(config["out_dir"], "results/hard_filter/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.SFS.png")
    shell:
        """
        python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {params.sfstxt} --ploidy diploid;
        Rscript {params.plotscriptdir}plot_sfs.R {params.sfstxt} {output.png}
        """

# X chromosome
rule SFSHardFilterXchr:
    input:
        vcf = os.path.join(config["out_dir"], "hard_filter/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"],
        plotscriptdir = config["popgen_plotting_scripts_dir"],
        sfstxt = os.path.join(config["out_dir"], "results/hard_filter/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.SFS.txt")
    output:
        sfstxt = os.path.join(config["out_dir"], "results/hard_filter/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.SFS.txt"),
        png = os.path.join(config["out_dir"], "results/hard_filter/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.SFS.png")
    run:
        if wildcards.vcf_options_x == "haploid":
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {params.sfstxt} --ploidy haploid")
            shell("Rscript {params.plotscriptdir}plot_sfs.R {params.sfstxt} {output.png}")
        else:
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {params.sfstxt} --ploidy diploid")
            shell("Rscript {params.plotscriptdir}plot_sfs.R {params.sfstxt} {output.png}")

#-------------------------------------------------------------------------------
# ALL VARIANT SITES #
# Plot SFSs in groups: raw, vqsr, and hard filter SFSs on 1 plot
# Autosomes
# NOTE: For R script order of input is important - Raw first, than VQSR, than hard filter
rule SFSAllAutos:
    input:
        rawfn = os.path.join(config["out_dir"], "results/pre_filter/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.SFS.txt"),
        vqsrfn = os.path.join(config["out_dir"], "results/vqsr/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.SFS.txt"),
        hfilterfn = os.path.join(config["out_dir"], "results/hard_filter/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.SFS.txt")
    params:
        scrpt = config["scripts_dir"]
    output:
        os.path.join(config["out_dir"], "results/merged/all_variable_sites/sfs/autosomes", "{chrs}", "{chrs}.{vcf_options_auto}.raw.vqsr.{filtering_options}.SFS.png")
    shell:
        "Rscript {params.scrpt}sfs_plotting.R {input.rawfn} {input.vqsrfn} {input.hfilterfn} {output}"

# X chromosome
rule SFSAllXchr:
    input:
        rawfn = os.path.join(config["out_dir"], "results/pre_filter/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.SFS.txt"),
        vqsrfn = os.path.join(config["out_dir"], "results/vqsr/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.SFS.txt"),
        hfilterfn = os.path.join(config["out_dir"], "results/hard_filter/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.SFS.txt")
    params:
        scrpt = config["scripts_dir"]
    output:
        os.path.join(config["out_dir"], "results/merged/all_variable_sites/sfs/sex_chromosomes", "{chrs}", "{chrs}.{vcf_options_x}.raw.vqsr.{filtering_options}.SFS.png")
    shell:
        "Rscript {params.scrpt}sfs_plotting.R {input.rawfn} {input.vqsrfn} {input.hfilterfn} {output}"

#-------------------------------------------------------------------------------
# ALL VARIANT SITES #
# Calculate pi per site

# Pre filter #
# Autosomes
rule piPreFilterAutos:
    input:
        vcf = os.path.join(config["in_vcf_dir"], "{chrs}.gatk.called.raw.{vcf_options_auto}.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"],
        beddirfn = os.path.join(config["hi_qual_neu_dir"], "{chrs}_GRCh38_high_quality_neutral_sites.short.bed")
    output:
        o1 = os.path.join(config["out_dir"], "results/pre_filter/pi/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.pi.txt"),
        o2 = os.path.join(config["out_dir"], "results/pre_filter/pi/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.pi.per.site.txt")
    shell:
        "python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --pi --pi_target --ploidy diploid --target_bed {params.beddirfn} --pi_target_out {output.o1} --pi_target_per_site_out {output.o2}"

# X chromosome
rule piPreFilterXchr:
    input:
        vcf = os.path.join(config["in_vcf_dir"], "{chrs}.gatk.called.raw.{vcf_options_x}.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"],
        beddirfn = os.path.join(config["hi_qual_neu_dir"], "{chrs}_GRCh38_high_quality_neutral_sites.short.bed"),
        parxtrbed = config["x_filter"],
        fvcf = os.path.join(config["out_dir"], "raw/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.noPARs.noXTR.vcf.gz")
    output:
        o1 = os.path.join(config["out_dir"], "results/pre_filter/pi/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.noPARs.noXTR.pi.txt"),
        o2 = os.path.join(config["out_dir"], "results/pre_filter/pi/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.noPARs.noXTR.pi.per.site.txt"),
        fvcf = os.path.join(config["out_dir"], "raw/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.noPARs.noXTR.vcf.gz")
    run:
        if wildcards.vcf_options_x == "haploid":
            shell("subtractBed -header -a {input.vcf} -b {params.parxtrbed} | bgzip -c > {params.fvcf}")
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {params.fvcf} --pi --pi_target --ploidy haploid --target_bed {params.beddirfn} --pi_target_out {output.o1} --pi_target_per_site_out {output.o2}")
        else:
            shell("subtractBed -header -a {input.vcf} -b {params.parxtrbed} | bgzip -c > {params.fvcf}")
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {params.fvcf} --pi --pi_target --ploidy diploid --target_bed {params.beddirfn} --pi_target_out {output.o1} --pi_target_per_site_out {output.o2}")

# Post VQSR #
# Autosomes
rule piVQSRAutos:
    input:
        vcf = os.path.join(config["out_dir"], "vqsr/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"],
        beddirfn = os.path.join(config["hi_qual_neu_dir"], "{chrs}_GRCh38_high_quality_neutral_sites.short.bed")
    output:
        o1 = os.path.join(config["out_dir"], "results/vqsr/pi/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.pi.txt"),
        o2 = os.path.join(config["out_dir"], "results/vqsr/pi/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.pi.per.site.txt")
    shell:
        "python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --pi --pi_target --ploidy diploid --target_bed {params.beddirfn} --pi_target_out {output.o1} --pi_target_per_site_out {output.o2}"

# X chromosome
rule piVQSRXchr:
    input:
        vcf = os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"],
        beddirfn = os.path.join(config["hi_qual_neu_dir"], "{chrs}_GRCh38_high_quality_neutral_sites.short.bed"),
        parxtrbed = config["x_filter"],
        fvcf = os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.noPARs.noXTR.vcf.gz")
    output:
        o1 = os.path.join(config["out_dir"], "results/vqsr/pi/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.noPARs.noXTR.pi.txt"),
        o2 = os.path.join(config["out_dir"], "results/vqsr/pi/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.noPARs.noXTR.pi.per.site.txt"),
        fvcf = os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.noPARs.noXTR.vcf.gz")
    run:
        if wildcards.vcf_options_x == "haploid":
            shell("subtractBed -header -a {input.vcf} -b {params.parxtrbed} | bgzip -c > {params.fvcf}")
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {params.fvcf} --pi --pi_target --ploidy haploid --target_bed {params.beddirfn} --pi_target_out {output.o1} --pi_target_per_site_out {output.o2}")
        else:
            shell("subtractBed -header -a {input.vcf} -b {params.parxtrbed} | bgzip -c > {params.fvcf}")
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {params.fvcf} --pi --pi_target --ploidy diploid --target_bed {params.beddirfn} --pi_target_out {output.o1} --pi_target_per_site_out {output.o2}")

# Post Hard Filter #
# Autosomes
rule piHardFilterAutos:
    input:
        vcf = os.path.join(config["out_dir"], "hard_filter/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"],
        beddirfn = os.path.join(config["hi_qual_neu_dir"], "{chrs}_GRCh38_high_quality_neutral_sites.short.bed")
    output:
        o1 = os.path.join(config["out_dir"], "results/hard_filter/pi/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.pi.txt"),
        o2 = os.path.join(config["out_dir"], "results/hard_filter/pi/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.pi.per.site.txt")
    shell:
        "python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --pi --pi_target --ploidy diploid --target_bed {params.beddirfn} --pi_target_out {output.o1} --pi_target_per_site_out {output.o2}"

# X chromosome
rule piHardFilterXchr:
    input:
        vcf = os.path.join(config["out_dir"], "hard_filter/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"],
        beddirfn = os.path.join(config["hi_qual_neu_dir"], "{chrs}_GRCh38_high_quality_neutral_sites.short.bed"),
        parxtrbed = config["x_filter"],
        fvcf = os.path.join(config["out_dir"], "hard_filter/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.noPARs.noXTR.vcf.gz")
    output:
        o1 = os.path.join(config["out_dir"], "results/hard_filter/pi/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.noPARs.noXTR.pi.txt"),
        o2 = os.path.join(config["out_dir"], "results/hard_filter/pi/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.noPARs.noXTR.pi.per.site.txt"),
        fvcf = os.path.join(config["out_dir"], "hard_filter/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.noPARs.noXTR.vcf.gz")
    run:
        if wildcards.vcf_options_x == "haploid":
            shell("subtractBed -header -a {input.vcf} -b {params.parxtrbed} | bgzip -c > {params.fvcf}")
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {params.fvcf} --pi --pi_target --ploidy haploid --target_bed {params.beddirfn} --pi_target_out {output.o1} --pi_target_per_site_out {output.o2}")
        else:
            shell("subtractBed -header -a {input.vcf} -b {params.parxtrbed} | bgzip -c > {params.fvcf}")
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {params.fvcf} --pi --pi_target --ploidy diploid --target_bed {params.beddirfn} --pi_target_out {output.o1} --pi_target_per_site_out {output.o2}")


#-------------------------------------------------------------------------------
# SITES RESTRICTED TO THE ARRAY #
# Calculate pi per site
# Pre filter #
# Autosomes
rule piPreFilterArraySitesAutos:
    input:
        vcf = os.path.join(config["out_dir"], "raw/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.array.sites.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"],
        beddirfn = os.path.join(config["hi_qual_neu_dir"], "{chrs}_GRCh38_high_quality_neutral_sites.short.bed")
    output:
        o1 = os.path.join(config["out_dir"], "results/pre_filter/pi/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.array.sites.pi.txt"),
        o2 = os.path.join(config["out_dir"], "results/pre_filter/pi/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.array.sites.pi.per.site.txt")
    shell:
        "python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --pi --pi_target --ploidy diploid --target_bed {params.beddirfn} --pi_target_out {output.o1} --pi_target_per_site_out {output.o2}"

# X chromosome
rule piPreFilterArraySitesXchr:
    input:
        vcf = os.path.join(config["out_dir"], "raw/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"],
        beddirfn = os.path.join(config["hi_qual_neu_dir"], "{chrs}_GRCh38_high_quality_neutral_sites.short.bed"),
        parxtrbed = config["x_filter"],
        fvcf = os.path.join(config["out_dir"], "raw/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites.noPARs.noXTR.vcf.gz")
    output:
        o1 = os.path.join(config["out_dir"], "results/pre_filter/pi/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites.noPARs.noXTR.pi.txt"),
        o2 = os.path.join(config["out_dir"], "results/pre_filter/pi/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites.noPARs.noXTR.pi.per.site.txt"),
        fvcf = os.path.join(config["out_dir"], "raw/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites.noPARs.noXTR.vcf.gz")
    run:
        if wildcards.vcf_options_x == "haploid":
            shell("subtractBed -header -a {input.vcf} -b {params.parxtrbed} | bgzip -c > {params.fvcf}")
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {params.fvcf} --pi --pi_target --ploidy haploid --target_bed {params.beddirfn} --pi_target_out {output.o1} --pi_target_per_site_out {output.o2}")
        else:
            shell("subtractBed -header -a {input.vcf} -b {params.parxtrbed} | bgzip -c > {params.fvcf}")
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {params.fvcf} --pi --pi_target --ploidy diploid --target_bed {params.beddirfn} --pi_target_out {output.o1} --pi_target_per_site_out {output.o2}")

# Post VQSR #
# Autosomes
rule piVQSRArraySitesAutos:
    input:
        vcf = os.path.join(config["out_dir"], "vqsr/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"],
        beddirfn = os.path.join(config["hi_qual_neu_dir"], "{chrs}_GRCh38_high_quality_neutral_sites.short.bed")
    output:
        o1 = os.path.join(config["out_dir"], "results/vqsr/pi/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.pi.txt"),
        o2 = os.path.join(config["out_dir"], "results/vqsr/pi/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.pi.per.site.txt")
    shell:
        "python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --pi --pi_target --ploidy diploid --target_bed {params.beddirfn} --pi_target_out {output.o1} --pi_target_per_site_out {output.o2}"

# X chromosome
rule piVQSRArraySitesXchr:
    input:
        vcf = os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"],
        beddirfn = os.path.join(config["hi_qual_neu_dir"], "{chrs}_GRCh38_high_quality_neutral_sites.short.bed"),
        parxtrbed = config["x_filter"],
        fvcf = os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.noPARs.noXTR.vcf.gz")
    output:
        o1 = os.path.join(config["out_dir"], "results/vqsr/pi/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.noPARs.noXTR.pi.txt"),
        o2 = os.path.join(config["out_dir"], "results/vqsr/pi/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.noPARs.noXTR.pi.per.site.txt"),
        fvcf = os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.noPARs.noXTR.vcf.gz")
    run:
        if wildcards.vcf_options_x == "haploid":
            shell("subtractBed -header -a {input.vcf} -b {params.parxtrbed} | bgzip -c > {params.fvcf}")
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {params.fvcf} --pi --pi_target --ploidy haploid --target_bed {params.beddirfn} --pi_target_out {output.o1} --pi_target_per_site_out {output.o2}")
        else:
            shell("subtractBed -header -a {input.vcf} -b {params.parxtrbed} | bgzip -c > {params.fvcf}")
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {params.fvcf} --pi --pi_target --ploidy diploid --target_bed {params.beddirfn} --pi_target_out {output.o1} --pi_target_per_site_out {output.o2}")

# Post Hard Filter #
# Autosomes
rule piHardFilterArraySitesAutos:
    input:
        vcf = os.path.join(config["out_dir"], "hard_filter/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"],
        beddirfn = os.path.join(config["hi_qual_neu_dir"], "{chrs}_GRCh38_high_quality_neutral_sites.short.bed")
    output:
        o1 = os.path.join(config["out_dir"], "results/hard_filter/pi/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.pi.txt"),
        o2 = os.path.join(config["out_dir"], "results/hard_filter/pi/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.pi.per.site.txt")
    shell:
        "python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --pi --pi_target --ploidy diploid --target_bed {params.beddirfn} --pi_target_out {output.o1} --pi_target_per_site_out {output.o2}"

# X chromosome
rule piHardFilterArraySitesXchr:
    input:
        vcf = os.path.join(config["out_dir"], "hard_filter/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"],
        beddirfn = os.path.join(config["hi_qual_neu_dir"], "{chrs}_GRCh38_high_quality_neutral_sites.short.bed"),
        parxtrbed = config["x_filter"],
        fvcf = os.path.join(config["out_dir"], "hard_filter/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.noPARs.noXTR.vcf.gz")
    output:
        o1 = os.path.join(config["out_dir"], "results/hard_filter/pi/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.noPARs.noXTR.pi.txt"),
        o2 = os.path.join(config["out_dir"], "results/hard_filter/pi/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.noPARs.noXTR.pi.per.site.txt"),
        fvcf = os.path.join(config["out_dir"], "hard_filter/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.noPARs.noXTR.vcf.gz")
    run:
        if wildcards.vcf_options_x == "haploid":
            shell("subtractBed -header -a {input.vcf} -b {params.parxtrbed} | bgzip -c > {params.fvcf}")
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {params.fvcf} --pi --pi_target --ploidy haploid --target_bed {params.beddirfn} --pi_target_out {output.o1} --pi_target_per_site_out {output.o2}")
        else:
            shell("subtractBed -header -a {input.vcf} -b {params.parxtrbed} | bgzip -c > {params.fvcf}")
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {params.fvcf} --pi --pi_target --ploidy diploid --target_bed {params.beddirfn} --pi_target_out {output.o1} --pi_target_per_site_out {output.o2}")


#-------------------------------------------------------------------------------
# ALL VARIANT SITES - NEUTRAL, HIGH QUALITY NO PARS AND XTR #
# Get number of sites, venn diagrams, and SFS
# 1) Filtering - Extract neutral and high quality sites, remove XTR and PARs #
# Pre filter
rule regionFiltersAllAutos:
    input:
        os.path.join(config["in_vcf_dir"], "{chrs}.gatk.called.raw.{vcf_options_auto}.vcf.gz")
    params:
        hiqneubed = os.path.join(config["hi_qual_neu_dir"], "{chrs}_GRCh38_high_quality_neutral_sites.bed")
    output:
        os.path.join(config["out_dir"], "raw/autosomes/filters", "{chrs}.gatk.called.raw.{vcf_options_auto}.high.qual.neutral.sites.vcf.gz")
    shell:
        "bedtools intersect -header -a {input} -b {params.hiqneubed} | bgzip -c > {output}"

rule regionFiltersAllXchr:
    input:
        os.path.join(config["in_vcf_dir"], "{chrs}.gatk.called.raw.{vcf_options_x}.vcf.gz")
    params:
        hiqneubed = os.path.join(config["hi_qual_neu_dir"], "{chrs}_GRCh38_high_quality_neutral_sites.bed"),
        parsxtrbed = config["x_filter"],
        filter1 = os.path.join(config["out_dir"], "raw/sex_chromosomes/filters", "{chrs}.gatk.called.raw.{vcf_options_x}.high.qual.neutral.sites.vcf.gz")
    output:
        os.path.join(config["out_dir"], "raw/sex_chromosomes/filters", "{chrs}.gatk.called.raw.{vcf_options_x}.high.qual.neutral.sites.noPARsXTR.vcf.gz")
    shell:
        """
        bedtools intersect -header -a {input} -b {params.hiqneubed} | bgzip -c > {params.filter1};
        subtractBed -header -a {params.filter1} -b {params.parsxtrbed} | bgzip -c > {output}
        """

# VQSR
rule regionFiltersVQSRAutos:
    input:
        os.path.join(config["out_dir"], "vqsr/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.vcf.gz")
    params:
        hiqneubed = os.path.join(config["hi_qual_neu_dir"], "{chrs}_GRCh38_high_quality_neutral_sites.bed")
    output:
        os.path.join(config["out_dir"], "vqsr/autosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.high.qual.neutral.sites.vcf.gz")
    shell:
        "bedtools intersect -header -a {input} -b {params.hiqneubed} | bgzip -c > {output}"

rule regionFiltersVQSRXchr:
    input:
        os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.vcf.gz")
    params:
        hiqneubed = os.path.join(config["hi_qual_neu_dir"], "{chrs}_GRCh38_high_quality_neutral_sites.bed"),
        parsxtrbed = config["x_filter"],
        filter1 = os.path.join(config["out_dir"], "vqsr/sex_chromosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.high.qual.neutral.sites.vcf.gz")
    output:
        os.path.join(config["out_dir"], "vqsr/sex_chromosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.high.qual.neutral.sites.noPARsXTR.vcf.gz")
    shell:
        """
        bedtools intersect -header -a {input} -b {params.hiqneubed} | bgzip -c > {params.filter1};
        subtractBed -header -a {params.filter1} -b {params.parsxtrbed} | bgzip -c > {output}
        """

# Hard Filter
rule regionFiltersHardFilterAutos:
    input:
        os.path.join(config["out_dir"], "hard_filter/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.vcf.gz")
    params:
        hiqneubed = os.path.join(config["hi_qual_neu_dir"], "{chrs}_GRCh38_high_quality_neutral_sites.bed")
    output:
        os.path.join(config["out_dir"], "hard_filter/autosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.high.qual.neutral.sites.vcf.gz")
    shell:
        "bedtools intersect -header -a {input} -b {params.hiqneubed} | bgzip -c > {output}"

rule regionFiltersHardFilterXchr:
    input:
        os.path.join(config["out_dir"], "hard_filter/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.vcf.gz")
    params:
        hiqneubed = os.path.join(config["hi_qual_neu_dir"], "{chrs}_GRCh38_high_quality_neutral_sites.bed"),
        parsxtrbed = config["x_filter"],
        filter1 = os.path.join(config["out_dir"], "hard_filter/sex_chromosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.high.qual.neutral.sites.vcf.gz")
    output:
        os.path.join(config["out_dir"], "hard_filter/sex_chromosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.high.qual.neutral.sites.noPARsXTR.vcf.gz")
    shell:
        """
        bedtools intersect -header -a {input} -b {params.hiqneubed} | bgzip -c > {params.filter1};
        subtractBed -header -a {params.filter1} -b {params.parsxtrbed} | bgzip -c > {output}
        """

# 2) Count number of sites for table #
# Raw
rule countSitesregionFiltersAllAutos:
    input:
        os.path.join(config["out_dir"], "raw/autosomes/filters", "{chrs}.gatk.called.raw.{vcf_options_auto}.high.qual.neutral.sites.vcf.gz")
    output:
        os.path.join(config["out_dir"], "stats_pre_filter/autosomes", "{chrs}.gatk.called.raw.{vcf_options_auto}.high.qual.neutral.sites", "{chrs}.gatk.called.raw.{vcf_options_auto}.high.qual.neutral.sites_num_sites.txt")
    shell:
        """
        zcat {input} | grep -v '#' | wc -l > {output}
        """

rule countSitesregionFiltersAllXchr:
    input:
        os.path.join(config["out_dir"], "raw/sex_chromosomes/filters", "{chrs}.gatk.called.raw.{vcf_options_x}.high.qual.neutral.sites.noPARsXTR.vcf.gz")
    output:
        os.path.join(config["out_dir"], "stats_pre_filter/sex_chromosomes", "{chrs}.gatk.called.raw.{vcf_options_x}.high.qual.neutral.sites.noPARsXTR", "{chrs}.gatk.called.raw.{vcf_options_x}.high.qual.neutral.noPARsXTR.sites_num_sites.txt")
    shell:
        """
        zcat {input} | grep -v '#' | wc -l > {output}
        """

# VQSR
rule countSitesregionFiltersVQSRAutos:
    input:
        os.path.join(config["out_dir"], "vqsr/autosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.high.qual.neutral.sites.vcf.gz")
    output:
        os.path.join(config["out_dir"], "stats_post_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.high.qual.neutral.sites", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.high.qual.neutral.sites_num_sites.txt")
    shell:
        """
        zcat {input} | grep -v '#' | wc -l > {output}
        """

rule countSitesregionFiltersVQSRXchr:
    input:
        os.path.join(config["out_dir"], "vqsr/sex_chromosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.high.qual.neutral.sites.noPARsXTR.vcf.gz")
    output:
        os.path.join(config["out_dir"], "stats_post_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.high.qual.neutral.sites.noPARsXTR", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.high.qual.neutral.sites.noPARsXTR_num_sites.txt")
    shell:
        """
        zcat {input} | grep -v '#' | wc -l > {output}
        """

# Hard Filter
rule countSitesregionFiltersHardFilterAutos:
    input:
        os.path.join(config["out_dir"], "hard_filter/autosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.high.qual.neutral.sites.vcf.gz")
    output:
        os.path.join(config["out_dir"], "stats_post_hard_filter/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.high.qual.neutral.sites", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.high.qual.neutral.sites_num_sites.txt")
    shell:
        """
        zcat {input} | grep -v '#' | wc -l > {output}
        """

rule countSitesregionFiltersHardFilterXchr:
    input:
        os.path.join(config["out_dir"], "hard_filter/sex_chromosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.high.qual.neutral.sites.noPARsXTR.vcf.gz")
    output:
        os.path.join(config["out_dir"], "stats_post_hard_filter/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.high.qual.neutral.sites.noPARsXTR", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.high.qual.neutral.sites.noPARsXTR_num_sites.txt")
    shell:
        """
        zcat {input} | grep -v '#' | wc -l > {output}
        """

# 3) Make venn diagrams comparing VQSR and hard filtering #
# 3.1. Make a list per vcf with chr_pos
# Autosomes
rule chrPosListsVQSRRegionFilterAutos:
    input:
        os.path.join(config["out_dir"], "vqsr/autosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.high.qual.neutral.sites.vcf.gz")
    output:
        os.path.join(config["out_dir"], "results/overlap_sites_filters/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.high.qual.neutral.sites.chr.pos.txt")
    shell:
        """
        zcat {input} | grep -v '#' | awk '{{ print $1'_'$2 }}' >  {output}
        """

rule chrPosListsHardFilterRegionFilterAutos:
    input:
        os.path.join(config["out_dir"], "hard_filter/autosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.high.qual.neutral.sites.vcf.gz")
    output:
        os.path.join(config["out_dir"], "results/overlap_sites_filters/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.high.qual.neutral.sites.chr.pos.txt")
    shell:
        """
        zcat {input} | grep -v '#' | awk '{{ print $1'_'$2 }}' >  {output}
        """

# X chromosome
rule chrPosListsVQSRRegionFilterXchr:
    input:
        os.path.join(config["out_dir"], "vqsr/sex_chromosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.high.qual.neutral.sites.noPARsXTR.vcf.gz")
    output:
        os.path.join(config["out_dir"], "results/overlap_sites_filters/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.high.qual.neutral.sites.noPARsXTR.chr.pos.txt")
    shell:
        """
        zcat {input} | grep -v '#' | awk '{{ print $1'_'$2 }}' >  {output}
        """

rule chrPosListsHardFilterRegionFilterXchr:
    input:
        os.path.join(config["out_dir"], "hard_filter/sex_chromosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.high.qual.neutral.sites.noPARsXTR.vcf.gz")
    output:
        os.path.join(config["out_dir"], "results/overlap_sites_filters/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.high.qual.neutral.sites.noPARsXTR.chr.pos.txt")
    shell:
        """
        zcat {input} | grep -v '#' | awk '{{ print $1'_'$2 }}' >  {output}
        """

# 3.2. Plot results
# Autosomes
rule plotOverlapRegionFilterAutos:
    input:
        vqsr = os.path.join(config["out_dir"], "results/overlap_sites_filters/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.high.qual.neutral.sites.chr.pos.txt"),
        hardf = os.path.join(config["out_dir"], "results/overlap_sites_filters/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.high.qual.neutral.sites.chr.pos.txt")
    params:
        srpt = config["scripts_dir"]
    output:
        os.path.join(config["out_dir"], "results/overlap_sites_filters/autosomes", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.high.qual.neutral.sites.venn.diagram.png")
    shell:
        "Rscript {params.srpt}overlappingSites.R {input.vqsr} {input.hardf} {output}"

# X chromosome
rule plotOverlapRegionFilterXchr:
    input:
        vqsr = os.path.join(config["out_dir"], "results/overlap_sites_filters/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.high.qual.neutral.sites.noPARsXTR.chr.pos.txt"),
        hardf = os.path.join(config["out_dir"], "results/overlap_sites_filters/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.high.qual.neutral.sites.noPARsXTR.chr.pos.txt")
    params:
        srpt = config["scripts_dir"]
    output:
        os.path.join(config["out_dir"], "results/overlap_sites_filters/sex_chromosomes", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.high.qual.neutral.sites.noPARsXTR.venn.diagram.png")
    shell:
        "Rscript {params.srpt}overlappingSites.R {input.vqsr} {input.hardf} {output}"

# 4) Generate SFSs #
# Raw #
# Autosomes
rule SFSPreFilterRegionFilterAutos:
    input:
        vcf = os.path.join(config["out_dir"], "raw/autosomes/filters", "{chrs}.gatk.called.raw.{vcf_options_auto}.high.qual.neutral.sites.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"]
    output:
        sfstxt = os.path.join(config["out_dir"], "results/pre_filter_region_filters/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.high.qual.neutral.sites.SFS.txt")
    shell:
        "python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {output.sfstxt} --ploidy diploid"

# X chromosome
rule SFSPreFilterRegionFilterXchr:
    input:
        vcf = os.path.join(config["out_dir"], "raw/sex_chromosomes/filters", "{chrs}.gatk.called.raw.{vcf_options_x}.high.qual.neutral.sites.noPARsXTR.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"]
    output:
        sfstxt = os.path.join(config["out_dir"], "results/pre_filter_region_filters/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.high.qual.neutral.sites.noPARsXTR.SFS.txt")
    run:
        if wildcards.vcf_options_x == "haploid":
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {output.sfstxt} --ploidy haploid")
        else:
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {output.sfstxt} --ploidy diploid")

# Post VQSR #
# Autosomes
rule SFSVQSRRegionFilterAutos:
    input:
        vcf = os.path.join(config["out_dir"], "vqsr/autosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.high.qual.neutral.sites.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"]
    output:
        sfstxt = os.path.join(config["out_dir"], "results/vqsr_region_filters/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.high.qual.neutral.sites.SFS.txt")
    shell:
        "python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {output.sfstxt} --ploidy diploid"

# X chromosome
rule SFSVQSRRegionFilterXchr:
    input:
        vcf = os.path.join(config["out_dir"], "vqsr/sex_chromosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.high.qual.neutral.sites.noPARsXTR.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"]
    output:
        sfstxt = os.path.join(config["out_dir"], "results/vqsr_region_filters/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.high.qual.neutral.sites.noPARsXTR.SFS.txt")
    run:
        if wildcards.vcf_options_x == "haploid":
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {output.sfstxt} --ploidy haploid")
        else:
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {output.sfstxt} --ploidy diploid")

# Post Hard Filter #
# Autosomes
rule SFSHardFilterRegionFilterAutos:
    input:
        vcf = os.path.join(config["out_dir"], "hard_filter/autosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.high.qual.neutral.sites.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"]
    output:
        sfstxt = os.path.join(config["out_dir"], "results/hard_filter_region_filters/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.high.qual.neutral.sites.SFS.txt")
    shell:
        "python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {output.sfstxt} --ploidy diploid"

# X chromosome
rule SFSHardFilterRegionFilterXchr:
    input:
        vcf = os.path.join(config["out_dir"], "hard_filter/sex_chromosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.high.qual.neutral.sites.noPARsXTR.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"]
    output:
        sfstxt = os.path.join(config["out_dir"], "results/hard_filter_region_filters/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.high.qual.neutral.sites.noPARsXTR.SFS.txt")
    run:
        if wildcards.vcf_options_x == "haploid":
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {output.sfstxt} --ploidy haploid")
        else:
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {output.sfstxt} --ploidy diploid")

# 5) SFSs plot combined #
# Autosomes
# NOTE: For R script order of input is important - Raw first, than VQSR, than hard filter
rule SFSAllRegionFilterAutos:
    input:
        rawfn = os.path.join(config["out_dir"], "results/pre_filter_region_filters/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.high.qual.neutral.sites.SFS.txt"),
        vqsrfn = os.path.join(config["out_dir"], "results/vqsr_region_filters/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.high.qual.neutral.sites.SFS.txt"),
        hfilterfn = os.path.join(config["out_dir"], "results/hard_filter_region_filters/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.high.qual.neutral.sites.SFS.txt")
    params:
        scrpt = config["scripts_dir"]
    output:
        os.path.join(config["out_dir"], "results/merged/all_variable_sites_region_filters/sfs/autosomes", "{chrs}", "{chrs}.{vcf_options_auto}.raw.vqsr.{filtering_options}.high.qual.neutral.sites.SFS.png")
    shell:
        "Rscript {params.scrpt}sfs_plotting.R {input.rawfn} {input.vqsrfn} {input.hfilterfn} {output}"

# X chromosome
rule SFSAllRegionFilterXchr:
    input:
        rawfn = os.path.join(config["out_dir"], "results/pre_filter_region_filters/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.high.qual.neutral.sites.noPARsXTR.SFS.txt"),
        vqsrfn = os.path.join(config["out_dir"], "results/vqsr_region_filters/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.high.qual.neutral.sites.noPARsXTR.SFS.txt"),
        hfilterfn = os.path.join(config["out_dir"], "results/hard_filter_region_filters/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.high.qual.neutral.sites.noPARsXTR.SFS.txt")
    params:
        scrpt = config["scripts_dir"]
    output:
        os.path.join(config["out_dir"], "results/merged/all_variable_sites_region_filters/sfs/sex_chromosomes", "{chrs}", "{chrs}.{vcf_options_x}.raw.vqsr.{filtering_options}.high.qual.neutral.sites.noPARs.noXTR.SFS.png")
    shell:
        "Rscript {params.scrpt}sfs_plotting.R {input.rawfn} {input.vqsrfn} {input.hfilterfn} {output}"

#-------------------------------------------------------------------------------
# Annotation plots for overlapping sites and sites uniquely retained post VQSR
#                                                           and hard filtering #
#-------------------------------------------------------------------------------
# All sites - neutral and high quality regions, no PARs and XTR
# Autosomes #
rule densityPlotsRegionFilterAutos:
    input:
        vqsr = os.path.join(config["out_dir"], "results/overlap_sites_filters/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.high.qual.neutral.sites.chr.pos.txt"),
        hardf = os.path.join(config["out_dir"], "results/overlap_sites_filters/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.high.qual.neutral.sites.chr.pos.txt"),
        vqsrvcf = os.path.join(config["out_dir"], "vqsr/autosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.high.qual.neutral.sites.vcf.gz"),
        hardfvcf = os.path.join(config["out_dir"], "hard_filter/autosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.high.qual.neutral.sites.vcf.gz")
    params:
        srpt = config["scripts_dir"],
        chrm = "{chrs}",
        out1 = os.path.join(config["out_dir"], "results/overlap_sites_filters/autosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_auto}.vqsr.high.qual.neutral.sites.unique.txt"),
        out2 = os.path.join(config["out_dir"], "results/overlap_sites_filters/autosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_auto}.{filtering_options}.high.qual.neutral.sites.unique.txt"),
        out3 = os.path.join(config["out_dir"], "results/overlap_sites_filters/autosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.high.qual.neutral.sites.overlap.txt"),
        vcf1 = os.path.join(config["out_dir"], "results/overlap_sites_filters/autosomes", "{chrs}", "density_plots", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.high.qual.neutral.sites.unique.vcf.gz"),
        vcf2 = os.path.join(config["out_dir"], "results/overlap_sites_filters/autosomes", "{chrs}", "density_plots", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.high.qual.neutral.sites.unique.vcf.gz"),
        vcf3 = os.path.join(config["out_dir"], "results/overlap_sites_filters/autosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.sv.high.qual.neutral.sites.overlap.vcf.gz"),
        annotation1 = os.path.join(config["out_dir"], "results/overlap_sites_filters/autosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.high.qual.neutral.sites.unique_annotations.txt"),
        annotation2 = os.path.join(config["out_dir"], "results/overlap_sites_filters/autosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_auto}.{filtering_options}.high.qual.neutral.sites.unique_annotations.txt"),
        annotation3 = os.path.join(config["out_dir"], "results/overlap_sites_filters/autosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.high.qual.neutral.sites.overlap_annotations.txt")
    output:
        outDP = os.path.join(config["out_dir"], "results/overlap_sites_filters/autosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.high.qual.neutral.sites.DP.png"),
        outQD = os.path.join(config["out_dir"], "results/overlap_sites_filters/autosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.high.qual.neutral.sites.QD.png")
    shell:
        """
        awk 'FNR==NR {{a[$0]++; next}} !a[$0]' {input.hardf} {input.vqsr} | sed -e 's/{params.chrm}/{params.chrm}\t/g' > {params.out1};
        awk 'FNR==NR {{a[$0]++; next}} !a[$0]' {input.vqsr} {input.hardf} | sed -e 's/{params.chrm}/{params.chrm}\t/g' > {params.out2};
        awk 'NR==FNR {{ lines[$0]=1; next }} $0 in lines' {input.hardf} {input.vqsr} | sed -e 's/{params.chrm}/{params.chrm}\t/g' > {params.out3};
        bcftools view -O z  {input.vqsrvcf} -T {params.out1} -o {params.vcf1};
        bcftools view -O z  {input.hardfvcf} -T {params.out2} -o {params.vcf2};
        bcftools view -O z  {input.vqsrvcf} -T {params.out3} -o {params.vcf3};
        python {params.srpt}extract_stats_from_vcf.py AN QD MQ DP --vcf {params.vcf1} --outfile {params.annotation1};
        python {params.srpt}extract_stats_from_vcf.py AN QD MQ DP --vcf {params.vcf2} --outfile {params.annotation2};
        python {params.srpt}extract_stats_from_vcf.py AN QD MQ DP --vcf {params.vcf3} --outfile {params.annotation3};
        Rscript {params.srpt}annotation_density_plots.R {params.annotation1} {params.annotation2} {params.annotation3} {output.outDP} {output.outQD}
        """

# X chromosome #
rule densityPlotRegionFiltersXchr:
    input:
        vqsr = os.path.join(config["out_dir"], "results/overlap_sites_filters/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.high.qual.neutral.sites.noPARsXTR.chr.pos.txt"),
        hardf = os.path.join(config["out_dir"], "results/overlap_sites_filters/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.high.qual.neutral.sites.noPARsXTR.chr.pos.txt"),
        vqsrvcf = os.path.join(config["out_dir"], "vqsr/sex_chromosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.high.qual.neutral.sites.noPARsXTR.vcf.gz"),
        hardfvcf = os.path.join(config["out_dir"], "hard_filter/sex_chromosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.high.qual.neutral.sites.noPARsXTR.vcf.gz")
    params:
        srpt = config["scripts_dir"],
        chrm = "{chrs}",
        out1 = os.path.join(config["out_dir"], "results/overlap_sites_filters/sex_chromosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_x}.vqsr.high.qual.neutral.sites.noPARsXTR.unique.txt"),
        out2 = os.path.join(config["out_dir"], "results/overlap_sites_filters/sex_chromosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_x}.{filtering_options}.high.qual.neutral.sites.noPARsXTR.unique.txt"),
        out3 = os.path.join(config["out_dir"], "results/overlap_sites_filters/sex_chromosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.high.qual.neutral.sites.noPARsXTR.overlap.txt"),
        vcf1 = os.path.join(config["out_dir"], "results/overlap_sites_filters/sex_chromosomes", "{chrs}", "density_plots", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.high.qual.neutral.sites.noPARsXTR.unique.vcf.gz"),
        vcf2 = os.path.join(config["out_dir"], "results/overlap_sites_filters/sex_chromosomes", "{chrs}", "density_plots", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.high.qual.neutral.sites.noPARsXTR.sv.unique.vcf.gz"),
        vcf3 = os.path.join(config["out_dir"], "results/overlap_sites_filters/sex_chromosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.sv.overlap.vcf.gz"),
        annotation1 = os.path.join(config["out_dir"], "results/overlap_sites_filters/sex_chromosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.high.qual.neutral.sites.noPARsXTR.unique_annotations.txt"),
        annotation2 = os.path.join(config["out_dir"], "results/overlap_sites_filters/sex_chromosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_x}.{filtering_options}.high.qual.neutral.sites.noPARsXTR.unique_annotations.txt"),
        annotation3 = os.path.join(config["out_dir"], "results/overlap_sites_filters/sex_chromosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.high.qual.neutral.sites.noPARsXTR.overlap_annotations.txt")
    output:
        outDP = os.path.join(config["out_dir"], "results/overlap_sites_filters/sex_chromosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.high.qual.neutral.sites.noPARsXTR.DP.png"),
        outQD = os.path.join(config["out_dir"], "results/overlap_sites_filters/sex_chromosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.high.qual.neutral.sites.noPARsXTR.QD.png")
    shell:
        """
        awk 'FNR==NR {{a[$0]++; next}} !a[$0]' {input.hardf} {input.vqsr} | sed -e 's/{params.chrm}/{params.chrm}\t/g' > {params.out1};
        awk 'FNR==NR {{a[$0]++; next}} !a[$0]' {input.vqsr} {input.hardf} | sed -e 's/{params.chrm}/{params.chrm}\t/g' > {params.out2};
        awk 'NR==FNR {{ lines[$0]=1; next }} $0 in lines' {input.hardf} {input.vqsr} | sed -e 's/{params.chrm}/{params.chrm}\t/g' > {params.out3};
        bcftools view -O z  {input.vqsrvcf} -T {params.out1} -o {params.vcf1};
        bcftools view -O z  {input.hardfvcf} -T {params.out2} -o {params.vcf2};
        bcftools view -O z  {input.vqsrvcf} -T {params.out3} -o {params.vcf3};
        python {params.srpt}extract_stats_from_vcf.py AN QD MQ DP --vcf {params.vcf1} --outfile {params.annotation1};
        python {params.srpt}extract_stats_from_vcf.py AN QD MQ DP --vcf {params.vcf2} --outfile {params.annotation2};
        python {params.srpt}extract_stats_from_vcf.py AN QD MQ DP --vcf {params.vcf3} --outfile {params.annotation3};
        Rscript {params.srpt}annotation_density_plots.R {params.annotation1} {params.annotation2} {params.annotation3} {output.outDP} {output.outQD}
        """

#-------------------------------------------------------------------------------
# SITES RESTRICTED TO THE ARRAY - NEUTRAL, HIGH QUALITY NO PARS AND XTR #
# Get number of sites, venn diagrams, and SFS

# 1) Filtering - Extract neutral and high quality sites, remove XTR and PARs #
# Pre filter
rule regionFiltersAllArraySitesAutos:
    input:
        os.path.join(config["out_dir"], "raw/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.array.sites.vcf.gz")
    params:
        hiqneubed = os.path.join(config["hi_qual_neu_dir"], "{chrs}_GRCh38_high_quality_neutral_sites.bed")
    output:
        os.path.join(config["out_dir"], "raw/autosomes/filters", "{chrs}.gatk.called.raw.{vcf_options_auto}.array.sites.high.qual.neutral.sites.vcf.gz")
    shell:
        "bedtools intersect -header -a {input} -b {params.hiqneubed} | bgzip -c > {output}"

rule regionFiltersAllArraySitesXchr:
    input:
        os.path.join(config["out_dir"], "raw/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites.vcf.gz")
    params:
        hiqneubed = os.path.join(config["hi_qual_neu_dir"], "{chrs}_GRCh38_high_quality_neutral_sites.bed"),
        parsxtrbed = config["x_filter"],
        filter1 = os.path.join(config["out_dir"], "raw/sex_chromosomes/filters", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites.high.qual.neutral.sites.vcf.gz")
    output:
        os.path.join(config["out_dir"], "raw/sex_chromosomes/filters", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites.high.qual.neutral.sites.noPARsXTR.vcf.gz")
    shell:
        """
        bedtools intersect -header -a {input} -b {params.hiqneubed} | bgzip -c > {params.filter1};
        subtractBed -header -a {params.filter1} -b {params.parsxtrbed} | bgzip -c > {output}
        """

# VQSR
rule regionFiltersVQSRArraySitesAutos:
    input:
        os.path.join(config["out_dir"], "vqsr/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.vcf.gz")
    params:
        hiqneubed = os.path.join(config["hi_qual_neu_dir"], "{chrs}_GRCh38_high_quality_neutral_sites.bed")
    output:
        os.path.join(config["out_dir"], "vqsr/autosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.high.qual.neutral.sites.vcf.gz")
    shell:
        "bedtools intersect -header -a {input} -b {params.hiqneubed} | bgzip -c > {output}"

rule regionFiltersVQSRArraySitesXchr:
    input:
        os.path.join(config["out_dir"], "vqsr/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.vcf.gz")
    params:
        hiqneubed = os.path.join(config["hi_qual_neu_dir"], "{chrs}_GRCh38_high_quality_neutral_sites.bed"),
        parsxtrbed = config["x_filter"],
        filter1 = os.path.join(config["out_dir"], "vqsr/sex_chromosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.high.qual.neutral.sites.vcf.gz")
    output:
        os.path.join(config["out_dir"], "vqsr/sex_chromosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.high.qual.neutral.sites.noPARsXTR.vcf.gz")
    shell:
        """
        bedtools intersect -header -a {input} -b {params.hiqneubed} | bgzip -c > {params.filter1};
        subtractBed -header -a {params.filter1} -b {params.parsxtrbed} | bgzip -c > {output}
        """

# Hard Filter
rule regionFiltersHardFilterArraySitesAutos:
    input:
        os.path.join(config["out_dir"], "hard_filter/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.vcf.gz")
    params:
        hiqneubed = os.path.join(config["hi_qual_neu_dir"], "{chrs}_GRCh38_high_quality_neutral_sites.bed")
    output:
        os.path.join(config["out_dir"], "hard_filter/autosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.vcf.gz")
    shell:
        "bedtools intersect -header -a {input} -b {params.hiqneubed} | bgzip -c > {output}"

rule regionFiltersHardFilterArraySitesXchr:
    input:
        os.path.join(config["out_dir"], "hard_filter/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.vcf.gz")
    params:
        hiqneubed = os.path.join(config["hi_qual_neu_dir"], "{chrs}_GRCh38_high_quality_neutral_sites.bed"),
        parsxtrbed = config["x_filter"],
        filter1 = os.path.join(config["out_dir"], "hard_filter/sex_chromosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.vcf.gz")
    output:
        os.path.join(config["out_dir"], "hard_filter/sex_chromosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.noPARsXTR.vcf.gz")
    shell:
        """
        bedtools intersect -header -a {input} -b {params.hiqneubed} | bgzip -c > {params.filter1};
        subtractBed -header -a {params.filter1} -b {params.parsxtrbed} | bgzip -c > {output}
        """

# 2) Count number of sites for table #
# Raw
rule countSitesregionFiltersAllArraySitesAutos:
    input:
        os.path.join(config["out_dir"], "raw/autosomes/filters", "{chrs}.gatk.called.raw.{vcf_options_auto}.array.sites.high.qual.neutral.sites.vcf.gz")
    output:
        os.path.join(config["out_dir"], "stats_pre_filter_array_sites/autosomes", "{chrs}.gatk.called.raw.{vcf_options_auto}.array.sites.high.qual.neutral.sites", "{chrs}.gatk.called.raw.{vcf_options_auto}.array.sites.high.qual.neutral.sites_num_sites.txt")
    shell:
        """
        zcat {input} | grep -v '#' | wc -l > {output}
        """

rule countSitesregionFiltersAllArraySitesXchr:
    input:
        os.path.join(config["out_dir"], "raw/sex_chromosomes/filters", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites.high.qual.neutral.sites.noPARsXTR.vcf.gz")
    output:
        os.path.join(config["out_dir"], "stats_pre_filter_array_sites/sex_chromosomes", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites.high.qual.neutral.sites.noPARsXTR", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites.high.qual.neutral.noPARsXTR.sites_num_sites.txt")
    shell:
        """
        zcat {input} | grep -v '#' | wc -l > {output}
        """

# VQSR
rule countSitesregionFiltersVQSRArraySitesAutos:
    input:
        os.path.join(config["out_dir"], "vqsr/autosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.high.qual.neutral.sites.vcf.gz")
    output:
        os.path.join(config["out_dir"], "stats_post_filter_array_sites/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.high.qual.neutral.sites", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.high.qual.neutral.sites_num_sites.txt")
    shell:
        """
        zcat {input} | grep -v '#' | wc -l > {output}
        """

rule countSitesregionFiltersVQSRArraySitesXchr:
    input:
        os.path.join(config["out_dir"], "vqsr/sex_chromosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.high.qual.neutral.sites.noPARsXTR.vcf.gz")
    output:
        os.path.join(config["out_dir"], "stats_post_filter_array_sites/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.high.qual.neutral.sites.noPARsXTR", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.high.qual.neutral.sites.noPARsXTR_num_sites.txt")
    shell:
        """
        zcat {input} | grep -v '#' | wc -l > {output}
        """

# Hard Filter
rule countSitesregionFiltersHardFilterArraySitesAutos:
    input:
        os.path.join(config["out_dir"], "hard_filter/autosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.vcf.gz")
    output:
        os.path.join(config["out_dir"], "stats_post_hard_filter_array_sites/autosomes", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites_num_sites.txt")
    shell:
        """
        zcat {input} | grep -v '#' | wc -l > {output}
        """

rule countSitesregionFiltersHardFilterArraySitesXchr:
    input:
        os.path.join(config["out_dir"], "hard_filter/sex_chromosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.noPARsXTR.vcf.gz")
    output:
        os.path.join(config["out_dir"], "stats_post_hard_filter_array_sites/sex_chromosomes", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.noPARsXTR", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.noPARsXTR_num_sites.txt")
    shell:
        """
        zcat {input} | grep -v '#' | wc -l > {output}
        """

# 3) Make venn diagrams comparing VQSR and hard filtering #
# 3.1. Make a list per vcf with chr_pos
# Autosomes
rule chrPosListsVQSRRegionFilterArraySitesAutos:
    input:
        os.path.join(config["out_dir"], "vqsr/autosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.high.qual.neutral.sites.vcf.gz")
    output:
        os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.high.qual.neutral.sites.chr.pos.txt")
    shell:
        """
        zcat {input} | grep -v '#' | awk '{{ print $1'_'$2 }}' >  {output}
        """

rule chrPosListsHardFilterRegionFilterArraySitesAutos:
    input:
        os.path.join(config["out_dir"], "hard_filter/autosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.vcf.gz")
    output:
        os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.chr.pos.txt")
    shell:
        """
        zcat {input} | grep -v '#' | awk '{{ print $1'_'$2 }}' >  {output}
        """

# X chromosome
rule chrPosListsVQSRRegionFilterArraySitesXchr:
    input:
        os.path.join(config["out_dir"], "vqsr/sex_chromosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.high.qual.neutral.sites.noPARsXTR.vcf.gz")
    output:
        os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.high.qual.neutral.sites.noPARsXTR.chr.pos.txt")
    shell:
        """
        zcat {input} | grep -v '#' | awk '{{ print $1'_'$2 }}' >  {output}
        """

rule chrPosListsHardFilterRegionFilterArraySitesXchr:
    input:
        os.path.join(config["out_dir"], "hard_filter/sex_chromosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.noPARsXTR.vcf.gz")
    output:
        os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.noPARsXTR.chr.pos.txt")
    shell:
        """
        zcat {input} | grep -v '#' | awk '{{ print $1'_'$2 }}' >  {output}
        """

# 3.2. Plot results
# Autosomes
rule plotOverlapRegionFilterArraySitesAutos:
    input:
        vqsr = os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.high.qual.neutral.sites.chr.pos.txt"),
        hardf = os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.chr.pos.txt")
    params:
        srpt = config["scripts_dir"]
    output:
        os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/autosomes", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.array.sites.high.qual.neutral.sites.venn.diagram.png")
    shell:
        "Rscript {params.srpt}overlappingSites.R {input.vqsr} {input.hardf} {output}"

# X chromosome
rule plotOverlapRegionFilterArraySitesXchr:
    input:
        vqsr = os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.high.qual.neutral.sites.noPARsXTR.chr.pos.txt"),
        hardf = os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.noPARsXTR.chr.pos.txt")
    params:
        srpt = config["scripts_dir"]
    output:
        os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/sex_chromosomes", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.array.sites.high.qual.neutral.sites.noPARsXTR.venn.diagram.png")
    shell:
        "Rscript {params.srpt}overlappingSites.R {input.vqsr} {input.hardf} {output}"


# 4) Generate SFSs #
# Pre filter #
# Autosomes
rule SFSPreFilterArraySitesRegionFilterAutos:
    input:
        vcf = os.path.join(config["out_dir"], "raw/autosomes/filters", "{chrs}.gatk.called.raw.{vcf_options_auto}.array.sites.high.qual.neutral.sites.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"]
    output:
        sfstxt = os.path.join(config["out_dir"], "results/pre_filter_array_sites_region_filters/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.array.sites.high.qual.neutral.sites.SFS.txt")
    shell:
        "python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {output.sfstxt} --ploidy diploid"

# X chromosome
rule SFSPreFilterArraySitesRegionFilterXchr:
    input:
        vcf = os.path.join(config["out_dir"], "raw/sex_chromosomes/filters", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites.high.qual.neutral.sites.noPARsXTR.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"]
    output:
        sfstxt = os.path.join(config["out_dir"], "results/pre_filter_array_sites_region_filters/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites.high.qual.neutral.sites.noPARsXTR.SFS.txt")
    run:
        if wildcards.vcf_options_x == "haploid":
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {output.sfstxt} --ploidy haploid")
        else:
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {output.sfstxt} --ploidy diploid")

# Post VQSR #
# Autosomes
rule SFSVQSRArraySitesRegionFilterAutos:
    input:
        vcf = os.path.join(config["out_dir"], "vqsr/autosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.high.qual.neutral.sites.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"]
    output:
        sfstxt = os.path.join(config["out_dir"], "results/vqsr_array_sites_region_filters/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.high.qual.neutral.sites.SFS.txt")
    shell:
        "python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {output.sfstxt} --ploidy diploid"


# X chromosome
rule SFSVQSRArraySitesRegionFilterXchr:
    input:
        vcf = os.path.join(config["out_dir"], "vqsr/sex_chromosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.high.qual.neutral.sites.noPARsXTR.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"]
    output:
        sfstxt = os.path.join(config["out_dir"], "results/vqsr_array_sites_region_filters/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.high.qual.neutral.sites.noPARsXTR.SFS.txt")
    run:
        if wildcards.vcf_options_x == "haploid":
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {output.sfstxt} --ploidy haploid")
        else:
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {output.sfstxt} --ploidy diploid")

# Post Hard Filter #
# Autosomes
rule SFSHardFilterArraySitesRegionFilterAutos:
    input:
        vcf = os.path.join(config["out_dir"], "hard_filter/autosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"]
    output:
        sfstxt = os.path.join(config["out_dir"], "results/hard_filter_array_sites_region_filters/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.SFS.txt")
    shell:
        "python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {output.sfstxt} --ploidy diploid"

# X chromosome
rule SFSHardFilterArraySitesRegionFilterXchr:
    input:
        vcf = os.path.join(config["out_dir"], "hard_filter/sex_chromosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.noPARsXTR.vcf.gz")
    params:
        scriptdir = config["popgen_scripts_dir"]
    output:
        sfstxt = os.path.join(config["out_dir"], "results/hard_filter_array_sites_region_filters/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.noPARsXTR.SFS.txt")
    run:
        if wildcards.vcf_options_x == "haploid":
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {output.sfstxt} --ploidy haploid")
        else:
            shell("python {params.scriptdir}popgen_tools.py --vcf_file {input.vcf} --sfs_all --sfs_all_out {output.sfstxt} --ploidy diploid")

# 5) SFSs plot combined #
# Autosomes
# NOTE: For R script order of input is important - Raw first, than VQSR, than hard filter
rule SFSAllArraySitesRegionFilterAutos:
    input:
        rawfn = os.path.join(config["out_dir"], "results/pre_filter_array_sites_region_filters/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_auto}.array.sites.high.qual.neutral.sites.SFS.txt"),
        vqsrfn = os.path.join(config["out_dir"], "results/vqsr_array_sites_region_filters/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.high.qual.neutral.sites.SFS.txt"),
        hfilterfn = os.path.join(config["out_dir"], "results/hard_filter_array_sites_region_filters/sfs/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.SFS.txt")
    params:
        scrpt = config["scripts_dir"]
    output:
        os.path.join(config["out_dir"], "results/merged/array_sites_region_filters/sfs/autosomes", "{chrs}", "{chrs}.{vcf_options_auto}.raw.vqsr.{filtering_options}.array.sites.high.qual.neutral.sites.SFS.png")
    shell:
        "Rscript {params.scrpt}sfs_plotting.R {input.rawfn} {input.vqsrfn} {input.hfilterfn} {output}"

# X chromosome
rule SFSAllArraySitesRegionFilterXchr:
    input:
        rawfn =  os.path.join(config["out_dir"], "results/pre_filter_array_sites_region_filters/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.raw.{vcf_options_x}.array.sites.high.qual.neutral.sites.noPARsXTR.SFS.txt"),
        vqsrfn = os.path.join(config["out_dir"], "results/vqsr_array_sites_region_filters/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.high.qual.neutral.sites.noPARsXTR.SFS.txt"),
        hfilterfn = os.path.join(config["out_dir"], "results/hard_filter_array_sites_region_filters/sfs/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.noPARsXTR.SFS.txt")
    params:
        scrpt = config["scripts_dir"]
    output:
        os.path.join(config["out_dir"], "results/merged/array_sites_region_filters/sfs/sex_chromosomes", "{chrs}", "{chrs}.{vcf_options_x}.raw.vqsr.{filtering_options}.array.sites.high.qual.neutral.sites.noPARsXTR.SFS.png")
    shell:
        "Rscript {params.scrpt}sfs_plotting.R {input.rawfn} {input.vqsrfn} {input.hfilterfn} {output}"


#-------------------------------------------------------------------------------
# Annotation plots for overlapping sites and sites uniquely retained post VQSR
#                                                           and hard filtering #
#-------------------------------------------------------------------------------
# Array sites - no additional filters
# Autosomes #
rule densityPlotsRegionFilterArraySitesAutos:
    input:
        vqsr = os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.high.qual.neutral.sites.chr.pos.txt"),
        hardf = os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/autosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.chr.pos.txt"),
        vqsrvcf = os.path.join(config["out_dir"], "vqsr/autosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.high.qual.neutral.sites.vcf.gz"),
        hardfvcf = os.path.join(config["out_dir"], "hard_filter/autosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.vcf.gz"),
        orgvqsrvcf = os.path.join(config["out_dir"], "vqsr/autosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.high.qual.neutral.sites.vcf.gz"),
        orghardfvcf = os.path.join(config["out_dir"], "hard_filter/autosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.high.qual.neutral.sites.vcf.gz")
    params:
        srpt = config["scripts_dir"],
        chrm = "{chrs}",
        out1 = os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/autosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_auto}.vqsr.array.sites.high.qual.neutral.sites.unique.txt"),
        out2 = os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/autosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_auto}.{filtering_options}.array.sites.high.qual.neutral.sites.unique.txt"),
        out3 = os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/autosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.array.sites.high.qual.neutral.sites.overlap.txt"),
        vcf1 = os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/autosomes", "{chrs}", "density_plots", "{chrs}.gatk.called.{vcf_options_auto}.vqsr.sv.array.sites.high.qual.neutral.sites.unique.vcf.gz"),
        vcf2 = os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/autosomes", "{chrs}", "density_plots", "{chrs}.gatk.called.{vcf_options_auto}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.unique.vcf.gz"),
        vcf3 = os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/autosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.sv.array.sites.high.qual.neutral.sites.overlap.vcf.gz"),
        annotation1 = os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/autosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.array.sites.high.qual.neutral.sites.unique_annotations.txt"),
        annotation2 = os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/autosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_auto}.{filtering_options}.array.sites.high.qual.neutral.sites.unique_annotations.txt"),
        annotation3 = os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/autosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.array.sites.high.qual.neutral.sites.overlap_annotations.txt")
    output:
        outDP = os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/autosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.array.sites.high.qual.neutral.sites.DP.png"),
        outQD = os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/autosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_auto}.vqsr.{filtering_options}.array.sites.high.qual.neutral.sites.QD.png")
    shell:
        """
        awk 'FNR==NR {{a[$0]++; next}} !a[$0]' {input.hardf} {input.vqsr} | sed -e 's/{params.chrm}/{params.chrm}\t/g' > {params.out1};
        awk 'FNR==NR {{a[$0]++; next}} !a[$0]' {input.vqsr} {input.hardf} | sed -e 's/{params.chrm}/{params.chrm}\t/g' > {params.out2};
        awk 'NR==FNR {{ lines[$0]=1; next }} $0 in lines' {input.hardf} {input.vqsr} | sed -e 's/{params.chrm}/{params.chrm}\t/g' > {params.out3};
        bcftools view -O z  {input.orgvqsrvcf} -T {params.out1} -o {params.vcf1};
        bcftools view -O z  {input.orghardfvcf} -T {params.out2} -o {params.vcf2};
        bcftools view -O z  {input.orgvqsrvcf} -T {params.out3} -o {params.vcf3};
        python {params.srpt}extract_stats_from_vcf.py AN QD MQ DP --vcf {params.vcf1} --outfile {params.annotation1};
        python {params.srpt}extract_stats_from_vcf.py AN QD MQ DP --vcf {params.vcf2} --outfile {params.annotation2};
        python {params.srpt}extract_stats_from_vcf.py AN QD MQ DP --vcf {params.vcf3} --outfile {params.annotation3};
        Rscript {params.srpt}annotation_density_plots.R {params.annotation1} {params.annotation2} {params.annotation3} {output.outDP} {output.outQD}
        """

# X chromosome #
rule densityPlotsRegionFilterArraySitesXchr:
    input:
        vqsr = os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.high.qual.neutral.sites.noPARsXTR.chr.pos.txt"),
        hardf = os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/sex_chromosomes", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.noPARsXTR.chr.pos.txt"),
        vqsrvcf = os.path.join(config["out_dir"], "vqsr/sex_chromosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.high.qual.neutral.sites.noPARsXTR.vcf.gz"),
        hardfvcf = os.path.join(config["out_dir"], "hard_filter/sex_chromosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.noPARsXTR.vcf.gz"),
        orgvqsrvcf = os.path.join(config["out_dir"], "vqsr/sex_chromosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.high.qual.neutral.sites.noPARsXTR.vcf.gz"),
        orghardfvcf = os.path.join(config["out_dir"], "hard_filter/sex_chromosomes/filters", "{chrs}", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.high.qual.neutral.sites.noPARsXTR.vcf.gz")
    params:
        srpt = config["scripts_dir"],
        chrm = "{chrs}",
        out1 = os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/sex_chromosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_x}.vqsr.array.sites.high.qual.neutral.sites.noPARsXTR.unique.txt"),
        out2 = os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/sex_chromosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_x}.{filtering_options}.array.sites.high.qual.neutral.sites.noPARsXTR.unique.txt"),
        out3 = os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/sex_chromosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.array.sites.high.qual.neutral.sites.noPARsXTR.overlap.txt"),
        vcf1 = os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/sex_chromosomes", "{chrs}", "density_plots", "{chrs}.gatk.called.{vcf_options_x}.vqsr.sv.array.sites.high.qual.neutral.sites.noPARsXTR.unique.vcf.gz"),
        vcf2 = os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/sex_chromosomes", "{chrs}", "density_plots", "{chrs}.gatk.called.{vcf_options_x}.snps.{filtering_options}.sv.array.sites.high.qual.neutral.sites.noPARsXTR.unique.vcf.gz"),
        vcf3 = os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/sex_chromosomes", "{chrs}", "density_plots", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.sv.array.sites.high.qual.neutral.sites.noPARsXTR.overlap.vcf.gz"),
        annotation1 = os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/sex_chromosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.array.sites.high.qual.neutral.sites.noPARsXTR.unique_annotations.txt"),
        annotation2 = os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/sex_chromosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_x}.{filtering_options}.array.sites.high.qual.neutral.sites.noPARsXTR.unique_annotations.txt"),
        annotation3 = os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/sex_chromosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.array.sites.high.qual.neutral.sites.noPARsXTR.overlap_annotations.txt")
    output:
        outDP = os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/sex_chromosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.array.sites.high.qual.neutral.sites.noPARsXTR.DP.png"),
        outQD = os.path.join(config["out_dir"], "results/overlap_sites_filters_array_sites/sex_chromosomes", "density_plots", "{chrs}", "{chrs}.{vcf_options_x}.vqsr.{filtering_options}.array.sites.high.qual.neutral.sites.noPARsXTR.QD.png")
    shell:
        """
        awk 'FNR==NR {{a[$0]++; next}} !a[$0]' {input.hardf} {input.vqsr} | sed -e 's/{params.chrm}/{params.chrm}\t/g' > {params.out1};
        awk 'FNR==NR {{a[$0]++; next}} !a[$0]' {input.vqsr} {input.hardf} | sed -e 's/{params.chrm}/{params.chrm}\t/g' > {params.out2};
        awk 'NR==FNR {{ lines[$0]=1; next }} $0 in lines' {input.hardf} {input.vqsr} | sed -e 's/{params.chrm}/{params.chrm}\t/g' > {params.out3};
        bcftools view -O z  {input.orgvqsrvcf} -T {params.out1} -o {params.vcf1};
        bcftools view -O z  {input.orghardfvcf} -T {params.out2} -o {params.vcf2};
        bcftools view -O z  {input.orgvqsrvcf} -T {params.out3} -o {params.vcf3};
        python {params.srpt}extract_stats_from_vcf.py AN QD MQ DP --vcf {params.vcf1} --outfile {params.annotation1};
        python {params.srpt}extract_stats_from_vcf.py AN QD MQ DP --vcf {params.vcf2} --outfile {params.annotation2};
        python {params.srpt}extract_stats_from_vcf.py AN QD MQ DP --vcf {params.vcf3} --outfile {params.annotation3};
        Rscript {params.srpt}annotation_density_plots.R {params.annotation1} {params.annotation2} {params.annotation3} {output.outDP} {output.outQD}
        """
