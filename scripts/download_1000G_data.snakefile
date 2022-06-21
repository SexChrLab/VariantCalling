import os

# Env: variant_calling_simulations_project

configfile: "download_1000G_data.config.json"

rule all:
    input:
        expand(os.path.join(config["scratchdir"], "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrm}.recalibrated_variants.vcf.gz"), chrm=config["chromosomes"]),
        expand(os.path.join(config["scratchdir"], "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrm}.recalibrated_variants.vcf.gz.tbi"), chrm=config["chromosomes"])


'''
#expand(os.path.join(config["scratchdir"], "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrm}.recalibrated_variants_subset_10_CEU_males.vcf.gz"), chrm=config["chromosomes"])

# rule: download and extract samples from vcfs (separated by chromosome)
rule downloadExtractSamples:
    params:
        inpath = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrm}.recalibrated_variants.vcf.gz",
        chrm = "{chrm}"
    output:
        os.path.join(config["scratchdir"], "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrm}.recalibrated_variants_subset_10_CEU_males.vcf.gz")
    shell:
        """
        PERL5LIB="";
        tabix -h {params.inpath} {params.chrm} | vcf-subset -c HG00096,HG00101,HG00103,HG00105,HG00107,HG00108,HG00109,HG00112,HG00113,HG00114 | bgzip -c > {output}
        """
'''

rule downloadvcf:
    params:
        inpath = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrm}.recalibrated_variants.vcf.gz",
    output:
        os.path.join(config["scratchdir"], "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrm}.recalibrated_variants.vcf.gz")
    shell:
        """
        wget {params.inpath}
        """

rule downloadtbi:
    params:
        inpathtbi = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrm}.recalibrated_variants.vcf.gz.tbi",
    output:
        os.path.join(config["scratchdir"], "20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrm}.recalibrated_variants.vcf.gz.tbi")
    shell:
        """
        wget {params.inpathtbi}
        """

# may want to try curl -C - -O {params.inpathtbi}
# wget -c {params.inpathtbi}
