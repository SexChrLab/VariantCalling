import os

# Environment: NEAT_env

configfile: "xy_var_calling_pilot_02.config.json"


rule all:
    input:
        expand(os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_read1.fastq.gz"), sample=config["males"]),
        expand(os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_read2.fastq.gz"), sample=config["males"]),

        #expand(os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf.gz"), sample=config["males"], chrms=config["chrY"]),
        #expand(os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_read1.fq.gz"), sample=config["males"], chrms=config["chrY"]),
        #expand(os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_read2.fq.gz"), sample=config["males"], chrms=config["chrY"]),

        #expand(os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf.gz"), sample=config["males"], chrms=config["chrX"]),
        #expand(os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_read1.fq.gz"), sample=config["males"], chrms=config["chrX"]),
        #expand(os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_read2.fq.gz"), sample=config["males"], chrms=config["chrX"]),

        #expand(os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf.gz"), sample=config["males"], chrms=config["chrX"]),
        #expand(os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_read1.fq.gz"), sample=config["males"], chrms=config["chrX"]),
        #expand(os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_read2.fq.gz"), sample=config["males"], chrms=config["chrX"]),

        #expand(os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf.gz"), sample=config["males"], chrms=config["autosomes"]),
        #expand(os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_read1.fq.gz"), sample=config["males"], chrms=config["autosomes"]),
        #expand(os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_read2.fq.gz"), sample=config["males"], chrms=config["autosomes"]),

        #expand(os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf.gz"), sample=config["males"], chrms=config["chrM"]),
        #expand(os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_read1.fq.gz"), sample=config["males"], chrms=config["chrM"]),
        #expand(os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_read2.fq.gz"), sample=config["males"], chrms=config["chrM"]),



        #expand(os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf.gz"), sample=config["females"], chrms=config["chrX"]),
        #expand(os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_read1.fq.gz"), sample=config["females"], chrms=config["chrX"]),
        #expand(os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_read2.fq.gz"), sample=config["females"], chrms=config["chrX"]),

        #expand(os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf.gz"), sample=config["females"], chrms=config["chrX"]),
        #expand(os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_read1.fq.gz"), sample=config["females"], chrms=config["chrX"]),
        #expand(os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_read2.fq.gz"), sample=config["females"], chrms=config["chrX"]),

        #expand(os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf.gz"), sample=config["females"], chrms=config["autosomes"]),
        #expand(os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_read1.fq.gz"), sample=config["females"], chrms=config["autosomes"]),
        #expand(os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_read2.fq.gz"), sample=config["females"], chrms=config["autosomes"]),

        #expand(os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf.gz"), sample=config["females"], chrms=config["chrM"]),
        #expand(os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_read1.fq.gz"), sample=config["females"], chrms=config["chrM"]),
        #expand(os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_read2.fq.gz"), sample=config["females"], chrms=config["chrM"]),

        expand(os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_read1.fastq.gz"), sample=config["females"]),
        expand(os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_read2.fastq.gz"), sample=config["females"])


#------------------------------------------------------------------------------#
# Step: Run NEAT
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#                                     MALES                                    #
#------------------------------------------------------------------------------#
# Y chromosome #
rule run_NEAT_Y:
    input:
        ref = os.path.join(config["tmp_scratch_path"], "refs/GRCh38_SCC/GRCh38_full_analysis_set_plus_decoy_hla_YPARmasked_{chrms}.fa"),
        vcf = os.path.join(config["tmp_scratch_path"], "vcfs/males/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_reformat.vcf"),
    params:
        neatpath = config["neat"],
        pypath = config["py_V3_8"],
        opath = os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}"),
        drbed = config["drbed"], # regions to remove
        seed = 4,
        coverage = 10
    priority: 50
    output:
        #bam = os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_nonPARs_golden.bam"),
        vcf = os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_golden.vcf.gz"),
        fq1 = os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_read1.fq.gz"),
        fq2 = os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_read2.fq.gz")
    shell:
        """
        {params.pypath} {params.neatpath}gen_reads.py -r {input.ref} -R 150 -o {params.opath} --bam --vcf --pe 300 30 -p 1 -v {input.vcf} -M 0 -c {params.coverage} --rng {params.seed}
        """

# X chromosome #
rule run_NEAT_X_PARs:
    input:
        ref = os.path.join(config["tmp_scratch_path"], "refs/GRCh38_SCC/GRCh38_full_analysis_set_plus_decoy_hla_YPARmasked_{chrms}.fa"),
        #vcf = os.path.join(config["proj_path"], "vcfs/EUR/males/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_reformat.vcf"),
        vcfPARs = os.path.join(config["tmp_scratch_path"], "vcfs/males/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_PARs_reformat.vcf"),
        #vcfnonPARs = os.path.join(config["tmp_scratch_path"], "vcfs/males/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_nonPARs_reformat.vcf")
    params:
        neatpath = config["neat"],
        pypath = config["py_V3_8"],
        opathPARs = os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs"),
        trPARs = config["trPARs"],
        trNonPARs = config["trNonPARs"],
        seed = 4,
        coverage = 20
    priority: 50
    output:
        vcfPARs = os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf.gz"),
        fq1PARs = os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_read1.fq.gz"),
        fq2PARs = os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_read2.fq.gz"),
    shell:
        """
        {params.pypath} {params.neatpath}gen_reads_edit.py -r {input.ref} -R 150 -o {params.opathPARs} --bam --vcf --pe 300 30 -p 2 -v {input.vcfPARs} -M 0 -c {params.coverage} -tr {params.trPARs} --rng {params.seed}
        """

rule run_NEAT_X_nonPARs:
    input:
        ref = os.path.join(config["tmp_scratch_path"], "refs/GRCh38_SCC/GRCh38_full_analysis_set_plus_decoy_hla_YPARmasked_{chrms}.fa"),
        #vcf = os.path.join(config["proj_path"], "vcfs/EUR/males/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_reformat.vcf"),
        #vcfPARs = os.path.join(config["tmp_scratch_path"], "vcfs/males/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_PARs.vcf"),
        vcfnonPARs = os.path.join(config["tmp_scratch_path"], "vcfs/males/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_nonPARs_reformat.vcf"),
    params:
        neatpath = config["neat"],
        pypath = config["py_V3_8"],
        opathNonPARs = os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs"),
        trPARs = config["trPARs"],
        trNonPARs = config["trNonPARs"],
        seed = 4,
        coverage = 10
    priority: 50
    output:
        vcfnonPARs = os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_golden.vcf.gz"),
        fq1nonPARs = os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_read1.fq.gz"),
        fq2nonPARs = os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_read2.fq.gz")
    shell:
        """
        {params.pypath} {params.neatpath}gen_reads_edit.py -r {input.ref} -R 150 -o {params.opathNonPARs} --bam --vcf --pe 300 30 -p 1 -v {input.vcfnonPARs} -M 0 -c {params.coverage} -tr {params.trNonPARs} --rng {params.seed}
        """


# Autosomes #
# I may add mtDNA here because there are het sites that we may want to simulate?
rule run_NEAT_autos:
    input:
        ref = os.path.join(config["tmp_scratch_path"], "refs/GRCh38_SCC/GRCh38_full_analysis_set_plus_decoy_hla_YPARmasked_{chrms}.fa"),
        vcf = os.path.join(config["tmp_scratch_path"], "vcfs/males/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_reformat.vcf"),
    params:
        neatpath = config["neat"],
        pypath = config["py_V3_8"],
        opath = os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}"),
        seed = 4,
        coverage = 20
    priority: 50
    output:
        #bam = os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.bam"),
        vcf = os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf.gz"),
        fq1 = os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_read1.fq.gz"),
        fq2 = os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_read2.fq.gz")
    shell:
        """
        {params.pypath} {params.neatpath}gen_reads.py -r {input.ref} -R 150 -o {params.opath} --bam --vcf --pe 300 30 -p 2 -v {input.vcf} -M 0 -c {params.coverage} --rng {params.seed}
        """


# mtDNA #
rule run_NEAT_chrM:
    input:
        ref = os.path.join(config["tmp_scratch_path"], "refs/GRCh38_SCC/GRCh38_full_analysis_set_plus_decoy_hla_YPARmasked_{chrms}.fa"),
        vcf = os.path.join(config["tmp_scratch_path"], "vcfs/males/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_reformat.vcf"),
    params:
        neatpath = config["neat"],
        pypath = config["py_V3_8"],
        opath = os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}"),
        seed = 4,
        coverage = 20
    output:
        vcf = os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf.gz"),
        fq1 = os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_read1.fq.gz"),
        fq2 = os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_read2.fq.gz")
    shell:
        """
        {params.pypath} {params.neatpath}gen_reads.py -r {input.ref} -R 150 -o {params.opath} --bam --vcf --pe 300 30 -p 1 -v {input.vcf} -M 0 -c {params.coverage} --rng {params.seed}
        """


#------------------------------------------------------------------------------#
# Step: Merge all per chromosome results
#------------------------------------------------------------------------------#
rule merge_fastqs_read1:
    input:
        fq1 = expand(os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_read1.fq.gz"), sample=config["males"], chrms=config["autosomes"]),
        fq1mtDNA = expand(os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_read1.fq.gz"), sample=config["males"], chrms=config["chrM"]),
        fq1nonPARsY = expand(os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_read1.fq.gz"), sample=config["males"], chrms=config["chrY"]),
        fq1nonPARs = expand(os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_read1.fq.gz"), sample=config["males"], chrms=config["chrX"]),
        fq1PARs = expand(os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_read1.fq.gz"), sample=config["males"], chrms=config["chrX"]),
    params:
        fqpth = os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT"),
        fq1 = os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_read1.fastq")
    priority: 1
    output:
        fq1 = os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_read1.fastq.gz"),
    shell:
        """
        zcat {params.fqpth}*_read1.fq.gz > {params.fq1};
        gzip {params.fq1}
        """

rule merge_fastqs_read2:
    input:
        fq2 = expand(os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_read2.fq.gz"), sample=config["males"], chrms=config["autosomes"]),
        fq1mtDNA = expand(os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_read2.fq.gz"), sample=config["males"], chrms=config["chrM"]),
        fq1nonPARsY = expand(os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_read2.fq.gz"), sample=config["males"], chrms=config["chrY"]),
        fq1nonPARs = expand(os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_10x_haploid_len150_PE_{chrms}_nonPARs_read2.fq.gz"), sample=config["males"], chrms=config["chrX"]),
        fq2PARs = expand(os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_read2.fq.gz"), sample=config["males"], chrms=config["chrX"])
    params:
        fqpth = os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT"),
        fq2 = os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_read2.fastq")
    priority: 1
    output:
        fq2 = os.path.join(config["proj_path"], "simulations/EUR/males/{sample}/{sample}_NEAT_read2.fastq.gz")
    shell:
        """
        zcat {params.fqpth}*_read2.fq.gz > {params.fq2};
        gzip {params.fq2}
        """



#------------------------------------------------------------------------------#
#                                   FEMALES                                    #
#------------------------------------------------------------------------------#
# X chromosome #
rule run_NEAT_X_PARs_females:
    input:
        ref = os.path.join(config["tmp_scratch_path"], "refs/GRCh38_SCC/GRCh38_full_analysis_set_plus_decoy_hla_YPARmasked_{chrms}.fa"),
        #vcf = os.path.join(config["proj_path"], "vcfs/EUR/females/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_reformat.vcf"),
        vcfPARs = os.path.join(config["tmp_scratch_path"], "vcfs/females/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_PARs_reformat.vcf"),
        #vcfnonPARs = os.path.join(config["tmp_scratch_path"], "vcfs/females/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_nonPARs_reformat.vcf")
    params:
        neatpath = config["neat"],
        pypath = config["py_V3_8"],
        opathPARs = os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs"),
        trPARs = config["trPARs"],
        trNonPARs = config["trNonPARs"],
        seed = 4,
        coverage = 20
    priority: 50
    output:
        vcfPARs = os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_golden.vcf.gz"),
        fq1PARs = os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_read1.fq.gz"),
        fq2PARs = os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_read2.fq.gz"),
    shell:
        """
        {params.pypath} {params.neatpath}gen_reads_edit.py -r {input.ref} -R 150 -o {params.opathPARs} --bam --vcf --pe 300 30 -p 2 -v {input.vcfPARs} -M 0 -c {params.coverage} -tr {params.trPARs} --rng {params.seed}
        """

rule run_NEAT_X_nonPARs_females_dip:
    input:
        ref = os.path.join(config["tmp_scratch_path"], "refs/GRCh38_SCC/GRCh38_full_analysis_set_plus_decoy_hla_YPARmasked_{chrms}.fa"),
        #vcf = os.path.join(config["proj_path"], "vcfs/EUR/females/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_reformat.vcf"),
        #vcfPARs = os.path.join(config["tmp_scratch_path"], "vcfs/females/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_PARs.vcf"),
        vcfnonPARs = os.path.join(config["tmp_scratch_path"], "vcfs/females/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_nonPARs_reformat.vcf"),
    params:
        neatpath = config["neat"],
        pypath = config["py_V3_8"],
        opathNonPARs = os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs"),
        trPARs = config["trPARs"],
        trNonPARs = config["trNonPARs"],
        seed = 4,
        coverage = 20
    priority: 50
    output:
        vcfnonPARs = os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_golden.vcf.gz"),
        fq1nonPARs = os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_read1.fq.gz"),
        fq2nonPARs = os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_read2.fq.gz")
    shell:
        """
        {params.pypath} {params.neatpath}gen_reads_edit.py -r {input.ref} -R 150 -o {params.opathNonPARs} --bam --vcf --pe 300 30 -p 2 -v {input.vcfnonPARs} -M 0 -c {params.coverage} -tr {params.trNonPARs} --rng {params.seed}
        """


# Autosomes #
# I may add mtDNA here because there are het sites that we may want to simulate?
rule run_NEAT_autos_females:
    input:
        ref = os.path.join(config["tmp_scratch_path"], "refs/GRCh38_SCC/GRCh38_full_analysis_set_plus_decoy_hla_YPARmasked_{chrms}.fa"),
        vcf = os.path.join(config["tmp_scratch_path"], "vcfs/females/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_reformat.vcf"),
    params:
        neatpath = config["neat"],
        pypath = config["py_V3_8"],
        opath = os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}"),
        seed = 4,
        coverage = 20
    priority: 50
    output:
        #bam = os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.bam"),
        vcf = os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.vcf.gz"),
        fq1 = os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_read1.fq.gz"),
        fq2 = os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_read2.fq.gz")
    shell:
        """
        {params.pypath} {params.neatpath}gen_reads.py -r {input.ref} -R 150 -o {params.opath} --bam --vcf --pe 300 30 -p 2 -v {input.vcf} -M 0 -c {params.coverage} --rng {params.seed}
        """


# mtDNA #
rule run_NEAT_chrM_females:
    input:
        ref = os.path.join(config["tmp_scratch_path"], "refs/GRCh38_SCC/GRCh38_full_analysis_set_plus_decoy_hla_YPARmasked_{chrms}.fa"),
        vcf = os.path.join(config["tmp_scratch_path"], "vcfs/females/{sample}/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{chrms}.recalibrated_variants_{sample}_biallelic_SNPs_reformat.vcf"),
    params:
        neatpath = config["neat"],
        pypath = config["py_V3_8"],
        opath = os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}"),
        seed = 4,
        coverage = 20
    output:
        #bam = os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_golden.bam"),
        vcf = os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_golden.vcf.gz"),
        fq1 = os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_read1.fq.gz"),
        fq2 = os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_read2.fq.gz")
    shell:
        """
        {params.pypath} {params.neatpath}gen_reads.py -r {input.ref} -R 150 -o {params.opath} --bam --vcf --pe 300 30 -p 1 -v {input.vcf} -M 0 -c {params.coverage} --rng {params.seed}
        """


#------------------------------------------------------------------------------#
# Step: Merge all per chromosome results
#------------------------------------------------------------------------------#
# Note - X was simulated as diploid i just accidentially named the file haploid
# Will haved to change naming of nonPARs to 20x_diploid from 10x_haploid
rule merge_fastqs_read1_females:
    input:
        fq1 = expand(os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_read1.fq.gz"), sample=config["females"], chrms=config["autosomes"]),
        fq1mtDNA = expand(os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_read1.fq.gz"), sample=config["females"], chrms=config["chrM"]),
        fq1nonPARs = expand(os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_read1.fq.gz"), sample=config["females"], chrms=config["chrX"]),
        fq1PARs = expand(os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_read1.fq.gz"), sample=config["females"], chrms=config["chrX"]),
    params:
        fqpth = os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT"),
        fq1 = os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_read1.fastq")
    priority: 1
    output:
        fq1 = os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_read1.fastq.gz"),
    shell:
        """
        zcat {params.fqpth}*_read1.fq.gz > {params.fq1};
        gzip {params.fq1}
        """

rule merge_fastqs_read2_females:
    input:
        fq2 = expand(os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_read2.fq.gz"), sample=config["females"], chrms=config["autosomes"]),
        fq1mtDNA = expand(os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_haploid_len150_PE_{chrms}_read2.fq.gz"), sample=config["females"], chrms=config["chrM"]),
        fq1nonPARs = expand(os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_nonPARs_read2.fq.gz"), sample=config["females"], chrms=config["chrX"]),
        fq2PARs = expand(os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_simulated_20x_diploid_len150_PE_{chrms}_PARs_read2.fq.gz"), sample=config["females"], chrms=config["chrX"])
    params:
        fqpth = os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT"),
        fq2 = os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_read2.fastq")
    priority: 1
    output:
        fq2 = os.path.join(config["proj_path"], "simulations/EUR/females/{sample}/{sample}_NEAT_read2.fastq.gz")
    shell:
        """
        zcat {params.fqpth}*_read2.fq.gz > {params.fq2};
        gzip {params.fq2}
        """
