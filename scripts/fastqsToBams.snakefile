import os

# Environment: kenya

configfile: "fastqsToBams.config.json"


## PATHS AND DIRECTORIES ##
# the fastq directory additionally has {sample_name}/Unaligned/
fastq_directory = "fastq_files/"

# Path to packages, user may change if they are not using conda installed packages
gatk_path_old = "/home/amtarave/packages/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar"
gatk_path = "/home/amtarave/packages/gatk-4.0.11.0/gatk-package-4.0.11.0-local.jar"
fastqc_path = "fastqc"
multiqc_path = "multiqc"
trimmomatic_path = "trimmomatic"
bbduksh_path = "bbduk.sh"
bbmerge_sh_path = "bbmerge.sh"
picard_path = "picard"
sambamba_path = "sambamba"
samtools_path = "samtools"
bwa_path = "bwa"


#this needs to be changed to reflect new symbolic file name
fastq_prefixes = [
	config[x]["fq1_sy"][:-9] for x in config["unique_identifier"]] + [
		config[x]["fq2_sy"][:-9] for x in config["unique_identifier"]]


rule all:
	input:
		expand("refs/{ref_type}.fa",
			ref_type=["Ref_GRCh38","Ref_GRCh38_Y_HardMasked","Ref_GRCh38_Y_PARsMasked"]),
		expand("refs/{ref_type}.fa.fai",
			ref_type=["Ref_GRCh38","Ref_GRCh38_Y_HardMasked","Ref_GRCh38_Y_PARsMasked"]),
		expand("refs/{ref_type}.dict",
			ref_type=["Ref_GRCh38","Ref_GRCh38_Y_HardMasked","Ref_GRCh38_Y_PARsMasked"]),
		expand("refs/{ref_type}.fa.amb",
			ref_type=["Ref_GRCh38","Ref_GRCh38_Y_HardMasked","Ref_GRCh38_Y_PARsMasked"]),
		expand(
			"fastq_files/{unique_identifier}_R1.fastq.gz",
			unique_identifier=config["unique_identifier"]),
		expand(
			"fastq_files/{unique_identifier}_R2.fastq.gz",
			unique_identifier=config["unique_identifier"]),
		#"multiqc/multiqc_report.html",
		expand(
			"/scratch/amtarave/Kenya_agave/whole_genome/trimmed_fastqs/{sample}_trimmed_R1.fastq.gz",
			sample=config["unique_identifier"]),
		expand(
			"/scratch/amtarave/Kenya_agave/whole_genome/trimmed_fastqs/{sample}_trimmed_R2.fastq.gz",
			sample=config["unique_identifier"]),
		#"multiqc_trimmed/multiqc_report.html",
		#expand("/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{unique_identifier}.GRCh38_minusYPARs.sorted.bam",
		#	unique_identifier=config["unique_identifier_males"]),
		#expand(
		#	"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{unique_identifier}.GRCh38_minusYPARs.sorted.bam.bai",
		#	unique_identifier=config["unique_identifier_males"]),
		#expand(
		#	"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{unique_identifier}.GRCh38_Ymasked.sorted.bam",
		#	unique_identifier=config["unique_identifier_females"]),
		#expand(
		#	"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{unique_identifier}.GRCh38_Ymasked.sorted.bam.bai",
		#	unique_identifier=config["unique_identifier_females"]),
		#expand(
		#	"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{sample_males}.GRCh38_minusYPARs.sorted.merged.bam",
		#	sample_males=config["males"]),
		#expand(
		#	"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{sample_females}.GRCh38_Ymasked.sorted.merged.bam",
		#	sample_females=config["females"]),
		#expand(
		#	"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{sample_males}.GRCh38_minusYPARs.sorted.merged.bam.bai",
		#	sample_males=config["males"]),
		#expand(
		#	"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{sample_females}.GRCh38_Ymasked.sorted.merged.bam.bai",
		#	sample_females=config["females"]),
		#expand(
		#	"/scratch/amtarave/Kenya_agave/whole_genome/stats/{sample}.{genome}.sorted.bam.stats".split(), zip,
		#	sample=config["samples_list_ordered"], genome=config["sex_specific_ref_list_ordered"]),
		#expand(
		#	"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam".split(), zip,
		#	sample=config["samples_list_ordered"], genome=config["sex_specific_ref_list_ordered"]),
		#expand(
		#	"/scratch/amtarave/Kenya_agave/whole_genome/stats/{sample}.{genome}.picard_mkdup_metrics.txt".split(), zip,
		#	sample=config["samples_list_ordered"], genome=config["sex_specific_ref_list_ordered"]),
		#expand(
		#	"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam.bai".split(), zip,
		#	sample=config["samples_list_ordered"], genome=config["sex_specific_ref_list_ordered"]),
		#expand(
		#	"/scratch/amtarave/Kenya_agave/whole_genome/stats/{sample}.{genome}.sorted.mkdup.bam.stats".split(), zip,
		#	sample=config["samples_list_ordered"], genome=config["sex_specific_ref_list_ordered"])
		expand("/scratch/amtarave/Kenya_agave/whole_genome/processed_bams_default_ref/{unique_identifier}.GRCh38_default.sorted.bam", unique_identifier=config["unique_identifier"]),
		expand("/scratch/amtarave/Kenya_agave/whole_genome/processed_bams_default_ref/{unique_identifier}.GRCh38_default.sorted.bam.bai", unique_identifier=config["unique_identifier"]),
		expand("/scratch/amtarave/Kenya_agave/whole_genome/processed_bams_default_ref/{samples}.GRCh38_default.sorted.merged.bam", samples=config["sample_names"]),
		expand("/scratch/amtarave/Kenya_agave/whole_genome/processed_bams_default_ref/{samples}.GRCh38_default.sorted.merged.bam.bai", samples=config["sample_names"]),
		expand("/scratch/amtarave/Kenya_agave/whole_genome/stats_bams_default_ref/{samples}.GRCh38_default.sorted.merged.bam.stats", samples=config["sample_names"]),
		expand("/scratch/amtarave/Kenya_agave/whole_genome/processed_bams_default_ref/{samples}.GRCh38_default.sorted.merged.mkdup.bam", samples=config["sample_names"]),
		expand("/scratch/amtarave/Kenya_agave/whole_genome/stats_bams_default_ref/{samples}.GRCh38_default.sorted.merged.picard_mkdup_metrics.txt", samples=config["sample_names"]),
		expand("/scratch/amtarave/Kenya_agave/whole_genome/processed_bams_default_ref/{samples}.GRCh38_default.sorted.merged.mkdup.bam.bai", samples=config["sample_names"]),
		expand("/scratch/amtarave/Kenya_agave/whole_genome/stats_bams_default_ref/{samples}.GRCh38_default.sorted.merged.mkdup.bam.stats", samples=config["sample_names"])


rule prep_refs_mk_sy_ln:
	input:
		ref = lambda wildcards: config["genome_paths"][wildcards.ref_type]
	output:
		ref_sy_ln = "refs/{ref_type}.fa",
	#priority:110
	shell:
		"""
		ln -s {input.ref} {output.ref_sy_ln}
		"""

rule prep_refs:
	input:
		"refs/{ref_type}.fa"
	output:
		fai = "refs/{ref_type}.fa.fai",
		dict = "refs/{ref_type}.dict",
		amb = "refs/{ref_type}.fa.amb"
	#priority:105
	params:
		samtools = samtools_path,
		bwa = bwa_path
	shell:
		"""
		{params.samtools} faidx {input};
		{params.samtools} dict -o {output.dict} {input};
		{params.bwa} index {input}
		"""

rule make_symbolic_link_for_fastqs:
	input:
		original_R1 = lambda wildcards: config[wildcards.unique_identifier]["fq_paths"] + config[wildcards.unique_identifier]["stem_name"] + "_R1_001.fastq.gz",
		original_R2 = lambda wildcards: config[wildcards.unique_identifier]["fq_paths"] + config[wildcards.unique_identifier]["stem_name"] + "_R2_001.fastq.gz"
	output:
		R1_out = "fastq_files/{unique_identifier}_R1.fastq.gz",
		R2_out = "fastq_files/{unique_identifier}_R2.fastq.gz"
	shell:
		"""
		ln -s {input.original_R1} {output.R1_out} && touch -h {output.R1_out};
		ln -s {input.original_R2} {output.R2_out} && touch -h {output.R2_out};
		"""


rule fastqc_analysis:
	input:
		os.path.join(fastq_directory, "{fq_prefix}.fastq.gz")
	output:
		"fastqc/{fq_prefix}_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o fastqc/ {input}"


rule multiqc_analysis:
	input:
		expand("fastqc/{fq_prefix}_fastqc.html", fq_prefix=fastq_prefixes)
	output:
		"multiqc/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} --interactive -o multiqc fastqc"


rule trim_adapters_paired_bbduk:
	input:
		fq1 = lambda wildcards: os.path.join(
			fastq_directory, config[wildcards.sample]["fq1_sy"]),
		fq2 = lambda wildcards: os.path.join(
			fastq_directory, config[wildcards.sample]["fq2_sy"])
	output:
		out_fq1 = "/scratch/amtarave/Kenya_agave/whole_genome/trimmed_fastqs/{sample}_trimmed_R1.fastq.gz",
		out_fq2 = "/scratch/amtarave/Kenya_agave/whole_genome/trimmed_fastqs/{sample}_trimmed_R2.fastq.gz"
	params:
		bbduksh = bbduksh_path
	shell:
		"{params.bbduksh} -Xmx3g in1={input.fq1} in2={input.fq2} "
		"out1={output.out_fq1} out2={output.out_fq2} "
		"ref=misc/adapter_sequence.fa ktrim=r k=21 mink=11 hdist=2 tbo tpe "
		"qtrim=rl trimq=15 minlen=75 maq=20"


rule fastqc_analysis_trimmed:
	input:
		fq1 = "/scratch/amtarave/Kenya_agave/whole_genome/trimmed_fastqs/{sample}_trimmed_R1.fastq.gz",
		fq2 = "/scratch/amtarave/Kenya_agave/whole_genome/trimmed_fastqs/{sample}_trimmed_R2.fastq.gz"
	output:
		html1 = "trimmed_fastqc/{sample}_trimmed_R1_fastqc.html",
		html2 = "trimmed_fastqc/{sample}_trimmed_R2_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o trimmed_fastqc/ {input.fq1} {input.fq2}"


rule multiqc_analysis_trimmed:
	input:
		expand(
			"trimmed_fastqc/{sample}_trimmed_{read}_fastqc.html",
			sample=config["unique_identifier"], read=["R1", "R2"])
	output:
		"multiqc_trimmed/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} --interactive -o multiqc_trimmed trimmed_fastqc"


### Before this step I have to prepare the ref files by using bwa index. I should probably
#	make a refs folder in the project directory.
rule map_to_sex_specific_refs_males:
	input:
		fq1 = "/scratch/amtarave/Kenya_agave/whole_genome/trimmed_fastqs/{unique_identifier}_trimmed_R1.fastq.gz",
		fq2 = "/scratch/amtarave/Kenya_agave/whole_genome/trimmed_fastqs/{unique_identifier}_trimmed_R2.fastq.gz",
		new = config["Ref_GRCh38_Y_PARsMasked"],
		fai = expand("refs/{ref_type}.fa.fai", ref_type=["Ref_GRCh38","Ref_GRCh38_Y_HardMasked","Ref_GRCh38_Y_PARsMasked"]),
		dict = expand("refs/{ref_type}.dict", ref_type=["Ref_GRCh38","Ref_GRCh38_Y_HardMasked","Ref_GRCh38_Y_PARsMasked"]),
		amb = expand("refs/{ref_type}.fa.amb", ref_type=["Ref_GRCh38","Ref_GRCh38_Y_HardMasked","Ref_GRCh38_Y_PARsMasked"])
	output:
		"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{unique_identifier}.GRCh38_minusYPARs.sorted.bam"
	params:
		id = lambda wildcards: config[wildcards.unique_identifier]["ID"],
		sm = lambda wildcards: config[wildcards.unique_identifier]["SM"],
		lb = lambda wildcards: config[wildcards.unique_identifier]["LB"],
		pu = lambda wildcards: config[wildcards.unique_identifier]["PU"],
		pl = lambda wildcards: config[wildcards.unique_identifier]["PL"],
		bwa = bwa_path,
		samtools = samtools_path,
		threads = 4
	threads: 4
	#priority: 100
	shell:
		" {params.bwa} mem -t {params.threads} -R "
	 	"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
		"{input.new} {input.fq1} {input.fq2}"
		" | {params.samtools} fixmate -O bam - - | {params.samtools} sort "
		"-O bam -o {output}"


rule map_to_sex_specific_refs_females:
	input:
		fq1 = "/scratch/amtarave/Kenya_agave/whole_genome/trimmed_fastqs/{unique_identifier}_trimmed_R1.fastq.gz",
		fq2 = "/scratch/amtarave/Kenya_agave/whole_genome/trimmed_fastqs/{unique_identifier}_trimmed_R2.fastq.gz",
		new = config["Ref_GRCh38_Y_HardMasked"],
		fai = expand("refs/{ref_type}.fa.fai", ref_type=["Ref_GRCh38","Ref_GRCh38_Y_HardMasked","Ref_GRCh38_Y_PARsMasked"]),
		dict = expand("refs/{ref_type}.dict", ref_type=["Ref_GRCh38","Ref_GRCh38_Y_HardMasked","Ref_GRCh38_Y_PARsMasked"]),
		amb = expand("refs/{ref_type}.fa.amb", ref_type=["Ref_GRCh38","Ref_GRCh38_Y_HardMasked","Ref_GRCh38_Y_PARsMasked"])
	output:
		"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{unique_identifier}.GRCh38_Ymasked.sorted.bam"
	params:
		id = lambda wildcards: config[wildcards.unique_identifier]["ID"],
		sm = lambda wildcards: config[wildcards.unique_identifier]["SM"],
		lb = lambda wildcards: config[wildcards.unique_identifier]["LB"],
		pu = lambda wildcards: config[wildcards.unique_identifier]["PU"],
		pl = lambda wildcards: config[wildcards.unique_identifier]["PL"],
		bwa = bwa_path,
		samtools = samtools_path,
		threads = 4
	threads: 4
	#priority: 95
	shell:
		" {params.bwa} mem -t {params.threads} -R "
	 	"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
		"{input.new} {input.fq1} {input.fq2}"
		" | {params.samtools} fixmate -O bam - - | {params.samtools} sort "
		"-O bam -o {output}"


rule index_bam_males:
	input:
		"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{sample}.GRCh38_minusYPARs.sorted.bam"
	output:
		"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{sample}.GRCh38_minusYPARs.sorted.bam.bai"
	params:
		samtools = samtools_path
	#priority: 90
	shell:
		"{params.samtools} index {input}"


rule index_bam_females:
	input:
		"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{sample}.GRCh38_Ymasked.sorted.bam"
	output:
		"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{sample}.GRCh38_Ymasked.sorted.bam.bai"
	params:
		samtools = samtools_path
	#priority: 85
	shell:
		"{params.samtools} index {input}"


rule merge_bams_males:
	input:
		bams = lambda wildcards: expand(
			"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{sample_males}.GRCh38_minusYPARs.sorted.bam",
			sample_males=config["samples_males"][wildcards.sample_males]),
		bais = lambda wildcards: expand(
			"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{sample_males}.GRCh38_minusYPARs.sorted.bam.bai",
			sample_males=config["samples_males"][wildcards.sample_males])
	output:
		"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{sample_males}.GRCh38_minusYPARs.sorted.merged.bam"
	threads: 4
	params:
		sambamba = sambamba_path,
		samtools = samtools_path,
		threads = 4
	#priority: 80
	shell:
		"{params.sambamba} merge -t {params.threads} {output} {input.bams}"


rule merge_bams_females:
	input:
		bams = lambda wildcards: expand(
			"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{sample_females}.GRCh38_Ymasked.sorted.bam",
			sample_females=config["samples_females"][wildcards.sample_females]),
		bais = lambda wildcards: expand(
			"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{sample_females}.GRCh38_Ymasked.sorted.bam.bai",
			sample_females=config["samples_females"][wildcards.sample_females])
	output:
		"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{sample_females}.GRCh38_Ymasked.sorted.merged.bam"
	threads: 4
	params:
		sambamba = sambamba_path,
		threads = 4
	#priority: 75
	shell:
		"{params.sambamba} merge -t {params.threads} {output} {input.bams}"



rule index_merged_bam_males:
	input:
		"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{sample_males}.GRCh38_minusYPARs.sorted.merged.bam"
	output:
		"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{sample_males}.GRCh38_minusYPARs.sorted.merged.bam.bai"
	params:
		samtools = samtools_path
	#priority: 70
	shell:
		"{params.samtools} index {input}"


rule index_merged_bam_females:
	input:
		"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{sample_females}.GRCh38_Ymasked.sorted.merged.bam"
	output:
		"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{sample_females}.GRCh38_Ymasked.sorted.merged.bam.bai"
	params:
		samtools = samtools_path
	#priority: 65
	shell:
		"{params.samtools} index {input}"


rule bam_stats:
	input:
		bam = "/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{sample}.{genome}.sorted.merged.bam",
		bai = "/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{sample}.{genome}.sorted.merged.bam.bai"
	output:
		"/scratch/amtarave/Kenya_agave/whole_genome/stats/{sample}.{genome}.sorted.bam.stats"
	params:
		samtools = samtools_path
	#priority: 60
	shell:
		"{params.samtools} stats {input.bam} | grep ^SN | cut -f 2- > {output}"


rule picard_mkdups:
	input:
		bam = "/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{sample}.{genome}.sorted.merged.bam",
		bai = "/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{sample}.{genome}.sorted.merged.bam.bai"
	output:
		bam = "/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam",
		metrics = "/scratch/amtarave/Kenya_agave/whole_genome/stats/{sample}.{genome}.picard_mkdup_metrics.txt"
	threads: 4
	params:
		picard = picard_path
	#priority: 55
	shell:
		"{params.picard} -Xmx14g MarkDuplicates I={input.bam} O={output.bam} "
		"M={output.metrics}"


rule index_mkdup_bam:
	input:
		"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam"
	output:
		"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam.bai"
	params:
		samtools = samtools_path
	#priority: 50
	shell:
		"{params.samtools} index {input}"


rule mkdup_bam_stats:
	input:
		bam = "/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam",
		bai = "/scratch/amtarave/Kenya_agave/whole_genome/processed_bams/{sample}.{genome}.sorted.merged.mkdup.bam.bai"
	output:
		"/scratch/amtarave/Kenya_agave/whole_genome/stats/{sample}.{genome}.sorted.mkdup.bam.stats"
	params:
		samtools = samtools_path
	#priority: 45
	shell:
		"{params.samtools} stats {input.bam} | grep ^SN | cut -f 2- > {output}"


#-------------------------------------------------------------------------------
# New rules here for generating bams for all samples mapped to the default
# reference
# Map all trimmed paired end fastq file to the default reference genome
rule map_to_default_ref:
	input:
		fq1 = "/scratch/amtarave/Kenya_agave/whole_genome/trimmed_fastqs/{unique_identifier}_trimmed_R1.fastq.gz",
		fq2 = "/scratch/amtarave/Kenya_agave/whole_genome/trimmed_fastqs/{unique_identifier}_trimmed_R2.fastq.gz",
		new = config["Ref_GRCh38"],
		fai = expand("refs/{ref_type}.fa.fai", ref_type=["Ref_GRCh38","Ref_GRCh38_Y_HardMasked","Ref_GRCh38_Y_PARsMasked"]),
		dict = expand("refs/{ref_type}.dict", ref_type=["Ref_GRCh38","Ref_GRCh38_Y_HardMasked","Ref_GRCh38_Y_PARsMasked"]),
		amb = expand("refs/{ref_type}.fa.amb", ref_type=["Ref_GRCh38","Ref_GRCh38_Y_HardMasked","Ref_GRCh38_Y_PARsMasked"])
	output:
		"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams_default_ref/{unique_identifier}.GRCh38_default.sorted.bam"
	params:
		id = lambda wildcards: config[wildcards.unique_identifier]["ID"],
		sm = lambda wildcards: config[wildcards.unique_identifier]["SM"],
		lb = lambda wildcards: config[wildcards.unique_identifier]["LB"],
		pu = lambda wildcards: config[wildcards.unique_identifier]["PU"],
		pl = lambda wildcards: config[wildcards.unique_identifier]["PL"],
		bwa = bwa_path,
		samtools = samtools_path,
		threads = 4
	threads: 4
	shell:
		" {params.bwa} mem -t {params.threads} -R "
	 	"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
		"{input.new} {input.fq1} {input.fq2}"
		" | {params.samtools} fixmate -O bam - - | {params.samtools} sort "
		"-O bam -o {output}"

# Index these bam files
rule index_default_mapped_bams:
	input:
		"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams_default_ref/{unique_identifier}.GRCh38_default.sorted.bam"
	output:
		"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams_default_ref/{unique_identifier}.GRCh38_default.sorted.bam.bai"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} index {input}"

# Merge bam files where needed. Some samples have multiple bam files and they
# were sequenced on different lanes, so at this point they need to be merged.
rule merge_default_mapped_bams:
	input:
		bams = lambda wildcards: expand(
			"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams_default_ref/{samples}.GRCh38_default.sorted.bam",
			samples=config["samples"][wildcards.samples]),
		bais = lambda wildcards: expand(
			"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams_default_ref/{samples}.GRCh38_default.sorted.bam.bai",
			samples=config["samples"][wildcards.samples])
	output:
		"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams_default_ref/{samples}.GRCh38_default.sorted.merged.bam"
	threads: 4
	params:
		sambamba = sambamba_path,
		samtools = samtools_path,
		threads = 4
	shell:
		"{params.sambamba} merge -t {params.threads} {output} {input.bams}"

# Index merged bam files
rule index_default_mapped_merged_bams:
	input:
		"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams_default_ref/{samples}.GRCh38_default.sorted.merged.bam"
	output:
		"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams_default_ref/{samples}.GRCh38_default.sorted.merged.bam.bai"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} index {input}"

# Get bam stats
rule bam_stats_default_mapped:
	input:
		bam = "/scratch/amtarave/Kenya_agave/whole_genome/processed_bams_default_ref/{samples}.GRCh38_default.sorted.merged.bam",
		bai = "/scratch/amtarave/Kenya_agave/whole_genome/processed_bams_default_ref/{samples}.GRCh38_default.sorted.merged.bam.bai"
	output:
		"/scratch/amtarave/Kenya_agave/whole_genome/stats_bams_default_ref/{samples}.GRCh38_default.sorted.merged.bam.stats"
	params:
		samtools = samtools_path
	#priority: 60
	shell:
		"{params.samtools} stats {input.bam} | grep ^SN | cut -f 2- > {output}"

# mark duplicates
rule picard_mkdups_default_mapped:
	input:
		bam = "/scratch/amtarave/Kenya_agave/whole_genome/processed_bams_default_ref/{samples}.GRCh38_default.sorted.merged.bam",
		bai = "/scratch/amtarave/Kenya_agave/whole_genome/processed_bams_default_ref/{samples}.GRCh38_default.sorted.merged.bam.bai"
	output:
		bam = "/scratch/amtarave/Kenya_agave/whole_genome/processed_bams_default_ref/{samples}.GRCh38_default.sorted.merged.mkdup.bam",
		metrics = "/scratch/amtarave/Kenya_agave/whole_genome/stats_bams_default_ref/{samples}.GRCh38_default.sorted.merged.picard_mkdup_metrics.txt"
	threads: 4
	params:
		picard = picard_path
	shell:
		"{params.picard} -Xmx14g MarkDuplicates I={input.bam} O={output.bam} "
		"M={output.metrics}"

# Get index
rule index_mkdup_bam_default_mapped:
	input:
		"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams_default_ref/{samples}.GRCh38_default.sorted.merged.mkdup.bam"
	output:
		"/scratch/amtarave/Kenya_agave/whole_genome/processed_bams_default_ref/{samples}.GRCh38_default.sorted.merged.mkdup.bam.bai"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} index {input}"

# get stats
rule mkdup_bam_stats_default_mapped:
	input:
		bam = "/scratch/amtarave/Kenya_agave/whole_genome/processed_bams_default_ref/{samples}.GRCh38_default.sorted.merged.mkdup.bam",
		bai = "/scratch/amtarave/Kenya_agave/whole_genome/processed_bams_default_ref/{samples}.GRCh38_default.sorted.merged.mkdup.bam.bai"
	output:
		"/scratch/amtarave/Kenya_agave/whole_genome/stats_bams_default_ref/{samples}.GRCh38_default.sorted.merged.mkdup.bam.stats"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} stats {input.bam} | grep ^SN | cut -f 2- > {output}"
