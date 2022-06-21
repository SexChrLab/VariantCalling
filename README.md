# VariantCalling

Assessment of mapping, variants calling and filtering strategies on human sex chromosomes through simulations


## Create conda environmnet

All necessary software for this environment are located here: `variant_calling_simulations_project.yml`.

To create the environment:
```
conda env create -n variant_calling_simulations_project --file variant_calling_simulations_project.yml


# To activate this environment, use
#
#     $ conda activate variant_calling_simulations_project
#
# To deactivate an active environment, use
#
#     $ conda deactivate

```

## Description of analyses and accompanying folders, files, and scripts

### 1) Coordinates of genomic features on X and Y
The folder `resources/genomic_features` contains coordinates for the different genomic features on X and Y (PARs, XTR, and amplicons). Readme can be found there: `resources/genomic_features/README.md`


### 2) Simulate sequence data
We simulated paired-end sequence reads using NExt-generation sequencing Analysis Toolkit (NEAT) software. Variants were inserted from 20 males and 20 females each from European, African, and East Asian populations (total of 120 individuals) from high coverage variant calls from 1000 Genomes.

#### 2.1) Prep 1000 genomes VCFs
Prior to running NEAT we downloaded high coverage variant calls from 1000 Genomes. Scripts for downloading and prepping individual VCFs can be found here: `scripts/`
1. `download_1000G_data.snakefile`: Snakefile to download 1000 Genomes high coverage VCFs per chromosome.
2. `download_1000G_data.config.json`: Accompanying json file
3. `download_1000G_data.sh`: Accompanying shell file to submit snakejobs to the cluster.

We next extracted the samples that were used to insert variants during the simulation step.
1. `xy_var_calling_pilot_01.snakefile`: Snakefile to extract samples.
2. `xy_var_calling_pilot_01.config.json`: Accompanying json file
3. `xy_var_calling_pilot_01.sh`: Accompanying shell file to submit snakejobs to the cluster.

These scripts above were for the European samples. This was also done with Asian and African samples and similar script can be found here: `scripts/ASN/` and `scripts/AFR/`.

#### 2.2) Simulate sequence reads
For simulated females the X chromosome was simulated to 20x coverage, while for simulated male samples the X and Y chromosomes were simulated to 10x coverage each. For each sample, NEAT outputs FASTQs along with a truth VCF. All chromomes were simulated separately and then merged into one FASTQ file prior to data processing.
1. `xy_var_calling_pilot_02.snakefile`: Snakefile to simulate sequence reads.
2. `xy_var_calling_pilot_02.config.json`: Accompanying json file
3. `xy_var_calling_pilot_02.sh`: Accompanying shell file to submit snakejobs to the cluster.

These scripts above were for the European samples. This was also done with Asian and African samples and similar script can be found here: `scripts/ASN/` and `scripts/AFR/`.


### 3) Process data
This includes trimming, alignment, variant calling, and variant filtering.
These scripts include alignment to the default and sex chromosome complement reference genomes. For males, the sex chromosomes were also called as haploid and diploid. For filtering, GATK recommended hard filters were implemented along with testing different filtering thresholds.
1. `xy_var_calling_pilot_03.snakefile`/ `xy_var_calling_pilot_03_haploid.snakefile`: Snakefile to trim, align to both default and sex chromosome complement informed reference genomes, calling variants, and filtering variants.
2. `xy_var_calling_pilot_03.config.json`/ `xy_var_calling_pilot_03_haploid.config.json`: Accompanying json file
3. `xy_var_calling_pilot_03.sh`/ `xy_var_calling_pilot_03_haploid.sh`: Accompanying shell file to submit snakejobs to the cluster.

These scripts above were for the European 10 simulated males and females. This was also done for 20 males and 20 female European, Asian and African samples. For the European samples scripts have `*_20_samples*` in the file names and these scripts for the Asian and African samples can be found here: `scripts/ASN/` and `scripts/AFR/`.

For a subset of male samples we aligned data masking the Y-linked XTR. Scripts have the prefix: `xy_var_calling_pilot_03_XTR.*`.

### 4) Analyze data
#### 4.1) Calculate performance metrics
For each sample and approach we compared the truth VCFs to the called VCFs. A true positive (TP) was defined as a SNP that was both called and simulated and where the reference and alternative allele matched between the golden and called VCFs, a false positive (FP) was defined as a SNP called but not simulated, and a false negative (FN) was defined as a SNP that was simulated but not called. This was done using the following script: `scripts/compare_VCFs_fix.py`.

1. `xy_var_calling_pilot_04.snakefile` / `xy_var_calling_pilot_04_haploid_20_samples.snakefile`: Snakefile to calculate performance metrics.
2. `xy_var_calling_pilot_04.config.json`: Accompanying json file
3. `xy_var_calling_pilot_04.sh`: Accompanying shell file to submit snakejobs to the cluster.

These scripts above were for the European 10 simulated males and females. This was also done for 20 males and 20 female European, Asian and African samples. For the European samples scripts have `*_20_samples*` in the file names and these scripts for the Asian and African samples can be found here: `scripts/ASN/` and `scripts/AFR/`.

For a subset of male samples we aligned data masking the Y-linked XTR. Scripts have the prefix: `xy_var_calling_pilot_03_XTR.*`.

For some of the analyses we compared the proportions of FPs and FNs to the total number of simulated variants. Script to get the total number of simulated variants per sample and for different genomic regions: `scripts/simulated_variants_by_region.sh`.

Tables for plotting were generated with the following scripts: `scripts/make_results_tables.sh` and `scripts/make_results_tables_filters.sh`.

#### 4.2) Window analysis
For some of the analyses, we looked at the performance metrics in 50 kb windows across the sex chromosomes. Script to do this can be found here: `scripts/metrics_per_window.sh`.

#### 4.3) By region
The number of TPs, FPs, and FNs were calculated across different regions of X and Y (PARs, XTR, amplicons, non-PARs minus amplicons and XTR). Script to get these values can be found here: `scripts/performance_metrics_by_region.sh`

#### 4.4) Calculate uncorrected diversity

##### Get putatively neutral bed files
For this I need to make a bed file with neutral and high quality sites. Neutral sites were obtained from UCSC table browser:

1. **Genic regions:** clade: Mammal; genome: Human; assembly: Dec. 2013 (GRCh38/hg38); group:
Genes and Gene Predictions; track: GENCODE v39; table: knownGene; region: genome; output
format: BED - browser extensible data; output file: genome.human.GRCh38.GenesandGenePredictions.GENCODEv39.knownGene.bed;
file type returned: gzip compressed; Create one BED record per: Whole Gene

2. **Conserved regions:** clade: Mammal; genome: Human; assembly: Dec. 2013 (GRCh38/hg38);
group: Comparative Genomics; track: Conservation; table: 100 Vert. El (phastConsElements100way); region: genome; output format: BED - browser extensible data; output file:
genome.human.GRCh38.ComparativeGenomics.Conservation.100VertEl.bed; file type returned: gzip compressed

3. **Repeats:** clade: Mammal; genome: Human; assembly: Dec. 2013 (GRCh38/hg38); group: Repeats; track: RepeatMasker; table: msk; region: genome; output format: BED - browser extensible data; output file: genome.human.GRCh38.Repeats.RepearMasker.msk.bed; file type returned: gzip compressed

4. **CpG islands:** clade: Mammal; genome: Human; assembly: Dec. 2013 (GRCh38/hg38); group:
regulation; track: CpG islands; table: cpgIslandExt; region: genome; output format: BED -
browser extensible data; output file: genome.human.GRCh38.Regulation.CpGIslands.cpgIslandExt.bed; file type returned: gzip compressed


```
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/resources/UCSC_genome_browser/for_diversity_calculation/

# Unzip
gunzip genome.human.GRCh38.*

# Sort
sort -k1,1 -k2,2n genome.human.GRCh38.ComparativeGenomics.Conservation.100VertEl.bed > genome.human.GRCh38.ComparativeGenomics.Conservation.100VertEl.sorted.bed
sort -k1,1 -k2,2n genome.human.GRCh38.GenesandGenePredictions.GENCODEv39.knownGene.bed > genome.human.GRCh38.GenesandGenePredictions.GENCODEv39.knownGene.sorted.bed
sort -k1,1 -k2,2n genome.human.GRCh38.Regulation.CpGIslands.cpgIslandExt.bed > genome.human.GRCh38.Regulation.CpGIslands.cpgIslandExt.sorted.bed
sort -k1,1 -k2,2n genome.human.GRCh38.Repeats.RepearMasker.msk.bed > genome.human.GRCh38.Repeats.RepearMasker.msk.bed.sorted.bed

# Make bed file with the chromomes for analysis (all regions)
echo -e "chr1\t0\t248956422" >> all.chromosomes.bed
echo -e "chr2\t0\t242193529" >> all.chromosomes.bed
echo -e "chr3\t0\t198295559" >> all.chromosomes.bed
echo -e "chr4\t0\t190214555" >> all.chromosomes.bed
echo -e "chr5\t0\t181538259" >> all.chromosomes.bed
echo -e "chr6\t0\t170805979" >> all.chromosomes.bed
echo -e "chr7\t0\t159345973" >> all.chromosomes.bed
echo -e "chr8\t0\t145138636" >> all.chromosomes.bed
echo -e "chr9\t0\t138394717" >> all.chromosomes.bed
echo -e "chr10\t0\t133797422" >> all.chromosomes.bed
echo -e "chr11\t0\t135086622" >> all.chromosomes.bed
echo -e "chr12\t0\t133275309" >> all.chromosomes.bed
echo -e "chr13\t0\t114364328" >> all.chromosomes.bed
echo -e "chr14\t0\t107043718" >> all.chromosomes.bed
echo -e "chr15\t0\t101991189" >> all.chromosomes.bed
echo -e "chr16\t0\t90338345" >> all.chromosomes.bed
echo -e "chr17\t0\t83257441" >> all.chromosomes.bed
echo -e "chr18\t0\t80373285" >> all.chromosomes.bed
echo -e "chr19\t0\t58617616" >> all.chromosomes.bed
echo -e "chr20\t0\t64444167" >> all.chromosomes.bed
echo -e "chr21\t0\t46709983" >> all.chromosomes.bed
echo -e "chr22\t0\t50818468" >> all.chromosomes.bed
echo -e "chrX\t0\t156040895" >> all.chromosomes.bed
echo -e "chrY\t0\t57227415" >> all.chromosomes.bed

# Merge (non-neutral)
bedtools merge -i genome.human.GRCh38.ComparativeGenomics.Conservation.100VertEl.sorted.bed > genome.human.GRCh38.ComparativeGenomics.Conservation.100VertEl.sorted.merge.bed
bedtools merge -i genome.human.GRCh38.GenesandGenePredictions.GENCODEv39.knownGene.sorted.bed > genome.human.GRCh38.GenesandGenePredictions.GENCODEv39.knownGene.sorted.merge.bed
bedtools merge -i genome.human.GRCh38.Regulation.CpGIslands.cpgIslandExt.sorted.bed > genome.human.GRCh38.Regulation.CpGIslands.cpgIslandExt.sorted.merge.bed
bedtools merge -i genome.human.GRCh38.Repeats.RepearMasker.msk.bed.sorted.bed > genome.human.GRCh38.Repeats.RepearMasker.msk.bed.sorted.merge.bed


# Subtract (all minus non-neutral to get putatively neutral)
bedtools subtract -a all.chromosomes.bed -b genome.human.GRCh38.ComparativeGenomics.Conservation.100VertEl.sorted.merge.bed genome.human.GRCh38.GenesandGenePredictions.GENCODEv39.knownGene.sorted.merge.bed genome.human.GRCh38.Regulation.CpGIslands.cpgIslandExt.sorted.merge.bed genome.human.GRCh38.Repeats.RepearMasker.msk.bed.sorted.merge.bed > all.chromosomes.putatively.neutral.bed
```

##### Get VCF statistics
To define high quality sites I will get info from the VCFs

Script to extract stats from VCFs from: https://github.com/tanyaphung/vcfhelper

Script to run analysis: `scripts/get_vcf_stats.sh`

##### Get high quality and neutral sites bed
High quality sites will be defined for autosomes as:
- AN at least 20 (total number of genotypes for 10 diploid samples)
- DP between 67 and 201 (50% and 150% mean DP across all autosome sites. Mean DP across all autosome sites is 134)

and for X chromosome as:
- AN at least 10 (total number of genotypes for 10 haploid samples)
- DP between 23 and 69 (50% and 150% mean DP across all X chromosome (nonPAR) sites. Mean DP across all X sites is 46)

This was done for default and sex chromosome complement aligned data.

Scripts for this analysis:
```
# To get high quality sites for each set of vcfs
scripts/get_high_quality_sites_all_sites_02.sh

# To make high quality and neutral sites bed file
scripts/get_high_quality_neutral_beds_all_sites_02.sh
```

##### Calculate pi per site
Script to calculate pi per site from: https://github.com/tanyaphung/popgen_tools
See: `scripts/calculate_pi_all_sites_02.sh`


### 5) Visualize results
Scripts to make the following figures can be found here: `figures/`

1. **Figure 1:** Comparison of true positives across X and Y aligning sequences to the default and sex chromosome complement reference genomes.
2. **Figure 2:** Number of called variants to number of simulated variants in XTR when aligning sequences to the sex chromosome complement reference genomes and to a reference genome with Y-link XTR masked.
3. **Figure 3:** Comparison of false positives and false negatives on X and Y in males with diploid and haploid calling.
4. **Figure 4:** Proportion of false positives and false negatives across different genomic features on X and Y.
5. **Figure 5:** Number of true positives and false positives for different filtering thresholds in males.
6. **Figure 6:** Proportion of false positives and false negatives when joint genotyping 10 and 20 samples.
7. **Figure 7:** Proportion of false positives and false negatives when across African, Asian, and European samples.
