# Best practices for variant calling on X and Y chromosomes

This readme describes best practices for calling variants on the X and Y chromosome as outlined in Taravella Oill et al. 2025.


After FASTQ visualization and filtering the following steps should be taken.

1. Alignment
2. Call variants 
3. Filter

## Alignment 
Alignment should be performed using a reference genome informed on the sex chromosome complement (SCC) of the individual. For XY individuals, this means the pseudoautosomal regions (PARs) on the Y chromosome are hard masked in the reference genome fasta file used for alignment. For XX individuals, this means the entire Y chromosome being hard masked in the reference genome fasta. Hard masked means the sequence is replaced with `N`s.

Here is an example for how to generate these reference genome fasta files.

```
Add code here
```

Then you may proceed to alignment. Here we used `bwa mem`. Some parameters may need to be tuned to your analysis (like threads or what you want to name your read groups).

For an XY individual use: `SCC_ref_XY.fa`
```
bwa mem -t 4 -R '@RG\tID:ID\tSM:SM\tLB:LB\tPL:PL' SCC_ref_XY.fa XY_sample.fq1 XY_sample.fq2 | samtools fixmate -O bam - - | samtools sort -O bam -o XY_sample_sorted.bam

# dont forget to index
samtools index XY_sample_sorted.bam
# and mark duplicates
picard -Xmx14g MarkDuplicates I=XY_sample_sorted.bam O=XY_sample_sorted_mkdups.bam M=XY_sample_sorted_mkdups_metrics.txt VALIDATION_STRINGENCY=LENIENT
# and index again
samtools index XY_sample_sorted_mkdups.bam
```

For an XX individual use: `SCC_ref_XX.fa`
```
bwa mem -t 4 -R '@RG\tID:ID\tSM:SM\tLB:LB\tPL:PL' SCC_ref_XX.fa XX_sample.fq1 XX_sample.fq2 | samtools fixmate -O bam - - | samtools sort -O bam -o XX_sample_sorted.bam

# dont forget to index
samtools index XX_sample_sorted.bam
# and mark duplicates
picard -Xmx14g MarkDuplicates I=XX_sample_sorted.bam O=XX_sample_sorted_mkdups.bam M=XY_sample_sorted_mkdups_metrics.txt VALIDATION_STRINGENCY=LENIENT
# and index again
samtools index XX_sample_sorted_mkdups.bam
```


## Call variants
We used GATK HaplotypeCaller to call variants and GenotypeGVCFs to joint genotype. We reccomend calling variants on the PARs and nonPARs separately. For XY samples, use ploidy of 1 on the nonPARs for X and Y. For XX samples, use ploidy of 2 on the nonPARs. PARs will always be called as diploid. 

For XY individuals, here is example code:
```
# 1. Non-PARs chromosome X
# Call variants per sample with HaplotypeCaller
gatk --java-options '-Xmx4g' HaplotypeCaller -R SCC_ref_XX.fa -I XY_sample_sorted_mkdups.bam -L chrY -XL YPARs.interval.list -ploidy 1 -O XY_sample_chrY_nonPARs_haploid.g.vcf.gz -ERC GVCF

# After all sample's gvcfs have been generated combine all samples

# Then joint genotype
```

For XX individuals, here is example code:
```

```

## Filter
