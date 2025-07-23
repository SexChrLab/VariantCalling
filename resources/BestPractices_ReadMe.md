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
```

For an XX individual use: `SCC_ref_XX.fa`
```
bwa mem -t 4 -R '@RG\tID:ID\tSM:SM\tLB:LB\tPL:PL' SCC_ref_XX.fa XX_sample.fq1 XX_sample.fq2 | samtools fixmate -O bam - - | samtools sort -O bam -o XX_sample_sorted.bam
```


## Call variants


## Filter
