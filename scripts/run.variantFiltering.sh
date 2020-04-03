#!/bin/bash
#SBATCH --job-name=run.variantFiltering # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH -n 1
#SBATCH -t 96:00:00
#SBATCH --mail-user=amtarave@asu.edu # send-to address

source activate varCalling

cd /home/amtarave/projects/variant_calling/scripts

snakemake -s variantFiltering.snakefile --dag | dot -Tpdf > variantFiltering.snakefile.pdf
snakemake -s variantFiltering.snakefile -j 20 --cluster "sbatch -n 2 -t 96:00:00 --mail-type=END,FAIL --mail-user=amtarave@asu.edu"

mv slurm* logs/

source deactivate varCalling
