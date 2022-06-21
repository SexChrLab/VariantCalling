#!/bin/bash
#SBATCH --job-name=download_1000G_data # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH -n 1
#SBATCH -q tempboost
#SBATCH -t 7-00:00:00
#SBATCH --mail-user=amtarave@asu.edu # send-to address

conda activate variant_calling_simulations_project

#PERL5LIB=""

cd /scratch/amtarave/variant_calling_simulations_project/pilot/data/all_chrs/

snakemake -s download_1000G_data.snakefile --rerun-incomplete -j 48 --cluster "sbatch -n 2 -t 7-00:00:00 -q tempboost --mail-type=END,FAIL --mail-user=amtarave@asu.edu"

source deactivate variant_calling_simulations_project
