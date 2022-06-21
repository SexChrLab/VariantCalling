#!/bin/bash
#SBATCH --job-name=ASN_02_xy_var_calling_pilot # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH -n 1
#SBATCH -t 7-00:00:00
#SBATCH -q tempboost
#SBATCH --mail-user=amtarave@asu.edu # send-to address
#SBATCH --mail-type=ALL

conda activate NEAT_env
cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/scripts/ASN/

snakemake -s ASN_xy_var_calling_pilot_02.snakefile --rerun-incomplete -j 41 --cluster "sbatch -n 2 -t 5-00:00:00 --mem-per-cpu=8000 -q tempboost --mail-type=ALL --mail-user=amtarave@asu.edu"
#snakemake -s ASN_xy_var_calling_pilot_02.snakefile --rerun-incomplete -j 55 --cluster "sbatch -n 2 --mem-per-cpu=8000 -t 4-00:00:00 -q tempboost --mail-type=ALL --mail-user=amtarave@asu.edu"

mv slurm* logs/

source deactivate
