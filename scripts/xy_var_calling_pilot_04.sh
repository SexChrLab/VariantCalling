#!/bin/bash
#SBATCH --job-name=04_xy_var_calling_pilot # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH -n 1
#SBATCH -t 1-00:00:00
#SBATCH -q tempboost
#SBATCH --mail-user=amtarave@asu.edu # send-to address
#SBATCH --mail-type=ALL

conda activate variant_calling_simulations_project

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/scripts/

#snakemake -s xy_var_calling_pilot_04.snakefile --rerun-incomplete -j 31 --cluster "sbatch -n 1 -t 7-00:00:00 -q tempboost --mail-type=END,FAIL --mail-user=amtarave@asu.edu"
snakemake -s xy_var_calling_pilot_04.snakefile --rerun-incomplete -j 1000 --cluster "sbatch -n 1 -t 1-00:00:00 -q tempboost --mail-type=ALL --mail-user=amtarave@asu.edu"

#mv slurm* ../logs/

source deactivate
