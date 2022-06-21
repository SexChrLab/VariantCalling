#!/bin/bash
#SBATCH --job-name=03_xy_var_calling_pilot # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH -n 1
#SBATCH -t 7-00:00:00
#SBATCH -q tempboost
#SBATCH --mail-user=amtarave@asu.edu # send-to address
#SBATCH --mail-type=ALL
conda activate variant_calling_simulations_project

cd /data/CEM/wilsonlab/projects/variant_calling_simulations_project/scripts/

#snakemake -s xy_var_calling_pilot_03.snakefile --rerun-incomplete -j 21 --cluster "sbatch -n 1 -t 7-00:00:00 --mem-per-cpu=16000 -q tempboost --mail-type=END,FAIL --mail-user=amtarave@asu.edu"
#snakemake -s xy_var_calling_pilot_03.snakefile --rerun-incomplete -j 50 --cluster "sbatch -n 2 -t 2-00:00:00 --mem-per-cpu=8000 --mail-type=END,FAIL --mail-user=amtarave@asu.edu"
snakemake -s xy_var_calling_pilot_03.snakefile --rerun-incomplete -j 422 --cluster "sbatch -n 2 -t 7-00:00:00 --mem-per-cpu=16000 -q tempboost --mail-type=ALL --mail-user=amtarave@asu.edu"
#snakemake -s xy_var_calling_pilot_03.snakefile --rerun-incomplete -j 21 --cluster "sbatch -n 2 -t 4-00:00:00 -q tempboost --mail-type=ALL --mail-user=amtarave@asu.edu"

#mv slurm* ../logs/

source deactivate variant_calling_simulations_project
