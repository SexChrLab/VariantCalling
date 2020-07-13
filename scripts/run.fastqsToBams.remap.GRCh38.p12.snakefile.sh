#!/bin/bash
#SBATCH --job-name=fastqsToBamsRemapGRCh38p12 # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH --mail-type=END,FAIL           # notifications for job done & fail
#SBATCH --mail-user=amtarave@asu.edu # send-to address
#SBATCH -n 2
#SBATCH -q tempboost
#SBATCH -t 96:00:00
cd /home/amtarave/projects/Kenya_sequencing/whole_genome/
conda activate kenya
snakemake -s fastqsToBams.remap.GRCh38.p12.snakefile -j 50 --latency-wait 60 --cluster-config cluster.json --cluster "sbatch -n {cluster.n} --nodes 1 -t {cluster.t}:00:00 --mail-type=END,FAIL --mail-user=amtarave@asu.edu "

mv slurm* logs_slurm/

conda deactivate
