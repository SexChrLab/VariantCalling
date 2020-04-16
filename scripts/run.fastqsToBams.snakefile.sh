#!/bin/bash
#SBATCH --job-name=fastqsToBams # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH --mail-type=END,FAIL           # notifications for job done & fail
#SBATCH --mail-user=amtarave@asu.edu # send-to address
#SBATCH -n 2
#SBATCH -t 96:00:00
cd /home/amtarave/projects/Kenya_sequencing/whole_genome/
source activate kenya
snakemake -s fastqsToBams.snakefile -j 30 --latency-wait 60 --cluster-config fastqsToBams.cluster.json --cluster "sbatch -n {fastqsToBams.cluster.n} --nodes 1 -t {fastqsToBams.cluster.t}:00:00 --mail-type=END,FAIL --mail-user=amtarave@asu.edu "

mv slurm* logs_slurm/

source deactivate kenya
