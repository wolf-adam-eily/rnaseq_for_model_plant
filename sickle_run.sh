#!/bin/bash
#SBATCH --job-name=sickle_run
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -o sickle_run_%j.out
#SBATCH -e sickle_run_%j.err
#SBATCH --partition=general

export TMPDIR=/home/CAM/$USER/tmp/

module load sickle

sickle se -f SRR3498212.fastq -t sanger -o trimmed_SRR3498212.fastq -q 35 -l 45

sickle se -f SRR3498213.fastq -t sanger -o trimmed_SRR3498213.fastq -q 35 -l 45

sickle se -f SRR3498215.fastq -t sanger -o trimmed_SRR3498215.fastq -q 35 -l 45

sickle se -f SRR3498216.fastq -t sanger -o trimmed_SRR3498216.fastq -q 35 -l 45
