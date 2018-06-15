#!/bin/bash
#SBATCH --job-name=data_dump
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -o data_dump_%j.out
#SBATCH -e data_dump_%j.err
#SBATCH --partition=general

export TMPDIR=/home/CAM/$USER/tmp/

module load sratoolkit

fastq-dump SRR3498212
fastq-dump SRR3498213
fastq-dump SRR3498215
fastq-dump SRR3498216

