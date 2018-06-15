#!/bin/bash
#SBATCH --job-name=quality_control
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -o quality_control_%j.out
#SBATCH -e quality_control_%j.err
#SBATCH --partition=general

export TMPDIR=/home/CAM/$USER/tmp/

module load fastqc
module load MultiQC

fastqc trimmed_SRR3498212.fastq

fastqc trimmed_SRR3498213.fastq

fastqc trimmed_SRR3498215.fastq

fastqc trimmed_SRR3498216.fastq

multiqc -f -n trimmed trimmed*
