#!/bin/bash
#SBATCH --job-name=sam_sort_bam
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -o sam_sort_bam_%j.out
#SBATCH -e sam_sort_bam_%j.err
#SBATCH --partition=general

export TMPDIR=/home/CAM/$USER/tmp/

module load samtools

samtools sort -@ 16 -o athaliana_root_1.bam athaliana_root_1.sam

samtools sort -@ 16 -o athaliana_root_2.bam athaliana_root_2.sam

samtools sort -@ 16 -o athaliana_shoot_1.bam athaliana_shoot_1.sam

samtools sort -@ 16 -o athaliana_shoot_2.bam athaliana_shoot_2.sam
