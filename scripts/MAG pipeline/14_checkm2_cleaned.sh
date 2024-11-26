#!/bin/bash -l
#
#SBATCH -n 24 #number cores
#SBATCH -p intel
#SBATCH -o logs/14_checkm.log
#SBATCH -e logs/14_checkm.log
#SBATCH -J sgmag_checkM
#SBATCH --mem 96G #memory in Gb
#SBATCH -t 96:00:00 #time in hours:min:sec


conda activate checkm2

BINFOLDER=magpurify_results/cleaned_mags
OUTPUT=magpurify_results/checkm2_cleaned
CPU=24

checkm2 predict --threads $CPU --input $BINFOLDER --output-directory $OUTPUT -x .fa
