#!/bin/bash -l
#
#SBATCH -n 24 #number cores
#SBATCH -p intel
#SBATCH -o logs/16_checkm.log
#SBATCH -e logs/16_checkm.log
#SBATCH -J sgmag_checkM
#SBATCH --mem 96G #memory in Gb
#SBATCH -t 96:00:00 #time in hours:min:sec

module load checkm


BINFOLDER=magpurify_results/cleaned_mags
OUTPUT=magpurify_results/checkm_original_cleaned
CPU=24

checkm lineage_wf -t $CPU -x fa $BINFOLDER $OUTPUT

checkm tree $BINFOLDER -x .fa -t $CPU $OUTPUT/tree

checkm tree_qa $OUTPUT/tree -f $OUTPUT/checkm_tree.txt
