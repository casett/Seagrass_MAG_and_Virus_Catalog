#!/bin/bash -l
#
#SBATCH -n 16 #number cores
#SBATCH -e logs/coverm.log
#SBATCH -o logs/coverm.log
#SBATCH --mem 50G #memory per node in Gb
#SBATCH -p intel
#SBATCH -J sg_virus_coverm


DIR=clustered_vOTUs
BAMDIR=$DIR/bowtie
conda activate sourmash_to_phylo

coverm contig -m trimmed_mean --min-covered-fraction 0.75 -b $BAMDIR/*.bam > $DIR/votu.tmean.tsv
coverm contig -m count -b $BAMDIR/*.bam > $DIR/votu.count.tsv
