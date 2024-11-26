#!/bin/bash -l
#
#SBATCH -n 24 #number cores
#SBATCH -e logs/drep.log
#SBATCH -o logs/drep.log
#SBATCH --mem 100G #memory per node in Gb
#SBATCH -p intel
#SBATCH -J hotvirus_drep



DIR=clustered_vOTUs
VIRUS=$DIR/data/split

conda activate anvio-7.1

dRep dereplicate dRep \
 -g $VIRUS/*.fa \
 --S_algorithm ANImf \
 -sa 0.95 \
 -nc 0.85 \
 -l 10000 \
 -N50W 0 \
 -sizeW 1 \
 --ignoreGenomeQuality \
 --clusterAlg single 
 
