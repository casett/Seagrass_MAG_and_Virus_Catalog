#!/bin/bash -l
#
#SBATCH -n 16 #number cores
#SBATCH -e logs/genomad.dRep.log
#SBATCH -o logs/genomad.dRep.log
#SBATCH --mem 50G #memory per node in Gb
#SBATCH -p intel
#SBATCH -J sg_virus_genomad

conda activate genomad

IN=clustered_vOTUs
#genomad download-database .
genomad annotate $IN/dRep_clustered.fa genomad_output genomad_db
