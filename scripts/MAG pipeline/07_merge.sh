#!/bin/bash -l
#
#SBATCH --ntasks 24 #number cores
#SBATCH --mem=650G #memory
#SBATCH -p highmem
#SBATCH -o logs/07_merge.log
#SBATCH -e logs/07_merge.log
#SBATCH -J sgmag_merge_anvio


SAMPFILE=samples.csv
COVDIR=coverage
CPU=24
MIN=1000

PROFDIR=anvio_profiles_1000minlen

module unload miniconda2
module unload anaconda3
module load miniconda3

source activate anvio-7

PREFIX=SG_Chytid

anvi-merge ${PROFDIR}/*profile/PROFILE.db -o ${PROFDIR}/$PREFIX'_SAMPLES_MERGED_1000min' -c $PREFIX.db
	
