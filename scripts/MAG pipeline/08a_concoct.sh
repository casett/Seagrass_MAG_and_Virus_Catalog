#!/bin/bash -l
#
#SBATCH --ntasks 32 #number cores
#SBATCH --mem=250G #memory
#SBATCH -p intel,batch
#SBATCH -o logs/08a_concoct.log
#SBATCH -e logs/08a_concoct.log
#SBATCH -J sgmag_concoct

SAMPFILE=samples.csv
COVDIR=coverage
CPU=32
MIN=1000

PROFDIR=anvio_profiles_1000minlen
PREFIX=SG_Chytid


module unload miniconda2
module unload anaconda3
module load miniconda3
module load concoct/1.1.0

source activate anvio-7
#must also install concoct to conda anvio env

anvi-cluster-contigs -p ${PROFDIR}/$PREFIX'_SAMPLES_MERGED_1000min'/PROFILE.db -c $PREFIX.db -C CONCOCT --driver concoct -T $CPU --just-do-it
	
