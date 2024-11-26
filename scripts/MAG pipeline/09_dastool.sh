#!/bin/bash -l
#
#SBATCH --ntasks 24 #number cores
#SBATCH --mem=98G #memory
#SBATCH -p intel,batch
#SBATCH -o logs/09_dastool.log
#SBATCH -e logs/09_dastool.log
#SBATCH -J sgmag_das_anvio


SAMPFILE=samples.csv
COVDIR=coverage
CPU=24
MIN=1000

PROFDIR=anvio_profiles_1000minlen
PREFIX=SG_Chytid


module unload miniconda2
module unload anaconda3
module load miniconda3
#module load concoct/1.1.0
#module load metabat/0.32.4
#module load maxbin/2.2.1 

module load diamond

source activate anvio-7
#made sure to install dastool via conda to anvio-7 env

#MAXBIN missing some dependencies so failed
#BINSanity also fail


anvi-cluster-contigs -p ${PROFDIR}/$PREFIX'_SAMPLES_MERGED_1000min'/PROFILE.db -c $PREFIX.db  -C DASTOOL --driver dastool -S CONCOCT,METABAT -T $CPU --search-engine diamond --just-do-it
	
