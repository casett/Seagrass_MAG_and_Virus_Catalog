#!/bin/bash -l
#
#SBATCH --ntasks 32 #number cores
#SBATCH --mem=700G #memory
#SBATCH -p highmem
#SBATCH -o logs/08c_maxbin.log
#SBATCH -e logs/08c_maxbin.log
#SBATCH -J sgmag_maxbin

SAMPFILE=samples.csv
COVDIR=coverage
CPU=48
MIN=1000

PROFDIR=anvio_profiles_1000minlen
PREFIX=SG_Chytid


module unload miniconda2
module unload anaconda3
module load miniconda3
#module load concoct/1.1.0
#module load metabat/0.32.4
module load maxbin/2.2.1 
#module load diamond

source activate anvio-7
#installed maxbin to anvio env


anvi-cluster-contigs -p ${PROFDIR}/$PREFIX'_SAMPLES_MERGED_1000min'/PROFILE.db -c $PREFIX.db -C MAXBIN --driver maxbin2 -T $CPU --just-do-it
	
