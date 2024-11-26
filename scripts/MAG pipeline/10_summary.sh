#!/bin/bash -l
#
#SBATCH --ntasks 24 #number cores
#SBATCH --mem=98G #memory
#SBATCH -p intel,batch
#SBATCH -o logs/10_summary.log
#SBATCH -e logs/10_summary.log
#SBATCH -J sgmag_summary

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


source activate anvio-7


#anvi-summarize -p ${PROFDIR}/$PREFIX'_SAMPLES_MERGED_1000min'/PROFILE.db -c $PREFIX.db -o ${PROFDIR}/$PREFIX'_SAMPLES_MERGED'/sample_summary_DASTOOL -C DASTOOL

#anvi-summarize -p ${PROFDIR}/$PREFIX'_SAMPLES_MERGED_1000min'/PROFILE.db -c $PREFIX.db -o ${PROFDIR}/$PREFIX'_SAMPLES_MERGED'/sample_summary_METABAT -C METABAT

anvi-summarize -p ${PROFDIR}/$PREFIX'_SAMPLES_MERGED_1000min'/PROFILE.db -c $PREFIX.db -o ${PROFDIR}/$PREFIX'_SAMPLES_MERGED'/sample_summary_CONCOCT -C CONCOCT

