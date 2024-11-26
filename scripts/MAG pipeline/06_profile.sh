#!/bin/bash -l
#
#SBATCH --ntasks 16 #number cores
#SBATCH --mem=250G #memory
#SBATCH -p intel,batch
#SBATCH -o logs/06_profile.log
#SBATCH -e logs/06_profile.log
#SBATCH -J sgmag_profile_anvio


COVDIR=coverage
CPU=16
MIN=1000
DIR=data/filt_seqs
SAMPFILE=$DIR/samples.csv
PREFIX=SG_Chytid
ASSEMDIR=SG_Chytid_CoAssem



module unload miniconda2
module unload anaconda3
module load miniconda3

OUT=anvio_profiles_1000minlen

mkdir $OUT

source activate anvio-7
IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read SAMPLE
do 
	anvi-profile -i ${COVDIR}/$SAMPLE'.bam' -c $PREFIX.db --num-threads $CPU --min-contig-length $MIN --cluster-contigs --output-dir ${OUT}/$SAMPLE'_profile'
done	

