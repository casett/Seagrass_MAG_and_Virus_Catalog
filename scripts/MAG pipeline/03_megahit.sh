#!/bin/bash -l
#
#SBATCH -e logs/megahit.log
#SBATCH -o logs/megahit.log
#SBATCH --nodes 1 --ntasks 24 -p highmem
#SBATCH -J SGMAG_megahit --mem 600G


module load megahit
CPU=24
DIR=data/filt_seqs
SAMPLEFILE=$DIR/samples.csv
OUT=SG_Chytid_CoAssem

megahit -1 $DIR/CE_108A_S70.ZM_removed_R1.fastq.gz,$DIR/CE_109A_S71.ZM_removed_R1.fastq.gz,$DIR/CE_140A_S72.ZM_removed_R1.fastq.gz -2 $DIR/CE_108A_S70.ZM_removed_R2.fastq.gz,$DIR/CE_109A_S71.ZM_removed_R2.fastq.gz,$DIR/CE_140A_S72.ZM_removed_R2.fastq.gz -o $OUT -t $CPU
