#!/bin/bash -l
#
#SBATCH --ntasks 16 #number cores
#SBATCH -e logs/04_cov.log
#SBATCH -o logs/04_cov.log
#SBATCH -J sgmag_cov_anvio
#SBATCH --mem=148G #memory
#SBATCH -p intel,batch

module unload miniconda2
module unload anaconda3

module load miniconda3

source activate anvio-7

COVDIR=coverage
CPU=16
DIR=data/filt_seqs
SAMPFILE=$DIR/samples.csv
PREFIX=SG_Chytid
ASSEMDIR=SG_Chytid_CoAssem

mkdir $COVDIR

anvi-script-reformat-fasta $ASSEMDIR/final.contigs.fa -o $ASSEMDIR/$PREFIX.fixed.fa -l 0 --simplify-names
 
bowtie2-build $ASSEMDIR/$PREFIX.fixed.fa ${COVDIR}/$PREFIX


IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read SAMPLE
do 
	bowtie2 --threads $CPU -x  ${COVDIR}/$PREFIX -1 ${DIR}/$SAMPLE.ZM_removed_R1.fastq.gz -2 ${DIR}/$SAMPLE.ZM_removed_R2.fastq.gz -S ${COVDIR}/$SAMPLE'.sam'
	samtools view -F 4 -bS ${COVDIR}/$SAMPLE'.sam' > ${COVDIR}/$SAMPLE'-RAW.bam'
	anvi-init-bam ${COVDIR}/$SAMPLE'-RAW.bam' -o ${COVDIR}/$SAMPLE'.bam'

done

anvi-gen-contigs-database -f $ASSEMDIR/$PREFIX.fixed.fa -o $PREFIX.db
anvi-run-hmms -c $PREFIX.db --num-threads $CPU
anvi-get-sequences-for-gene-calls -c $PREFIX.db -o $PREFIX.gene.calls.fa
anvi-get-sequences-for-gene-calls -c $PREFIX.db --get-aa-sequences -o $PREFIX.amino.acid.sequences.fa


