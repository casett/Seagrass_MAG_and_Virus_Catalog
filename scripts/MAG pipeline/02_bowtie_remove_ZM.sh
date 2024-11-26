#!/bin/bash
#SBATCH --nodes 1 --ntasks 16 --mem 90G --out logs/bowtie.%a.log -p intel,batch 
#SBATCH -J SGMAG_remove_ZM

module load bowtie2
module load samtools


DB=ZM_ref
ZM_ASM=Zostera_marina.combined.fasta
DIR=data/filt_seqs
SAMPLEFILE=$DIR/samples.csv
CPU=16

bowtie2-build $DB/$ZM_ASM $DB


IFS=,
tail -n +2 $SAMPLEFILE | sed -n ${N}p | while read SAMPLE 
do

	bowtie2 -p $CPU -x $DB -1 $DIR/${SAMPLE}_R1_filtered.fastq.gz -2 $DIR/${SAMPLE}_R2_filtered.fastq.gz -S $DIR/$SAMPLE.mapped_and_unmapped.sam

	samtools view -bS $DIR/$SAMPLE.mapped_and_unmapped.sam > $DIR/$SAMPLE.mapped_and_unmapped.bam

	samtools view -b -f 12 -F 256 $DIR/$SAMPLE.mapped_and_unmapped.bam > $DIR/$SAMPLE.bothReadsUnmapped.bam 

	samtools sort -n -@ $CPU $DIR/$SAMPLE.bothReadsUnmapped.bam -o $DIR/$SAMPLE.bothReadsUnmapped_sorted.bam

	samtools fastq -@ $CPU $DIR/$SAMPLE.bothReadsUnmapped_sorted.bam -1 $DIR/$SAMPLE.ZM_removed_R1.fastq.gz -2 $DIR/$SAMPLE.ZM_removed_R2.fastq.gz -0 /dev/null -s /dev/null -n

done

