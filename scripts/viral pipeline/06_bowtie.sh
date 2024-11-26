#!/bin/bash -l
#
#SBATCH -n 16 #number cores
#SBATCH -e logs/bowtie.log
#SBATCH -o logs/bowtie.log
#SBATCH --mem 50G #memory per node in Gb
#SBATCH -p intel
#SBATCH -J sg_virus_map


CPU=16
SAMPFILE=data/sra_all.csv
IN=data/filtered
OUT=clustered_vOTUs/bowtie
DIR=clustered_vOTUs

mkdir $OUT

module load bowtie2
module load samtools

gunzip $IN/*/*

# build reference DB 
bowtie2-build $DIR/dRep_clustered.fa drep_votu_db

IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read PREFIX 
do

	if [ -f $IN/$PREFIX/$PREFIX'_1_filtered.fastq' ]; then
		echo "$PREFIX has split fastqs"
		bowtie2 -x drep_votu_db -1 $IN/$PREFIX/$PREFIX'_1_filtered.fastq' -2 $IN/$PREFIX/$PREFIX'_2_filtered.fastq' -S $OUT/$PREFIX.sam --sensitive


	elif [ -f $IN/$PREFIX/$PREFIX'_interleaved_filtered.fastq' ]; then
		#echo "$PREFIX has interleaved fastq"
		#echo "$PREFIX has single end reads"
		bowtie2 -x drep_votu_db  --interleaved $IN/$PREFIX/$PREFIX'_interleaved_filtered.fastq' -S $OUT/$PREFIX.sam --sensitive

	else
		echo "Something went wrong with $PREFIX"
	fi

# make bam files from sam files
samtools view -F 4 -bS $OUT/$PREFIX.sam | samtools sort > $OUT/$PREFIX'_sortedIndexed.bam'

# index the bam file
samtools index $OUT/$PREFIX'_sortedIndexed.bam'

done
