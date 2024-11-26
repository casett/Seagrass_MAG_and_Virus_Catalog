#!/bin/bash -l
#
#SBATCH -n 16 #number cores
#SBATCH -e logs/bbduk.more.log
#SBATCH -o logs/bbduk.more.log
#SBATCH --mem 50G #memory per node in Gb
#SBATCH -p intel
#SBATCH -J virus_qc

IN=data
SAMPFILE=data/sra_samples.csv

OUT=data/filtered

mkdir $OUT

IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read PREFIX 
do
	mkdir $OUT/$PREFIX
	
	if [ -f $IN/$PREFIX/$PREFIX'_1.fastq.gz' ]; then
		echo "$PREFIX has split fastqs"
		/rhome/cassande/bigdata/software/bbmap/bbduk.sh in1=$IN/$PREFIX/$PREFIX'_1.fastq.gz' in2=$IN/$PREFIX/$PREFIX'_2.fastq.gz' out1=$OUT/$PREFIX/$PREFIX'_1_filtered.fastq' out2=$OUT/$PREFIX/$PREFIX'_2_filtered.fastq' ref=/rhome/cassande/bigdata/software/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=10 maq=10

	elif [ -f $IN/$PREFIX/$PREFIX'.fastq.gz' ]; then
		echo "$PREFIX has interleaved fastq"
		/rhome/cassande/bigdata/software/bbmap/bbduk.sh in=$IN/$PREFIX/$PREFIX'.fastq.gz' out=$OUT/$PREFIX/$PREFIX'_interleaved_filtered.fastq' ref=/rhome/cassande/bigdata/software/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=10 maq=10

	fi 

done

