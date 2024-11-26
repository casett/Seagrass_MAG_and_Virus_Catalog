#!/bin/bash -l
#
#SBATCH -n 16 #number cores
#SBATCH -e logs/bbrepair.log
#SBATCH -o logs/bbrepair.log
#SBATCH --mem 50G #memory per node in Gb
#SBATCH -p intel,batch
#SBATCH -J virus_repair

IN=data
SAMPFILE=data/sra_samples.csv
OUT=data/filtered

#mkdir $OUT

IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read PREFIX 
do
	mkdir $OUT/$PREFIX
	
	if [ -f $OUT/$PREFIX/$PREFIX'_interleaved_filtered.fastq' ]; then
		echo "$PREFIX has interleaved fastq - fixing"
		/rhome/cassande/bigdata/software/bbmap/bbsplitpairs.sh in=$OUT/$PREFIX/$PREFIX'_interleaved_filtered.fastq' out=$OUT/$PREFIX/$PREFIX'_interleaved_filtered_fixed.fastq' outs=$OUT/$PREFIX/$PREFIX'_interleaved_filtered_singeltons.fastq' fint

	fi 

done

