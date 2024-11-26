#!/bin/bash -l
#
#SBATCH -n 16 #number cores
#SBATCH -e logs/vcontact2.dRep.log
#SBATCH -o logs/vcontact2.dRep.log
#SBATCH --mem 50G #memory per node in Gb
#SBATCH -p intel
#SBATCH -J sg_virus_vcontact2


IN=dRep_clustered.fa
DIR=clustered_vOTUs
OUT=vcontact2_results_dRep_clustered_plus1Aug2023inphared
#PREFIX=dRep_clustered
#OUTFILE=dRep_clustered_prot.fa

PREFIX=dRep_clustered_plus1Aug2023inphared
OUTFILE=dRep_clustered_plus1Aug2023inphared_prot.faa

module load centos
centos.sh

#vcontact2 on vOTUs

conda activate vContact2

#prodigal -i $DIR/$PREFIX.fa -o $DIR/$PREFIX.coords.gbk -a $DIR/$OUTFILE -p meta

#vcontact2_gene2genome -p $DIR/$OUTFILE -o $DIR/$PREFIX.gene2genome.csv -s Prodigal-FAA

#combined with ICTV files
#https://github.com/RyanCook94/inphared#supplementing-and-annotating-vcontact2-clusters
#combine the date_vConTACT2_proteins.faa with your own fasta of file of translated ORFs, and combine date_vConTACT2_gene_to_genome.csv with your own mapping file (watch out for duplicated headers in the gene_to_genome.csv file if your file already has headers). Then run vConTACT2 as normal using the --db 'None' option, as this will avoid RefSeq duplicates.

vcontact2 --raw-proteins $DIR/$OUTFILE --rel-mode Diamond --proteins-fp $DIR/$PREFIX.gene2genome.csv --db 'None' --pcs-mode MCL --vcs-mode ClusterONE  --output-dir $DIR/$OUT

