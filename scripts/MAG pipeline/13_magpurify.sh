#!/bin/bash -l
##
#SBATCH -o logs/13_MAGpurify.log
#SBATCH -e logs/13_MAGpurify.log
#SBATCH --nodes=1
#SBATCH --ntasks=12 # Number of cores
#SBATCH --mem=60G # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -J sg_chytid_magpurify
#SBATCH -p batch

conda activate /rhome/cassande/.conda/envs/magpurify

export MAGPURIFYDB=/rhome/cassande/bigdata/software/MAGpurify-db-v1.0


INPUTDIR=dastool_bac_bins
OUT=magpurify_results

for MAG in $(ls $INPUTDIR);
do 
	BASE=$(basename $MAG .fa)
	OUTPUTDIR=$OUT/$BASE.magpurify
	OUTPUTMAG=$OUTPUTDIR/$BASE.cleaned.fa 
	
	magpurify phylo-markers $INPUTDIR/$MAG $OUTPUTDIR
	magpurify clade-markers $INPUTDIR/$MAG $OUTPUTDIR
	magpurify tetra-freq $INPUTDIR/$MAG $OUTPUTDIR
	magpurify gc-content $INPUTDIR/$MAG $OUTPUTDIR
	magpurify known-contam $INPUTDIR/$MAG $OUTPUTDIR

	magpurify clean-bin $INPUTDIR/$MAG $OUTPUTDIR $OUTPUTMAG
	
done
