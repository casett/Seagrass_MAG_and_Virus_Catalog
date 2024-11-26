#!/bin/bash -l
#SBATCH --ntasks=24 # Number of cores
#SBATCH --mem=200G # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -p intel,batch
#SBATCH -o logs/12_gtdbtk.log
#SBATCH -e logs/12_gtdbtk.log
#SBATCH -J sgmag_gtdbtk


module unload miniconda2
module load miniconda3

conda activate gtdbtk-1.3.0


INPUT=bin_fasta

OUTPUT=gtbdk_results
CPU=24
PREFIX=SG_Chytid

gtdbtk classify_wf --genome_dir $INPUT --out_dir $OUTPUT -x .fa --cpus $CPU --prefix $PREFIX.gtbdk
