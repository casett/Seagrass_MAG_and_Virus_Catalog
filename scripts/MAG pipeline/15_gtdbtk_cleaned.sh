#!/bin/bash -l
#SBATCH --ntasks=16 # Number of cores
#SBATCH --mem=400G # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -p intel
#SBATCH -o logs/15_gtdbtk.log
#SBATCH -e logs/15_gtdbtk.log
#SBATCH -J sgmag_gtdbtk

conda activate gtdbtk-2.2.6

INPUT=magpurify_results/cleaned_mags
OUTPUT=magpurify_results/gtbdk_results_v2_2_6
CPU=16
PREFIX=SG_Chytid_magpurify

GTDBTK_DATA_PATH=/rhome/cassande/bigdata/.conda/envs/gtdbtk-2.2.6/share/gtdbtk-2.2.6/db

gtdbtk classify_wf --genome_dir $INPUT --out_dir $OUTPUT -x .fa --cpus $CPU --prefix $PREFIX.gtbdk --mash_db $OUTPUT/mash_db

