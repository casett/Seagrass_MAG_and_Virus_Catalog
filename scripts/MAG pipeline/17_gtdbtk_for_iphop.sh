#!/bin/bash -l
#SBATCH --ntasks=16 # Number of cores
#SBATCH --mem=400G # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -p intel
#SBATCH -o logs/17_gtdbtk_for_iphop.log
#SBATCH -e logs/17_gtdbtk_for_iphop.log
#SBATCH -J sgmag_gtdbtk_for_iphop



conda activate gtdbtk-2.2.6

INPUT=magpurify_results/cleaned_mags
OUTPUT=magpurify_results/gtbdk_results_v2_2_6_for_iphop
CPU=16
PREFIX=SG_Chytid_magpurify

GTDBTK_DATA_PATH=/rhome/cassande/bigdata/.conda/envs/gtdbtk-2.2.6/share/gtdbtk-2.2.6/db


gtdbtk de_novo_wf --genome_dir $INPUT --bacteria --outgroup_taxon p__Patescibacteria --out_dir $OUTPUT --cpus $CPU --force --extension fa
