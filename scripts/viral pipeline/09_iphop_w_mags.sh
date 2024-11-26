#!/bin/bash -l
#
#SBATCH -n 16 #number cores
#SBATCH -e logs/iphop.dRep.mags.log
#SBATCH -o logs/iphop.dRep.mags.log
#SBATCH --mem 50G #memory per node in Gb
#SBATCH -p intel
#SBATCH -J sg_virus_iphop


DIR=clustered_vOTUs
VOTU=dRep_clustered.fa
OUTDIR=iphop_mags
DB=/rhome/cassande/bigdata/software/Sept_2021_pub_rw_w_SGMAG

conda activate iphop_env


#iphop add_to_db --fna_dir /rhome/cassande/bigdata/eisenlab/sg_chytrid/stajichlab/magpurify_results/cleaned_mags --gtdb_dir /rhome/cassande/bigdata/eisenlab/sg_chytrid/stajichlab/magpurify_results/gtbdk_results_v2_2_6_for_iphop/ --out_dir Sept_2021_pub_rw_w_SGMAG --db_dir /path/to/iphop_db/Sept_2021_pub_rw/
 
iphop predict --fa_file $DIR/$VOTU --out_dir $OUTDIR --db_dir $DB
