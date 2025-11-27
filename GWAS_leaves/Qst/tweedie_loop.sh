#!/bin/bash -l

set -euo pipefail

# Usage: bash tweedie_loop.sh 
# (OBS! Adjust Qst_tweedie_pipeline.R according to your data and give the correct path to the R file before you sbatch!)

# Set arguments
Tweedie_file="/mnt/picea/projects/aspseq/nstreet/swasp/Sara/GWAS_leaves/Pop_genetics/Qst/tweedie_extras/SPGs_buds_Outlier_removed_tweedie_extras_data.3.reps_pop_lat.txt"
OUT_DIR="/mnt/picea/projects/aspseq/nstreet/swasp/Sara/GWAS_leaves/Pop_genetics/Qst/tweedie_extras/res"

# Create a column file 
awk 'BEGIN{ FS="\t" }
       { for(fn=1;fn<=NF;fn++) {print fn" = "$fn}; exit; }
      ' $Tweedie_file > "${OUT_DIR}/tweedie_names.tsv"

# Number of phenotypes
N=$(cat "${OUT_DIR}/tweedie_names.tsv" | wc -l)

# Start at 5 to skip the meta columns
for i in $(seq 5 ${N}); do
COLNR=$(cat "${OUT_DIR}/tweedie_names.tsv" | awk -v lnr="$i" 'NR==lnr''{ print $1 }')
TRAIT=$(cat "${OUT_DIR}/tweedie_names.tsv" | awk  -v lnr="$i" 'NR==lnr''{ print $3 }')
echo "Submitting trait ${TRAIT}"
sbatch /mnt/picea/projects/aspseq/nstreet/swasp/Sara/GWAS_leaves/scripts/R/Pop_genetics/Qst/lmer/functions/runR.sh /mnt/picea/projects/aspseq/nstreet/swasp/Sara/GWAS_leaves/scripts/R/Pop_genetics/Qst/tweedie_extras/Qst_tweedie_pipeline.R ${COLNR} ${TRAIT} ${OUT_DIR}
done


