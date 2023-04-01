#!/usr/bin/env bash
#SBATCH -p kellis
#SBATCH --array=1-7
#SBATCH --nice
source $HOME/src/conda_init.sh
source $HOME/.compbiorc
conda activate scanpy
export CELLTYPES_CSV="${HOME}/data-group/Benjamin/AD_ATAC/TSS6/celltypes.csv"
export ATAC=$(tail -n+2 ${CELLTYPES_CSV} | awk -F, -v TASK=${TASK} 'NR == TASK { print $1 }')
echo "celltype: ${ATAC}"
cd ~/data-group/Benjamin/AD_ATAC/linking/pseudobulk
`type -P time` -v python3 ~/prj/AD_ATAC/epiclust_predict_links.py \
	       --interactions "interactions_${ATAC}.tsv.gz" \
	       -i links_"${ATAC}"_*_nonsquared.tsv.gz \
	       --activity activity_per_celltype.tsv.gz \
	       --celltype "${ATAC}" \
	       --activity-columns mean_counts log1p_total_counts \
	       -o "predicted_${ATAC}.tsv.gz"
