#!/usr/bin/env bash
#SBATCH -p kellis
#SBATCH --nice
#SBATCH --array=1-7
source $HOME/src/conda_init.sh
source $HOME/.compbiorc
conda activate scanpy
export CELLTYPES_CSV="${HOME}/data-group/Benjamin/AD_ATAC/TSS6/celltypes.csv"
export ATAC_CT=$(tail -n+2 ${CELLTYPES_CSV} | awk -F, -v TASK=${TASK} 'NR == TASK { print $1 }')

export PEAKS="/home/benjames/data/AD_ATAC/fragments/Matrices/${ATAC_CT}.h5ad"
export OUTPUT="/net/bmc-lab5/data/kellis/group/Benjamin/AD_ATAC/linking/pseudobulk/${ATAC_CT}.h5ad"
python3 ~/prj/AD_ATAC/coaccessibility.py --peaks "${PEAKS}" --output "${OUTPUT}"
