#!/usr/bin/env bash
#SBATCH -p kellis
#SBATCH --array=1-7
# module load miniconda3
# source /home/software/conda/miniconda3/bin/condainit
source $HOME/src/conda_init.sh
source $HOME/.compbiorc
cd ${HOME}/data-group/Benjamin/AD_ATAC/TSS6
conda activate scanpy
export ATAC_CT=$(tail -n+2 "celltypes.csv" | awk -F, -v TASK=${TASK} 'NR == TASK { print $1 }')
export RNA_CT=$(tail -n+2 "celltypes.csv" | awk -F, -v TASK=${TASK} 'NR == TASK { print $2 }')

export ITER_PREV=${1:-0}
export ITER_CUR=$((ITER_PREV + 1))
export ATAC_GM="GeneMatrix_estimated_${ITER_PREV}.h5ad"
export RNA_GM="../RNA.h5ad"
echo "Using ATAC Gene matrix: ${ATAC_GM}"
echo "Using ATAC Celltype: ${ATAC_CT}"
echo "Using RNA Gene matrix: ${RNA_GM}"
echo "Using RNA Celltype: ${RNA_CT}"
export IPREFIX="integration/${RNA_CT}_${ITER_CUR}"
`type -P time` -v python3 $HOME/prj/AD_ATAC/integrate.py \
	       --subset-atac CellType="${ATAC_CT}" \
	       --subset-rna CellType="${RNA_CT}" \
	       --atac "${ATAC_GM}" --rna "${RNA_GM}" \
	       --hvg 10000 --leiden 1 2 \
	       -o "${IPREFIX}.h5ad" -t "${IPREFIX}.tsv.gz" \
