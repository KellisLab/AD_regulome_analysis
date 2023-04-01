#!/usr/bin/env bash
#SBATCH -p kellis
#SBATCH --array=1-140%12
#SBATCH --nice
source $HOME/src/conda_init.sh
source $HOME/.compbiorc
conda activate scanpy
export PARAMS="/net/bmc-lab5/data/kellis/group/Benjamin/AD_ATAC/linking/pseudobulk/params.txt"

export POWER=$(awk -F, < ${PARAMS} -v TASK=${TASK} 'NR == TASK { print $1 }')
export ASSAY=$(awk -F, < ${PARAMS} -v TASK=${TASK} 'NR == TASK { print $2 }')
export CT=$(awk -F, < ${PARAMS} -v TASK=${TASK} 'NR == TASK { print $3 }')
export SQ=$(awk -F, < ${PARAMS} -v TASK=${TASK} 'NR == TASK { print $4 }')
export ADATA="${HOME}/data-group/Benjamin/AD_ATAC/linking/pseudobulk/${ASSAY}/${CT}.h5ad"
export INTER="${HOME}/data-group/Benjamin/AD_ATAC/linking/pseudobulk/interactions_${CT}.tsv.gz"
export OUT_PREFIX="${HOME}/data-group/Benjamin/AD_ATAC/linking/pseudobulk/links_${CT}_${ASSAY}_${POWER}_${SQ}"
`type -P time` -v python3 ~/prj/AD_ATAC/epiclust_linking.py -i "${ADATA}" --interactions "${INTER}" -o "${OUT_PREFIX}.tsv.gz" -p "${POWER}" "--${SQ}-correlation"
