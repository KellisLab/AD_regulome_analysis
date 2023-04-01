#!/usr/bin/env bash
#SBATCH -p kellis
#SBATCH --array=1-93%5
cd /home/benjames/tmp/tmp/
source ~/src/conda_init.sh
conda activate snapatac
export FRAG=$(< "/home/benjames/data-group/Benjamin/AD_ATAC/pfc_fragment_files.txt" awk -v TASK=${SLURM_ARRAY_TASK_ID} 'NR == TASK { print $1 }')
export SAMPLE=$(< "/home/benjames/data-group/Benjamin/AD_ATAC/pfc_fragment_files.txt" awk -v TASK=${SLURM_ARRAY_TASK_ID} 'NR == TASK { print $2 }')
echo "Computing fragment $FRAG for sample $SAMPLE"

`type -P time` -v python3 ~/prj/AD_ATAC/load_fragments_from_annot.py \
	-i "${FRAG}" \
	-o "${SAMPLE}.h5ad" \
	-s "${SAMPLE}" \
	-a "/net/bmc-lab5/data/kellis/group/Benjamin/AD_ATAC/TSS6/PeakMatrix/barcodes.tsv.gz"
