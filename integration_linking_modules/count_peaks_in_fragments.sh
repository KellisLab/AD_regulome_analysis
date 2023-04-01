#!/usr/bin/env bash
#SBATCH -p kellis
#SBATCH --nice
#SBATCH --array=1-93%15

module load samtools
module load bedtools
export SAMPLE=$(awk -v TASK="$SLURM_ARRAY_TASK_ID" 'NR == TASK' < /net/bmc-lab5/data/kellis/users/benjames/AD_ATAC/fragments/samples.txt)

export PEAKDIR="/net/bmc-lab5/data/kellis/users/benjames/AD_ATAC/fragments/Peaks/"
export MATDIR="/net/bmc-lab5/data/kellis/users/benjames/AD_ATAC/fragments/Matrices/"
export FRAG="/net/bmc-lab5/data/kellis/users/benjames/AD_ATAC/fragments/${SAMPLE}_fragments.bed.gz"
while read celltype; do
    echo "Celltype: \"${celltype}\""
    bedtools intersect -a "${FRAG}" -b "${PEAKDIR}/${celltype}.bed.gz" -wa -wb \
	| awk -v SAMPLE="${SAMPLE}" '{ print SAMPLE"#"$4,$9,$5 }' OFS="\t" \
	| sort -k1 -k2 \
	| gzip -c > "${MATDIR}/${SAMPLE}_${celltype}.tsv.gz"
done < /net/bmc-lab5/data/kellis/users/benjames/AD_ATAC/fragments/celltypes.txt
