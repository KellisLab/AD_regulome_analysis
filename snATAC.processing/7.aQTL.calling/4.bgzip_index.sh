for ct in Ast Ex In Microglia OPC Oligo  PerEndo
do
awk 'NR == 1; NR > 1 {print $0 | "sort -k1,1 -k2n,2"}' $ct.TMM.voom.normed.mtx.bed |sed -e 's/chr//' > $ct.TMM.voom.normed.mtx.sort.bed
bgzip $ct.TMM.voom.normed.mtx.sort.bed  && tabix -p bed $ct.TMM.voom.normed.mtx.sort.bed.gz
done

