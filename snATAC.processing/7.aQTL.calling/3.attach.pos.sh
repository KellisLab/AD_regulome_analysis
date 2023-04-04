for ct in Ast Ex In Microglia OPC Oligo  PerEndo
do
perl attach.pos.pl $ct.peaks.hg19.bed $ct.TMM.voom.normed.mtx $ct.TMM.voom.normed.mtx.bed 
done
