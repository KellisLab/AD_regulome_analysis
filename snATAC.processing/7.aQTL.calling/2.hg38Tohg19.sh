for ct in Ast Ex In Microglia OPC Oligo  PerEndo
do
sed -i '1d' $ct.peaks.hg38.bed
liftOver $ct.peaks.hg38.bed /net/bmc-lab5/data/kellis/users/xiongxs/Reference/human/hg38ToHg19.over.chain.gz $ct.peaks.hg19.bed unmap
done
