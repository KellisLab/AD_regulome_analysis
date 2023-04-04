#! /usr/bin/perl -w
use strict;

my $ct = $ARGV[0];
my $num=$ARGV[1];
my $type="excl_AD";

open OUT2,">2.$ct.pf2-$num.$type.runFastqtl.sh" or die $!;
print OUT2 "#!/bin/bash
#SBATCH -N 1                                # Request one node
#SBATCH -n 4                               # Request n cpu (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-020:00                          # Runtime in D-HH:MM format
#SBATCH -p kellis
#SBATCH -o $ct.pf2-$num.$type.fastqtl.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e $ct.pf2-$num.$type.fastqtl.err                 # File to which STDERR will be written, including job ID
#SBATCH --mem=20000                         ## member used
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=xiongxushen19s2\@gmail.com    # Email to which notifications will be sent

";

	for(my $i=1;;$i++){
		last if $i==23;
		print OUT2 "/net/bmc-lab5/data/kellis/users/xiongxs/software/fastqtl/FastQTL/bin/fastQTL.static --vcf /home/xiongxs/data-lab5/Proj_snATAC_AD/Genotypes/rosmap1709-chr$i.xs.more.vcf.gz --bed ../$ct.TMM.voom.normed.mtx.sort.bed.gz --region $i --out output/$ct.chr$i.out.txt.gz --cov ../7.PEER.cor/$ct.PEER_covariates.pf2-$num.$type.txt --exclude-sites /net/bmc-lab5/data/kellis/users/xiongxs/Proj_snATAC_AD/Processing_0414_AllPFC_ATAC/Use_ArchR/Try.atacQTL.GTavailable.2/FastQTL.test/excl/PFC.ATAC.chr$i.excl.snp --window 1e6 --permute 1000\n";
	}
	close OUT2;
	




