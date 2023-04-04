#! /usr/bin/perl -w
use strict;

open IN,"/net/bmc-pub14/data/kellis/users/xiongxs/Proj_snATAC_AD/CellMarker/FromCarles/PEC_markers_2.tsv" or die $!;
my %ct;

<IN>;
while(<IN>){
	chomp;
	my @sp=split /\t/;
	$ct{$sp[0]}=$sp[1];
}
close IN;

my %mat;
open IN2,"GeneMat_by_clusters.Feat50k.txt" or die $!;

chomp(my $head=<IN2>);
my @head=split /\t/,$head;

while(<IN2>){
	chomp;
	my @sp=split /\t/;
	next if !exists $ct{$sp[0]};
	for(my $i=1;;$i++){
		last if !defined $sp[$i];
		$mat{$ct{$sp[0]}}{$head[$i]}+=$sp[$i];
	}
}
close IN2;

open OUT,">PEC_markers_2.cluster.Feat50k.mat.txt" or die $!;
print OUT "Celltype";
for(my $i=1;;$i++){
	last if !defined $head[$i];	
	print OUT "\t$head[$i]";
}
print OUT  "\n";

foreach my $ct(sort keys %mat){
	print OUT "$ct";
	for(my $i=1;;$i++){
		last if !defined $head[$i];
		print OUT "\t$mat{$ct}{$head[$i]}";
	}
	print OUT  "\n";
}







