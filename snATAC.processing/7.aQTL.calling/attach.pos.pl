#! /usr/bin/perl -w
use strict;

open IN,$ARGV[0] or die $!;
my %pos;

while(<IN>){
	chomp;
	my @sp=split /\t/;
	$pos{$sp[3]}="$sp[0]\t$sp[1]\t$sp[2]";   ## attach hg19 positions
}
close IN;

open IN2,$ARGV[1] or die $!;
open OUT,">$ARGV[2]" or die $!;

chomp(my $head=<IN2>);
my @head=split /\s+/,$head;
print OUT "#Chr\tstart\tend\tID";

for(my $i=1;;$i++){
	last if !defined $head[$i];
	$head[$i]=~ s/X//;
	print OUT "\t$head[$i]";  ## to match the vcf file
}
print OUT "\n";

while(<IN2>){
	chomp;
        my @sp=split /\t/;
	next if !exists $pos{$sp[0]};
	print OUT "$pos{$sp[0]}\t$_\n";
}

