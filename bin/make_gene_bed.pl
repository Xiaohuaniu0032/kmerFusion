use strict;
use warnings;

my ($cosmicGeneFile,$ensGeneAnnotFile,$chrname,$outdir) = @ARGV;

my $gene_bed_file = "$outdir/all.gene.bed";

# link ENST and gene name
my %enst2gene;
open IN, "$cosmicGeneFile" or die;
while (<IN>){
	chomp;
	my @arr = split /\t/;
	$enst2gene{$arr[1]} = $arr[0]; # ENST=>gene
}
close IN;

open O, ">$gene_bed_file" or die;
# use ENST to get gene name, and its region
open IN, "$ensGeneAnnotFile" or die;
while (<IN>){
	chomp;
	my @arr = split /\t/;
	my $enst = $arr[1];

	my $chr = $arr[2];
	if ($chrname ne "with_chr_prefix"){
		$chr =~ s/^chr//;
	}

	if (exists $enst2gene{$$enst}){
		my $gene_name = $enst2gene{$enst};
		my $gene_start = $arr[4];
		my $gene_end = $arr[5];
		
		print O "$chr\t$gene_start\t$gene_end\t$gene_name\n";
	}
}

close IN;
close O;
