use strict;
use warnings;

my ($cosmicGeneFile,$ensGeneAnnotFile,$softclipPosPASSAnnotGeneFile,$chrname,$outdir) = @ARGV;

# link enst with gene name
my %gene2enst;
open IN, "$cosmicGeneFile" or die;
while (<IN>){
	chomp;
	my @arr = split /\t/;
	$gene2enst{$arr[0]} = $arr[1]; # gene=>enst
}
close IN;

# get enst exon region info
my %enstExonAnnotInfo;
my %enst_chr;

open IN, "$ensGeneAnnotFile" or die;
while (<IN>){
	chomp;
	my @arr = split /\t/;
	my $enst = $arr[1];
	my $strand = $arr[3];
	
	my $chr = $arr[2];
	if ($chrname eq "no_chr_prefix"){
		$chr =~ s/^chr//;
	}

	$enst_chr{$enst} = $chr; # chr for fusion csv file
	
	my $start_exon = $arr[-7];
	$start_exon =~ s/\,//;

	my $end_exon = $arr[-6];
	$end_exon =~ s/\,//;
	
	my @start_exon = split /\,/, $start_exon;
	my @end_exon = split /\,/, $end_exon;
	
	my $exon_num = @start_exon;
	
	my $exon_idx = 0;
	
	if ($strand eq "+"){
		for (my $i=1;$i<=$exon_num;$i++){
			my $start = $start_exon[$i-1];
			my $end = $end_exon[$i-1];
			$exon_idx += 1;
			$enstExonAnnotInfo{$enst}{$exon_idx} = "$start\t$end"; # enst=>exon=>region
		}
	}else{
		for (my $i=1;$i<=$exon_num;$i++){
			my $start = $start_exon[$i-1];
			my $end = $end_exon[$i-1];
			$exon_idx = $exon_num - $i + 1;
			$enstExonAnnotInfo{$enst}{$exon_idx} = "$start\t$end";
		}
	}
}
close IN;

# get each gene's soft clip region
my %soft_clip_pos;
my %genes;

open IN, "$softclipPosPASSAnnotGeneFile" or die;
while (<IN>){
	chomp;
	my @arr = split /\t/;
	my $pos = $arr[1] + 1;
	$soft_clip_pos{$arr[-2]}{$pos} = 1; # gene=>pos=>1
	$genes{$arr[-2]} = 1;
}
close IN;

# make final fusion csv list
my $fusion_list = "$outdir/fusion.csv";
open O, ">$fusion_list" or die;

foreach my $gene (keys %genes){
	my $enst = $gene2enst{$gene}; # get enst
	my $chr = $enst_chr{$enst};
	my @pos_sort = sort {$a <=> $b} keys %{$soft_clip_pos{$gene}};
	
	my $pos_sort_left = $pos_sort[0] - $padding_len; # extent soft clip pos to left 1k bp
	my $pos_sort_right = $pos_sort[-1] + $padding_len; # extent soft clip pos to right 1k bp
	
	# print sv region for genefuse to search
	print O "\>$gene\,$chr\:$pos_sort_left\-$pos_sort_right\n";
	
	# print all exon info
	my @exon_sort = sort {$a <=> $b} keys %{$enstExonAnnotInfo{$enst}};
	for my $e (@exon_sort){
		my @r = split /\t/, $enstExonAnnotInfo{$enst}{$e};
		print O "$e\,$r[0]\,$r[1]\n";
	}
	
	print O "\n";
}
close O;



