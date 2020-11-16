use strict;
use warnings;
use File::Basename;

my $dir = "/home/fulongfei/workdir/fusion_work/kmerFusion_test_20201113_new";

my @res = glob "$dir/*/*.xls";

for my $res (@res){
	my $dir = dirname($res);	
	my $fcsv = "$dir/fusion.csv";
	my $fgene_n = `less $fcsv | grep ">" | wc -l`;
	chomp $fgene_n;
	
	my @she = glob "$dir/*sh.e*";
	for my $f (@she){
		my $flag = 0;
		open IN, "$f" or die;
		while (<IN>){
			chomp;
			if (/^terminate/){
					
