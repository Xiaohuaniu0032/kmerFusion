use strict;
use warnings;

# 1.统计fusion.list.pass.filter.csv文件信息
# 2.统计genefuse结果文件*.json中检测到的fusion信息


my ($flist,$fjson) = @ARGV;

my ($fgene_num,$fgene_search_len_bp) = (0,0);
open IN, "$flist" or die;
while (<IN>){
	chomp;
	if (/^>/){
		my $reg = (split /\:/, $_)[1];
		my @reg = split /\-/, $reg;
		$fgene_num += 1;
		my $len = $reg[1] - $reg[0];
		$fgene_search_len_bp += $len;
	}
}
close IN;

print "This sample have $fgene_num likely fusion gene(s)\nThe genefuse search len is $fgene_search_len_bp bp\n";

print "The detected fusion gene(s) are:\n";
my $cmd = "less $fjson | grep "fusion"";
system($cmd);
print "Stat Over\nGood Luck\n";

