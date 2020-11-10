use strict;
use warnings;

my $refFlat = "/home/fulongfei/workdir/git_repo/GeExCNV/public_db/refFlat.txt";
open IN, "$refFlat" or die;
while (<IN>){
	chomp;
	my @arr = split /\t/;
	my $chr = $arr[2];
	$chr =~ s/^chr//;
	print "$chr\t$arr[4]\t$arr[5]\t$arr[0]\n";
}
close IN;	
