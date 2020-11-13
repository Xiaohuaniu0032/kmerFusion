use strict;
use warnings;

my ($json,$name,$outdir) = @ARGV;

my $fusion_result_tmp = "$outdir/fusion.result.tmp\.$name\.xls";
my $cmd = "less $json | grep \"Fusion\" | grep -v \"command\" >$fusion_result_tmp";
system($cmd);

#my $fusion_result_tmp = "$outdir/fusion.result.tmp\.$name\.xls";
my $final_result = "$outdir/$name\.fusion.xls";
open O, ">$final_result" or die;
#print O "gene1\tgene1_chr\tgene1_breakpoint\tgene1_annot_region\tgene1_ENST\tgene2\tgene2_chr\tgene2_breakpoint\tgene2_annot_region\tgene2_ENST\ttotal_support_reads\tuniq_support_reads\n";
print O "fusion_info\ttotal_support_reads\tuniq_support_reads\n";

open IN, "$fusion_result_tmp" or die;
while (<IN>){
	chomp;
	#print "$_\n";
	if (/Fusion: (.+)\s+\(total: (\d+)\, unique:(\d+)\)/){
		print O "$1\t$2\t$3\n";
	}
}
close IN;
close O;
`rm $fusion_result_tmp`;
