use strict;
use warnings;

my ($json,$name,$outdir) = @ARGV;

my $cmd = "less $json | grep \"Fusion\" | grep -v \"command\" >$outdir/fusion.result.tmp\.$name";
system($cmd);

my $fusion_result_tmp = "$outdir/fusion.result.tmp\.$name";
my $final_result = "$outdir/fusion.result.final\.$name";
open O, ">$final_result" or die;
#print O "gene1\tgene1_chr\tgene1_breakpoint\tgene1_annot_region\tgene1_ENST\tgene2\tgene2_chr\tgene2_breakpoint\tgene2_annot_region\tgene2_ENST\ttotal_support_reads\tuniq_support_reads\n";
print O "fusion_info\ttotal_support_reads\tuniq_support_reads\n";

open IN, "$fusion_result_tmp" or die;
while (<IN>){
	chomp;
	my $fusion_info = (split /\s+/, $_)[1];
	if (/.*\(total: (\d+)\, unique: (\d+)\)/){
		my $total_n = $1;
		my $uniq_n = $2;
		print O "$fusion_info\t$total_n\t$uniq_n\n";
	}
}
close IN;
close O;
