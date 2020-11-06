use strict;
use warnings;
use File::Basename;

my ($bam,$sampleName,$fa_ref,$outdir);

GetOptions(
    "bam:s" => \$bam,
    "n:s" => \$sampleName,
    "fa:s" => \$fa_ref,
    "od:s" => \$outdir,
    ) or die "unknown args\n";


my $runsh = "$od/$sampleName\.genefuse.sh";
open O, ">$runsh" or die;

# convert bam into fq using bedtools bamtofastq
# sort by name, not by coord
my $bam_sort_by_name = "$outdir/$sampleName\.sort_by_name.bam";