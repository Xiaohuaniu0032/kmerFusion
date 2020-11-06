use strict;
use warnings;
use File::Basename;
use FindBin qw/$Bin/;
use Getopt::Long;

my ($bam,$sampleName,$fa_ref,$annot_file,$main_NM,$chrNaming,$outdir);

GetOptions(
    "bam:s" => \$bam,                  # Need
    "n:s" => \$sampleName,             # Need
    "fa:s" => \$fa_ref,                # Default: /data1/database/b37/human_g1k_v37.fasta
    "annot:s" => \$annot_file,         # Default: /home/fulongfei/workdir/git_repo/GeExCNV/public_db/refFlat.txt
    "NM:s" => \$main_NM,               # Default: /home/fulongfei/workdir/git_repo/GeExCNV/public_db/MainNM_yrt_20191111.txt
    "chrName:s" => \$chrNaming,        # Default: no_chr_prefix
    "od:s" => \$outdir,                # Need
    ) or die "unknown args\n";


# default value
if (not defined $fa_ref){
    $fa_ref = "/data1/database/b37/human_g1k_v37.fasta";
}

if (not defined $annot_file){
    $annot_file = "/home/fulongfei/workdir/git_repo/GeExCNV/public_db/refFlat.txt";
}

if (not defined $main_NM){
    $main_NM = "/home/fulongfei/workdir/git_repo/GeExCNV/public_db/MainNM_yrt_20191111.txt";
}

if (not defined $chrNaming){
    $chrNaming = "no_chr_prefix";
}


# check outdir
if (!-d $outdir){
    `mkdir -p $outdir`;
}

#my $runsh = "$outdir/$sampleName\.genefuse.sh";
my $runsh = "$outdir/genefuse\.$sampleName\.sh";
open O, ">$runsh" or die;

# convert bam into fq using bedtools bamtofastq
# sort by name, not by coord
my $bam_sort_by_name = "$outdir/$sampleName\.sort_by_name.bam";
my $cmd = "samtools sort -n -o $bam_sort_by_name $bam";
print O "$cmd\n";

# bamtofastq
my $fq1 = "$outdir/$sampleName\.R1.fastq";
my $fq2 = "$outdir/$sampleName\.R2.fastq";

$cmd = "bedtools bamtofastq -i $bam_sort_by_name -fq $fq1 -fq2 $fq2";
print O "$cmd\n";

# make csv list
my $bed = "$Bin/bin/all_genes.bed";
$cmd = "bedtools intersect -abam $bam -b $bed -bed -wo >$outdir/covered_genes.txt";
print O "$cmd\n";

my $fusion_list = "$outdir/fusion_list.csv";
$cmd = "perl $Bin/bin/make_genefuse_csv_list.pl -gene $outdir/covered_genes.txt -annot $annot_file -NM $main_NM -chrname $chrNaming -of $fusion_list";
print O "$cmd\n";

# run genefuse
$cmd = "$Bin/genefuse -r $fa_ref -f $fusion_list -1 $fq1 -2 $fq2 -h $outdir/$sampleName\.genefuse.html -j $outdir/$sampleName\.genefuse.json";
print O "$cmd\n";

close O;