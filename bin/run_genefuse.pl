use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my ($sampleName,$sampleType,$bam,$fq1,$fq2,$fasta,$fusion_list,$outdir);
GetOptions(
    "n:s" => \$sampleName,     # NEED
    "t:s" => \$sampleType,     # NEED [bam|fq]
    "bam:s" => \$bam,          # Optional
    "fq1:s" => \$fq1,          # Optional
    "fq2:s" => \$fq2,          # Optional
    "fa:s" => \$fasta,         # Optional </data1/database/b37/human_g1k_v37.fasta>
    "flist:s" => \$fusion_list,# NEED
    "od:s" => \$outdir,        # NEED
    ) or die;

# you need to add bedtools,samtools,genefuse into your $PATH env

if (not defined $sampleName || not defined $sampleType || not defined $fusion_list || not defined $outdir){
    die "you must specify -n/-t/-od args\n";
}

if ($sampleType ne "bam" and $sampleType ne "fq"){
    die "-t must be <bam> or <fq>\n";
}


if (not defined $fasta){
    $fasta = "/data1/database/b37/human_g1k_v37.fasta";
}

open O, ">$outdir/genefuse\.$sampleName.sh";


my $html = "$outdir/$sampleName\.report.html";
my $json = "$outdir/$sampleName\.report.json";

if ($sampleType eq "bam"){
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

    $cmd = "genefuse -r $fasta -f $fusion_list -1 $fq1 -2 $fq2 -h $html -j $json";
    print O "$cmd\n";
}else{
    my $cmd = "genefuse -r $fasta -f $fusion_list -1 $fq1 -2 $fq2 -h $html -j $json";
    print O "$cmd\n";
}
close O;
