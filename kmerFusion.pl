use strict;
use warnings;
use File::Basename;
use FindBin qw/$Bin/;
use Getopt::Long;

my ($bam,$sampleName,$fa_ref,$padding_len,$outdir);

GetOptions(
    "bam:s"         => \$bam,                   # Need
    "n:s"           => \$sampleName,            # Need
    "fa:s"          => \$fa_ref,                # Default: </data1/database/b37/human_g1k_v37.fasta>
    "Plen:i"        => \$padding_len,           # Default: 1000
    "od:s"          => \$outdir,                # Need
    ) or die;

# default value
if (not defined $fa_ref){
    $fa_ref = "/data1/database/b37/human_g1k_v37.fasta";
}

if (not defined $padding_len){
    $padding_len = 1000;
}

my $soft_clip_count = 2;

my $runsh = "$outdir/gfuse\.$sampleName\.sh";
open RUNSH, ">$runsh" or die;

# bam2fq
my $cmd = "perl $Bin/bin/bam2fq.pl $bam $sampleName $outdir";
#system($cmd);
print RUNSH "$cmd\n";

# get soft-clip region info
my $soft_clip_pos_file = "$outdir/$sampleName\.softclip.pos.txt";
$cmd = "samtools view $bam | awk '\$6~/S/ {print \$3\"\\t\"\$4}\' >$soft_clip_pos_file";
#system($cmd);
print RUNSH "$cmd\n";

# filter soft clip pos
# if a pos has <= 2 clip reads, then this pos will be filtered out
$cmd = "perl $Bin/bin/filter_sv_pos.pl $soft_clip_pos_file $sampleName $outdir";
#system($cmd);
print RUNSH "$cmd\n";

# make gene bed
my $cosmicGene = "$Bin/public_db/cosmic_transcript.txt.gz";
my $ensGeneAnnot = "$Bin/public_db/hg19_ensGene.txt.gz";
my $chr_naming = &check_chr_naming($fa_ref);
$cmd = "perl $Bin/bin/make_gene_bed.pl $cosmicGene $ensGeneAnnot $chr_naming $outdir";
#system($cmd);
print RUNSH "$cmd\n";

# annot soft clip pos with gene name using bedtools
my $all_gene_bed = "$outdir/all.gene.bed";
my $pass_filter_pos_file = "$outdir/$sampleName\.softclip.pos.PASS_FILTER.txt";
$cmd = "bedtools intersect -a $pass_filter_pos_file -b $all_gene_bed -wo -bed >$outdir/$sampleName\.softclip.pos.PASS.annot.Gene.txt";
#system($cmd);
print RUNSH "$cmd\n";

# make fusion csv file
$cmd = "perl $Bin/bin/make_fusion_list.pl $cosmicGene $ensGeneAnnot $outdir/$sampleName\.softclip.pos.PASS.annot.Gene.txt $chr_naming $padding_len $outdir";
#system($cmd);
print RUNSH "$cmd\n";

# call genefuse to detect fusion gene
my $flist = "$outdir/fusion.csv";
my $fq1 = "$outdir/$sampleName\.R1.fastq";
my $fq2 = "$outdir/$sampleName\.R2.fastq";

$cmd = "$Bin/genefuse -t 12 -r $fa_ref -f $flist -1 $fq1 -2 $fq2 -h $outdir/$sampleName\.genefuse.html -j $outdir/$sampleName\.genefuse.json";
#system($cmd);
print RUNSH "$cmd\n";

# get fusion result from json file
my $json = "$outdir/$sampleName\.genefuse.json";
$cmd = "perl $Bin/bin/get_final_fusion_result.pl $json $sampleName $outdir";
#system($cmd);
print RUNSH "$cmd\n";

# stat genefuse result info


close RUNSH;




sub check_chr_naming{
    my ($fa) = @_;
    open IN, "$fa" or die;
    my $first_line = <IN>;
    close IN;
    
    my $naming;
    if ($first_line =~ /^>chr/){
        $naming = "with_chr_prefix";
    }else{
        $naming = "no_chr_prefix";
    }

    return($naming);
}


