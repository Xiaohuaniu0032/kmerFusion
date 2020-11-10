use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my ($sv_read_name,$fq1,$fq2,$name,$outdir);
GetOptions(
    "list:s" => \$sv_read_name,        # NEED
    "fq1:s" => \$fq1,                  # NEED
    "fq2:s" => \$fq2,                  # NEED
    "n:s" => \$name,                   # NEED
    "od:s" => \$outdir,                # NEED
    ) or die;


my $if_fq_gz;
if ($fq1 =~ /\.gz$/){
    $if_fq_gz = "yes";
    print "INFO:input fq is .gz\n";
}else{
    $if_fq_gz = "no";
    print "INFO:input fq is not .gz\n";
}

my %read_name;
open IN, "$sv_read_name" or die;
while (<IN>){
    chomp;
    $read_name{$_} = 1;
}
close IN;

my $sv_fq1 = "$outdir/$name\.R1.fastq.fusionmap.gz";
my $sv_fq2 = "$outdir/$name\.R2.fastq.fusionmap.gz";

open FQ1, "| gzip >$sv_fq1" or die;
open FQ2, "| gzip >$sv_fq2" or die;

if ($if_fq_gz eq "yes"){
    # .gz
    open IN1, "gunzip -dc $fq1 |" or die;
    while (<IN1>){
        chomp;
        my $seq = <IN1>;
        <IN1>;
        my $qual = <IN1>;
        my $name = (split /\s/, $_)[0];
        $name =~ s/^\@//;
        if (exists $read_name{$name}){
            print FQ1 "$_\n";
            print FQ1 "$seq";
            print FQ1 "+\n";
            print FQ1 "$qual";
        }
    }
    close IN1;

    open IN2, "gunzip -dc $fq2 |" or die;
    while (<IN2>){
        chomp;
        my $seq = <IN2>;
        <IN2>;
        my $qual = <IN2>;
        my $name = (split /\s/, $_)[0];
        $name =~ s/^\@//;
        if (exists $read_name{$name}){
            print FQ2 "$_\n";
            print FQ2 "$seq";
            print FQ2 "+\n";
            print FQ2 "$qual";
        }
    }
    close IN2;
}else{
    open IN1, "$fq1" or die;
    while (<IN1>){
        chomp;
        my $seq = <IN1>;
        <IN1>;
        my $qual = <IN1>;
        my $name = (split /\s/, $_)[0];
        $name =~ s/^\@//;
        if (exists $read_name{$name}){
            print FQ1 "$_\n";
            print FQ1 "$seq";
            print FQ1 "+\n";
            print FQ1 "$qual";
        }
    }
    close IN1;

    open IN2, "$fq2" or die;
    while (<IN2>){
        chomp;
        my $seq = <IN2>;
        <IN2>;
        my $qual = <IN2>;
        my $name = (split /\s/, $_)[0];
        $name =~ s/^\@//;
        if (exists $read_name{$name}){
            print FQ2 "$_\n";
            print FQ2 "$seq";
            print FQ2 "+\n";
            print FQ2 "$qual";
        }
    }
    close IN2;
}

close FQ1;
close FQ2;

# check lines
my $fq1_line = (split /\s/, `gunzip -dc $sv_fq1 | wc -l`)[0]/4;
my $fq2_line = (split /\s/, `gunzip -dc $sv_fq2 | wc -l`)[0]/4;

if ($fq1_line == $fq2_line){
    print "INFO:find $fq1_line likely SV reads\n";
}else{
    die "ERROR:find $fq1_line likely SV reads for FQ1, and find $fq2_line likely SV reads for FQ2\n";
}


