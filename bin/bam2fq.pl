use strict;
use warnings;

my ($bam,$name,$outdir) = @ARGV;

my $t1 = time();

my $sort_by_name_bam = "$outdir/$name\.sort_by_name.bam";
my $cmd = "samtools sort -n -o $sort_by_name_bam $bam";
my $start_time = localtime();
print "sort bam by seq name START at $start_time\n";
system($cmd);
my $end_time = localtime();
print "sort bam by seq name END at $end_time\n";

my $fq1 = "$outdir/$name\.R1.fastq";
my $fq2 = "$outdir/$name\.R2.fastq";

$cmd = "bedtools bamtofastq -i $sort_by_name_bam -fq $fq1 -fq2 $fq2";
$start_time = localtime();
print "convert bam to fastq START at $start_time\n";
system($cmd);
$end_time = localtime();
print "convert bam to fastq END at $end_time\n";

my $t2 = time();

my $time_used_minutes = int(($t2-$t1)/60);
print "bam2fq step used $time_used_minutes Minutes\n";
