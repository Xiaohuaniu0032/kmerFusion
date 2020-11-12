use strict;
use warnings;


my ($infile,$name,$outdir) = @ARGV;

my $fcount = 2; # default value

my %info;
open IN, "$infile" or die;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    push @{$info{$arr[0]}}, $arr[1];
}
close IN;

my $outfile = "$outdir/$name\.softclip.pos.PASS_FILTER.txt";
open O, ">$outfile" or die;

foreach my $chr (keys %info){
	print "processing chr $chr\.\.\.\n";
	my @pos_sort = sort {$a <=> $b} @{$info{$chr}}; # may contain repeat pos

	my %rep_pos;
	for my $pos (@pos_sort){
		$rep_pos{$pos} += 1;
	}

	my $flag = 0;
	foreach my $pos (sort {$a <=> $b} keys %rep_pos){
		my $n = $rep_pos{$pos};
		if ($n > $fcount){
			# effective pos
			$flag += 1;
			print "POS $pos is a OK sv pos\n";
			my $start_pos = $pos - 1; # 0-based
			my $end_pos = $pos;
			my $flag_new = "POS".$flag;
			print O "$chr\t$start_pos\t$end_pos\t$flag_new\n";
		}
	}
}
close O;






