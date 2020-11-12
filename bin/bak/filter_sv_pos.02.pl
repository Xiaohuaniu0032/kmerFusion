use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my ($file,$fcount,$outdir);

GetOptions(
    "f:s" => \$file,                # NEED <sv.pos.ALL.txt>
    "c:i" => \$fcount,              # Default: 2
    "od:s" => \$outdir              # NEED
    ) or die "unknown args\n";

# default value
if (not defined $fcount){
    $fcount = 2;
}

my %info;
open IN, "$file" or die;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    push @{$info{$arr[0]}}, $arr[1];
}
close IN;


my $outfile = "$outdir/sv.pos.PASS_FILTER.txt";
print "$outfile is created\n";
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
            $flag += 1;
            # pass filter
            print "POS $pos is a OK sv pos\n";
            #print O "$chr\t$pos\n";
            my $start_pos = $pos - 1; # 0-based
            my $end_pos = $pos;
            my $flag_new = "POS".$flag;
            print O "$chr\t$start_pos\t$end_pos\t$flag_new\n";
        }
    }
}

close O;
