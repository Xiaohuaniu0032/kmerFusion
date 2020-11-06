use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my ($sv_read_pos,$gene_annot_file,$main_NM,$chrNaming,$outfile);

GetOptions(
    "svRead:s" => \$sv_read_pos,
    "annot:s" => \$gene_annot_file,
    "NM:s" => \$main_NM,
    "chrname:s" => \chrNaming,
    "of:s" => \$outfile,
    ) or die;


# get each gene's start and end pos
my %gene_pos;
my %gene_chr;
open IN, "$gene_annot_file" or die;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    push @{$gene_pos{$arr[0]}}, $arr[4]; # a gene may have many transcripts, we just use the max region of this gene
    push @{$gene_pos{$arr[0]}}, $arr[5];
    my $chr;
    if ($chrNaming eq "with_chr_prefix"){
        $chr = $arr[2];
    }else{
        $chr1 = $arr[2];
        $chr1 =~ s/^chr//;
        $chr = $chr1;
    }
    $gene_chr{$arr[0]} = $chr;
}
close IN;


my %sv_pos;
open IN, "$sv_read_pos" or die;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    $sv_pos{$arr[0]}{$arr[1]} = 1; # chr/pos
}
close IN;


my %occur_gene;
foreach my $gene (keys %gene_pos){
    my $chr = $gene_chr{$gene};
    my @pos_sort = sort {$a <=> $b} @{$gene_pos{$gene}};
    for my $pos ($pos_sort[0]..$pos_sort[-1]){
        if (exists $sv_pos{$chr}{$pos}){
            $occur_gene{$gene} = 1;
        }
    }
}


my %gene_NM_exon_info;
open IN, "$gene_annot_file" or die;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    
    my $chr;
    if ($chrNaming eq "with_chr_prefix"){
        $chr = $arr[2];
    }else{
        $chr1 = $arr[2];
        $chr1 =~ s/^chr//;
        $chr = $chr1;
    }

    my $gene_start = $arr[4];
    my $gene_end = $arr[5];

    my $strand = $arr[3];
    
    my $start_exon = $arr[-2];
    $start_exon =~ s/\,$//;

    my $end_exon = $arr[-1];
    $end_exon =~ s/\,$//;

    my @start_exon = split /\,/, $start_exon;
    my @end_exon = split /\,/, $end_exon;

    my $exon_num = $arr[-3];
    my $exon_idx = 0;

    my $nm = $arr[1]; # NM

    if ($strand eq "+"){
        for (my $i=1;$i<=$exon_num;$i++){
            my $start = $start_exon[$i-1];
            my $end = $end_exon[$i-1];
            $exon_idx += 1;
            $gene_NM_exon_info{$arr[0]}{$chr}{$nm}{"pos"} = "$gene_start\-$gene_end";
            $gene_NM_exon_info{$arr[0]}{$chr}{$nm}{"exon"}{$exon_idx} = "$start\t$end";
        }
    }else{
        for (my $i=1;$i<=$exon_num;$i++){
            my $start = $start_exon[$i-1];
            my $end = $end_exon[$i-1];
            $exon_idx = $exon_num - $i + 1;
            $gene_NM_exon_info{$arr[0]}{$chr}{$nm}{"pos"} = "$gene_start\-$gene_end";
            $gene_NM_exon_info{$arr[0]}{$chr}{$nm}{"exon"}{$exon_idx} = "$start\t$end";
    }
}
close IN;


open O, ">$outfile" or die;


# NM
my %gene_NM;
open IN, "$main_NM" or die;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    $gene_NM{$arr[0]} = $arr[1]; # chr=>NM;
}
close IN;


# if main NM is in annot file's NM col? this may be a problem
# need to validate


# make csv list
for my $gene (@genes){
    if (exists $gene_NM{$gene}){
        # exists NM
        my $nm = $gene_NM{$gene};
        if (exists $gene_NM_exon_info{$gene}){
            # exists annot info
            my $chr = $gene_chr{$gene};
            my $pos = $gene_NM_exon_info{$gene}{$chr}{$nm}{"pos"};
            my @pos = split /\-/, $pos;
            print O "\>$gene\,$chr\:$pos[0]\-$pos[1]\n";
            my @exon_sort = sort {$a <=> $b} @{$gene_NM_exon_info{$gene}{$chr}{$nm}{"exon"}};
            for my $ex (@exon_sort){
                my $e = $gene_NM_exon_info{$gene}{$chr}{$nm}{"exon"}{$ex};
                my @e = split /\t/, $e;
                print O "$ex\,$e[0]\,$e[1]\n";
            }
        }
    }
}

close O;








