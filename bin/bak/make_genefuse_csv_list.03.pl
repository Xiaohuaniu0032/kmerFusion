use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Data::Dumper;

my ($covered_gene_file,$gene_annot_file,$main_NM,$chrNaming,$outdir);

GetOptions(
    "gene:s" => \$covered_gene_file,            # NEED <sv.pos.PASS_FILTER.txt>
    "annot:s" => \$gene_annot_file,             # Default: </home/fulongfei/workdir/git_repo/GeExCNV/public_db/refFlat.txt>
    "NM:s" => \$main_NM,                        # Default: </home/fulongfei/workdir/git_repo/GeExCNV/public_db/MainNM_yrt_20191111.txt>
    "chrname:s" => \$chrNaming,                 # Default: no_chr_prefix
    "od:s" => \$outdir,                         # NEED
    ) or die;



# Default value
if (not defined $gene_annot_file){
    $gene_annot_file = "/home/fulongfei/workdir/git_repo/GeExCNV/public_db/refFlat.txt";
}

if (not defined $main_NM){
    $main_NM = "/home/fulongfei/workdir/git_repo/GeExCNV/public_db/MainNM_yrt_20191111.txt";
}

if (not defined $chrNaming){
    $chrNaming = "no_chr_prefix";
}


my $fusion_csv = "$outdir/fusion.list";

open O, ">$fusion_csv" or die;

my %gene;
open IN, "$covered_gene_file" or die;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    $gene{$arr[-2]} = 1;
}
close IN;

print "Finished reading covered_gene_file\n";

my @genes = keys %gene; # all genes, both genes in BED and genes related SV region
my $n_gene = scalar(@genes);
print "find $n_gene genes according the BAM file\n";

# get each gene's exon info
my %exon_info;
my %region_info;
my %gene2chr;

open IN, "$gene_annot_file" or die;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    my $chr;
    if ($chrNaming eq "with_chr_prefix"){
        $chr = $arr[2];
    }else{
        my $chr_temp = $arr[2];
        $chr_temp =~ s/^chr//;
        $chr = $chr_temp;
    }

    $gene2chr{$arr[0]} = $chr; # chr info

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
            $region_info{$arr[0]}{$arr[1]}{"region"} = "$gene_start\-$gene_end"; # gene=>NM=>
            $exon_info{$arr[0]}{$arr[1]}{$exon_idx} = "$start\t$end";
        }
    }else{
        for (my $i=1;$i<=$exon_num;$i++){
            my $start = $start_exon[$i-1];
            my $end = $end_exon[$i-1];
            $exon_idx = $exon_num - $i + 1;

            $region_info{$arr[0]}{$arr[1]}{"region"} = "$gene_start\-$gene_end"; # gene=>NM=>
            $exon_info{$arr[0]}{$arr[1]}{$exon_idx} = "$start\t$end";
        }
    }
}
close IN;




## NM
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
    my $nm;
    if (exists $gene_NM{$gene}){
        $nm = $gene_NM{$gene};
    }else{
        $nm = "NA";
        print "gene does not have NM in main_NM file, will be skipped\n";
        next;
    }

    # check if refFlat exists this gene and its NM
    if (exists $region_info{$gene}){
        if (exists $region_info{$gene}{$nm}){
            # get this gene's chr/pos/exon info
            print "make csv for $gene\n";
            my $chr = $gene2chr{$gene};
            my $region = $region_info{$gene}{$nm}{"region"};
            print O "\>$gene\,$chr\:$region\n";
            
            foreach my $exon (sort {$a <=> $b} keys %{$exon_info{$gene}{$nm}}){
                my $ex = $exon_info{$gene}{$nm}{$exon};
                my @ex = split /\t/, $ex;
                print O "$exon\,$ex[0]\,$ex[1]\n";
            }
            print O "\n";
        }else{
            print "$gene NM $nm does not exists in refFlat DB, will be skipped\n";
            next;
        }
    }else{
        print "$gene does not exists in refFlat DB, will be skipped\n";
        next;
    }
}

close O;

