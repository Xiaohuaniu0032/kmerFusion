use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

my ($fusion_list,$sv_gene_bedtools,$padding_len,$outdir);

GetOptions(
    "flist:s" => \$fusion_list,           # MUST GIVEN
    "slpos:s" => \$sv_gene_bedtools,      # MUST GIVEN
    "len:i"   => \$padding_len,           # Default: 500
    "od:s"    => \$outdir,                # MUST GIVEN
    ) or die "unknown args\n";


# default value
if (not defined $padding_len){
    $padding_len = 500;
}



# first, let us record soft clipped position
my %soft_clip_pos;
open IN, "$sv_gene_bedtools" or die;
while (<IN>){
    chomp;
    my @arr = split /\t/;
    my $pos = $arr[1] + 1;
    $soft_clip_pos{$arr[-2]}{$arr[0]}{$pos} = 1; # gene=>chr=>pos=>1
}
close IN;
#print(Dumper(\%soft_clip_pos));


# then, parse fusion.list file to get its exon annot info
my ($gene,$chr,$region);

#my %gene_info; # a gene's chr/region info
my %fusion_list_exon; # 记录待过滤gene的exon信息

open IN, "$fusion_list" or die;
while (<IN>){
    chomp;
    next if /^$/;
    if (/^>/){
        $gene = (split /\,/, $_)[0];
        $gene =~ s/^>//;
        $chr = (split /\:/, $_)[0];
        $chr =~ s/(.*)\,//;
        $region = (split /\:/, $_)[1];

        #$gene_info{$gene}{$chr} = $region;
    }else{
        my @arr = split /\,/, $_;
        $fusion_list_exon{$gene}{$chr}{$arr[0]} = "$arr[1]\t$arr[2]"; # gene=>chr=>exon=>"start/end"
    }
}
close IN;
#print(Dumper(\%fusion_list_exon));


#my %fusion_list_intron; # 记录待过滤gene的intron信息

#foreach my $gene (keys %fusion_list_exon){
#    my @chrs = keys %{$fusion_list_exon{$gene}};
#    
#    if (@chrs > 1){
#        print "Warning: $gene exists on >=2 chroms, please pay attention\n";
#    }
#
#    for my $chr (@chrs){
#        my @exon_sort = sort {$a <=> $b} keys %{$fusion_list_exon{$gene}{$chr}};
#        
#        # if +/- strand?
#        my $ex_first = $exon_sort[0];
#        my $ex_last = $exon_sort[-1];
#
#        my $first_exon_region = $fusion_list_exon{$gene}{$chr}{$ex_first};
#        my $last_exon_region = $fusion_list_exon{$gene}{$chr}{$ex_last};
#
#        my @first_exon_region = split /\t/, $first_exon_region;
#        my @last_exon_region = split /\t/, $last_exon_region;
#
#        my $intron_num = $exon_sort[-1] - 1;
#        #print "gene $gene has $intron_num introns\n";
#        # determine +/- strand
#
#        if ($first_exon_region[0] > $last_exon_region[-1]){
#            # - strand
#            for my $e (1..$intron_num){
#                #print "processing Gene:$gene Chr:$chr Intron:$e\n";
#                my $e_later = $e + 1;
#                my $e_later_region = $fusion_list_exon{$gene}{$chr}{$e_later};
#
#                my @e = split /\t/, $fusion_list_exon{$gene}{$chr}{$e};
#                my @e_later = split /\t/, $fusion_list_exon{$gene}{$chr}{$e_later};
#
#                my $reg_start = $e_later[1] + 1;
#                my $reg_end = $e[0] - 1;
#
#                my $reg = "$reg_start\t$reg_end";
#                $fusion_list_intron{$gene}{$chr}{$e} = $reg;
#            }
#        }else{
#            # + strand
#            for my $e (1..$intron_num){
#                #print "$e\n";
#                my $e_later = $e + 1;
#                my $e_later_region = $fusion_list_exon{$gene}{$chr}{$e_later};
#
#                my @e = split /\t/, $fusion_list_exon{$gene}{$chr}{$e};
#                my @e_later = split /\t/, $fusion_list_exon{$gene}{$chr}{$e_later};
#
#                my $reg_start = $e[1] + 1;
#                my $reg_end = $e_later[0] - 1;
#
#                my $reg = "$reg_start\t$reg_end";
#                $fusion_list_intron{$gene}{$chr}{$e} = $reg;
#            }
#        }
#    }
#}



# 循环每个基因,判断是exon还是intron.
# 如果是intron,则需输出相邻exon信息.
# 如果是exon,则输出该exon信息.

# finally, for each soft clipped pos, output an entire item (pos + all exon annot)
my $final_fusion_list = "$outdir/fusion.list.pass.filter.csv";
open O, ">$final_fusion_list" or die;

my $merge_log = "$outdir/pos_merge.log";
open LOG, ">$merge_log" or die;

my $pos_flag = 0;
foreach my $gene (keys %soft_clip_pos){
    #next if ($gene ne "ALK" and $gene ne "EML4");

    my @chrs = keys %{$soft_clip_pos{$gene}};
    for my $chr (@chrs){
        if (exists $fusion_list_exon{$gene}{$chr}){
            my @pos_sort = sort {$a <=> $b} keys %{$soft_clip_pos{$gene}{$chr}};
            my $pos_sort_left = $pos_sort[0];
            my $pos_sort_right = $pos_sort[-1];
            my $P = "$pos_sort_left\_$pos_sort_right";
            my @PP;
            push @PP, $P;
            #my $merged_pos_aref = &merge_pos(\@pos_sort);
            #for my $pos (@{$merged_pos_aref}){
            for my $pos (@PP){
                #$pos_flag += 1;
                print LOG "processing Gene:$gene Chr:$chr Pos:$pos\n";
                if ($pos =~ /\_/){
                    my @seg = split /\_/, $pos;
                    my $pos_left = $seg[0] - $padding_len;
                    my $pos_right = $seg[-1] + $padding_len;

                    if ($pos_left <= 0){
                        $pos_left = 1; # pos can not be negtive value
                    }

                    # print sv region for genefuse to search
                    print O "\>$gene\,$chr\:$pos_left\-$pos_right\n";
                    #print O "\>$gene\_$pos_flag\,$chr\:$pos_left\-$pos_right\n";
                    # print all exon info
                    my @exon_sort = sort {$a <=> $b} keys %{$fusion_list_exon{$gene}{$chr}};
                    for my $e (@exon_sort){
                        my @r = split /\t/, $fusion_list_exon{$gene}{$chr}{$e};
                        print O "$e\,$r[0]\,$r[1]\n";
                    }

                    print O "\n";
                }else{
                    my $pos_left = $pos - $padding_len;
                    my $pos_right = $pos + $padding_len;

                    if ($pos_left <= 0){
                        $pos_left = 1;
                    }

                    print O "\>$gene\,$chr\:$pos_left\-$pos_right\n";
                    #print O "\>$gene\_$pos_flag\,$chr\:$pos_left\-$pos_right\n";

                    my @exon_sort = sort {$a <=> $b} keys %{$fusion_list_exon{$gene}{$chr}};
                    for my $e (@exon_sort){
                        my @r = split /\t/, $fusion_list_exon{$gene}{$chr}{$e};
                        print O "$e\,$r[0]\,$r[1]\n";
                    }

                    print O "\n";
                }
            }
        }
    }
}
close O;
close LOG;



sub merge_pos{
    my ($aref) = @_;
    my @merged_pos;
    my $first_pos = shift @{$aref};
    push @merged_pos, $first_pos;
    for my $pos (@{$aref}){
        my $former_pos = $merged_pos[-1];
        if ($former_pos =~ /\_/){
            my @former_pos = split /\_/, $former_pos;
            my $last_pos = $former_pos[-1];
            if ($pos - $last_pos <= 50){
                # 如果当前POS与前一个POS的距离小于50bp,那将这个POS合并为一个大的segment,格式如:100_104_132,256,321_322
                # OK to merge
                my $new_seg = $former_pos.'_'.$pos;
                pop @merged_pos; # remove last old segment
                push @merged_pos, $new_seg; # use new to replace old segment
            }else{
                # can not merge
                push @merged_pos, $pos;
            }
        }else{
            if ($pos - $former_pos <= 50){
                # OK to merge
                my $new_seg = $former_pos.'_'.$pos;
                pop @merged_pos;
                push @merged_pos, $new_seg;
            }
        }
    }

    return(\@merged_pos);
}



#my %sv_gene_info;

#foreach my $gene (keys %soft_clip_pos){
#    #print "$gene\n";
#    my @chrs = keys %{$soft_clip_pos{$gene}};
#
#    if (@chrs > 1){
#        my $n = @chrs;
#        print "Warning: gene $gene occured on $n chroms, please pay attention\n";
#    }
#
#    for my $chr (@chrs){
#        my @pos = sort {$a <=> $b} keys %{$soft_clip_pos{$gene}{$chr}}; # get all soft clip pos
#        for my $pos (@pos){
#            print "processing Gene:$gene Chr:$chr Pos:$pos\n";
#            # 循环该基因的exon/intron区域,判断该pos是哪个区域
#            
#            # for exon region
#            my $where_pos; # format: exon/5;intron/6
#            my @exon_sort = sort {$a <=> $b} keys %{$fusion_list_exon{$gene}{$chr}};
#            for my $exon (@exon_sort){
#                # get region
#                my @exon_region = split /\t/, $fusion_list_exon{$gene}{$chr}{$exon};
#
#                if ($pos >= $exon_region[0] and $pos <= $exon_region[1]){
#                    # pos in this region
#                    $sv_gene_info{$gene}{$chr}{"exon"}{$exon} = 1;
#                    print "Gene:$gene Chr:$chr Pos:$pos FIND in Exon($exon) <$exon_region[0]\-$exon_region[1]>\n"
#                }else{
#                    print "Gene:$gene Chr:$chr Pos:$pos NOT FIND in Exon($exon) <$exon_region[0]\-$exon_region[1]>\n";
#                }
#            }
#
#            my @intron_sort = sort {$a <=> $b} keys %{$fusion_list_intron{$gene}{$chr}};
#            for my $i (@intron_sort){
#                # get region
#                my @intron_region = split /\t/, $fusion_list_intron{$gene}{$chr}{$i};
#
#                if ($pos >= $intron_region[0] and $pos <= $intron_region[1]){
#                    # pos in this region
#                    $sv_gene_info{$gene}{$chr}{"intron"}{$i} = 1;
#                    print "Gene:$gene Chr:$chr Pos:$pos FIND in INTRON($i) <$intron_region[0]\-$intron_region[1]>\n";
#                }else{
#                    print "Gene:$gene Chr:$chr Pos:$pos NOT FIND in INTRON($i) <$intron_region[0]\-$intron_region[1]>\n#";#
#                }#
#            }#
#        }#
#    }
#}

#print(Dumper(\%sv_gene_info));


#my %final_sv_region_to_search;

#foreach my $gene (keys %sv_gene_info){
#    my @chrs = keys %{$sv_gene_info{$gene}};
#
#    for my $chr (@chrs){
#        if (exists $gene_info{$gene}{$chr}){
#            # check exon
#            if (exists $sv_gene_info{$gene}{$chr}{"exon"}){
#                # just output directly
#                my @exon_sort = sort {$a <=> $b} keys %{$sv_gene_info{$gene}{$chr}{"exon"}};
#                for my $e (@exon_sort){
#                    my @exon_region = split /\t/, $fusion_list_exon{$gene}{$chr}{$e};
#                    my $reg = "$exon_region[0]\-$exon_region[1]";
#                    #$final_sv_region_to_search{$gene}{$chr} = "$exon_region[0]\-$exon_region[1]";
#                    $final_sv_region_to_search{$gene}{$chr}{$reg} = 1;
#                }
#            }
#
#            # check intron
#            if (exists $sv_gene_info{$gene}{$chr}{"intron"}){
#                # 输出该intron相邻的exon
#                my @intron_sort = sort {$a <=> $b} keys %{$sv_gene_info{$gene}{$chr}{"intron"}};
#                for my $intron (@intron_sort){
#                    my @intron_region = split /\t/, $fusion_list_intron{$gene}{$chr}{$intron};
#                    my $reg = "$intron_region[0]\-$intron_region[1]";
#                    $final_sv_region_to_search{$gene}{$chr}{$reg} = 1;
#                }
#            }
#        }
#    }
#}


# output final fusion csv list
#my $final_fusion_list = "$outdir/fusion.list.pass.filter.csv";
#open O, ">$final_fusion_list" or die;


#foreach my $gene (keys %final_sv_region_to_search){
#    my @chrs = keys %{$final_sv_region_to_search{$gene}};
#    for my $chr (@chrs){
#        my @regs = keys %{$final_sv_region_to_search{$gene}{$chr}};
#        for my $r (@regs){
#            # 对每个可能的soft-clip POS,左右各扩500bp,输出完整的一个fusion csv item
#            print O "\>$gene\,$chr\:$r\n";
#
#            # print all exon info
#            my @exon_sort = sort {$a <=> $b} keys %{$fusion_list_exon{$gene}{$chr}};
#            for my $e (@exon_sort){
#                my @exon_reg = split /\t/, $fusion_list_exon{$gene}{$chr}{$e};
#                print O "$e\,$exon_reg[0]\,$exon_reg[1]\n";
#            }
#            print O "\n";
#        }
#    }
#}
#close O;

