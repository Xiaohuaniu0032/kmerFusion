use strict;
use warnings;

my ($cosmicGeneFile,$ensGeneAnnotFile,$softclipPosPASSAnnotGeneFile,$chrname,$padding_len,$outdir) = @ARGV;

# link enst with gene name
my %gene2enst;
open IN, "gunzip -dc $cosmicGeneFile |" or die;
while (<IN>){
	chomp;
	my @arr = split /\t/;
	$gene2enst{$arr[0]} = $arr[1]; # gene=>enst
}
close IN;

# get enst exon region info
my %enstExonAnnotInfo;
my %enst_chr;

open IN, "gunzip -dc $ensGeneAnnotFile |" or die;
while (<IN>){
	chomp;
	my @arr = split /\t/;
	my $enst = $arr[1];
	my $strand = $arr[3];
	
	my $chr = $arr[2];
	if ($chrname eq "no_chr_prefix"){
		$chr =~ s/^chr//;
	}

	$enst_chr{$enst} = $chr; # chr for fusion csv file
	
	my $start_exon = $arr[-7];
	$start_exon =~ s/\,$//;

	my $end_exon = $arr[-6];
	$end_exon =~ s/\,$//;
	
	my @start_exon = split /\,/, $start_exon;
	my @end_exon = split /\,/, $end_exon;
	
	my $exon_num = @start_exon;
	
	my $exon_idx = 0;
	
	if ($strand eq "+"){
		for (my $i=1;$i<=$exon_num;$i++){
			my $start = $start_exon[$i-1];
			my $end = $end_exon[$i-1];
			$exon_idx += 1;
			$enstExonAnnotInfo{$enst}{$exon_idx} = "$start\t$end"; # enst=>exon=>region
		}
	}else{
		for (my $i=1;$i<=$exon_num;$i++){
			my $start = $start_exon[$i-1];
			my $end = $end_exon[$i-1];
			$exon_idx = $exon_num - $i + 1;
			$enstExonAnnotInfo{$enst}{$exon_idx} = "$start\t$end";
		}
	}
}
close IN;

# get each gene's soft clip region
my %soft_clip_pos;
my %genes;

open IN, "$softclipPosPASSAnnotGeneFile" or die;
while (<IN>){
	chomp;
	my @arr = split /\t/;
	my $pos = $arr[1] + 1;
	$soft_clip_pos{$arr[-2]}{$pos} = 1; # gene=>pos=>1
	$genes{$arr[-2]} = 1;
}
close IN;

# make final fusion csv list
my $fusion_list = "$outdir/fusion.csv";
open O, ">$fusion_list" or die;

foreach my $gene (keys %genes){
	my $enst = $gene2enst{$gene}; # get enst
	my $chr = $enst_chr{$enst};
	my @pos_sort = sort {$a <=> $b} keys %{$soft_clip_pos{$gene}};

	# if pos2 - pos1 > 2 * padding_len, then these two pos can not merge into one, otherwise, can merge into one new pos
	my $merged_pos_aref = &merge_pos(\@pos_sort, $padding_len);

	for my $pos (@{$merged_pos_aref}){
		my ($pos_sort_left,$pos_sort_right);

		if ($pos =~ /\_/){
			# get first and last item
			my @pos = split /\_/, $pos;
			my $first_item = $pos[0];
			my $last_item = $pos[-1];

			$pos_sort_left = $first_item - $padding_len;
			$pos_sort_right = $last_item + $padding_len;

			if ($pos_sort_left <= 0){
				$pos_sort_left = 1;
			}
		}else{
			# single pos
			$pos_sort_left = $pos - $padding_len;
			$pos_sort_right = $pos + $padding_len;
		}

		if ($pos_sort_left <= 0){
			$pos_sort_left = 1;
		}

		# print sv region for genefuse to search
		print O "\>$gene\,$chr\:$pos_sort_left\-$pos_sort_right\n";

		# print all exon info
		my @exon_sort = sort {$a <=> $b} keys %{$enstExonAnnotInfo{$enst}};
		for my $e (@exon_sort){
			my @r = split /\t/, $enstExonAnnotInfo{$enst}{$e};
			print O "$e\,$r[0]\,$r[1]\n";
		}

		print O "\n";
	}
}

close O;


sub merge_pos{
	my ($aref, $gap_len) = @_;
	my $new_gap_len = $gap_len + 100; # add a small value to avoid
	my @new_pos;
	my $first_pos = shift @{$aref};
	push @new_pos, $first_pos;

	for my $pos (@{$aref}){
		my $former_pos = $new_pos[-1];
		if ($former_pos =~ /\_/){
			my $former_item = (split /\_/, $former_pos)[-1];
			if ($pos - $former_item > 2 * $new_gap_len){
				# can not merge
				push @new_pos, $pos;
			}else{
				# can merge
				pop @new_pos; # remove last seg
				my $new_seg = $former_pos.'_'.$pos;
				push @new_pos, $new_seg; # replace last seg with new seg
			}
		}else{
			if ($pos - $former_pos > 2 * $new_gap_len){
				push @new_pos, $pos;
			}else{
				pop @new_pos;
				my $new_seg = $former_pos.'_'.$pos;
				push @new_pos, $new_seg;
			}
		}
	}

	return(\@new_pos);
}

