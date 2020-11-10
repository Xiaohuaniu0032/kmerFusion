#zcat genefuse.R1.fq.gz|head -n 10000 >test.1000.R1.fq
#zcat genefuse.R2.fq.gz|head -n 10000 >test.1000.R2.fq

perl refFlat2geneBed.pl >all_genes.bed
