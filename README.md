# kmerFusion
detect DNA fusion using **`genefuse`** software, which is a DNA fusion caller based on kmer method.

## Usage
`perl /path/kmerFusion/kmerFusion.pl -bam <bam> -n <name> -od <outdir>`

### parameter spec

`-bam`: tumor bam file

`-n`: tumor name

`-od`: output dir

`-fa`: fasta ref

`-Plen`: padding lenght on soft-clip pos

## Testing
1. `cd /path/kmerFusion/test/`
2. `sh run.sh`
3. `nohup sh gfuse.ZZ20052708_CL01156.sh >gfuse.log 2>&1 &`

### output files
```
-rw-rw-r-- 1 fulongfei fulongfei     477849 11月 13 09:37 all.gene.bed
-rw-rw-r-- 1 fulongfei fulongfei    3234591 11月 13 15:53 fusion.csv
-rw-rw-r-- 1 fulongfei fulongfei       2328 11月 13 15:53 gfuse.ZZ20052708_CL01156.sh
-rw-rw-r-- 1 fulongfei fulongfei        176 11月 13 09:17 run.sh
-rw-rw-r-- 1 fulongfei fulongfei   18980591 11月 13 15:51 softclip.pos.PASS_FILTER.Annot_Genes.txt
-rw-rw-r-- 1 fulongfei fulongfei   13887873 11月 13 15:48 softclip.pos.PASS_FILTER.txt
-rw-rw-r-- 1 fulongfei fulongfei   28015211 11月 13 15:47 softclip.pos.txt
-rw-rw-r-- 1 fulongfei fulongfei       8673 11月 13 10:43 ZZ20052708_CL01156.fusion.xls
-rw-rw-r-- 1 fulongfei fulongfei    3503730 11月 13 10:14 ZZ20052708_CL01156.genefuse.html
-rw-rw-r-- 1 fulongfei fulongfei     279516 11月 13 10:14 ZZ20052708_CL01156.genefuse.json
-rw-rw-r-- 1 fulongfei fulongfei 3008045788 11月 13 09:32 ZZ20052708_CL01156.R1.fastq
-rw-rw-r-- 1 fulongfei fulongfei 3015106236 11月 13 09:32 ZZ20052708_CL01156.R2.fastq
-rw-rw-r-- 1 fulongfei fulongfei 1631348801 11月 13 09:29 ZZ20052708_CL01156.sort_by_name.bam
```

### file spec

`*.sort_by_name.bam`: bam file sort by read name

`*.R1.fastq`: fastq R1 converted from `*.sort_by_name.bam`

`*.R2.fastq`: fastq R2 converted from `*.sort_by_name.bam`

`softclip.pos.txt`: soft clip pos

`softclip.pos.PASS_FILTER.txt`: soft clip pos after filter by support count

`softclip.pos.PASS_FILTER.Annot_Genes.txt`: soft clip pos annot with gene name by bedtools intersect

`fusion.csv`: fusion gene list used by `genefuse` to search

`*.genefuse.html`: fusion result HTML file

`*.genefuse.json`: fusion result JSON file

`*.fusion.xls`: final fusion result XLS file

## Method
1. converted bam file into fq1 & fq2
2. get soft clip pos
3. filter soft clip pos if a pos has <=2 SoftClip support
4. annot soft clip pos with its gene name
5. make a dynamic fusion list
6. call `genefuse` software to detect DNA fusion
7. get final fusion result from JSON file created by `genefuse`

> *You can see `genefuse` papre for detailed method*

`https://github.com/OpenGene/GeneFuse`

## Why not use `genefuse` directly?
* `genefuse` need a `fusion.csv` file as input file, the genes in this file will be searched for fusion. 
* However, for example, if one sample has a rare `ALK` fusion in which the partner fusion gene of `ALK` is not in `fusion.csv` file, then `genefuse` will failed to detect this rare ALK fusion. 
* This pipeline just solved one problem: dynamically creat `fusion.csv` from a specified tumor BAM file.


## Why use kmer-based software to detect DNA fusion?
* `kmer` method is very **sensitive** to detect DNA fusion
* `genefuse` is a good kmer-based DNA fusion caller, and its visualizing function is very powerful and useful

## Software Needed
*you need to add `samtools` and `bedtools` binary path in your Linux PATH env*


