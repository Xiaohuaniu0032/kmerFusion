#sh get_sv_pos.01.sh /data3/Projects/panel338/200530_A00682_0367_AH3N5WDSXY/02_aln/ZZ20052708_CL01156.sv_rmdup.bam $PWD
#perl filter_sv_pos.02.pl -f sv.pos.ALL.txt -od $PWD
#bedtools intersect -a sv.pos.PASS_FILTER.txt -b /home/fulongfei/workdir/git_repo/kmerFusion/bin/bak/all_genes.bed -wo -bed >sv.genes.txt
#perl make_genefuse_csv_list.03.pl -gene sv.genes.txt -od $PWD
#perl filter_csv_list.04.pl fusion.list sv.genes.txt $PWD
#perl filter_csv_list.04.pl -flist fusion.list -slpos sv.genes.txt -len 1000 -od $PWD

name='ZZ20052708_CL01156'
flist='/home/fulongfei/workdir/git_repo/kmerFusion/bin/fusion.list.pass.filter.csv'
#flist='/home/fulongfei/workdir/git_repo/kmerFusion/data/fusions.csv'
#flist='/home/fulongfei/workdir/git_repo/kmerFusion/bin/fusions.csv'
#flist='/home/fulongfei/workdir/git_repo/kmerFusion/bin/fusions.csv.ALK.EML4'

/data1/workdir/fulongfei/git_repo/kmerFusion/genefuse -t 8 -r /data1/database/b37/human_g1k_v37.fasta -f $flist -1 /home/fulongfei/workdir/fusion_work/kmerFusion_test_20201106/ZZ20052708_CL01156/ZZ20052708_CL01156.R1.fastq -2 /home/fulongfei/workdir/fusion_work/kmerFusion_test_20201106/ZZ20052708_CL01156/ZZ20052708_CL01156.R2.fastq -h $PWD/${name}.genefuse.html -j $PWD/${name}.genefuse.json
