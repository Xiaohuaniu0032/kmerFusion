### bwa align
# data download from: http://opengene.org/dataset.html
fq1='/data1/workdir/fulongfei/git_repo/kmerFusion/data/genefuse.R1.fq.gz'
fq2='/data1/workdir/fulongfei/git_repo/kmerFusion/data/genefuse.R2.fq.gz'

#bwa mem -M -R "@RG\tID:GeneFuse\tSM:test_sample\tPL:illumina" /data1/database/b37/human_g1k_v37.fasta $fq1 $fq2 | samtools view -b -o genefuse.test.bam -

#samtools sort genefuse.test.bam -o genefuse.test.sort.bam
#samtools index genefuse.test.sort.bam

# genefuse
# please note the chr naming in your bam and in *csv file
flist='/home/fulongfei/workdir/git_repo/kmerFusion/data/fusions.csv'
genefuse -r /data1/database/b37/human_g1k_v37.fasta -f $flist -1 $fq1 -2 $fq2 -h test.html -j test.json
