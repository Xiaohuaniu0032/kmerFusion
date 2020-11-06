bwa mem -M -R "@RG\tID:GeneFuse\tSM:test_sample\tPL:illumina" /data1/database/b37/human_g1k_v37.fasta /home/fulongfei/workdir/FusionScan/test/data/R1.fq /home/fulongfei/workdir/FusionScan/test/data/R2.fq | samtools view -b -o genefuse.test.bam -

samtools sort genefuse.test.bam -o genefuse.test.sort.bam
samtools index genefuse.test.sort.bam
