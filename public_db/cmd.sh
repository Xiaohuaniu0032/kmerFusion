cp /data1/workdir/wangce/database/humandb/cosmic/cosmic_transcript.txt ./
cp /data1/software/annovar/humandb/hg19_ensGene.txt ./

chmod 444 cosmic_transcript.txt
chmod 444 hg19_ensGene.txt

gzip cosmic_transcript.txt
gzip hg19_ensGene.txt
