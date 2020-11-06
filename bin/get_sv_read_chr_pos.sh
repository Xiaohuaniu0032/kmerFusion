bam=$1
od=$2

samtools view $bam | awk '$6~/S/' | awk '{print $3"\t"$4}' >$od/soft_clip_read.pos
samtools view $bam | awk '$6~/H/' | awk '{print $3"\t"$4}' >$od/hard_clip_read.pos

cat $od/soft_clip_read.pos $od/hard_clip_read.pos | sort | uniq -c | awk '{print $2"\t"$3}' >$od/sv_read.pos
