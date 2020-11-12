bam=$1
outdir=$2

samtools view $bam | awk '$6~/S/ {print $3"\t"$4}' >$outdir/sv.pos.ALL.txt

