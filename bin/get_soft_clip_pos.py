import pysam
import sys
import argparse
import re

def parse_args():
    AP = argparse.ArgumentParser("extract soft-clip pos from bam file")
    AP.add_argument('-bam',help='tumor bam file',dest='bam')
    AP.add_argument('-od',help='output dir',dest='outdir')

    return AP.parse_args()

def count_cigar_S(cigar_list):
    s_num = 0
    for i in cigar_list:
        if i[0] == 4:
            s_num += 1

    return s_num


def main():
    args = parse_args()
    samfile = pysam.AlignmentFile(args.bam,'rb')

    outfile = args.outdir + '/softclip.pos.txt'
    of = open(outfile,'w')

    for read in samfile.fetch():
        if read.is_duplicate or read.is_qcfail or read.is_secondary or read.is_supplementary:
            continue

        if read.mapping_quality < 10:
            continue

        chrom = read.reference_name
        left_aln_pos = read.reference_start + 1

        cigar_string = read.cigarstring
        cigar_tuple = read.cigartuples

        '''
        how many S in cigar?

        1. only one S, can be left or right
        2. two S, both left and right

        '''

        s_num = count_cigar_S(cigar_tuple)

        if s_num == 0:
            continue
        else:
            if s_num == 1:
                # check left or right
                pos = ''

                first_item = cigar_tuple[0]
                last_item = cigar_tuple[-1]

                if cigar_tuple[0][0] == 4:
                    # pos is left
                    pos = left_aln_pos
                else:
                    # pos is right
                    # count how many ref base(s) consumed by CIGAR
                    ref_consume_len = 0
                    cigar_tuple.pop() # remove last item
                    for i in cigar_tuple:
                        '''
                        M   BAM_CMATCH      0
                        I   BAM_CINS        1
                        D   BAM_CDEL        2
                        N   BAM_CREF_SKIP   3
                        S   BAM_CSOFT_CLIP  4
                        H   BAM_CHARD_CLIP  5
                        P   BAM_CPAD        6
                        =   BAM_CEQUAL      7
                        X   BAM_CDIFF       8
                        B   BAM_CBACK       9

                        M/D/N/=/X will consume ref
                        '''
                        if i[0] == 0 or i[0] == 2 or i[0] == 3 or i[0] == 7 or i[0] == 8:
                            # will consume ref
                            ref_consume_len += i[1]

                    pos = left_aln_pos + ref_consume_len - 1

                val = "%s\t%s" % (chrom,pos)
                of.write(val+'\n')
            
            else:
                # left and right both contains S
                
                # for left S
                pos = left_aln_pos
                val = "%s\t%s" % (chrom,pos)
                of.write(val+'\n')

                # for right S
                ref_consume_len = 0
                cigar_tuple.pop() # remove last item
                for i in cigar_tuple:
                    # S do not consume ref
                    if i[0] == 0 or i[0] == 2 or i[0] == 3 or i[0] == 7 or i[0] == 8:
                        ref_consume_len += i[1]

                pos = left_aln_pos + ref_consume_len - 1
                val = "%s\t%s" % (chrom,pos)
                of.write(val+'\n')

    of.close()


if __name__ == "__main__":
    main()











