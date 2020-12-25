import pysam
import sys
import argparse
import re

def parse_args():
    AP = argparse.ArgumentParser("check the uniqness of alignment of fusion support split reads")
    AP.add_argument('-bam',help='tumor bam',dest='bam')
    AP.add_argument('-f',help='*.fusion.xls file',dest='fsFile')
    AP.add_argument('-mpq',help='map quality cutoff',dest='mpq',default=30)
    AP.add_argument('-refFlag',help='refFlat to annot gene strand info',dest='refFlat',default='/home/fulongfei/workdir/git_repo/GeExCNV/public_db/refFlat.txt')
    AP.add_argument('-od',help='output dir',dest='outdir')

    return AP.parse_args()


def main():
    args = parse_args()
    fsFH = open(args.fsFile,'r')


    # get gene strand info
    gene_strand = {}
    refFlatFH = open(args.refFlat,'r')
    for line in refFlatFH:
        arr = line.strip().split('\t')
        gene_strand[arr[0]] = arr[3]

    refFlatFH.close()


    for line in fsFH:
        arr = line.strip().split('\t')
        if arr[0] == 'SampleID':
            # skip header
            continue

        if arr[9] != 'YES':
            # skip not report fusion
            continue

        print("[INFO]: check %s fusion info..." % (arr[1]))

        hgene = arr[1].split('->')[0]
        tgene = arr[1].split('->')[1]

        # fusion gene pos info
        fs_info = {}
        fs_info[arr[10]] = int(arr[12]) # gene -> pos
        fs_info[arr[16]] = int(arr[18])

        samfile = pysam.AlignmentFile(args.bam,'rb')
        
        bp1_left = int(arr[12]) - 10
        bp1_right = int(arr[12]) + 10

        
        hgene_strand = gene_strand[hgene]
        print("[INFO]: check Hgene %s (%s) info..." % (hgene,hgene_strand))
        
        hgene_bp = fs_info[hgene]

        hgene_left = hgene_bp - 10
        hgene_right = hgene_bp + 10

        hgene_strand = gene_strand[hgene]
        tgene_strand = gene_strand[tgene]

        for read in samfile.fetch(arr[11],hgene_left,hgene_right):
            if read.is_read1:
                r1r2 = 'READ1'
            else:
                r1r2 = 'READ2'

            if read.is_secondary or read.is_supplementary or read.is_duplicate:
                # skip
                continue

            cigarS = read.cigarstring

            if read.is_reverse:
                read_rev = 'REVERSE'
            else:
                read_rev = 'FORWARD'


            if read.has_tag('SA'):
                sa_tag = read.get_tag('SA')

                # count how many SA
                sa_count = sa_tag.count(';')
                sa_strand = sa_tag.split(',')[2] # +/-
                if sa_count > 1:
                    # skip
                    continue

                # skip mapQ < 10
                mapQ = int(sa_tag.split(',')[4])
                if mapQ < 10:
                    continue

                # check Hgene's sup read's ori
                if hgene_strand == '+' and tgene_strand == '-':
                    # EML4->ALK
                    # https://cancer.sanger.ac.uk/cosmic/fusion/overview?fid=180308&gid=150139
                    
                    if cigarS[-1] != 'S':
                        continue

                    # for hgene
                    if read_rev == 'FORWARD':
                        # split read need on '-'
                        if sa_strand == '-':
                            pass
                        else:
                            continue
                    else:
                        # need on +
                        if sa_strand == '+':
                            pass
                        else:
                            continue


                if hgene_strand == '+' and tgene_strand == '+':
                    # FGFR3->TACC3
                    # https://cancer.sanger.ac.uk/cosmic/fusion/overview?fid=90231&gid=9677

                    if cigarS[-1] != 'S':
                        continue

                    # for hgene
                    if read_rev == 'FORWARD':
                        if sa_strand == '+':
                            pass
                        else:
                            continue
                    else:
                        if sa_strand == '-':
                            pass
                        else:
                            continue


                if hgene_strand == '-' and tgene_strand == '-':
                    # CD74->ROS1
                    # https://cancer.sanger.ac.uk/cosmic/fusion/overview?fid=141941&gid=72432

                    if cigarS[-1] != 'M':
                        continue

                    # for hgene
                    if read_rev == 'FORWARD':
                        # split read need on '+'
                        if sa_strand == '+':
                            pass
                        else:
                            continue
                    else:
                        # split read need on '-'
                        if sa_strand == '-':
                            pass
                        else:
                            continue

                if hgene_strand == '-' and tgene_strand == '+':
                    # KIF5B->RET
                    # https://cancer.sanger.ac.uk/cosmic/fusion/overview?fid=9508&gid=190031

                    if cigarS[-1] != 'M':
                        continue

                    # for hgene
                    if read_rev == 'FORWARD':
                        # need on '-'
                        if sa_strand == '-':
                            pass
                        else:
                            continue
                    else:
                        # need on '+'
                        if sa_strand == '+':
                            pass
                        else:
                            continue




                    

                # infer breakpoint pos
                left_pos = read.reference_start + 1 # 0-based leftmost coordinate
                cigar_tuple = read.cigartuples
                
                ref_consume_len = 0
                for i in cigar_tuple:
                    if i[0] == 0 or i[0] == 2 or i[0] == 3 or i[0] == 7 or i[0] == 8:
                        ref_consume_len += i[1]

                bnd_pos = left_pos + ref_consume_len - 1 + 1

                print(read.query_name,r1r2,read.reference_name,bnd_pos,cigarS,read_rev,sa_tag)

        






        
        print("[INFO]: check Tgene %s info..." % (tgene))
        
        tgene_bp = fs_info[tgene]
        tgene_left = tgene_bp - 10
        tgene_right = tgene_bp + 10


        for read in samfile.fetch(arr[11],tgene_left,tgene_right):
            if read.is_read1:
                r1r2 = 'READ1'
            else:
                r1r2 = 'READ2'

            if read.is_secondary or read.is_supplementary or read.is_duplicate:
                continue

            cigarS = read.cigarstring

            if read.is_reverse:
                read_rev = 'REVERSE'
            else:
                read_rev = 'FORWARD'

            if read.has_tag('SA'):
                sa_tag = read.get_tag('SA')
                sa_strand = sa_tag.split(',')[2] # +/-

                # count how many SA
                sa_count = sa_tag.count(';')

                if sa_count > 1:
                    # skip
                    continue

                # skip mapQ < 10
                mapQ = int(sa_tag.split(',')[4])
                if mapQ < 10:
                    continue


                # check Tgene's sup read's ori
                if hgene_strand == '+' and tgene_strand == '-':
                    # EML4->ALK
                    # https://cancer.sanger.ac.uk/cosmic/fusion/overview?fid=180308&gid=150139
                    
                    # for tgene
                    if cigarS[-1] != 'S':
                        # skip
                        continue

                    # for Tgene
                    if read_rev == 'FORWARD':
                        # split read need on '-'
                        if sa_strand == '-':
                            pass
                        else:
                            continue
                    else:
                        if sa_strand == '+':
                            pass
                        else:
                            continue


                if hgene_strand == '+' and tgene_strand == '+':
                    # FGFR3->TACC3
                    # https://cancer.sanger.ac.uk/cosmic/fusion/overview?fid=90231&gid=9677

                    # for tgene
                    if cigarS[-1] != 'M':
                        # skip
                        continue

                    # for Tgene
                    if read_rev == 'FORWARD':
                        if sa_strand == '+':
                            pass
                        else:
                            continue
                    else:
                        if sa_strand == '-':
                            pass
                        else:
                            continue


                if hgene_strand == '-' and tgene_strand == '-':
                    # CD74->ROS1
                    # https://cancer.sanger.ac.uk/cosmic/fusion/overview?fid=141941&gid=72432

                    # for tgene
                    if cigarS[-1] != 'S':
                        continue

                    # for Tgene
                    if read_rev == 'FORWARD':
                        # split read need on '+'
                        if sa_strand == '+':
                            pass
                        else:
                            continue
                    else:
                        # split read need on '-'
                        if sa_strand == '-':
                            pass
                        else:
                            continue

                if hgene_strand == '-' and tgene_strand == '+':
                    # KIF5B->RET
                    # https://cancer.sanger.ac.uk/cosmic/fusion/overview?fid=9508&gid=190031

                    # for tgene
                    if cigarS[-1] != 'M':
                        continue

                    # for Tgene
                    if read_rev == 'FORWARD':
                        # need on '-'
                        if sa_strand == '-':
                            pass
                        else:
                            continue
                    else:
                        # need on '+'
                        if sa_strand == '+':
                            pass
                        else:
                            continue

                

                # infer breakpoint pos
                left_pos = read.reference_start + 1 # 0-based leftmost coordinate
                cigar_tuple = read.cigartuples
                ref_consume_len = 0
                for i in cigar_tuple:
                    if i[0] == 0 or i[0] == 2 or i[0] == 3 or i[0] == 7 or i[0] == 8:
                        ref_consume_len += i[1]

                bnd_pos = left_pos + ref_consume_len - 1 + 1
                print(read.query_name,r1r2,read.reference_name,bnd_pos,cigarS,read_rev,sa_tag)





        samfile.close()
    fsFH.close()


if __name__ == "__main__":
    main()



