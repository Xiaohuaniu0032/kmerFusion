import pysam
import sys
import argparse
import re
import os

def parse_args():
    AP = argparse.ArgumentParser("check the uniqness of alignment of fusion support split reads")
    AP.add_argument('-bam',help='tumor bam',dest='bam')
    #AP.add_argument('-n',help='sample name',dest='name')
    AP.add_argument('-f',help='*.fusion.xls file',dest='fsFile')
    AP.add_argument('-mq',help='map quality cutoff',dest='mapq',default=60)
    AP.add_argument('-refFlag',help='refFlat to annot gene strand info',dest='refFlat',default='/home/fulongfei/workdir/git_repo/GeExCNV/public_db/refFlat.txt')
    AP.add_argument('-od',help='output dir',dest='outdir')

    return AP.parse_args()


def main():
    args = parse_args()
    fsFH = open(args.fsFile,'r')
    name = os.path.basename(args.bam).split('.')[0]
    outfile = '%s/%s.sr_mapq_uniqness_check.txt' % (args.outdir,name)

    ofFH = open(outfile,'w')
    # fusion gene | gene | gene_type | read name | READ1/2 | chr | sv_pos | split_read_pos | pos_check | SR MAPQ | mapQ check | final result (check pos and mapq)
    outfile_header = ['fusionGene','Gene','GeneType(Hgene or Tgene)','read_name','READ1/2','chr','sv_pos','sr_pos','pos_check','sr_mapq','sr_mapq_check','final_check']
    ofFH.write('\t'.join(outfile_header)+'\n') # write outfile header

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

        hgene = arr[1].split('-')[0] # not use ->
        tgene = arr[1].split('-')[1]

        # fs chrom info
        fs_chr = {}
        fs_chr[arr[11]] = arr[12] # gene => chr
        fs_chr[arr[17]] = arr[18]


        # fusion gene pos info
        fs_info = {}
        fs_info[arr[11]] = int(arr[13]) # gene -> pos
        fs_info[arr[17]] = int(arr[19])

        samfile = pysam.AlignmentFile(args.bam,'rb')
        
        bp1_left = int(arr[13]) - 10
        bp1_right = int(arr[13]) + 10


        # 统计hgene/tgene支持split read的SA MAPQ
        # sv断点在5bp以内的primary read (with SA tag)，算作一条有效fusion-support read
        # 统计有效read的MAPQ，>=60 MAPQ占比需>=90%

        hgene_mapQ = [] # 
        
        hgene_strand = gene_strand[hgene]
        print("[INFO]: check Hgene %s (%s) info..." % (hgene,hgene_strand))
        
        hgene_bp = fs_info[hgene]

        hgene_left = hgene_bp - 10
        hgene_right = hgene_bp + 10

        #print(hgene,hgene_left,hgene_right)

        hgene_strand = gene_strand[hgene]
        tgene_strand = gene_strand[tgene]

        partner_chr = fs_chr[tgene]

        for read in samfile.fetch(arr[12],hgene_left,hgene_right):
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
                sr_mapq = int(sa_tag.split(',')[4])
                if sr_mapq < 10:
                    continue

                # check SA's chr
                sa_chr = sa_tag.split(',')[0]
                if sa_chr != partner_chr:
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
                            #print(read.query_name)
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
                
                if cigar_tuple[0][0] == 4:
                    # 50S100M
                    bnd_pos = left_pos
                elif cigar_tuple[-1][0] == 4:
                    # 100M50S
                    # cal ref consume
                    ref_consume_len = 0
                    for i in cigar_tuple:
                        if i[0] == 0 or i[0] == 2 or i[0] == 3 or i[0] == 7 or i[0] == 8:
                            ref_consume_len += i[1]

                    bnd_pos = left_pos + ref_consume_len - 1 + 1
                else:
                    # skip
                    continue

                # 比较sv给出的断点位置与split read计算的断点位置
                if abs(bnd_pos -hgene_bp) <= 10:
                    # 断点位置符合
                    bnd_pos_diff = 'PASS'
                else:
                    bnd_pos_diff = 'NOPASS'




                if sr_mapq >= args.mapq:
                    sr_mapq_check = 'PASS'
                else:
                    sr_mapq_check = 'NOPASS'

                if bnd_pos_diff == 'PASS' and sr_mapq_check == 'PASS':
                    final_check = 'PASS'
                else:
                    final_check = 'NOPAS'

                # fusion gene | gene | gene_type | read name | READ1/2 | chr | sv_pos | split_read_pos | pos_check | SR MAPQ | mapQ check | final result (check pos and mapq)
                ofFH.write('\t'.join([arr[1],hgene,'Hgene',read.query_name,r1r2,str(read.reference_name),str(hgene_bp),str(bnd_pos),bnd_pos_diff,str(read.mapping_quality),str(sr_mapq),sr_mapq_check,final_check])+'\n')

        


        tgene_strand = gene_strand[tgene]
        print("[INFO]: check Tgene %s (%s) info..." % (tgene,tgene_strand))
        
        
        tgene_bp = fs_info[tgene]
        tgene_left = tgene_bp - 10
        tgene_right = tgene_bp + 10
        #print(tgene,tgene_left,tgene_right)
        partner_chr = fs_chr[hgene]

        for read in samfile.fetch(arr[18],tgene_left,tgene_right):
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
                sr_mapq = int(sa_tag.split(',')[4])
                if sr_mapq < 10:
                    continue

                # check SA's chr
                sa_chr = sa_tag.split(',')[0]
                if sa_chr != partner_chr:
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
                #print(cigar_tuple)
                if cigar_tuple[0][0] == 4:
                    # 50S100M
                    bnd_pos = left_pos
                elif cigar_tuple[-1][0] == 4:
                    # 100M50S
                    # cal ref consume
                    ref_consume_len = 0
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
                    '''
                    for i in cigar_tuple:
                        if i[0] == 0 or i[0] == 2 or i[0] == 3 or i[0] == 7 or i[0] == 8:
                            ref_consume_len += i[1]

                    bnd_pos = left_pos + ref_consume_len - 1 + 1
                else:
                    # skip
                    continue


                # 比较sv给出的断点位置与split read计算的断点位置
                if abs(bnd_pos -tgene_bp) <= 10:
                    # 断点位置符合
                    bnd_pos_diff = 'PASS'
                else:
                    bnd_pos_diff = 'NOPASS'



                if sr_mapq >= args.mapq:
                    sr_mapq_check = 'PASS'
                else:
                    sr_mapq_check = 'NOPASS'

                if bnd_pos_diff == 'PASS' and sr_mapq_check == 'PASS':
                    final_check = 'PASS'
                else:
                    final_check = 'NOPAS'

                # fusion gene | gene | gene_type | read name | READ1/2 | chr | sv_pos | split_read_pos | pos_check | SR MAPQ | mapQ check | final result (check pos and mapq)
                ofFH.write('\t'.join([arr[1],tgene,'Tgene',read.query_name,r1r2,str(read.reference_name),str(tgene_bp),str(bnd_pos),bnd_pos_diff,str(read.mapping_quality),str(sr_mapq),sr_mapq_check,final_check])+'\n')




        samfile.close()
    fsFH.close()


if __name__ == "__main__":
    main()



