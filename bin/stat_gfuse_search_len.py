import os
import sys
from collections import defaultdict

fusion_csv = sys.argv[1]

fcsv = open(fusion_csv,'r')
gene2len = defaultdict(int)

total_len = 0
for i in fcsv:
	if '>' in i:
		reg = i.split(':')[1]
		gene = i.split(',')[0][1:]	
		start = int(reg.split('-')[0])
		end = int(reg.split('-')[1])
		reg_len = end - start + 1
		gene2len[gene] += reg_len
		#total_len += reg_len
		#val = "gene %s len is:%s" % (gene,reg_len)
		#print(val)

total_len = 0;
for gene,n in gene2len.items():
	#print(gene,n)
	val = "gene %s len is: %s" % (gene,n)
	print(val)
	total_len += n

val = "total search len is: %s(bp)" % (total_len)
print(val)
