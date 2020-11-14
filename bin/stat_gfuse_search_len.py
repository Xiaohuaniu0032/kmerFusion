import os
import sys

fusion_csv = sys.argv[1]

fcsv = open(fusion_csv,'r')
total_len = 0
for i in fcsv:
	if '>' in i:
		reg = i.split(':')[1]
		gene = i.split(',')[0][1:]	
		start = int(reg.split('-')[0])
		end = int(reg.split('-')[1])
		reg_len = end - start + 1
		total_len += reg_len
		val = "gene %s len is:%s" % (gene,reg_len)
		print(val)

val = "total search len is: %s(bp)" % (total_len)
print(val)
