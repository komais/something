#! /usr/bin/env python3
import argparse
import sys
import os
import re
import random
import bisect

__author__='Liu Tao'
__mail__= 'taoliu@annoroad.com'

pat1=re.compile('^\s+$')

def find_overlap(start,end,point):
	region = list(zip(start,end))
	index_start = bisect.bisect_left(sorted(start),point[1])  ### find the position of end 
	index_end   = bisect.bisect_left(sorted(end),point[0])    ## find the position of start 
	record = []
	if index_start < len(region) - index_end:   ## find all correct position, and choose a smaller way
		for i in range(index_start):
			sorted_region = sorted(region)
			end_point = sorted_region[i][1]
			if end_point > point[0] :
				record.append(sorted_region[i])
	else :
		for i in range(index_end,len(end)):
			sorted_region = sorted(region,key=lambda r :r[1])
			end_start = sorted_region[i][0]
			if end_start < point[1]:
				record.append(sorted_region[i])
	return(record)

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--input',help='input file',dest='input',type=open,required=True)
	parser.add_argument('-o','--output',help='output file',dest='output',type=argparse.FileType('w'),required=True)
	args=parser.parse_args()

	start = random.sample(list(range(10000000)),100000)
	end = random.sample(list(range(10000000)),100000)
	point = [10,100]
	print(find_overlap(start,end,point))

if __name__ == '__main__':
	main()
