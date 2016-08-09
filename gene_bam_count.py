#! /usr/bin/env python2.7
'''
output single base depth, coverage, initial read count  per transcript 
'''
import argparse
import sys
import os
import re
import pysam
import pprint
sys.path.append("/annoroad/bioinfo/PMO/liutao/config/bin/common_python")
import gtf_reader

__author__='Liu Tao'
__mail__= 'taoliu@annoroad.com'

pat1=re.compile('^\s+$')
pat2=re.compile('^chr',re.I)
pat3=re.compile('scaffold')
samtools = '/annoroad/bioinfo/PMO/share/software/samtools-0.1.19/samtools'

def bam_reader(bam,quality=0):
	'''
	return a dict 
	dict -> chr -> position = depth
	'''
	r_dict={}
	print(bam)
	samfile = pysam.Samfile(bam,'rb')
	for chr in samfile.references:
		nchr = chr
		if not pat2.search(chr) : nchr = 'chr'+chr
		if pat3.search(chr) : continue
		if not nchr in r_dict:
			r_dict[nchr]={}
		for read in samfile.fetch(chr):
			if read.mapq < quality : continue
			strand = '+'
			if read.is_reverse: strand = '-'
			left = read.pos  ## 0 base
			pos = read.positions
			if not left in r_dict[nchr]:
				r_dict[nchr][left] = 0 
			r_dict[nchr][left] += 1
	samfile.close()
	return r_dict

class Counter:
	count = 0
	def __call__(self,alignment):
		self.count += 1

def head(iter):
	for a_alignment in iter.pileups:
		if a_alignment.is_head == 1 :
			return True
	return False

def tr(record,g_len):
	tt = ''
	for a_pos in range(g_len):
		if not a_pos in record:
			record[a_pos] = 0
		tt += '{0} '.format(record[a_pos])
	return tt.rstrip()

def output(bam, gene,  output):
	samfile = pysam.Samfile(bam,'rb')
	output.write('gene\tversion\tread_count\tgene_length\tcover_length\tcover_length/gene_length\tfirst_cover_length\tfirst_cover_length/gene_length\tper_base_depth\n')
	for a_gen in gene:
	#for a_gen in ['NM_014851']:
		for a_version in gene[a_gen]:
			a_transcript = gene[a_gen][a_version]
			gene_length , first_cover_length ,read_count,cover_length = 0 , 0 , 0 ,0
			chr = a_transcript.chr
			c = Counter()
			record = {}
			header = [] 
			for a_region in a_transcript.region:
				tt = list(a_region)
				start,end = tt[0],tt[-1]
				samfile.fetch(chr,start,end,callback=c)
				read_count += c.count
				
				for iter in samfile.pileup(chr,start,end):
					if iter.pos+1  < start  or iter.pos >= end : continue 
					pos = iter.pos - start + gene_length + 1
					#print(start,end,iter.pos,gene_length,pos)
					count = iter.n 
					record[pos] = count 
					if head(iter) == True:
						header.append(pos)
				gene_length += end - start + 1
			
				#print(record)

			cover_length = len(record.keys())
			first_cover_length = len(header)
			
			per_base_depth=''
			if not len(record.keys())== 0 :
				per_base_depth = tr(record,gene_length)
			else:
				per_base_depth = '0 ' * gene_length
			
			output.write('{0}\n'.format("\t".join([str(i) for i in [a_gen,a_version,read_count,gene_length,cover_length,cover_length/gene_length,first_cover_length,first_cover_length/gene_length,per_base_depth]])))
		#sys.exit()

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--input',help='input file',dest='input',required=True)
	#parser.add_argument('-q','--quality',help='read mapping quality ',dest='quality',type=int , default = 10)
	parser.add_argument('-g','--gtf',help='gtf file',dest='gtf',type=argparse.FileType('r'),required=True)
	parser.add_argument('-o','--output',help='output file',dest='output',type=argparse.FileType('w'),required=True)
	args=parser.parse_args()
	
	#site_dict = bam_reader(args.input,args.quality)

	#mapping_rate(args.input,args.output)
	exon_description = gtf_reader.gtf_reader(args.gtf,type='exon')
	output(args.input,exon_description,args.output)
	#tt = bam_reader(args.input)
	#pprint.pprint(tt)

if __name__ == '__main__':
	main()
