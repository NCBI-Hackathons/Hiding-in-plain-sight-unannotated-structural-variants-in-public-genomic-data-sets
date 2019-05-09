#!/usr/bin/env python

## import packages
import argparse
from pybedtools import BedTool

def get_args():
	parser = argparse.ArgumentParser(description='Overlap merged bed',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input', 
		help='the input bed files for overlap', 
		required = False, type = str,
		default = '../tests/chr1_deletions_merged.bed')
	parser.add_argument('-r','--ref', 
		help='the bed file for overlap reference', 
		required = False, type = str,
		default = '../resources/RefSeq_hg19_introns.bed')
	parser.add_argument('-p','--prop',
		help='the proportion for overlap',
		required = False,default=0.9)
	parser.add_argument('-o','--output',
		help='the ouput bed files for overlap', 
		required = False, type=str,
		default = '../tests/chr1_detect_psuedogene.bed')
	
	
	args = parser.parse_args()
	return(args)

def get_header(filename):
	f = open(filename)
	headerline = f.readline()
	ret = headerline.rstrip().split('\t')
	return(ret)


def overlap_intron(args):
	dat = BedTool(args.input)
	headerdat = get_header(args.input)
	ref = BedTool(args.ref)
	headerref = get_header(args.ref)
	headerref = [s + '_introns' for s in headerref]

	intersectdat = dat.intersect(ref, wo = True, f=args.prop)
	datf = intersectdat.to_dataframe(header=None)
	datf.columns = headerdat + headerref + ['overlap']

	seq_columns = headerdat + [headerref[3]]
	datf = datf.reindex(columns = seq_columns)

	datf.to_csv(args.output,header=True,index=False,sep='\t')
	return


def main():
	args = get_args()
	overlap_intron(args)

if __name__=='__main__':
	main()



