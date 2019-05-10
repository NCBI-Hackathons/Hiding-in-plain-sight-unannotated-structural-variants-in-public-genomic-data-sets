#!/usr/bin/env python

## import packages
import argparse
from pybedtools import BedTool

def get_args():
	parser = argparse.ArgumentParser(description='Merge the MEIs',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input', 
		help='the input bed files for MEIs', 
		required = False, type = str,
		default = '../tests/chr1_MEIs.bed')
	parser.add_argument('-o','--output', 
		help='the output bed files for MEIs', 
		required = False, type = str,
		default = '../tests/chr1_MEIs_merged.bed')
	parser.add_argument('-b','--basepair', 
		help='the length of segment to expand', 
		required = False, type = int,
		default = 10)
	
	args = parser.parse_args()
	return(args)

def changechr(line):
	IDlist = line[3].split('_')
	line[0] = '_'.join([IDlist[0],IDlist[1]])
	return(line)

def changebase(line,startbase=10,endbase=10):
	line.start = line.start-startbase
	line.end = line.end+endbase
	return(line)

'''
def reversechangebase(line,startbase=10,endbase=10):
	
	#to_dataframe won't work after it.
	
	line.start = line.start+startbase
	line.end = line.end-endbase
	return(line)
'''

def reversechangebase(datf,startbase=10,endbase=10):
	datf['start'] = datf['start']+startbase
	datf['end'] = datf['end']-endbase
	return(datf)

def bedtoframe(beddat,chrom):
	ret = beddat.to_dataframe()
	ret = ret.drop(columns=['start','end'])
	ret.columns=['#chr','start','end','merged_IDs','mean_bitscore','SV_type']
	ret['SV_type'] = 'MET'
	ret['uniqueID'] = [item.split('|')[0] for item in list(ret['merged_IDs'])]
	ret['start'] = ret['start'].astype(int)
	ret['end'] = ret['end'].astype(int)
	ret['#chr'] = str(chrom)

	seq_columns = ['#chr','start','end','uniqueID','SV_type','mean_bitscore','merged_IDs']
	
	ret = ret.reindex(columns = seq_columns)
	return(ret)


def merge_MEI(args):
	dat = BedTool(args.input)
	chrom = dat[0][0]
	dat = dat.each(changechr)
	dat = dat.sort()
	changebase_dat = dat.each(changebase,startbase=args.basepair,endbase=args.basepair)
	ret = changebase_dat.merge(c=(2,3,4,5,6), o=('mean','mean','collapse','mean','collapse'), delim='|')
	ret = bedtoframe(ret,chrom)
	ret = reversechangebase(ret,startbase=args.basepair,endbase=args.basepair)
	ret.to_csv(args.output,header=True,index=False,sep='\t')
	return()

def main():
	args = get_args()
	merge_MEI(args)

if __name__=='__main__':
	main()










