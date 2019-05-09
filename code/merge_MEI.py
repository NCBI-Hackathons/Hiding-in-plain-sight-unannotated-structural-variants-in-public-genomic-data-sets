#!/usr/bin/env python

## import packages
import argparse
from pybedtools import BedTool

def get_args():
	parser = argparse.ArgumentParser(description='Merge the MEIs',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input', 
		help='the input bed files for MEIs', 
		required = False, type = str,
		default = '../test/chr1_MEIs.bed')
	parser.add_argument('-o','--output', 
		help='the output bed files for MEIs', 
		required = False, type = str,
		default = '../out/chr1_MEIs_merged.bed')
	
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

def bedtoframe(beddat):
	ret = beddat.to_dataframe()
	ret = ret.drop(columns=['start','end'])
	ret.columns=['#chr','start','end','uniqueID','score','SV_type']
	ret['start'] = ret['start'].astype(int)
	ret['end'] = ret['end'].astype(int)
	return(ret)


def merge_MEI(args):
	dat = BedTool(args.input)
	dat = dat.each(changechr)
	dat = dat.sort()
	changebase_dat = dat.each(changebase)
	ret = changebase_dat.merge(c=(2,3,4,5,6), o=('mean','mean','collapse','mean','collapse'), delim='|')
	ret = bedtoframe(ret)
	ret.to_csv(args.output,header=True,index=False,sep='\t')
	return()

def main():
	args = get_args()
	merge_MEI(args)

if __name__=='__main__':
	main()









