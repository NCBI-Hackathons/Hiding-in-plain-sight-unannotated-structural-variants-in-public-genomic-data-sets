#!/usr/bin/env python

## import packages
import argparse
from pybedtools import BedTool

def get_args():
	parser = argparse.ArgumentParser(description='Overlap merged bed',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input', 
		help='the input bed files for overlap', 
		required = False, type = str,
		default = '../tests/chr1_MEIs_merged.bed')
	parser.add_argument('-r','--ref', 
		help='the bed file for overlap reference', 
		required = False, type = str,
		default = '../resources/1kg_MEI.bed')
	parser.add_argument('-p','--prop',
		help='the proportion for overlap',
		required = False,default=0.9)
	parser.add_argument('-o','--output',
		help='the ouput bed files for overlap')
	parser.add_argument('-s','--overlapstyle', 
		help='the overlapstyle', 
		required = False, choices = ['r','n'], default='r')
	parser.add_argument("-m", "--MEI", 
		help="turn on if input is MEI",action="store_true")
	parser.add_argument('-b','--basepair', 
		help='the length of segment to expand', 
		required = False, type = int,
		default = 10)
	
	args = parser.parse_args()
	return(args)
'''
def get_overlap(dat,prop_thresh):
	overlaplist = []
	for index, row in dat.iterrows():
		prop = (float(row[2])-float(row[1]))/float(row[len(row)-1])
		if prop>=prop_thresh:
			overlaplist.append(1)
		else:
			overlaplist.append(0)
	return(overlaplist)
'''

def get_overlap(dat,prop_thresh):
	overlaplist = []
	for index, row in dat.iterrows():
		if int(row[len(row)-1])>0:
			overlaplist.append(1)
		else:
			overlaplist.append(0)
	return(overlaplist)

def get_header(filename):
	f = open(filename)
	headerline = f.readline()
	ret = headerline.rstrip().split('\t')
	return(ret)

def changebase(line,startbase=10,endbase=10):
	line.start = line.start-startbase
	line.end = line.end+endbase
	return(line)

def overlap_bed(args):
	dat = BedTool(args.input)
	headerdat = get_header(args.input)
	ref = BedTool(args.ref)
	#headerref = get_header(args.ref)
	#headerref = [s + '_ref' for s in headerref]

	if args.MEI:
		ref = ref.each(changebase,startbase=args.basepair,endbase=args.basepair)
		intersectdat = dat.intersect(ref, c=True)
	else:
		if args.overlapstyle=='r':
			intersectdat = dat.intersect(ref, c=True,f=args.prop,r=True)
		elif args.overlapstyle=='n':
			intersectdat = dat.intersect(ref, c=True)

	datf = intersectdat.to_dataframe(header=None)
	overlap = get_overlap(datf,prop_thresh=args.prop)
	datf['overlap'] = overlap
	#datf.columns=headerdat+headerref+['length','overlap']
	datf.columns=headerdat+['count','overlap']
	datf.to_csv(args.output,header=True,index=False,sep='\t')

	return

def main():
	args = get_args()
	overlap_bed(args)

if __name__=='__main__':
	main()







