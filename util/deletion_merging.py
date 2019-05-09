import pybedtools
import argparse
from statistics import mean

class SV:
    
    def __init__(self, chr, start, stop, id_num, score, sv_type):
        self.chr = chr
        self.start = int(start)
        self.stop = int(stop)
        self.id_num = id_num
        self.score = float(score)
        self.sv_type = sv_type
        self.inter_sv = []

    def equal(self, new_sv):
        ''' check if two SVs are the same based on start and stop coordinates '''
        return(self.start == new_sv.start and self.stop == new_sv.stop)

    def add_sv(self, new_sv):
        ''' add an overlapping SV '''
        self.inter_sv.append(new_sv)

    def merge(self):
        ''' find the mean start and stop coordinates from all overlapping SVs '''

        start_coords = []
        stop_coords = []
        align_scores = []
        ids = []
        ## collect attributes across overlapping SVs
        for sv in self.inter_sv:
            start_coords.append(sv.start)
            stop_coords.append(sv.stop)
            align_scores.append(sv.score)
            ids.append(sv.id_num)

        ## merge attributes across overlapping SVs by mean
        self.start = round(mean(start_coords))
        self.stop  = round(mean(stop_coords))
        self.score = round(mean(align_scores),1)
        ## combine sv IDs
        self.id_num = ','.join(ids)

    def print(self):
        ''' print SV object for testing purposes '''
        print(f'{self.chr}\t{self.start}\t{self.stop}\t{self.id_num}\t{self.score}\t{self.sv_type}')

    def to_string(self):
        ''' print SV object for testing purposes '''
        return(f'{self.chr}\t{self.start}\t{self.stop}\t{self.id_num}\t{self.score}\t{self.sv_type}\n')

parser = argparse.ArgumentParser()
parser.add_argument('bedfile', help='path to the bedfile to merge')
parser.add_argument('outfile', help='merged bedfile path and name')
parser.add_argument('-v', '--overlap', default=0.9, type=float, help='amount of overlap required to merge two events') 
args = parser.parse_args()

## Import bedfile to merge
a = pybedtools.BedTool(args.bedfile)

## Bedtools intersect of the same file is used here to take
## advatage of the reciprocal overlap option
ab_intersect = a.intersect(a, wao=True, f=args.overlap, r=True)

## List to store the merged SV objects
merged_svs = []

## prev_sv holds A file SV object that B file SVs are currently being added to
prev_sv = None

for item in ab_intersect:
    chr1,start1, end1, id1, score1, sv_type1 = item[0:6]
    chr2, start2, end2, id2, score2, sv_type2 = item[6:12]
    
    curr_sv = SV(chr1, start1, end1, id1, score1, sv_type1)
    next_sv = SV(chr2, start2, end2, id2, score2, sv_type2)

    ## First SV in the file
    if prev_sv == None:
        curr_sv.add_sv(next_sv)
        prev_sv = curr_sv
    ## A file SV has not changed
    elif curr_sv.equal(prev_sv):
        prev_sv.add_sv(next_sv)
    ## A file SV has changed
    else:
        ## merge and save the previous A file SV
        prev_sv.merge()
        merged_svs.append(prev_sv)
        ## Change prev_sv to new A file SV
        curr_sv.add_sv(next_sv)
        prev_sv = curr_sv

## Write merged SVs to the output file
with open(args.outfile, 'w') as f:
    for items in merged_svs:
        f.write(items.to_string())

