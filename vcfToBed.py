import sys
import vcf
import time
import logging
import argparse
import os
import re
import subprocess
import pandas as pd
import numpy as np


class MakeFastaFile:

    def __init__(self,inVCFFile,outFastaFile,chr,minInsertLength,logger):
        self.inVCFFile = inVCFFile
        self.outFile = outFastaFile
        self.chr = chr
        self.minInsertLength = minInsertLength
        self.logger = logger
        self.convertToFasta()

    def convertToFasta(self):
        start = time.time()
        if self.inVCFFile.endswith('.gz'):

            progressCounter = 0
            # vcf_reader = vcf.Reader(open(vcf_file, 'r'))
            vcf_reader = vcf.Reader(filename=self.inVCFFile)
            record = next(vcf_reader)
            all_coord = {}
            for record in vcf_reader.fetch(self.chr):

                # for record in vcf_reader:
                if (record.var_type != "snp") and (record.var_subtype == "ins"):
                    alt_seq = record.ALT
                    ref_seq = record.REF
                    chrom = record.CHROM
                    pos = record.POS
                    ref_seq = str(ref_seq).replace('[', '').replace(']', '')
                    alt_seq = str(alt_seq).replace('[', '').replace(']', '')
                    alt_seq_length = len(alt_seq) - len(ref_seq)
                    if alt_seq_length > self.minInsertLength:
                        coord = (">chr" + str(chrom) + ":" + str(pos))
                        if coord in all_coord:
                            all_coord[coord] += 1
                        else:
                            all_coord[coord] = 1
                        text = (coord + "_" + str(all_coord[coord]) + "\n" + alt_seq)
                        # print(text)
                        # add status print statements every 1e7 records
                        if progressCounter == 10000:
                            print ("10000 records reached")
                            progressCounter = 0
                        progressCounter += 1
                        with open(self.outFile, 'a') as fasta:
                            fasta.write(text + '\n')

            # print("done")
            end = time.time()
            program_time = (end - start)
            self.logger.info("The VCF file for chr" + self.chr + " has been converted to a fastA file")
            print (program_time)


class MakeBlastFile:

    def __init__(self,inFastFile,outBlastFile,blastDir,outMeiFile,chr,meiOnly,logger):
        self.inFastFile = inFastFile
        self.outBlastFile = outBlastFile
        self.outMeiFile = outMeiFile
        self.blastDir = blastDir
        self.chr = chr
        self.blastDb = self.blastDir + "/chr" + str(chr) + ".fa"
        self.blastMEI = self.blastDir + "/MEI_refs.fa"
        self.meiOnly = meiOnly
        self.logger = logger
        self.createBlastFile()

    def createBlastFile(self):
        # TODO: make blast dictionary configurable
        if self.meiOnly is False:    
            print("Using genome BLAST database: " + self.blastDb)

        print("Using MEI BLAST database: " + self.blastMEI)

        if self.meiOnly is False:
            blastCommand = 'blastn -task megablast -query ' + self.inFastFile \
                            + ' -db ' + self.blastDb + ' -best_hit_score_edge 0.1 -gapopen 5 -gapextend 2 -reward 1 -penalty -1 -word_size 11 -perc_identity 90 -out ' + self.outBlastFile \
                            + ' -outfmt "6 qseqid sseqid qcovs stitle qlen slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos sstrand" '
            blast = os.system(blastCommand)

        meiBlastCommand = 'blastn -task megablast -query ' + self.inFastFile \
                        + ' -db ' + self.blastMEI + ' -best_hit_score_edge 0.1 -gapopen 5 -gapextend 2 -reward 1 -penalty -1 -word_size 11 -perc_identity 90 -out ' + self.outMeiFile \
                        + ' -outfmt "6 qseqid sseqid qcovs stitle qlen slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos sstrand" '

        # print(blastCommand)
        meiBlast = os.system(meiBlastCommand)
        self.logger.info("The fastA file for chr" + self.chr + " has been BLASTed")
        # print (program_time)
# 

class ApplyFilters:

    def __init__(self,chr,inBlastFile,inMeiFile,outBEDFile,meiOnly,logger):

        self.inBlastFile = inBlastFile
        self.inMeiFile = inMeiFile
        self.outBEDFile = outBEDFile
        self.meiOnly = meiOnly
        self.logger = logger
        self.chr = 'chr' + str(chr)

        self.filter()
        # filter(self)
# 


    def filter(self):

        for i in range(2):
            if i==0:
                if self.meiOnly is True:
                    continue
                else:
                    df = pd.read_csv(self.inBlastFile, sep='\t', header=None)
                    mei = False
            elif i==1:
                df = pd.read_csv(self.inMeiFile, sep='\t', header=None)
                mei = True 
            else:
                break

            df.columns = ["qseqid", "sseqid", "qcovs", "stitle", "qlen" ,"slen",
             "qstart", "qend", "sstart", "send", "qseq", "sseq", "evalue",
              "bitscore", "score", "length", "pident", "nident", "mismatch",
              "positive", "gapopen", "gaps", "ppos", "sstrand"]

            df['original_pos'] = df.qseqid.str.extract(r'^chr.*:(\d*)').astype(int)
            df['bp_length'] = df.sstart.astype(int) - df.original_pos.astype(int)
            df['blast_length'] = (df.length.astype(int)) / (df.qlen.astype(int))

            ## GET GOOD ALIGNMENTS
            good_alignments = df[(df.bitscore > 50) & ( (df.length / df.qlen) > 0.9) & (df.pident > 90) & (df.gaps < 3) ].copy()  

            if mei is False:
                self.deletions(good_alignments,self.chr)

                self.inversions(good_alignments,self.chr)

                self.tandumDup(good_alignments,self.chr)
            else:
                self.meis(good_alignments,self.chr)
            #print qseqid, "\t", bitscore, "\t", blast_length, "\t", pident, "\t", gaps, "\t", sstrand, "\t", send, "\t", sstart

        self.logger.info("The BLAST results have been filtered out output as BED")
            # print (program_time)


    def deletions(self,good_alignments,chr):
        #require that alignment is on plus strand nad downstream from org site
        dels =  good_alignments[( good_alignments['bp_length'] > 50 ) & ( good_alignments['bp_length'] < 1e5 ) & (good_alignments['sstrand']=='plus')]

        ## keep smallest call with best bitscore per start position
        dels = dels.sort_values(by=['bitscore', 'bp_length'], ascending=[False, True])
        dels = dels.drop_duplicates(subset='qseqid')
        ## CONVERT TO BED
        dels_bed = dels[['sseqid', 'original_pos', 'sstart', 'bitscore']].copy()
        dels_bed['SV_type'] = 'Deletion'
        lis = range(len(dels_bed))
        lis =  ["{:02d}".format(x) for x in lis]
        dels_bed['uniqueID'] = dels_bed.SV_type.str.cat(lis, sep='_')
        dels_bed.columns = ['chr', 'start', 'end', 'score', 'SV_type', 'uniqueID']
        dels_bed = dels_bed[['chr', 'start', 'end', 'uniqueID', 'score', 'SV_type']]
        # print(dels_bed)
        dels_bed.to_csv(path_or_buf=(chr + '_deletions.bed'), sep='\t', index=False)


    def meis(self,good_alignments,chr):
        #stop is org position of inserted sequence that is in query id
        #start is stop -1
        good_alignments['start'] = good_alignments.original_pos.astype(int) - 1
        good_alignments = good_alignments.sort_values(by='bitscore', ascending=False)
        meis = good_alignments.drop_duplicates(subset='qseqid').copy()
        meis['SV_type'] = 'MEI'
        lis = range(len(meis))
        lis =  ["{:02d}".format(x) for x in lis]
        meis['uniqueID'] = meis.sseqid.str.cat(meis.sstrand, sep="_").str.cat(lis, sep="_")
        meis_bed = meis[['qseqid', "start", "original_pos", "uniqueID", 'bitscore', 'SV_type']].copy()
        meis_bed.columns = ['chr', 'start', 'end', 'uniqueID', 'score', 'SV_type']
        meis_bed['chr'] = chr
        # print(meis_bed)
        meis_bed.to_csv(path_or_buf=(chr + '_MEIs.bed'), sep='\t', index=False)

    def inversions(self,good_alignments,chr):
        good_alignments['abs_length'] = abs(good_alignments.bp_length)
        invs =  good_alignments[( abs(good_alignments['abs_length']) > 50 ) & ( abs(good_alignments['abs_length']) < 1e5 ) & (good_alignments['sstrand']=='minus')]

        ## keep smallest call with best bitscore per start position
        invs = invs.sort_values(by=['bitscore', 'abs_length'], ascending=[False, True])
        invs = invs.drop_duplicates(subset='qseqid')
        ## CONVERT TO BED
        invs['SV_type'] = 'Inversion'
        lis = range(len(invs))
        lis =  ["{:02d}".format(x) for x in lis]
        invs['uniqueID'] = invs.SV_type.str.cat(lis, sep='_')
        invs['start'] = invs.loc[:, ['original_pos', 'sstart']].min(axis=1) 
        invs['end']  = invs.loc[:, ['original_pos', 'sstart']].max(axis=1)
        invs_bed = invs[['sseqid', 'start', 'end', 'uniqueID', 'bitscore', 'SV_type']].copy()
        invs_bed.columns = ['chr', 'start', 'end', 'uniqueID', 'bitscore', 'SV_type']
        # print(invs_bed)
        invs_bed.to_csv(path_or_buf=(chr + '_inversions.bed'), sep='\t', index=False)


    def tandumDup(self,good_alignments,chr):
        good_alignments['abs_length'] = abs(good_alignments.bp_length)
        tdups =  good_alignments[( abs(good_alignments['abs_length']) > 50 ) & ( abs(good_alignments['abs_length']) < 1e5 ) & (good_alignments['sstrand']=='plus') & (good_alignments.sstart.astype(int) < good_alignments.original_pos.astype(int))]

        ## keep smallest call with best bitscore per start position
        tdups = tdups.sort_values(by=['bitscore', 'abs_length'], ascending=[False, True])
        tdups = tdups.drop_duplicates(subset='qseqid')
        ## CONVERT TO BED
        tdups['SV_type'] = 'TandemDup'
        lis = range(len(tdups))
        lis =  ["{:02d}".format(x) for x in lis]
        tdups['uniqueID'] = tdups.SV_type.str.cat(lis, sep='_')
        tdups_bed = tdups[['sseqid', 'sstart', 'original_pos', 'uniqueID', 'bitscore', 'SV_type']].copy()
        tdups_bed.columns = ['chr', 'start', 'end', 'uniqueID', 'bitscore', 'SV_type']
        # print(tdups_bed)
        tdups_bed.to_csv(path_or_buf=(chr + '_tandemDups.bed'), sep='\t', index=False)



def main():

    #inputfile = ''
    #outputfile = ''
    #minInsertLength = 25
    # chromosome='21'


    parser = argparse.ArgumentParser(prog='vcfToFasta.py',
                                     description='''Convert a vcf file to a fasta file''')

    parser.add_argument('-inFile', '--inputFile', type=str, required=False,
                        help='''VCF file to be parsed''')

    parser.add_argument('-outFile', '--outputFile', type=str, required=False,
                        help='''Name of fasta file that gets outputted''')

    parser.add_argument('-chr', '--chromosome', type=str, required=False,
                        help='''Chromosome''')

    parser.add_argument('-dir', '--blastDir', type=str, required=False,
                        help='''BLAST database directory''')

    parser.add_argument('-meiOnly', '--meiOnly', type=bool, required=False,
                        help='''Analyze MEIs only''')

    parser.add_argument('-v', '--version', action='version', version='0.0.1 ')

    args = parser.parse_args()

    meiOnly = False
    inVCFFile = args.inputFile
    outFastaFile = args.outputFile
    chr = args.chromosome
    minInsertLength = 25
    blastDir = args.blastDir
    outBlastFile = chr + 'Blast.txt'
    outMeiFile = chr + 'BlastMEI.txt'
    outBEDFile = chr + 'BedFile.bed'
    meiOnly = args.meiOnly


    # create logger
    logger = logging.getLogger('vcf_to_fasta')
    logger.setLevel(logging.DEBUG)

    # create console handler
    vcf2Fasta = logging.StreamHandler()
    vcf2Fasta.setLevel(logging.DEBUG)

    # create formatter and add it to the handler
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    vcf2Fasta.setFormatter(formatter)

    # add the handler to the logger
    logger.addHandler(vcf2Fasta)
    # print ("here")
    try:
        outputfile = chr + "FastA_AS.fa"
        # simple test to try and open file
        f = open(inVCFFile, 'r')
        print ('Input file is ', inVCFFile)
        print ('Output file is ', outFastaFile)
        print ('Minimum insert length is ', minInsertLength)
        print ('Chromosome is ', str(chr))

    except IOError as e:
        logger.critical('Cannot open file.\n%s' % e)
        exit(1)

    #convert VCF to Fasta file
    MakeFastaFile(inVCFFile,outFastaFile,chr,minInsertLength,logger)

    # take newly created fasta file and create a blast file and a mei blast file
    inFastFile = outFastaFile
    MakeBlastFile(inFastFile,outBlastFile,blastDir,outMeiFile,chr,meiOnly,logger)

    #take newly created blast file and mei blast file and apply filters to create bed file
    inBlastFile = outBlastFile
    inMeiFile = outMeiFile

    # blast_file = open('/Volumes/bioinfo/users/Alexandrea/hack-a-thon/blastchr21_output6.txt', 'r')
    ApplyFilters(chr,inBlastFile,inMeiFile,outBEDFile,meiOnly,logger)
    # ApplyFilters(blast_file) #inBlastFile, inMeiFile, outBEDFile)

if __name__ == "__main__":
    main()
