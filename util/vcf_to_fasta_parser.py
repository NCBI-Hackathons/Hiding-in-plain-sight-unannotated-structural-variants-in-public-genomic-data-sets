import sys, getopt
import vcf

#########
# This script has been modified with permission
# from an original script written by students at UMUC:
# Mike Fox, David Jalali, Eric Keller,
# Mary Schramko, Upasana Pandey, Albakry Koroma
#########

def main(argv):
    inputfile = ''
    outputfile = ''
    minInsertLength = 25
    chromosome='21'
    try:
        opts, args = getopt.getopt(argv,"hi:o:m:c:",["ifile=","ofile=", "minlen=", "chrom="])
    except getopt.GetoptError:
        print ('vcf_to_fasta_parser.py -i <inputfile> -o <outputfile> -m <minInsertLength> -c <chromosome>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('vcf_to_fasta_parser.py -i <inputfile> -o <outputfile> -m <minInsertLength> -c <chromosome>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-m", "--minlen"):
            minInsertLength = arg
        elif opt in ("-c", "--chrom"):
            chromosome = arg

    print ('Input file is ', inputfile)
    print ('Output file is ', outputfile)
    print ('Minimum insert length is ', minInsertLength)
    print ('Chromosome is ', str(chromosome))


    # TODO: add check for .gz extension and read in file appropriatey
    # TODO: check that file exists
    # TODO: add status print statements every 1e7 records
    # TODO: send error messages to log file
    # vcf_file = "/mnt/hnas/bioinfo/projects/VariantInfoDownloads/ExAC/gnomAD/Version2.1/genomes/gnomad.genomes.r2.1.sites.vcf.gz"
    vcf_file = inputfile
    fastA = outputfile
    minInsertLength = int(minInsertLength)
    chromosome = str(chromosome)


    # vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    vcf_reader = vcf.Reader(filename=vcf_file)
    record = next(vcf_reader)
    all_coord = {}
    for record in vcf_reader.fetch(chromosome):
    # for record in vcf_reader:
        if (record.var_type != "snp") and (record.var_subtype == "ins"):
            alt_seq = record.ALT
            ref_seq = record.REF
            chrom = record.CHROM
            pos = record.POS
            ref_seq = str(ref_seq).replace('[', '').replace(']', '')
            alt_seq = str(alt_seq).replace('[', '').replace(']', '')
            alt_seq_length = len(alt_seq) - len(ref_seq)
            if alt_seq_length > minInsertLength:
                coord = (">chr" + str(chrom) + ":" + str(pos))
                if coord in all_coord:
                    all_coord[coord]+=1
                else:
                    all_coord[coord]=1
                text = (coord + "_"  + str(all_coord[coord]) +  "\n" + alt_seq)
                # print(text)
                with open(fastA, 'a') as fasta:
                    fasta.write(text + '\n')

    print("done")

if __name__ == "__main__":
   main(sys.argv[1:])
