'''
This function takes a fastq file and strips the fastq formatting plus strips
the i7 index and anything after it, outputting the final strands to a .txt file
'''

i7_dict = {
    "D701": "ATTACTCG",
    "D702": "TCCGGAGA",
    "D703": "CGCTCATT",
    "D704": "GAGATTCC",
    "D705": "ATTCAGAA",
    "D706": "GAATTCGT",
    "D707": "CTGAAGCT",
    "D708": "TAATGCGC",
    "D709": "CGGCTATG",
    "D710": "TCCGCGAA",
    "D711": "TCTCGCGC",
    "D712": "AGCGATAG",
    }


def fastq_stripping(fileinput, strand_array):
    with open(fileinput, 'r') as r:
        strand = 0
        for line in r:
            if line[0] == '@':
                # Every strand has an @ starting in the line before
                strand = 1
            elif strand == 1:
                strand_array.append(line)
                strand = 0


def index_alignment(fileinput, strand_array, index, outputdirectory,skip_alignment):
    index_formula = ("AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" + i7_dict[index] +
                     "ATCTCGTATGCCGTCTTCTGCTTG")  # Adds formula & converts
    if outputdirectory[-1] is not "/":
        _outputdirectory=outputdirectory+"/"
    else:
        _outputdirectory=outputdirectory
    bad=0
    with open(_outputdirectory + fileinput +"_Stripped.txt", 'w') as w:  # Output filename formatting
        for line in strand_array:
            if skip_alignment==True:
                w.write(line.rstrip()+'\n') #if we don't need to align to an index strand just strip
                continue
            for a in pairwise2.align.localms(line, index_formula,
                                             1, -10, -10, -10):
                align1 = a[0]  # Necessary for index finding equation
                begin = a[3]  # Necessary for index finding equation
                index_start = len(align1[:begin]) - align1[:begin].count('-')
                # This was taken from the source code of format_alignment
                if len(line[0:index_start])>0: w.write(line[0:index_start] + '\n')
                if len(line[0:index_start])==0: bad+=1
                # Write the final strand from 0 until the start of the index
                break
        #print("Done!")
        #print "% bad: {}".format(float(bad)/float(len(strand_array)))


if __name__ == "__main__":
    import argparse
    import json
    import os
    from Bio import pairwise2
    from Bio.pairwise2 import format_alignment

    parser = argparse.ArgumentParser(description="Parse file paths")

    parser.add_argument("--input", help="Input File") #path to raw fastq file
    parser.add_argument("--config", help="Config File")
    parser.add_argument("--output", help="File Output Directory")
    parser.add_argument("--skip_alignment",action="store_true")

   
    args = parser.parse_args()

    strand_array = []

    with open(args.config, 'r') as config:
        for line in config:
            config = line.split()  # config[0] = filename, config[1] = index
            if config[0] in args.input: break # look for the file of interest
        fastq_stripping(args.input, strand_array)
        index_alignment(config[0], strand_array, config[1], args.output,args.skip_alignment)
