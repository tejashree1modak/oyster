#This script was written to compare degs from different comparisons to get the common degs 
#For examples in design5 I compared each control to the other and each treatment to the other
#If I wanted to know which genes are common degs from comparison of con 12d vs 5d and con 12d vs 16d
#Here is how I would run the script: python compare.py Gene_dfSig_des5_con_trans_12vs5_0.1.csv Gene_dfSig_des5_con_trans_12vs16_0.1.csv

import argparse
import csv

#This function enumerates rows and we select r[0] meaning first element which is the gene name.
#Using the {} makes a set that eliminates duplicates.
def get_genes(infile):
    return { r[0] for i, r in enumerate(csv.reader(infile)) if i > 1 }

p = argparse.ArgumentParser()
p.add_argument('files', nargs='+', type=argparse.FileType('r'), help='the input files')
#nargs='+' meaning can give more than 1 input files
args = p.parse_args()

genes= None
for f in args.files:
    if genes is None:
        genes = get_genes(f) #Making a set of all genes in file1 in first loop
    else:
        genes = genes.intersection(get_genes(f)) #second time loop runs genes obj is already full so get intersection 

if genes:
    print "\n".join(genes)
    #output on stdout one gene per line,
