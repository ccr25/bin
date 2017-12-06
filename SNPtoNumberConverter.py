#/usr/bin/env python3
#This scritpt changes all bases identical to the reference base to '0' and all bases different from the reference base to '1'. 

import sys
import os
import csv

__version__ = "0.0.1"
__author__ = "Chandler Roe <chandler.c.roe@nau.edu"

#Usage: SNPtoNumberConverter.py bestsnp.tsv outfileName
SNP_file = sys.argv[1]
output = sys.argv[2]
with open(SNP_file) as handle, open(output, 'w') as outfile:
    for line in handle:
        if line.startswith('LocusID'):
            header = line.split('\t')
            locus = header[header.index('LocusID')]
            outfile.write(line)
            writer = csv.DictWriter(outfile, fieldnames=header, delimiter='\t')
            reference_state = header[header.index('LocusID') + 1]
            #print (reference_state)
            next(handle)
            for row in csv.DictReader(handle, fieldnames=header, delimiter='\t'):
                outfile.write((row[locus] + '\n'))
            
