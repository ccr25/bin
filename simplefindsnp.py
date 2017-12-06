#! /usr/bin/env python3

__author__ = 'croe and jtravis'
import csv
import sys


# Usage: simplefindsnp input.vcf output.vcf min_DP min_ratio max_ratio

input = sys.argv[1]
output = sys.argv[2]
min_dp = sys.argv[3]
min_ratio = sys.argv[4]
max_ratio = sys.argv[5]

with open(input) as handle, open(output, 'w') as outfile:
# TODO: skip metadata

    for line in handle:
        if line.startswith('#CHROM'):
            header = line.split('\t')
            #break
            # Uncomment this line to include metadata in the outfile
            #outfile.write(line)
            writer = csv.DictWriter(outfile, fieldnames=header, delimiter='\t')
            
            sample_name = header[header.index('FORMAT') + 1]

            for row in csv.DictReader(handle, fieldnames=header, delimiter='\t'):
                row[sample_name] = dict(zip(row['FORMAT'].split(':'), row[sample_name].split(':')))


                if row['ALT'] != '.' and float(row[sample_name]['DP']) >= float(min_dp):
                    # "63,0" => ['63', '0']
                    ad = row[sample_name]['AD'].split(',')
                    # 63
                    allele1 = float(ad[0])
                    # 0
                    allelel2 = float(ad[1])

                    ratio = allele1 / float(row[sample_name]['DP'])

            # ratio should be a value between 0 and 1	
                    if ratio < float(min_ratio) or ratio > float(max_ratio):
                        continue

                    writer.writerow(row)
