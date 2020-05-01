#! /bin/python
from sys import argv
#filename = "illumina.somatic.adjusted.vcf"
filename = argv[1]
with open(filename, 'r') as infile:
    l = []
    for line in infile:
        line = line.rstrip()
        if line.startswith("#"):
            print line
            continue
        ident = line.split('\t')[2][:-1]
        if ident in l:
            continue
        else:
            print line
            l.append(ident)
        
