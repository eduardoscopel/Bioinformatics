#!/usr/bin/python

import sys
infile = open(sys.argv[1])
outfile = []

for line in infile:
    if line.startswith(">"):
        if(outfile != []): outfile.close()
        strain_name = line.split("chr",1)[0]
        filename = line[1:-1]+".fasta"
        outfile = open(filename, 'w')
        outfile.write(strain_name + '\n')
    else:
        outfile.write(line)


outfile.close()
