# Using Python 3.4.3 on Sapelo1 import pandas package to deal with large tables

import pandas as pd
import csv
import re
import sys

infile = sys.argv[1]
outfile = sys.argv[2]

def ploidytab(infile):
    chr_cov = [[],[],[]]
    df = pd.read_csv(infile, sep ='\t', header = None, names = ['Chromosome', 'Position', 'Depth'])
    strain = re.sub('.depth', '', str(infile))
    cats = df['Chromosome'].unique()
    chr_cov[0]=strain
    for cat in cats:
        if re.search("chrM", cat):
            continue
        chr_cov[1].append(cat)
        chr_cov[2].append(df[(df['Chromosome'] == cat)].loc[:,'Depth'].median())


    with open(outfile, 'w') as f:
        w = csv.writer(f,delimiter="\t")
        f.write('%s\n' % (str(chr_cov[0])))
        w.writerows(chr_cov[1:3])

#if __name__ == '__main__':
ploidytab()
