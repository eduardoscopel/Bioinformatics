# Using Python 3.4.3 on Sapelo1 import pandas package to deal with large tables

import timeit

start = timeit.default_timer()

import os
import glob
import pandas as pd
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import re
import csv

def makehist():
# Open read depth tables
    files = glob.glob('*.coverage') # create a list of filemanes in wd, ending in .coverage
    for filename in files:
        depthtab = pd.read_csv(filename, sep ='\t', header = None, names = ['Chromosome', 'Position', 'Depth'])
        a = re.sub('.coverage', '', filename)
        b = re.sub('Depth_', '', a)
        med = 3*np.median(depthtab['Depth'])
        name = b + '_hist_py.pdf'
        with PdfPages(name) as pdf:
            plt.figure()
            plt.hist(depthtab['Depth'], 5000, range=[0,med], color = 'k')
            plt.suptitle('\n'.join([b, str(int(med/3)) + 'x']))
            plt.ylabel('Frequency')
            plt.xlabel('Read Depth')
            plt.subplots_adjust(left=0.2)
            plt.axvline(x=med/3, color = 'b')
            pdf.savefig()
            plt.close()

# subsetting the table by chromosome
        cats = df['Chromosome'].unique()
        for cat in cats:
            chrm = df[(df['Chromosome'] == cat)]
            chrm = chrm.iloc[:,1:3]
            chrm.head
# Create sliding window
            W = list(range(0, int(max(chrm["Position"])/1000)*1000,1000))
            W.append(max(W)+1000)


stop = timeit.default_timer()

print stop - start
