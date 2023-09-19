#!/Users/es47540/opt/anaconda3/bin/python3

#import numpy as np
import sys, getopt
import pandas as pd
import jenkspy

def main(argv):
    infile = ''
    strain1 = ''
    strain2 = ''
    ref = ''
    k = ''
    seed = 1
    arg_list = sys.argv[1:]
    short_options = 'i:a:b:r:k:s:h:'
    long_options = ['infile=', 'strain1=','strain2=','ref=','clusters=','seed=','help']
    usage = 'getGenDist.py\n \
    Usage:\n -h, --help \n \
    -i, --infile, genetic distance matrix by Phylip [req] \n \
    -a, --strain1, name of first strain of the pair \n \
    -b, --strain2, name of second strain of the pair \n \
    -r, --ref, name of the reference strain \n \
    -k, --clusters, number of bins to group strains by \n \
    -s, --seed, random seed number'

    try:
        options, args = getopt.getopt(arg_list,short_options,long_options)
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)

    if len(arg_list)<1 or not options:
        print(usage)
        sys.exit(2)

    for arg, val in options:
        if arg in ('-i', '--infile'):
            print('Using genetic distance matrix in', val)
            infile = val
        elif arg in ('-a', '--strain1'):
            print('Using strain1:', val)
            strain1 = val
        elif arg in ('-b', '--strain2'):
            print('Using strain2:', val)
            strain2 = val
        elif arg in ('-r', '--ref'):
            print('Using ref:', val)
            ref = val
        elif arg in ('-k', '--clusters'):
            print('Using clusters:', val)
            k = int(val)
        elif arg in ('-s', '--seed'):
            print('Using seed:', val)
            seed = int(val)
        elif arg in ('-h', '--help'):
            print('This program returns the genetic distance between 2 given strains.\n')
            print(usage)
            sys.exit()
        else:
            print(usage)
            sys.exit()

    Nstrains, distmat = openFile(infile)
    strainlist, newdistmat = makeStrainMatrix(distmat, Nstrains)
    print(strainlist)
    index1, index2 = getStrainIndex(strainlist, strain1, strain2)
# minlist, meanlist, maxlist, avgD = getGenDistSummary(newdistmat)
    D = getPairwiseD(newdistmat, index1, index2)
    print('The genetic distance between', strain1, 'and', strain2, 'is',D)
# get reference vector (second part)
    refVec = dist2one(ref, strainlist, newdistmat)
    refDict = orderDistVec(refVec, strainlist)
    sampledStrains = jenksNatBreaks(refDict, k, seed)
    print(sampledStrains)

def openFile(infile):
    Nstrains = int(open(infile).read().splitlines()[0])
    distmat = open(infile).read().splitlines()[1:] # add [1:] if first line has the number of strains
    return(Nstrains, distmat)

def makeStrainMatrix(distmat, Nstrains):
    # create new distance matrix with a a list of lists of strains
    Nlines = round(Nstrains/8 + 0.5)
    temp=str()
    newdistmat = []
    for i in enumerate(distmat):
        if i[0] % Nlines != 0 and i[0] !=0:
            temp = temp + i[1]
        else:
            temp = temp.split()
            newdistmat.append(temp)
            temp = str()
            temp = temp + i[1]
    newdistmat.append(temp.split())
    newdistmat = newdistmat[1:]
    # create list of strains with same indices as distance matrix
    strainlist = [i[0] for i in newdistmat]
    # make numbers float for each strain
    for list in newdistmat:
        for item in list:
            curr_index = list.index(item)
            if curr_index != 0:
                list[curr_index] = float(item)
    return(strainlist, newdistmat)

# get index for user-defined pair of strain
def getStrainIndex(strainlist, strain1, strain2):
    index1 = strainlist.index(strain1)
    index2 = strainlist.index(strain2)
#    for list in newdistmat:
#        if(strain1 in list):
#            index1 = newdistmat.index(list)
#        if(strain2 in list):
#            index2 = newdistmat.index(list)
    return(index1, index2)

# get pairwise genetic distance between strains 1 and 2
def getPairwiseD(newdistmat, index1, index2):
    print(index1, index2+1)
    print(newdistmat[index2][index1+1])
    return(newdistmat[index1][index2 + 1])

# get genetic distance summary (minimum, average, max)
def getGenDistSummary(newdistmat):
    maxlist = [max(list[1:]) for list in newdistmat]
    #meanlist = [np.mean(list[1:]) for list in newdistmat]
    minlist = [min(item for item in list[1:] if item != 0) for list in newdistmat]
    #avgD = np.mean([item for item in sum(newdistmat,[]) if (item !=0 and not isinstance(item, str))])
    return(minlist, maxlist)

# get distances from all strains to one user-defined strain
def dist2one(ref, strainlist, newdistmat):
    refIndex = strainlist.index(ref)
    refVec = newdistmat[refIndex]
    return(refVec)

def orderDistVec(refVec, strainlist):
    sortedRefVec = sorted(refVec[1:])
    sortedIndexList = []
    for item in sortedRefVec:
        if refVec.index(item) in sortedIndexList:
            sortedIndexList.append(refVec[1:].index(item, refVec[1:].index(item) + 1)+1)
        else:
            sortedIndexList.append(refVec[1:].index(item)+1)
    sortedIndexList = [item-1 for item in sortedIndexList]
    sortedStrainList = []
    sortedStrainList[:] = [strainlist[i] for i in sortedIndexList]
    refDict = {'strain': sortedStrainList, 'dist': sortedRefVec}
    return(refDict)

### Use 1d k-means clustering (Jenks Fisher Optimization) to group strains into k bins
def jenksNatBreaks(refDict, k, seed):
    df = pd.DataFrame(refDict)
    labels = range(k)
    #df['quantile'] = pd.qcut(df['dist'], q = 26, labels = k)
    #df['cut_bins'] = pd.cut(df['dist'], bins = 26, labels = k)
    breaks = jenkspy.jenks_breaks(df['dist'], nb_class=k+1)
    df['cuts_jenks'] = pd.cut(df['dist'], bins = breaks, labels = labels, include_lowest = True, duplicates = 'drop')
    #df.groupby(['cuts_jenks']).size().reset_index(name='counts')
    grouped = df.groupby('cuts_jenks')
    sampledStrains = grouped.apply(lambda x: x.sample(n=1, replace=True, random_state=seed))['strain'].tolist()
    return(sampledStrains)

if __name__ == "__main__":
    main(sys.argv[1:])
