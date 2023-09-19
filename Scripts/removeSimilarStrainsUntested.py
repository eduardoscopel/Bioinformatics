# awk 'NR%7 {printf("%s", $0); next}
# {print $0} ' WG.cont_hh0.dnadist > Sc.dnadist
import numpy as np

def openFile(infile):
    with open(infile, 'r') as f:
        distmat = f.readlines()[1:]
    return(distmat)

def makeStrainMatrix(distmat):
    # create list of strains with same indices as distance matrix
    strainlist = [distmat[i].split()[0] for i in range(len(distmat))]
    # replace tab after strain name with space for every line
    newdistmat = [sub.replace('       ',' ').split() for sub in distmat]
    # make numbers float for each strain
    for list in newdistmat:
        for item in list:
            curr_index = list.index(item)
            if curr_index != 0:
                list[curr_index] = float(item)
    return(strainlist, newdistmat)

strain1 = 'cont'
strain2 = 'hh0'
strain1 = '08-19-02-3'
strain2 = 'AF65'
# get index for user-defined pair of strain
def getStrainIndex(newdistmat, strain1, strain2):
    for list in newdistmat:
        if(strain1 in list):
            index1 = newdistmat.index(list)
        if(strain2 in list):
            index2 = newdistmat.index(list)
    return(newdistmat)

def getGenDistSummary(newdistmat):
maxlist = [max(list[1:]) for list in newdistmat]
meanlist = [np.mean(list[1:]) for list in newdistmat]
minlist = [min(item for item in list[1:] if item != 0) for list in newdistmat]
avgD = np.mean([item for item in sum(newdistmat,[]) if (item !=0 and not isinstance(item, str))])
    return(minlist, meanlist, maxlist, avgD)
