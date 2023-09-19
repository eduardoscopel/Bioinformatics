### First try with bins using mean/median/stdev/mad/mad2
from statistics import stdev
from scipy import stats
from numpy import mean, absolute
gaps = [y-x for x , y in zip(sortedRefVec[:-1], sortedRefVec[1:])]
sd = stdev(gaps)
mad = stats.median_abs_deviation(gaps)
mad2 = mean(absolute(gaps - mean(gaps)))
avg = np.mean(gaps)
med = np.median(gaps)
lists = [[sortedRefVec[0]]]
for x in sortedRefVec[1:]:
    #if(x - lists[-1][-1]) > avg:
    if(x - lists[-1][-1]) > med:
    #if(x - lists[-1][-1]) / sd > 1:
    #if(x - lists[-1][-1]) / mad > 2:
        lists.append([])
    lists[-1].append(x)

### Second try with histogram
mylist = np.array(sortedRefVec)
bins = np.arange(min(mylist),max(mylist), sd/10)
for i in range(1,52):
    mylist[np.digitize(mylist, bins) == i]
hist=np.histogram(mylist, bins = 26)
