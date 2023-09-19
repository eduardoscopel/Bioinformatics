import time
start = time.time()
import re
chr=[]
c=0
i=0
num_lines = sum(1 for line in open("19F.vcf", 'rb'))
with open("19F.vcf", 'rb') as f:
    for lines in f:
        if lines.startswith('#'):
            num_lines = num_lines-1
        else:
            #if re.search('AF1=0;', lines):
            #    continue
            #else:
            try:
                if c == 0:
                    pos=[0]*num_lines
                    af=[0]*num_lines
                    chr.append(re.search('(.+?)\t',lines).group(1))
                    print chr[c]
                    c+=1
                else:
                    if re.search('(.+?)\t',lines).group(1) != chr[c-1]:
                        chr.append(re.search('(.+?)\t',lines).group(1))
                        print chr[c]
                        c+=1
            except AttributeError:
                print 'chr not found'
            try:
                #pos.append(re.search('\t(.+?)\t',lines).group(1))
                #pos[k] = int(pos[k])
                pos[i] = int(re.search('\t(.+?)\t',lines).group(1))
                af[i] = re.split(r',',(re.search('DP4=(.+?);', lines).group(1)))
                for j in range(len(af[i])):
                    af[i][j] = int(af[i][j])
            except AttributeError:
                print 'pos/af not found'
            #try:
                #af.append(re.search('DP4=(.+?);', lines).group(1))
                #af[k] = re.split(r',',af[k])
            #except AttributeError:
            #    print 'af not found'
            i+=1


#for i in range(len(pos)):
#    if i == 0:
#        print chr[i],pos[i]
#        continue
#    else:
#        if pos[i]<pos[i-1]:
#            print chr[i], pos[i]



end = time.time()
print (end-start)
