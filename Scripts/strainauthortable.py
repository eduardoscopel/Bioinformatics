# loops through folders to create a file with the filename (strain) next to the folder name (author)

import os
import re

authorlist = list()

for root, dirs, files in os.walk(os.getcwd()):
    author = re.sub('/scratch/es47540/wkdir/pdfs/', '', root)
    with open('%s.csv' % author, 'w') as f:
        for filename in files:
            strain = re.sub('.raw.5000.pdf', '', filename)
            authorlist.append(author)
            for index in range(len(authorlist)):
                f.write(str(strain) + '\t' + str(author) + '\n')
            authorlist = list()
