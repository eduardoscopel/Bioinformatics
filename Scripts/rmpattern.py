### remove pattern (_merged in this case) from a filename in every file in a folder

from os import rename, listdir

files = listdir('.')
pattern = '_merged'

for file in files:
  rename(file, file.replace(pattern, '', 1))
  
