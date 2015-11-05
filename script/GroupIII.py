#!/usr/bin/python
import os
import sys
import glob
from string import atof

infile = open('data/GroupIIECor','r')
totalcount = 0
count = 0
dynamic = [7,9,11,12,14,16,33,37,38,40,54,56]
for line in infile.xreadlines():
  line = line.rstrip().rsplit(' ')
  mut = line[0]
  pos = int(mut[1:-1])+1
  totalcount += 1
  if pos in dynamic:
    count += 1
#    print mut[0]+str(pos)
infile.close()
print count, totalcount, atof(count)/atof(totalcount)
