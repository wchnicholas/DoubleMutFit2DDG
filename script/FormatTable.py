#!/usr/bin/python
import os
import sys
import glob

infile = open('data/Hit2result','r')
for line in infile.xreadlines():
  if 'Mut' in line: 
    print line.rstrip()
  else:
    line = line.rstrip().rsplit("\t")
    mut  = line[0]
    info = line[1::]
    print mut[0]+str(int(mut[1:-1])+1)+mut[-1]+"\t"+"\t".join(info)
infile.close()
