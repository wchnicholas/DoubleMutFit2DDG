#!/usr/bin/python
import os
import sys
import glob
from string import atof

def hashindist(filename):
  infile = open(filename,'r')
  H = {}
  countline = 0
  for line in infile.xreadlines():
    if countline == 0 or countline == 1: countline += 1; continue
    line  = line.rstrip().rsplit("\t")
    dists = line[8::]
    for i in range(len(dists)):
      pos1 = str(i+2)
      pos2 = line[1]
      H[pos1+'-'+pos2] = dists[i]
      H[pos2+'-'+pos1] = dists[i]
  infile.close()
  return H

def hashinhc(filename,outfile,disthash):
  infile  = open(filename,'r')
  outfile = open(outfile,'w')
  IDs = infile.readline().rstrip().rsplit(' ')
  infile.close()
  outfile.write(' '.join(IDs)+"\n")
  for ID1 in IDs[::-1]:
    outfile.write(ID1)
    p1 = str(int(ID1.replace('"','')[1:-1])+1)
    for ID2 in IDs:
      p2 = str(int(ID2.replace('"','')[1:-1])+1)
      dist = 0
      if p1 != p2: dist = disthash[str(p1)+'-'+str(p2)]
      outfile.write("\t"+str(dist))
    outfile.write("\n")
  outfile.close()
      
    


disthash = hashindist('Doc/distance')
hchash   = hashinhc('data/LandhcMatrix','data/LanddisMatrix',disthash)
#hchash   = hashinhc('data/MutPothcMatrix','data/MutPotdisMatrix',disthash)
