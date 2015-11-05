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

def hashin(filename):
  infile = open(filename,'r')
  H = {}
  C = {}
  for line in infile.xreadlines():
    if 'Mut' in line: continue
    line = line.rstrip().rsplit("\t")
    mut = line[0]
    DNAcount = atof(line[1])
    DNAfreq  = DNAcount/atof(line[7])
    Selfreq  = atof(line[6])/atof(line[12])
    fit = Selfreq/DNAfreq
    H[mut] = fit
    C[mut] = DNAcount
  return H, C

Shash, C_Shash = hashin('Doc/SMutList')
Dhash, C_Dhash = hashin('Doc/DMutList')
disthash = hashindist('Doc/distance')

outfile = open('pairwiseC','w')
header  = "\t".join(['Dmut','DNAcount','dist','mut1fit','mut2fit','dmutfit'])
outfile.write(header+"\n")
for dmut in Dhash.keys():
  m1 = dmut.rsplit('-')[0]
  m2 = dmut.rsplit('-')[1]
  p1 = str(int(m1[1:-1])+1)
  p2 = str(int(m2[1:-1])+1)
  if m1[-1] == 'C' and m2[-1] == 'C':
    outfile.write("\t".join([dmut,str(int(C_Dhash[dmut])),disthash[p1+'-'+p2],str(Shash[m1]),str(Shash[m2]),str(Dhash[dmut])])+"\n")
outfile.close()
