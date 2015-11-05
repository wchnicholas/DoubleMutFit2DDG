#!/usr/bin/python
import os
import sys
import glob
from string import atof

def hashin(filename):
  infile = open(filename,'r')
  H = {}
  for line in infile.xreadlines():
    if 'Mut' in line: continue
    line = line.rstrip().rsplit("\t")
    mut = line[0]
    DNAcount = atof(line[1])
    DNAfreq  = DNAcount/atof(line[7])
    Selfreq  = atof(line[6])/atof(line[12])
    fit = Selfreq/DNAfreq
    if DNAcount > 10: #TUNABLE
      H[mut] = fit
  return H

#MAIN#
Shash = hashin('../General/data/SMutList')
infile  = open('data/hcMatrix','r')
IDs = infile.readline().rstrip().rsplit(' ')
infile.close()
outfile = open('data/SFitMatrix','w')
for ID in IDs:
  ID = ID.replace('"','')
  outfile.write(ID+"\t"+str(Shash[ID])+"\n")
outfile.close()

