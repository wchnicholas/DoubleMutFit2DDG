#!/usr/bin/python
import os
import sys
import glob
import numpy as np
from string import atof
from math import log

def minfloor(param):
  if param < 1:
    return atof(1)
  else:
    return param

def hashin(filename):
  infile = open(filename,'r')
  H = {}
  for line in infile.xreadlines():
    if 'Mut' in line: continue
    line = line.rstrip().rsplit("\t")
    mut = line[0]
    DNAcount = atof(line[1])
    DNAfreq  = DNAcount/atof(line[7])
    Selfreq  = minfloor(atof(line[6]))/atof(line[12])
    fit = Selfreq/DNAfreq
    if DNAcount > 10: #TUNABLE
      H[mut] = fit
  return H

def Scale(Shash, HScale, AASizeH):
  resis = [x[0:-1] for x in Shash.keys()]
  fmap  = {}
  for resi in resis:
    fmap[resi] = HScale[resi[0]]*1
  for mut in Shash.keys():
    resi  = mut[0:-1]
    wtaa  = mut[0]
    mutaa = mut[-1]
    fit   = Shash[mut]
    fmap[resi] += (HScale[mutaa])*fit/(AASizeH[wtaa]-AASizeH[mutaa]+1)**2
  return fmap

#HASH IN HSCALE
infile = open('Doc/HydrophobicityScale','r')
HScale = {}
for line in infile.xreadlines():
  if 'AA' in line: continue
  line = line.rstrip().rsplit("\t")
  HScale[line[0]] = atof(line[2])
infile.close()

#HASH IN AA SIZE
infile  = open('Doc/AASize','r')
AASizeH = {}
for line in infile.xreadlines():
  if 'AA' in line: continue
  line = line.rstrip().rsplit("\t")
  AASizeH[line[0]] = atof(line[1])
infile.close()

#MAIN#
Shash = hashin('Doc/SMutList')
fmap  = Scale(Shash, HScale, AASizeH)
outfile=open('Doc/Scale1','w')
for resi in fmap.keys():
  outfile.write(resi+"\t"+str(fmap[resi])+"\n")
outfile.close()
