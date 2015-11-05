#!/usr/bin/python
import os
import sys
import glob
from string import atof

benchmarkpos = [8,2,4,16,55,7,53,18,15,17,19,42]
infile = open('data/Hit2result','r')
phash   = {}
phash_A = {}
for line in infile.xreadlines():
  if 'Mut' in line: continue
  line = line.rstrip().rsplit("\t")
  mut = line[0]
  fit = line[-1]
  pos = int(mut[1:-1])+1
  aa  = mut[-1]
  wta = mut[0]
  if pos in benchmarkpos: 
    phash[str(pos)+aa] = atof(fit)
    phash[str(pos)+wta] = atof(0)
    if aa == 'A':
      phash_A[str(pos)] = atof(fit)
infile.close()

aas = ['E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W']
outfile = open('data/BetaBenchMark','w')
outfile.write('aa'+"\t"+"\t".join(map(str,benchmarkpos))+"\n")
for aa in aas:
  outfile.write(aa)
  for pos in benchmarkpos:
    pos = str(pos)
    outfile.write("\t"+str(phash[pos+aa]-phash_A[pos]))
  outfile.write("\n")
outfile.close()
