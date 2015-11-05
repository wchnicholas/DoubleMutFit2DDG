#!/usr/bin/python
import os
import sys
import glob
import copy
import networkx as nx
import matplotlib.pyplot as plt
from string import atof
from scipy.stats.stats import pearsonr

def OptimizeCor(E1,E2):
  pvalue = pearsonr(E1,E2)[1]
  print pearsonr(E1,E2)
  for i in range(len(E1)):
    E1_tmp = copy.deepcopy(E1)
    E2_tmp = copy.deepcopy(E2)
    del E1_tmp[i]
    del E2_tmp[i]
    print len(E1_tmp), pearsonr(E1_tmp,E2_tmp)

aas = ['E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W']

#HASH IN ENERGY 
infile = open('data/Hit2result','r')
Ehash  = {}
for line in infile.xreadlines():
  if 'Median' in line: continue
  line = line.rstrip().rsplit("\t")
  mut  =  line[0]
  aa   =  mut[-1]
  pos  =  int(mut[1:-1])
  E    =  line[-1]
  if Ehash.has_key(pos): Ehash[pos][aa] = E
  else: Ehash[pos] = {aa:E}
infile.close()

#FILTERING POSITION WITH INSUFFICIENT INFO
for pos in Ehash.keys():
  if len(Ehash[pos].keys()) < 10:
    del Ehash[pos]

#NETWORK BUIDLILNG
poss = sorted(Ehash.keys())
outfile = open('data/PropCorMatrix','w')
header = 'pos'+"\t"+"\t".join(map(str,poss))
outfile.write(header+"\n")
for p1 in poss:
  outfile.write(str(p1))
  for p2 in poss:
    E1 = []
    E2 = []
    for aa in aas:
      if Ehash[p1].has_key(aa) and Ehash[p2].has_key(aa):
        e1 = Ehash[p1][aa]
        e2 = Ehash[p2][aa]
        if e1 != 'nan' and e2!= 'nan':
          E1.append(atof(e1))
          E2.append(atof(e2))
    outfile.write("\t"+str(pearsonr(E1,E2)[0]))
  outfile.write("\n")
outfile.close()
