#!/usr/bin/python
import os
import sys
import glob
import copy
import rpy
import operator
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import itertools
from pymodelfit.builtins import LinearModel
from scipy import stats
from string import atof
from scipy.stats.stats import pearsonr

def cal_ss_tot(l):
  ss_tot = 0
  for i in l:
    ss_tot += (i-np.mean(l))**2
  return ss_tot

def cal_ss_res(x,y,slope,itc):
  ss_res = 0
  for i in range(len(x)):
    ss_res += (y[i]-(slope*x[i]+itc))**2
  return ss_res
    
#LinearModel: count, return R square, intercept
#r means r-square
def OptimizeCor(E1,E2,AA):
  r = 0
  cutoff = 0.8
  assert(len(E1) == len(E2))
  count = len(E1)
  while r < cutoff:
    if count == 1: return count, 'bad', 'bad'
    indexs = itertools.combinations(range(len(E1)), count)
    for index in indexs:
      x  = np.array([E1[i] for i in index])
      y  = np.array([E2[i] for i in index])
      aa = [AA[i] for i in index]
      lm     = LinearModel.fitBasic(x,y,fixslope=1)[0]
      slope  = lm[0]
      itc    = lm[1]
      ss_res = cal_ss_res(x,y,slope,itc)
      ss_tot = cal_ss_tot(y)
      r      = 1-ss_res/ss_tot
      if r > cutoff: return count, r, slope
    count += -1

def numberOfNonNans(data):
    count = 0
    for i in data:
        if not np.isnan(i):
            count += 1
    return count 
  
############MAIN################
aas = ['E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W']

#HASH IN BETA-PROPENSITY SCALE
infile = open('data/BetaPropensity','r')
Phash  = {}
for line in infile.xreadlines():
  if 'aa' in line: continue
  line = line.rstrip().rsplit("\t")
  aa = line[0]
  if aa != 'P':
    Phash[aa] = atof(line[1])
infile.close()

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

for pos in Ehash.keys(): 
  if numberOfNonNans(map(atof,Ehash[pos].values())) < 10:
    del Ehash[pos]

#NETWORK BUIDLILNG
poss = sorted(Ehash.keys())
header = ['position','totalaa','remainaa','r2','slope']
outfile = open('data/BetaCorResult','w')
outfile.write("\t".join(header)+"\n")
for pos in poss:
  print 'working on position', pos
  E1 = []
  E2 = []
  AA = []
  for aa in aas:
    if Phash.has_key(aa) and Ehash[pos].has_key(aa):
      e1 = Phash[aa]
      e2 = Ehash[pos][aa]
      if e1 != 'nan' and e2!= 'nan':
        E1.append(atof(e1))
        E2.append(atof(e2))
        AA.append(aa)
  if len(E1) > 0 and len(E2) > 0:
    outfile.write(str(pos)+"\t"+str(len(E1))+"\t"+"\t".join(map(str,OptimizeCor(E1,E2,AA)))+"\n")
outfile.close()
