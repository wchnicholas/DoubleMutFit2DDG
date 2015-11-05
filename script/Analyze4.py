#!/usr/bin/python
import sys
import os
import glob
import Pycluster
import numpy as np
from string import atof
from scipy.stats.stats import pearsonr
from math import log, log10

def ceil(fit):
  if fit > 7.6: return 7.6
  else: return fit

def floor(fit):
  if fit < 0.01: return 0.01
  else: return fit

def F2E(fit):
  Bmax = 8
  return -(log(Bmax/floor(ceil(atof(fit)))-1)-log(Bmax-1))

def filterE(Es):
  newEs = []
  for E in Es:
    if atof(E) < 1: newEs.append(E)
  return newEs

def extractbybgsID(pot,bgID):
  newpot = []
  for b in bgID:
    newpot.append(pot[b])
  return newpot

def identifybgsID(muts,bgs):
  bgID = []
  for b in bgs:
    bgID.append(muts.index(b))
  return bgID

def hashinsinglemut(filename):
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

def filtertwolists(i,j):
  assert(len(i)==len(j))
  new_i = []
  new_j = []
  for n in range(len(i)):
    if np.isnan(i[n]) == False and np.isnan(j[n]) == False:
      new_i.append(i[n])
      new_j.append(j[n])
  return  new_i, new_j

def hscore(BGs):
  BGs = BGs.rsplit('-')
  l   = []
  for i in [x[0] for x in BGs]:
    l.append(HScale[i])
  return l
  
def scaling(BGs):
  BGs = BGs.rsplit('-')
  l   = []
  for i in [x[0:-1] for x in BGs]:
    l.append(Scale1[i])
  return l

#CALCULATE ENERGY
def energy(BGs,Sol):
  infile   = open('data/MutLandscapes','r')
  IDhash   = {}
  DDhash   = {}
  l_sol    = []
  allmuts  = infile.readline().rstrip().rsplit("\t")[2::]
  for i in range(len(allmuts)):
    if allmuts[i] in Ehash_L.keys():
      IDhash[allmuts[i]] = i
      DDhash[i] = []

  for BG in BGs:
    pos = BG[1:-1]
    l_sol.append(Sol[pos])

  for line in infile.xreadlines(): #HASH IN VALUES
    line = line.rstrip().rsplit("\t")
    mut  = line[0]
    fit  = atof(line[1])
    pot  = line[2::]
    if mut in BGs:
      for i in DDhash.keys():
        DDhash[i].append(pot[i])
  infile.close()
  
  l_DD = []
  l_bmark = []
  for mut in Ehash_L.keys(): #COMPILE THE MEAN/MEDIAN
    index = IDhash[mut]
    while 'NA' in DDhash[index]: DDhash[index].remove('NA')
    DD = np.median(map(atof,DDhash[index]))
    #DD = np.mean(map(atof,DDhash[index]))
    l_DD.append(DD)
    l_bmark.append(Ehash_L[mut])
    #print mut, DD, Ehash_L[mut]
  l_DD, l_bmark = filtertwolists(l_DD, l_bmark)
  BGs = '-'.join(BGs)
#  print pearsonr(l_DD,l_bmark)[0], np.mean(l_sol), np.median(l_sol), max(l_sol), min(l_sol),BGs#, F2E(Shash[BGs])
#  print BGs, pearsonr(l_DD,l_bmark)[0], np.mean(l_sol), Shash[BGs]
  return pearsonr(l_DD,l_bmark)[0], np.mean(l_sol), BGs, np.mean(hscore(BGs)), np.mean(scaling(BGs))
  
#FOR OUTPUT ENERGY
def allenergy(BGs):
  DDhash = {}
  infile   = open('data/MutLandscapes','r')
  allmuts  = infile.readline().rstrip().rsplit("\t")[2::]
  for i in range(len(allmuts)):
    DDhash[i] = []
  for line in infile.xreadlines(): #HASH IN VALUES
    line = line.rstrip().rsplit("\t")
    mut  = line[0]
    fit  = atof(line[1])
    pot  = line[2::]
    if mut in BGs:
      for i in range(len(pot)): 
        DDhash[i].append(pot[i])
  infile.close()
  print 'Mut'+"\t"+"\t".join(BGs)+"\t"+'Median'
  for i in range(len(allmuts)):
    V = "\t".join(DDhash[i])
    while 'NA' in DDhash[i]: DDhash[i].remove('NA')
    DD = np.median(map(atof,DDhash[i]))
    print allmuts[i]+"\t"+V+"\t"+str(DD)

#HASH IN BENCHMARK
Ehash_L = {}
infile = open('Doc/LiteratureE','r')
for line in infile.xreadlines():
  line = line.rstrip().rsplit("\t")
  while '' in line: line.remove('')
  mut  = line[0]
  mut  = mut[0].upper()+str(int(mut[1:-1])-1)+mut[-1].upper()
  if 'A43' in mut or 'A52' in mut: continue
  E    = line[1]
  Ehash_L[mut] = atof(E)
infile.close()

#HASH IN ROSETTA
infile = open('Doc/Rosetta1PGA','r')
Rhash  = {}
for line in infile.xreadlines():
  if 'description' in line: continue
  line = line.rstrip().rsplit(" ")
  while '' in line: line.remove('')
  mut   = line[1]
  score = line[2]
  mut  = mut[0].upper()+str(int(mut[1:-1])-1)+mut[-1].upper()
  Rhash[mut] = atof(score)
infile.close()

#HASH IN HSCALE
infile = open('Doc/HydrophobicityScale','r')
HScale = {}
for line in infile.xreadlines():
  if 'AA' in line: continue
  line = line.rstrip().rsplit("\t")
  HScale[line[0]] = atof(line[2])
infile.close()

#HASH IN SCALE1
infile = open('Doc/Scale1','r')
Scale1 = {}
for line in infile.xreadlines():
  line = line.rstrip().rsplit("\t")
  Scale1[line[0]] = atof(line[1])
infile.close()

#COMPUTE CORRELATION BETWEEN ROSETTA AND LITERATURE BENCHMARK
#for mut in Ehash_L.keys():
#  if Rhash.has_key(mut): print Ehash_L[mut], -Rhash[mut]
#sys.exit()

#Hash in single mut
Shash = hashinsinglemut('Doc/SMutList')

#COMPUTE CORRELATION BETWEEN WT BG AND LITERATURE BENCHMARK
#for mut in Ehash_L.keys():
#  if Shash.has_key(mut): print Ehash_L[mut], F2E(Shash[mut])
#sys.exit()

#Hash in solvation
Sol = {}
infile = open('Doc/1PGA.sol','r')
for line in infile.xreadlines():
  line = line.rstrip().rsplit("\t")
  Sol[line[1]] = atof(line[2])
infile.close()

#hash in landscape correlation
infile = open('data/LandCorMatrix','r')
Fhash  = {}
mutIDs  = []
mutcors = []
mutfits = []
countline = 0
minfit = 0
maxfit = 1
for line in infile.xreadlines():
  countline += 1
  if countline == 1: continue
  line = line.rstrip().rsplit("\t")
  mut  = line[0]
  pos  = mut[1:-1]
  fit  = Shash[mut]
  if fit < minfit or fit > maxfit: continue
  print mut+"\t"+str(fit)
  cor  = map(atof,line[1::])
  mutIDs.append(mut)
  mutfits.append(fit)
  mutcors.append(cor)
  Fhash[mut] = fit
infile.close()
sys.exit()#############
#for mutID in mutIDs:
#  mutID = [mutID]
#  energy(mutID,Sol)
#sys.exit()

#BGs  = ['A25T','L4N','L4S','L4T','Y2A','Y2C']
#BGs  = ['F29W','F29Y','L4N','L4S','L4T','V53S']
#BGs = ['A25T','L4S','L4T','Y2A','Y2C']
#energy(BGs,Sol)
#allenergy(BGs)
#sys.exit()

#Kmeans Clustering
k = int(sys.argv[1])
seed = int(sys.argv[2])
labels, error, nfound = Pycluster.kcluster(mutcors,k) #npass=1
H   = {}
for n in range(len(labels)):
  label = labels[n]
  mutID = mutIDs[n]
  if H.has_key(label): H[label].append(mutID)
  else: H[label] = [mutID]

lowests = 99
outputr = 0
finalg  = ''
for label in H.keys():
  BGs   = H[label]
  if len(BGs) > 3: #4
    r, s, g, h, f = energy(BGs,Sol)
    print k, r, s, '-'.join(BGs)
    if s < lowests: #Sol 
    #if h < lowests: #Hydrophobicity Scale 
    #if f < lowests: #Fitness Scale
      lowests = s #Sol
      #lowests = h #Hydro Scale
      #lowests = f
      outputr = r
      finalg  = g
#print k, outputr, lowests, finalg #NEED THIS FOR OUTPUT
