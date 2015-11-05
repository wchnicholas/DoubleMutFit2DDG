#!/usr/bin/python
import os
import sys
import glob
import numpy as np
from string import atof
from math import log

def ceil(fit):
  if fit > 7.6: return 7.6
  else: return fit

def floor(fit):
  if fit < 0.01: return 0.01
  else: return fit

def F2E(fit):
  Bmax = 8
  return -(log(Bmax/floor(ceil(atof(fit)))-1)-log(Bmax-1))

def formatbgs(bgs):
  newbgs = []
  for b in bgs:
    newbgs.append(b[1:-1])
  return newbgs

def identifybgsID(muts,bgs):
  bgID = []
  for b in bgs:
    bgID.append(muts.index(b))
  return bgID

def extractbybgsID(pot,bgID):
  newpot = []
  for b in bgID:
    newpot.append(pot[b])
  return newpot

#HASH IN BENCHMARK
infile = open('../Epistasis/data/SMutFE','r')
Ehash  = {}
for line in infile.xreadlines():
  line = line.rstrip().rsplit("\t")
  Ehash[line[0]] = atof(line[1])
infile.close()

#SEARCH CORRELATION
infile = open('data/hcMatrix','r')
bgs    = infile.readline().rstrip().rsplit(' ')
infile.close()
bgs = formatbgs(bgs[bgs.index('"N34Y"')::])

#bgs = ['Y2A','Y2S','L4S','L4T','A25T'] #BENCHMARK BACKGROUND

#MAIN#
infile = open('data/MutPotentiation','r')
bgID   = identifybgsID(infile.readline().rstrip().rsplit("\t")[2::],bgs)
for line in infile.xreadlines():
  line = line.rstrip().rsplit("\t")
  mut  = line[0]
  fit  = atof(line[1])
  pot  = extractbybgsID(line[2::],bgID)
  while 'NA' in pot: pot.remove('NA')
  pot  = sorted(map(atof,pot))
  print mut,F2E(ceil(floor(fit))),F2E(np.mean(pot)),F2E(np.median(pot)),Ehash[mut]
infile.close()
