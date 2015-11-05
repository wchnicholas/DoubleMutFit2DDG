#!/usr/bin/python
import os
import sys
import glob
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
      #print "\t".join(map(str,[mut, DNAcount, fit]))
  return H

def floor(param):
  if param < 0.01:
    return atof(0.01)
  else:
    return param

def cap(param):
  if param > 7.6:
    return atof(7.6)
  else:
    return param

def fit2energy(param):
  #param = -(log(8/atof(param)-1)-log(7))
  param = 7/(8/atof(param)-1)
  return param

def extractDM(Shash,Dhash,outfile):
  outfile = open(outfile,'w')
  Smuts = sorted(Shash.keys())
  outfile.write('Background'+"\t"+'Fit'+"\t"+"\t".join(Smuts)+"\n")
  countben = 0
  for mut1 in Smuts:
    if Shash[mut1] > 1: countben += 1; continue
    outfile.write(mut1+"\t"+str(Shash[mut1]))
    for mut2 in Smuts:
      Dmut= ''
      fit = '' #fit actually means fraction folded
      Smut1Fit = fit2energy(cap(floor(atof(Shash[mut1]))))
      Smut2Fit = fit2energy(cap(floor(atof(Shash[mut2]))))
      if Dhash.has_key(mut1+'-'+mut2): Dmut = mut1+'-'+mut2
      elif Dhash.has_key(mut2+'-'+mut1): Dmut = mut2+'-'+mut1
      if Dmut != '':
        #fit = fit2energy(cap(floor(atof(Dhash[Dmut]))))/Smut1Fit/Smut2Fit*Shash[mut1] #
        fit = Dhash[Dmut]/Shash[mut2] #Another way of calculating, this way without capping and flooring matches anders
        if fit >= 1 or fit <= 0: fit = 'NA'
        else: fit = -1.9858775*296/1000*log((1-fit)/fit) - (-1.9858775*296/1000*log((1-Shash[mut1])/Shash[mut1])) #
      else: fit = 'NA'
      outfile.write("\t"+str(fit))
    outfile.write("\n")
  outfile.close()
  print "There is a total of %d beneficial mutations" % countben

#MAIN#
Dhash = hashin('Doc/DMutList')
Shash = hashin('Doc/SMutList')
extractDM(Shash,Dhash,'data/MutLandscapes')
