#!/usr/bin/python
import os
import sys
import glob
from string import atof

badmut = ['Q1A','Q1F','Q1I','Q1L','Q1M','Q1P','Q1T','Q1V']

def double(Shash):
  infile  = open('Doc/DMutList','r')
  outfile = open('result/DoubleSub.txt','w')
  header  = "\t".join(['Substitution1-WTaa','Substitution1-Pos','Substitution1-Mutaa','Substitution2-WTaa','Substitution2-Pos','Substitution2-Mutaa','InputCount','SelectionCount(SumOfTriplicates)','ExpectedFit','Sub1fit','Sub2fit'])
  outfile.write(header+"\n")
  for line in infile.xreadlines():
    if 'Mut' in line: continue
    line   = line.rstrip().rsplit("\t")
    mut1   = line[0].rsplit('-')[0] 
    mut2   = line[0].rsplit('-')[1]
    DNA    = line[1]
    sel1   = int(line[3])
    sel2   = int(line[4])
    sel3   = int(line[5])
    SEL    = line[6] 
    wtaa1  = mut1[0]
    pos1   = str(int(mut1[1:-1])+1)
    mutaa1 = mut1[-1]
    wtaa2  = mut2[0]
    pos2   = str(int(mut2[1:-1])+1)
    mutaa2 = mut2[-1]
    mut1   = wtaa1+pos1+mutaa1
    mut2   = wtaa2+pos2+mutaa2
    expfit = str("%.3f" % (Shash[mut1]*Shash[mut2]))
    m1fit  = str("%.3f" % Shash[mut1])
    m2fit  = str("%.3f" % Shash[mut2])
    out    = "\t".join([wtaa1,pos1,mutaa1,wtaa2,pos2,mutaa2,DNA,SEL,expfit,m1fit,m2fit])
    outfile.write(out+"\n")
  infile.close()
  outfile.close()

def single():
  infile = open('Doc/SMutList','r')
  outfile = open('result/SingleSub.txt','w')
  header  = "\t".join(['Substitution-WTaa','Substitution-Pos','Substitution-Mutaa','InputCount','SelectionCount(SumOfTriplicates)'])
  outfile.write(header+"\n")
  Shash = {}
  for line in infile.xreadlines():
    if 'Mut' in line: continue
    line   = line.rstrip().rsplit("\t")
    mut    = line[0]
    DNA    = line[1]
    sel1   = int(line[3])
    sel2   = int(line[4])
    sel3   = int(line[5])
    SEL    = line[6]
    wtaa   = mut[0]
    pos    = str(int(mut[1:-1])+1)
    mutaa  = mut[-1]
    out    = "\t".join([wtaa,pos,mutaa,DNA,SEL])
    outfile.write(out+"\n")
    DNAfreq  = atof(DNA)/atof(line[7])
    Selfreq  = atof(SEL)/atof(line[12])
    fit = Selfreq/DNAfreq
    Shash[wtaa+pos+mutaa] = fit
  infile.close()
  outfile.close()
  return Shash

Shash = single()
Dhash = double(Shash)
