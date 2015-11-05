#!/usr/bin/python
import os
import sys
import glob

infile = open('data/PropCorhc','r')
IDs    = infile.readline().rstrip().replace('"','').rsplit(' ')[::-1]
infile.close()

sshash  = {}
infile  = open('Doc/1PGA.SS','r')
for line in infile.xreadlines():
  line = line.rstrip().rsplit(' ')
  while '' in line: line.remove('')
  pos = str(int(line[3])-1)
  ss  = line[5]
  sshash[pos] = ss
infile.close()

solhash = {}
infile = open('Doc/1PGA.sol','r')
for line in infile.xreadlines():
  line = line.rstrip().rsplit("\t")
  pos = line[1]
  sol = line[2]
  solhash[pos] = sol
infile.close()

outfile = open('data/PropInfo','w')
header = "\t".join(['Pos','SS','Sol'])
outfile.write(header+"\n")
for ID in IDs:
  outfile.write(ID+"\t"+sshash[ID]+"\t"+solhash[ID]+"\n")
outfile.close()
