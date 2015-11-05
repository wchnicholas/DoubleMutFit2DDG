#!/usr/bin/python
import os
import sys
import glob

def main():
  filenames = glob.glob('data/Group*ECor')
  for infile in filenames:
    outfile = infile+'.adjpos'
    infile  = open(infile,'r')
    outfile = open(outfile,'w')
    for line in infile.xreadlines():
      if 'Mut' in line: outfile.write(line); continue
      line = line.rstrip().rsplit(' ')
      Mut  = line[0]
      R    = line[1]
      RSA  = line[2]
      Fit  = line[3]
      outfile.write(' '.join([Mut[0]+str(int(Mut[1:-1])+1)+Mut[-1],R,RSA,Fit])+"\n")
    infile.close()
    outfile.close()
 
if __name__ == '__main__':
  main()
