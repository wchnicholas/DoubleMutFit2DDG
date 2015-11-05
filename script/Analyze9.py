#hash in landscape correlation
from string import atof

infile = open('data/LandCorMatrix','r')
row    = []
col    = []
table  = []
countline = 0
for line in infile.xreadlines():
  countline += 1
  line = line.rstrip().rsplit("\t")
  if countline == 1: row = line[1::]
  else:
    col.append(line[0])
    table.append(line[1::])
infile.close()

#hash in rsa
infile = open('Doc/SMutProperty','r')
RSAhash = {}
for line in infile.xreadlines():
  if 'Mut' in line: continue
  line = line.rstrip().rsplit(' ')
  RSAhash[line[0]] = atof(line[2])
infile.close()

#MAIN#
outfile = open('data/PairwiseCorvsRSA','w')
outfile.write("\t".join(['Mut','Cor','RSA'])+"\n")
for c in range(0,len(col)):
  for r in range(0,len(row)):
    m1  = col[c]
    m2  = row[r]
    if m1 != m2 and c > r:
      cor = table[c][r]
      mut = m1+'-'+m2
      rsa = abs(RSAhash[m1]-RSAhash[m2])
      outfile.write("\t".join(map(str,[mut,cor,rsa]))+"\n")
outfile.close()

