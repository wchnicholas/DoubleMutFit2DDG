#!/usr/bin/python
import os
import sys
import glob
import numpy as np
#import networkx as nx
#from networkx.utils import open_file, make_str
from string import atof
from scipy.stats.stats import pearsonr

def connected_components(G):
    if G.is_directed():
        raise nx.NetworkXError("""Not allowed for directed graph G.
              Use UG=G.to_undirected() to create an undirected graph.""")
    seen={}
    components=[]
    for v in G:      
        if v not in seen:
            c=nx.single_source_shortest_path_length(G,v)
            components.append(list(c.keys()))
            seen.update(c)
    components.sort(key=len,reverse=True)            
    return components  

def connected_component_subgraphs(G):
    cc=connected_components(G)
    graph_list=[]
    for c in cc:
        graph_list.append(G.subgraph(c).copy())
    return graph_list

def formatpot(l):
  newl = []
  for i in l:
    if i == 'NA': newl.append(-9)
    else: newl.append(atof(i))
  return newl

def cleanpot(l):
  newl = []
  for i in l:
    if i == 'NA': continue
    else: newl.append(atof(i))
  return newl

def cor(l1,l2):
  assert(len(l1)==len(l2))
  l1_clean = []
  l2_clean = []
  for i in range(len(l1)):
    if l1[i] != 'NA' and l2[i] != 'NA':
      l1_clean.append(atof(l1[i]))
      l2_clean.append(atof(l2[i]))
  return pearsonr(l1_clean,l2_clean)

#READ FILE
infile = open('data/MutLandscapes','r')
Phash  = {}
Fhash  = {}
Dhash  = {}
Mhash  = {}
for line in infile.xreadlines():
  if 'Fit' in line: continue
  else: 
    line = line.rstrip().rsplit("\t")
    mut  = line[0]
    fit  = line[1]
    pot  = line[2::]
#    cpot = cleanpot(pot)
#    Phash[mut] = formatpot(pot)
    Phash[mut] = pot
    Fhash[mut] = atof(fit)
#    Dhash[mut] = np.std(cpot)
#    Mhash[mut] = np.mean(cpot)
  #  print mut,fit,np.mean(cpot),np.std(cpot)
infile.close()

#COR/NETWORK
#g=nx.Graph()
#fitmin = 0
#fitmax = 99
#corcut = 0.6
count = 0
muts = sorted(Phash.keys())
#outfile = open('data/MutPotCorMatrix','w')
outfile = open('data/LandCorMatrix','w')
outfile.write('Mut'+"\t"+"\t".join(muts)+"\n")
for i in range(len(muts)):
  outfile.write(muts[i])
  for j in range(len(muts)):
    m1 = muts[i]
    m2 = muts[j]
    p1 = m1[1:-1]
    p2 = m2[1:-1]
    r  = 0
    if m1 != m2: 
      r = cor(Phash[m1],Phash[m2])[0]
      if r == 0: print m1,m2, 'cor = 0'
    outfile.write("\t"+str(r))
#    if Fhash[m1] < fitmin or Fhash[m2] < fitmin: continue
#    if Fhash[m1] > fitmax or Fhash[m2] > fitmax: continue
#    if not g.has_node(m1): g.add_node(m1)
#    if not g.has_node(m2): g.add_node(m2)
#    if cor(Phash[m1],Phash[m2])[0] > corcut: g.add_edge(m1,m2)
  outfile.write("\n")
outfile.close()
#subgraphs = connected_component_subgraphs(g)
#nx.write_dot(subgraphs[0],'out.dot')
