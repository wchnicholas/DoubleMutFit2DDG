#!/usr/bin/python
import os
import sys
import glob
import numpy
import numpy as np
import networkx as nx
from scipy.stats.stats import pearsonr
from string import atof
from math import log

def min_weighted_vertex_cover(graph, weight=None):
  weight_func = lambda nd: nd.get(weight, 1)
  cost = dict((n, weight_func(nd)) for n, nd in graph.nodes(data=True))

  for u,v in graph.edges_iter():
      min_cost = min([cost[u], cost[v]])
      cost[u] -= min_cost
      cost[v] -= min_cost
  return set(u for u in cost if cost[u] == 0)

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


def hashin(filename):
  H = {}
  infile  = open(filename,'r')
  IDs = infile.readline().rstrip().rsplit("\t")
  IDs = IDs[1::]
  for line in infile.xreadlines():
    line = line.rstrip().rsplit("\t")
    mutA = line[0]
    cors = map(atof,line[1::])
    for n in range(len(IDs)):
      mutB = IDs[n]
      H[mutA+'-'+mutB] = cors[n]
  infile.close()
  return H

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

#HASHIN SINGLE MUT FITNESS
Shash = hashinsinglemut('../General/data/SMutList')


#READ IN ID
infile  = open('data/MutPotCorMatrix','r')
IDs = infile.readline().rstrip().rsplit("\t")
IDs = IDs[1::]
infile.close()

#HASH IN MUT POT AND LANSCAPES
Pothash = hashin('data/MutPotCorMatrix')
Lanhash = hashin('data/LandCorMatrix')
BGs     = []
g       = nx.Graph()
lowcut  = 0.4
highcut = 0.8
for i in range(len(IDs)):
  for j in range(len(IDs)):
    if i < j:
      ID1 = IDs[i]
      ID2 = IDs[j]
      ID  = ID1+'-'+ID2
      if Shash[ID1] < lowcut or Shash[ID2] < lowcut: continue
      if Shash[ID1] > highcut or Shash[ID2] > highcut: continue
      if Lanhash[ID] > 0.95:
        BGs.extend([ID1,ID2])
        g.add_node(ID1)
        g.add_node(ID2)
        g.add_edge(ID1,ID2)

#BGs = sorted(list(set(BGs)))

#BGs = list(min_weighted_vertex_cover(g))

#BGs = []
#for node in g.nodes():
#  if g.degree(node) < 30: BGs.append(node)

subgraphs = connected_component_subgraphs(g)
BGs = subgraphs[2].nodes()
#print BGs
#for subgraph in subgraphs:
#  print len(subgraph.nodes())
#  print subgraph.nodes()
#sys.exit()

#BGs = ['Y2A','Y2S','L4S','L4T','A25T']
#HASH IN BENCHMARK
infile = open('../Epistasis/data/SMutFE','r')
Ehash  = {}
for line in infile.xreadlines():
  line = line.rstrip().rsplit("\t")
  Ehash[line[0]] = atof(line[1])
infile.close()

#HASH IN BENCHMARK
Ehash_L = {}
infile = open('Doc/LiteratureE','r')
for line in infile.xreadlines():
  line = line.rstrip().rsplit("\t")
  while '' in line: line.remove('')
  mut  = line[0]
  mut  = mut[0].upper()+str(int(mut[1:-1])-1)+mut[-1].upper()
  E    = line[1]
  Ehash_L[mut] = atof(E)
infile.close()


#CALCULATE ENERGY
infile = open('data/MutLandscapes','r')
bgID   = identifybgsID(infile.readline().rstrip().rsplit("\t")[2::],BGs)
l_mean   = []
l_median = []
l_bmark  = []
for line in infile.xreadlines():
  line = line.rstrip().rsplit("\t")
  mut  = line[0]
  fit  = atof(line[1])
  pot  = extractbybgsID(line[2::],bgID)
  while 'NA' in pot: pot.remove('NA')
#  if fit < 0.24: continue
  if not Ehash_L.has_key(mut): continue
  pot  = sorted(map(atof,pot))
  l_mean.append(np.mean(pot))
  l_median.append(np.median(pot))
  l_bmark.append(Ehash_L[mut])
#  print mut,F2E(ceil(floor(fit))),F2E(np.mean(pot)),F2E(np.median(pot)),Ehash[mut]
infile.close()

print 'Mean vs Bmark:', pearsonr(l_mean,l_bmark)[0]
print 'Median vs Bmark:', pearsonr(l_median,l_bmark)[0]
print 'Bmark vs Bmark:', pearsonr(l_bmark,l_bmark)[0]


