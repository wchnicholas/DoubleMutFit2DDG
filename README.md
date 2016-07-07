THIS README DESCRIBES THE SCRIPTS (SEE "STEP BY STEP PROTOCOL") USED FOR THE ANALYSES IN:  
[High-throughput identification of protein mutant stability computed from a double mutant fitness landscape](http://onlinelibrary.wiley.com/doi/10.1002/pro.2840/full)

#EXPLANATION OF RELEVANT FILES FOR THE ANALYSES IN THE PAPER ARE PROVIDED (SEE "FILE DESCRIPTION")
#OTHER SCRIPTS/FILES ARE NOT USED FOR THE ANALYSES IN THE PAPER. THEY WERE GENERATED DURING "IDEA TESTING" STAGE OF THE PROJECT. SOME OF THEM ARE DESCRIBED (SEE "OTHER SCRIPTS"). 

#################################################################
#####################STEP BY STEP PROTOCOL#######################
#################################################################
#fitness of mutation A or B (Potentiation or Landscape) on background mutation; output file = data/MutPotentiation or data/MutLandscapes
#Input file: Doc/DMutList and Doc/SMutList
#Output file: data/MutLandscapes
python script/Format1.py

#Produce a correlation matrix based on single substitution fitness landscape among all mutation backgrounds; output file = data/CorMatrix
#Input file: data/MutLandscapes
#Output file: data/LandCorMatrix
python script/Analyze1.py

#Plot a heatmap by clustering mutation background based on fitness landscape similarity; output file = data/hcMatrix
#Input file: data/LandCorMatrix
#Output file: data/LandhcMatrix, graph/LandCorHeatmap.png
R CMD BATCH script/Heatmap1.R

#Simulation to generate correlating cluster, and some other functions to generate hit result
#This script will produce a standard output arranged as: [k]\t[Rliterature]\t[meanRSA]\t[S_BGgroup]
#Fitness range filtering of S_BG can be controled by varaibles minfit and maxfit in line 216 and 217. 
#script/Sim1.sh allows iterate run with different k and random seeding
#script/tmp.sh allows parsing the output of Sim1.sh to search for S_BGgroup with highest meanRSA within a given k among different random seeding
#Input file: Doc/SMutList, data/MutLandscapes, Doc/LiteratureE, Doc/Rosetta1PGA, Doc/HydrophobicityScale, Doc/Scale1, Doc/1PGA.sol, data/LandCorMatrix
python script/Analyze4.py [k] [RANDOM SEED]

#-----------------END OF STEP BY STEP PROTOCOL------------------#


#################################################################
########################FILE DESCRIPTION#########################
#################################################################
data/SingleECor: All correlation (R_literature, SASA, fitness) for background mutation with non-zero and < 1 fitness
data/n4to8Sim: Simulation result using minfit = 0.4, maxfit = 0.8, k ranges from 1 to 50, 100 different seedings
Doc/SMutList: Single substitution data from (Olson et al. PMID: 25455030)
Doc/DMutList: Double substitution data from (Olson et al. PMID: 25455030)
data/MutLandscapes: DDG of unfolding of mutation B (column name) based on S_BG (i.e. mutation A, row name). Also see equation 11 in the paper
data/LandCorMatrix: A symmetric matrix that states the Pearson correlation coefficient between the mutational energy profile computed based on different S_BG
data/LandhcMatrix: Hierachy clustering of the S_BG based on the correlation matrix (LandCorMatrix). Values are the Pearson correlation coefficient. This file is basically a rearrangement of data/LandCorMatrix
Doc/LiteratureE: 84 DDG from Literature as benchmark. Two of them (A44R, and A53E) are not used in this study due to the mismatch of WT amino acid. The first column is the substitution identity and the second column is the DDG value from Literature. 
Doc/Rosetta1PGA: Output from Rosetta DDG_monomer prediction using 1PGA as input
Doc/HydrophobicityScale: Hydrophobicity scales for different amino acids. Hopp-Woods scale (Third column) was used for this study.
Doc/Scale1: A proposed scale as a function of amino acid size and Hydrophobicity scales. This was generated during idea testing stage but were not used for the final version of analysis in this study (variable f in line 263 of Analysis4.py).
Doc/1PGA.sol: RSA of each residue on Protein GB 1 (based on 1PGA).

#-------------------END OF FILE DESCRIPTION---------------------#

#################################################################
#########################OTHER SCRIPTS###########################
#################################################################

#Plot a heatmap by distance based on the clustering in data/hcMatrix
#Input file: data/LanddisMatrix
#Output file: graph/LandDisHeatMap.png
R CMD BATCH script/Heatmap2.R

#Format Distance matrix according to the clustering in data/hcMatrix; output file = data/disMatrix
#Input file: Doc/distance, data/LandhcMatrix
#Output file: data/LanddisMatrix
python script/Format2.py

#Generate position correlation matrix by propensity
#Inout file: data/Hit2result'
#Output file: data/PropCorMatrix
python script/Analyze5.py

#Plot a heatmap by distance based on the clustering in data/PropCorMatrix 
#Input file: data/PropCorMatrix
#Output file: data/PropCorhc, graph/PropCor.png
R CMD BATCH script/Heatmap3.R

#Generate a compile table that contains the secondary structure classification (SS), and RSA (based on 1PGA) infomration for each residue. 
#Input file: Doc/1PGA.SS, data/PropCorhc, Doc/1PGA.sol
#Output file: data/PropInfo
python script/Analyze6.py

#---------------------END OF OTHER SCRIPTS----------------------#
