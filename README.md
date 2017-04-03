This readme describes the scripts used for the analyses in:  
[High-throughput identification of protein mutant stability computed from a double mutant fitness landscape](http://onlinelibrary.wiley.com/doi/10.1002/pro.2840/full)
* Explanation of relevant files for the analyses in the paper is provided (see section "FILE DESCRIPTIONS")
* Protocol for running the analysis is described (see section "STEP BY STEP PROTOCOL")
* Other scripts/files that are not used for the analyses in the paper. They were generated during "idea testing" stage of the project. Some of them are described (see section "OTHER SCRIPTS")

## FILE DESCRIPTIONS
* doc/SMutList:  Single substitution data from [Olson et al. 2014](http://www.cell.com/current-biology/abstract/S0960-9822(14)01268-8)
* doc/DMutList:  Double substitution data from [Olson et al. 2014](http://www.cell.com/current-biology/abstract/S0960-9822(14)01268-8)
* data/SingleECor: All correlation (R\_literature, SASA, fitness) for background mutation with non-zero and < 1 fitness
* data/n4to8Sim: Simulation result using minfit = 0.4, maxfit = 0.8, k ranges from 1 to 50, 100 different seedings
* data/MutLandscapes: DDG of unfolding of mutation B (column name) based on S\_BG (i.e. mutation A, row name). Also see equation 11 in the paper
* data/LandCorMatrix: A symmetric matrix that states the Pearson correlation coefficient between the mutational energy profile computed based on different S\_BG
* data/LandhcMatrix: Hierachy clustering of the S\_BG based on the correlation matrix (LandCorMatrix). Values are the Pearson correlation coefficient. This file is basically a rearrangement of data/LandCorMatrix
* Doc/LiteratureE: 84 DDG from Literature as benchmark. Two of them (A44R, and A53E) are not used in this study due to the mismatch of WT amino acid. The first column is the substitution identity and the second column is the DDG value from Literature. This table is the same as Table S4 from [Olson et al. 2014](http://www.cell.com/current-biology/abstract/S0960-9822(14)01268-8)
* Doc/Rosetta1PGA: Output from Rosetta DDG\_monomer prediction using 1PGA as input
* Doc/HydrophobicityScale: Hydrophobicity scales for different amino acids. Hopp-Woods scale (Third column) was used for this study.
* Doc/Scale1: A proposed scale as a function of amino acid size and Hydrophobicity scales. This was generated during idea testing stage but were not used for the final version of analysis in this study (variable f in line 263 of Analysis4.py).
* Doc/1PGA.sol: RSA of each residue on Protein GB 1 (based on 1PGA).

## STEP BY STEP PROTOCOL
1. python script/Format1.py: fitness of mutation A or B (Potentiation or Landscape) on background mutation
  * Input file:
    * Doc/DMutList
    * Doc/SMutList
  * Output file: data/MutLandscapes

2. python script/Analyze1.py: Produce a correlation matrix based on single substitution fitness landscape among all mutation backgrounds
  * Input file: data/MutLandscapes
  * Output file: data/LandCorMatrix

3. R CMD BATCH script/Heatmap1.R: Plot a heatmap by clustering mutation background based on fitness landscape similarity
  * Input file: data/LandCorMatrix
  * Output file: 
    * data/LandhcMatrix
    * graph/LandCorHeatmap.png

4. python script/Analyze4.py [k] [RANDOM SEED]: Simulation to generate correlating cluster, and some other functions to generate hit result. This script will produce a standard output arranged as: [k]\t[Rliterature]\t[meanRSA]\t[S\_BGgroup]. Fitness range filtering of S\_BG can be controled by varaibles minfit and maxfit in line 216 and 217. 
  * Input file:
    * Doc/SMutList
    * data/MutLandscapes
    * Doc/LiteratureE
    * Doc/Rosetta1PGA
    * Doc/HydrophobicityScale
    * Doc/Scale1
    * Doc/1PGA.sol
    * data/LandCorMatrix
  * script/Sim1.sh allows iterate run with different k and random seeding
  * script/tmp.sh allows parsing the output of Sim1.sh to search for S\_BGgroup with highest meanRSA within a given k among different random seeding


## OTHER SCRIPTS
script/Heatmap2.R: Plot a heatmap by distance based on the clustering in data/hcMatrix
  * Input file: data/LanddisMatrix
  * Output file: graph/LandDisHeatMap.png

script/Format2.py: Format Distance matrix according to the clustering in data/hcMatrix
  * Input file:
    * Doc/distance
    * data/LandhcMatrix
  * Output file: data/LanddisMatrix

script/Analyze5.py: Generate position correlation matrix by propensity
  * Inout file: data/Hit2result
  * Output file: data/PropCorMatrix

script/Heatmap3.R: Plot a heatmap by distance based on the clustering in data/PropCorMatrix 
  * Input file: data/PropCorMatrix
  * Output file: 
    * data/PropCorhc 
    * graph/PropCor.png

script/Analyze6.py: Generate a compile table that contains the secondary structure classification (SS), and RSA (based on 1PGA) infomration for each residue. 
  * Input file: 
    * Doc/1PGA.SS
    * data/PropCorhc
    * Doc/1PGA.sol
  * Output file: data/PropInfo
