import pandas as pd
from sys import argv
import numpy as np

#get protein descriptions:
descript_df = pd.read_csv("D:\jmcmurry\Documents\Davis_data\id_dictionary.tsv", sep = "\t")
descript_dict = {}
for x in descript_df.index:
	descript_dict[descript_df.loc[x, "ID"]] = descript_df.loc[x, "description"]
	
#get output file name
output = argv[1]

###Read in list of peptide-spectrum-match files  Should be either the WT or KO files.
file = "all_PSMs.csv" #.tsv" #, "QEPlus2_02012016_9_96B.PSMs.tsv","QEPlus2_02012016_19_97A.PSMs.tsv", "QEPlus2_02012016_13_97b.PSMs.tsv"]

###Make a data frame containing the data specified in "files"
raw = pd.read_csv(file, sep = ",")

#filter to get PSMs with a low Q-value (<0.1%) and that originate from a real protein, not a decoy protein
psms = raw[(raw["Q-Value (%)"] < 1/1) & raw["Target?"] == True]
#reindex the dataframe.  Don't quite remember why this is necessary, but it most certainly is.
psms = psms.set_index(np.arange(psms.count()[0]))

#make a blank dictionary.  
#keys will be base peptide sequence.  
#values will be a dict with the # of single-T PSMs # of Single-T-coding-but-Fth-modd'd PSMs, also the peptide name.
seq_dict = {}
prot_dict = {}
#iterate over the unified data frame of PSMs
for i in psms.index:
   #get seq with mods, base seq, and the protein from whence it all came
   seq = psms.loc[i, "Peptide Sequence"]
   bseq = psms.loc[i, "Base Peptide Sequence"]
   protein = psms.loc[i, "Protein Description"]
   #count # of threonines in the base seq and # of fluorothreonines in the mod seq.
   Tmatch = bseq.count("T")
   Fmatch = seq.count("(Fluorothreonine)")
   #If we haven't seen the  peptide yet, make a dictionary for it !
   if (Fmatch > 0) :
	   if seq not in seq_dict.keys():
		   seq_dict[seq] = {"Fth":1, "Protein":protein}
	   else:
		   seq_dict[seq]["Fth"] += 1
	   if protein not in prot_dict.keys():
		   prot_dict[protein] = {"Fth":1}
	   else:
		   prot_dict[protein]["Fth"] += 1
		   
#Turn the seq_dictionary into a dataframe and flip it sideways.
seq_df = pd.DataFrame(seq_dict)
seq_df = seq_df.transpose()

#Sort by frequency.  Eliminate rows with no hits.  Print it out.
seq_df = seq_df.sort_values("Fth", ascending = False)
seq_df = seq_df[seq_df["Fth"] > 0 ]
seq_df.to_csv(("peptides" + output + ".tsv"), sep ="\t")


#Turn the prot_dictionary into a dataframe and flip it sideways.
prot_df = pd.DataFrame(prot_dict)
prot_df = prot_df.transpose()

#Sort by frequency.  Eliminate rows with no hits.  Print it out.
prot_df = prot_df.sort_values("Fth", ascending = False)
prot_df = prot_df[prot_df["Fth"] > 0 ]
prot_descript = pd.Series(" " * len(prot_df.index), index = prot_df.index, name = "Description")
print("made_series")
for i in prot_df.index :
	prot_df.loc[i,"Description"] = descript_dict[i]
	print(i)

prot_df.to_csv(("proteins" + output + ".tsv"), sep ="\t")
