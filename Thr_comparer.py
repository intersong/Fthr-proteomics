import pandas as pd
from sys import argv
import numpy as np

###Read in list of peptide-spectrum-match files  Should be either the WT or KO files.
files = ["QEPlus2_12212015_60_JM_96A.PSMs.tsv", "QEPlus2_02012016_9_96B.PSMs.tsv","QEPlus2_02012016_19_97A.PSMs.tsv", "QEPlus2_02012016_13_97b.PSMs.tsv"]
#files = argv[1]
###Make a list of data frames containing the data specified in "files"
csv_list = []
for f in files:
	pre_df = pd.read_csv(f, sep = "\t")
	csv_list.append(pre_df)

#concatenate the data frames in csv_list by stacking vertically
raw = pd.concat(csv_list, axis = 0)  #stack vertically
#filter to get PSMs with a low Q-value (<0.1%) and that originate from a real protein, not a decoy protein
psms = raw[(raw["Q-Value (%)"] < 1/1) & raw["Protein Description"].str.startswith("embl-cds")]
#reindex the dataframe.  Don't quite remember why this is necessary, but it most certainly is.
psms = psms.set_index(np.arange(psms.count()[0]))

#make a blank dictionary.  
#keys will be base peptide sequence.  
#values will be a dict with the # of single-T PSMs # of Single-T-coding-but-Fth-modd'd PSMs, also the peptide name.
seq_dict = {}

#iterate over the unified data frame of PSMs
for i in psms.index:
   #get seq with mods, base seq, and the protein from whence it all came
   seq = psms.loc[i, "Peptide Sequence"]
   bseq = psms.loc[i, "Base Peptide Sequence"]
   protein = psms.loc[i, "Protein Description"]
   #count # of threonines in the base seq and # of fluorothreonines in the mod seq.
   Tmatch = bseq.count("T")
   Fmatch = seq.count("(Fluorothreonine)")
   #If we haven't seen the base peptide yet, make a dictionary for it !
   if bseq not in seq_dict.keys():
      seq_dict[bseq] = {"Fth":0, "Thr":0, "Protein":protein}
   #Does it have one Thr and no Fluorothr?  Increment the Thr-counter by one.
   if (Tmatch == 1) and (Fmatch == 0):
      seq_dict[bseq]["Thr"] += 1
   #Does it code for one Thr but actually have one Fth? Increment the Fth-counter by 1
   if (Fmatch == 1) and (Tmatch == 1):							#this will give problems with peps with 1x F but 2x T
      seq_dict[bseq]["Fth"] += 1

#Turn the dictionary into a dataframe and flip it sideways.
seq_df = pd.DataFrame(seq_dict)
seq_df = seq_df.transpose()
#make two series of zeros of appropriate length and concatenate them with the seq_df.  Giving them the same idex as seq_df facilitates concatenation.
protcount_T = pd.Series(np.zeros(len(seq_df.index)), index = seq_df.index, name = "Thr PSMs Per Protein")
protcount_F = pd.Series(np.zeros(len(seq_df.index)), index = seq_df.index, name = "Fth PSMs Per Protein")
seq_df = pd.concat([seq_df, protcount_T, protcount_F], axis = 1)

#go through seq_df.  at each row, grab the protein ID corresponding to the PSM.  
#Count the # of times that protein ID occurred in single-T peptides and the # of times it occured in single-T, Fthr substituted peptides. 
for p in seq_df.index:
	protein = seq_df.loc[p, "Protein"]
	seq_df.loc[p, "Thr PSMs Per Protein"] = seq_df[seq_df['Protein'] == protein]["Thr"].sum()
	seq_df.loc[p, "Fth PSMs Per Protein"] = seq_df[seq_df['Protein'] == protein]["Fth"].sum()

#Sort by frequency.  Eliminate rows with no hits.  Print it out.
seq_df = seq_df.sort_values("Thr PSMs Per Protein", ascending = False)
seq_df = seq_df[seq_df["Thr"] > 0 ]
seq_df.to_csv("p0564_Fcounts_vs_Tcounts_20160916.tsv", sep ="\t")
