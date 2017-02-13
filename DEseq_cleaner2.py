###import libraries###
import pandas as pd
from sys import argv
import numpy as np

#get data
file = argv[1]
sample_id = file.split(".")[-3].split("\\")[-1].split("_")[-1]
raw = pd.read_csv(file, sep = "\t")

#filter to get PSMs with a low Q-value (<0.1%) and that originate from a real protein, not a decoy protein
psms = raw[(raw["Q-Value (%)"] < 1/1) & raw["Target?"] == True]
#reindex the dataframe.  Don't quite remember why this is necessary, but it most certainly is.
psms = psms.set_index(np.arange(psms.count()[0]))

#make a blank dictionary.  keys will be base peptide sequence.  values will be a dict with the # of single-T PSMs # of Single-T-coding-but-Fth-modd'd PSMs, also the peptide name.
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
   if protein not in prot_dict.keys():				#If we haven't seen the  peptide yet, make a dictionary for it !
      prot_dict[protein] = {"Fth_counts":0, "Thr_counts":0}
   if Fmatch > 0:									#If there is a fluorothreonine add it to the protdict
      prot_dict[protein]["Fth_counts"] += 1
   if Fmatch == 0 and Tmatch > 0:					#ditto if there is no fluorothreonine.
      prot_dict[protein]["Thr_counts"] += 1
      
#Turn the prot_dictionary into a dataframe and flip it sideways.  Cull rows w/o hits
prot_df = pd.DataFrame(prot_dict)
prot_df = prot_df.transpose()

#Sort by frequency.  Eliminate rows with no hits.  Print it out.
prot_df = prot_df[prot_df["Thr_counts"] > 0 ]
prot_df.to_csv(("proteins" + sample_id + ".tsv"), sep ="\t")
