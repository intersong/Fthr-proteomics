import pandas as pd
from sys import argv
import numpy as np

###Read in list of peptide-spectrum-match files  Should be either the WT or KO files.
#feed it the XXXPSMs.tsv file
#files = ["QEPlus2_12212015_60_JM_96A.PSMs.tsv"]#, "QEPlus2_02012016_9_96B.PSMs.tsv","QEPlus2_02012016_19_97A.PSMs.tsv", "QEPlus2_02012016_13_97b.PSMs.tsv"]
files = argv[1]
###Make a list of data frames containing the data specified in "files"

raw = pd.read_csv(files, sep = "\t")

#filter to get PSMs with a low Q-value (<0.1%) and that originate from a real protein, not a decoy protein
psms = raw[(raw["Q-Value (%)"] < 1/1) & raw["Protein Description"].str.startswith("embl-cds")]
#reindex the dataframe.  Don't quite remember why this is necessary, but it most certainly is.
psms = psms.set_index(np.arange(psms.count()[0]))

#make a blank dictionary.  
#keys will be base peptide sequence.  
#values will be a dict with the # of single-T PSMs # of Single-T-coding-but-Fth-modd'd PSMs, also the peptide name.
proteins_list_F = []
proteins_list = []

#iterate over the unified data frame of PSMs
for i in psms.index:
   seq = psms.loc[i, "Peptide Sequence"]
   Fmatch = seq.count("(Fluorothreonine)")
   proteins_list.append(psms.loc[i, "Protein Description"])
   if Fmatch > 0:
      proteins_list_F.append(psms.loc[i, "Protein Description"])


print("# of proteins with F")
print(len(set(proteins_list_F)))
print("# of proteins total")
print(len(set(proteins_list)))
