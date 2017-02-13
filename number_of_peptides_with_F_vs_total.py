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
psms = raw[(raw["Q-Value (%)"] < 1/1) & raw["Target?"] == 1]
#reindex the dataframe.  Don't quite remember why this is necessary, but it most certainly is.
psms = psms.set_index(np.arange(psms.count()[0]))

Total_psm_count = 0
F_psm_count = 0
#iterate over the unified data frame of PSMs
for i in psms.index:
   #increment total counter
   Total_psm_count += 1
   #get seq with mods, base seq, and the protein from whence it all came
   seq = psms.loc[i, "Peptide Sequence"]
   #count # of threonines in the base seq and # of fluorothreonines in the mod seq.
   Fmatch = seq.count("(Fluorothreonine)")
   #If we haven't seen the  peptide yet, make a dictionary for it !
   if (Fmatch > 0) :
	   F_psm_count += 1
   #print (psms.iloc[i,30])

print("# of PSMs with F")
print(F_psm_count)
print("# of PSMs total")
print(Total_psm_count)
