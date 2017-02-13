import pandas as pd
from glob import glob
import numpy as np


#files = glob("*.PSMs.tsv")
files = ["QEPlus2_02012016_11WxA.PSMs.tsv", "QEPlus2_02012016_15_Wxb.PSMs.tsv", "QEPlus2_12212015_62_JM_w2A.PSMs.tsv", "QEPlus2_02012016_17_Wzb.PSMs.tsv"]

csv_list = []
for f in files:
	csv_list.append(pd.read_csv(f, sep = "\t"))
	#raw = pd.read_csv(f, sep = "\t")

raw = pd.concat(csv_list, axis = 0)  #stack vertically

#bring me all your good hits,
psms = raw[(raw["Q-Value (%)"] < 1) & raw['Peptide Sequence'].str.contains("Fluorothreonine")  & raw['Protein Description'].str.startswith("embl-cds")] 
#				#good hits 				# just Fth containing 										#no decoys, decoys start with DECOY for morpheus

#make a new column to count the number of times each peptide sequence occurs across all PSMs being processed.
pepcount = pd.Series(np.zeros(psms.count()[0]), index = psms.index, name = "Number of PSMs")
psms = psms.join(pepcount)
print ("joined pepcount to psms")
reind_psms = psms.set_index(np.arange(psms.count()[0])) #somewhat incomrehensibly, you can have a non-ambiguous index and that happened in this case,so I reset the index.
#Just a dummy variable that will store whether this sequence has been seen before.
old_seq = []
#iterate thru the table of data.
for i in reind_psms.index:
	seq = reind_psms.loc[i,"Peptide Sequence"]
	#if we haven't seen this peptide sequence before....
	if seq not in old_seq:
		#set the value in the "Number of PSMs" column equal to the number of occurences of the sequence in the table of filtered data
		reind_psms.loc[i,"Number of PSMs"] = psms['Peptide Sequence'].value_counts()[seq]
	old_seq.append(seq)

#Sort the 
psms_sorted = reind_psms.sort_values(["Protein Description", "Peptide Sequence"])
#only peptides with atleast 2 hits over all datasetsse
psms_sorted_pared = psms_sorted[(psms_sorted["Number of PSMs"] > 1)]						
new_file = "WT_parsed.tsv"
psms_sorted_pared.to_csv(new_file, sep = "\t")
