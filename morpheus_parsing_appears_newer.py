import pandas as pd
from glob import glob
import numpy as np

#files = glob("*.PSMs.tsv")
files = ["QEPlus2_02012016_11WxA.PSMs.tsv", "QEPlus2_02012016_15_Wxb.PSMs.tsv", "QEPlus2_12212015_62_JM_w2A.PSMs.tsv", "QEPlus2_02012016_17_Wzb.PSMs.tsv"]
#files = ["QEPlus2_12212015_60_JM_96A.PSMs.tsv", "QEPlus2_02012016_9_96B.PSMs.tsv","QEPlus2_02012016_19_97A.PSMs.tsv", "QEPlus2_02012016_13_97b.PSMs.tsv"]

csv_list = []
for f in files:
	pre_df = pd.read_csv(f, sep = "\t")
	#psmcount = pd.Series(np.zeros(pre_df.count()[0]), index = pre_df.index, name = "Individual File PSMs")
	#psms = pre_df.join(psmcount)
	#old_seq = []
	#for i in pre_df.index:
	   #seq = pre_df.loc[i,"Peptide Sequence"]
	   #if seq not in old_seq:
	      #pre_df.loc[i,"Individual File PSMs] = pre_df['Peptide Sequence'].value_counts()[seq]
	csv_list.append(pre_df)

raw = pd.concat(csv_list, axis = 0)  #stack vertically

psms = raw[(raw["Q-Value (%)"] <  1/10) & raw['Peptide Sequence'].str.contains("Fluorothreonine")  & raw['Protein Description'].str.startswith("embl-cds")] #bring me all your good hits,
#				#good hits 				# just Fth containing 										#no decoys, decoys start with DECOY for morpheus
###ADD BLANK COLUMNS WITH ZEROS###
pepcount = pd.Series(np.zeros(len(psms.index)), index = psms.index, name = "Number of PSMs")
protcount = pd.Series(np.zeros(len(psms.index)), index = psms.index, name = "Protein Number of PSMs")
psms = pd.concat([psms, pepcount, protcount], axis = 1)
print ("joined pepcount to psms")

###De-redundify and set values####
reind_psms = psms.set_index(np.arange(psms.count()[0])) #somewhat incomrehensibly, you can have a non-ambiguous index and that happened in this case,so I reset the index.

old_seq = []
for i in reind_psms.index:
	seq = reind_psms.loc[i,"Peptide Sequence"]
	prot_des = reind_psms.loc[i,"Protein Description"]
	if seq not in old_seq:
		reind_psms.loc[i,"Number of PSMs"] = reind_psms['Peptide Sequence'].value_counts()[seq]
		reind_psms.loc[i,"Protein Number of PSMs"] = reind_psms['Protein Description'].value_counts()[prot_des]
	old_seq.append(seq)

###refilter to remove PSMS with fewere than 2x hits ###
#reind_psms = reind_psms[(reind_psms["Number of PSMs"] > 1)]		#only peptides with atleast 2 hits over all datasets-more trouble than it's worth.
#old_seq = []
#for i in reind_psms.index:
#	seq = reind_psms.loc[i,"Peptide Sequence"]
#	prot_des = reind_psms.loc[i,"Protein Description"]
#	if seq not in old_seq:
#		reind_psms.loc[i,"Protein Number of PSMs"] = ((reind_psms['Protein Description'] == prot_des) * reind_psms['Number of PSMs']).sum()
#	old_seq.append(seq)

###SOME ODD-JOB PROCESSING###
reind_psms = reind_psms[(reind_psms["Number of PSMs"] > 0)]
psms_sorted = reind_psms.sort_values(["Protein Number of PSMs", "Protein Description","Number of PSMs", "Peptide Sequence"], ascending = [False, True, False, True])

###PRINT THE FILE###
new_file = "WT_parsed.tsv"
#new_file = "delp0564_parsed.tsv"
psms_sorted.to_csv(new_file, sep = "\t")


###print protein-oriented file####

for p in psms_sorted.index:
	prot_des = psms_sorted.loc[p,"Protein Description"]
	pindex = (psms_sorted['Protein Description'] == prot_des).index[0]
	print( pindex)
	#psms_sorted.loc[pindex, "Precursor Charge"] = 999999

print( psms_sorted)



###TODO###
#retention time validation
#sort
