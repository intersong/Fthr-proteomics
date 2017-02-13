###USAGE###
#python codon_puller.py 
 
###IMPORT MODULES ####
from Bio import SeqIO
import re
from pandas import read_csv
import sys

###READ DATA###
genome  = SeqIO.read('D:\jmcmurry\Documents\cattleya-complete-chromosome.gb','genbank') #you MUST tell SeqIO what format is being read. fix XXXXX
#genome_dict = SeqIO.to_dict(genome)
plasmid = SeqIO.read('D:\jmcmurry\Documents\cattleya-plasmid.gb','genbank')
#plasmid_dict = SeqIO.to_dict(plasmid)
peptides = read_csv(sys.argv[1], sep = "\t")
#peptides = read_csv("D:\jmcmurry\Documents\Davis_data\morph_withF\delp0564_parsed.tsv", sep = "\t")
peptides = peptides[(peptides["Q-Value (%)"] <  1/10) & peptides['Protein Description'].str.startswith("embl-cds")] #

###MAKE LIST OF FTH CONTAINING PEPTIDE####
Fth_list = []

for i in peptides.index:
	if "(Fluorothreonine)" not in peptides.loc[i,"Peptide Sequence"] and peptides.loc[i, "Base Peptide Sequence"].count("T") == 1: #x is row with modified sequence; whilst (phos.. is for morpheus
		Fth_list.append(peptides.loc[i, "Base Peptide Sequence"]) #y is row witth base seq.  " . " corresponds to cut site.  Splitting just getts the base seq #peptides["Base Peptide Sequence"].split(".")[1])

Fth_list = list(set(Fth_list))   #prevent redundancy
print(Fth_list)
print("peptides indexed")
###VALIDATE RETENTION TIMES###
#for pep in Fth_list:
#    peptides[peptides

#df.ix[df["Base Peptide Sequence"] == "PVVVLLTGVAALAAIAVPAASLEMGLSGDGSK"]["Retention Time (minutes)"]
#df["Peptide Sequence"].ix[38785].replace("(Fluorothreonine)", "T")


###MAKE LIST OF CDS IN GENBANK FILE###
CDS_list = []
for feature in genome.features:
   if feature.type == 'CDS':
	   CDS_list.append(feature.extract(genome.seq))
for feature in plasmid.features:
   if feature.type == 'CDS':
	   CDS_list.append(feature.extract(plasmid.seq))

print("CDSs indexed")
###MAKE A DICT WITH PEPTIDES AS KEYS AND CODONS AS  
codon_dict = {}
for peptide in Fth_list:
	codon_dict[peptide] = []
	
for CDS in CDS_list:
   CDS_aa = str(CDS.translate())
   CDS_na = str(CDS)
   for peptide in Fth_list:
	   T_matches = re.finditer("T", peptide)
	   peptide_match = re.search(peptide, CDS_aa)#should use finditer or else a cds with the same peptide twice could get confused?
	   if peptide_match :  			#I mean to say if peptide_match != null
		   for Thr in T_matches:
			   codon_start = (peptide_match.start() + Thr.start()) * 3   # take the start posnof peptide in AA and add to start of T in AA.  Then X 3.
			   codon_dict[peptide].append(CDS_na[codon_start : codon_start + 3]) # pull the codon
			   print(peptide + "        " + "            " + str(codon_start) + "         " + CDS_na[codon_start : codon_start + 3])
			   
###GO BACK AND PRINT CODON USAGE ###
print("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")
for key in codon_dict.keys():
   codons = set(codon_dict[key])
   if len(codons) == 1:
      stupid_list = []
      stupid_list.append(key)         
      stupid_list.append("\t" + list(codons)[0])
      print (' '.join(stupid_list))

#         peptide_na = CDS.sequence()[peptide_match.start() * 3 : (peptide_match.end() + 1) *3 - 1]  #pull the nucleic acid sequence of peptide
