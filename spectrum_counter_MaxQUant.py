###USAGE###
#python codon_puller_MSstats.py D:\jmcmurry\Documents\Davis_data\combined\txt\MSMS.txt
###IMPORT MODULES ####
from Bio import SeqIO
import re
import pandas as pd

###READ DATA###
peptides = pd.read_csv('D:\\jmcmurry\\Documents\\Davis_data\\combined\\txt\\msms.txt', sep = "\t")
#Hi quality only
peptides = peptides[(peptides["PEP"] <  .01)]
#No hits to decoy/reversed
peptides = peptides[peptides["Reverse"] != "+"]
#Just a Boolean TRUE if it is a proteotyptic peptide (only one protein with this peptide)
#peptides = peptides[list(map(lambda x: len(x.split(";")) == 1, peptides["Proteins"]))]
#Somewhat cryptic, but just a Boolean TRUE if 1 T in seq, FALSE otherwise
#peptides = peptides[list(map(lambda x: x.count("T") == 1, peptides["Sequence"]))]
#Eliminate hits to contaminant database.
peptides = peptides[list(map(lambda x: str(x).startswith("embl-cds"), peptides['Proteins']))]

samples = set(peptides["Raw file"])

sampdict = {}
for s in samples :
    table = peptides[peptides["Raw file"] == s]
    #How many rows have a fluorothreonine?
    sampdict[s] = {}
    sampdict[s]["Fthr"] = sum(table["Fluorothreonine"] > 0)
    sampdict[s]["Fproteins"] = len(set(table[(table["Fluorothreonine"] > 0)]["Proteins"]))
    sampdict[s]["allproteins"] = len(set(table["Proteins"]))
    #How many rows are there, period?
    sampdict[s]["All"]  = len(table["Sequence"])

print(sampdict)
for s in sampdict.keys():
    print(sampdict[s])
