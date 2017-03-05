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
peptides = peptides[list(map(lambda x: len(x.split(";")) == 1, peptides["Proteins"]))]
#Somewhat cryptic, but just a Boolean TRUE if 1 T in seq, FALSE otherwise
peptides = peptides[list(map(lambda x: x.count("T") == 1, peptides["Sequence"]))]
#Eliminate hits to contaminant database.
peptides = peptides[list(map(lambda x: str(x).startswith("embl-cds"), peptides['Proteins']))]
#Only KO stuff:
peptides = peptides[list(map(lambda x: x[-3] == "9", peptides['Raw file']))]
#With list comprehensions instead of map!
#        = peptides[x[-3] == "9" for x in peptides["Raw file"]]

print("DOne filtering the MS  data!!!")
###MAKE LIST OF FTH CONTAINING PEPTIDE####
def get_spec_counts(table) :
    ppdict   = {}
    pset = set(table['Sequence'])
    for i in pset :
        subtab = peptides[peptides['Sequence'] == i]
        prot   = list(subtab['Proteins'])[0].split(":")[1]
        counts = len(table[table["Sequence"] == i])
        if not prot in ppdict :
            ppdict[prot] = {}
        if not i in ppdict[prot] :
            ppdict[prot][i] = counts
    return(ppdict)


Upeps = peptides[peptides["Modifications"] == "Unmodified"]
Fpeps = peptides[peptides["Modifications"] == "Fluorothreonine"]
#Only peps where location of Fthr is known
#Fpeps = Fpeps[Fpeps['Localization prob'] > 0.95]
#Fseq  = list(map(lambda x: x[1:-1].replace('T(fl)', 'X'), Fpeps['Modified sequence']))

print("Got the fluoro/threonine peptides in tables...")
Fth_dict = get_spec_counts(Fpeps)
Unm_dict = get_spec_counts(Upeps)
print("Got the peptides/counts in dictionaries")
###MAKE LIST OF CDS IN GENBANK FILE###
genplas   = SeqIO.parse('D:\jmcmurry\Documents\cattleya-cchromandplas.gb','genbank') #you MUST tell SeqIO what format is being read. fix XXXXX
CDS_dict = {}
for r in genplas :
    for f in r.features :
        if f.type == 'CDS' :
            if 'protein_id' in f.qualifiers:
                embl_id   = f.qualifiers['protein_id'][0].split(".")[0]
            else:
                embl_id   = f.qualifiers['db_xref'][0].split(":")[1].split(".")[0]
            naseq     = f.extract(r.seq)
            CDS_dict  [embl_id] = {
            "CDS_na" : str(naseq),
            "CDS_aa" : str(naseq.translate())
            }
print("CDSs indexed")
###MAKE A DICT WITH PEPTIDES AS KEYS AND CODONS AS values
def get_cod(pep_dict, CDS_dict) :
    codon_dict = {"ACA": 0,
                  "ACG": 0,
                  "ACC": 0,
                  "ACT": 0}
    counter = 0
    for CDS in CDS_dict.keys():
        if CDS in pep_dict :
            CDS_aa            = CDS_dict[CDS]['CDS_aa']
            CDS_na            = CDS_dict[CDS]['CDS_na']
            for peptide in  pep_dict[CDS].keys() :
                T_matches     = re.finditer("T", peptide)
                peptide_match = re.search(peptide, CDS_aa)
                if peptide_match :  			                                         #I mean to say if peptide_match != null
                    for Thr in T_matches:
                        codon_start        = (peptide_match.start() + Thr.start()) * 3   # take the start posnof peptide in AA and add to start of T in AA.  Then X 3.
                        codon              = CDS_na[codon_start : codon_start + 3]
                        codon_dict[codon]  += pep_dict[CDS][peptide]                                           # Increment the counter
                        if codon == "ACA" :
                            print(peptide)
                            print(pep_dict[CDS][peptide])
                            print(codon + "\n\n\n")
        #print("Done with CDS " + CDS)
    return(codon_dict)

###GO BACK AND PRINT CODON USAGE ###
#AACK !  Right now it adds one count per peptide.  I should probably have it add one count/PSM.
Fcods = get_cod(pep_dict = Fth_dict, CDS_dict = CDS_dict)
Ucods = get_cod(pep_dict = Unm_dict, CDS_dict = CDS_dict)
print("Fthr codons: \n")
print(Fcods)
print("Unmod codons: \n")
print(Ucods)
