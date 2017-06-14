This repository contains a series of scripts for extracting information from tabular output from proteomics database searches.  
1.  spectrum_counter_MaxQUant.py
This script parses MaxQuant output, and prints a table with counts for FThr-countaining PSMs, total PSMs, FThr-countaining proteins, and total proteins.  It was used to generate table 1a.

2.  codon_puller_MaxQuant.py
This script parses MaxQuant output, and prints a table with the incidence of PSMs corresponding to single-Thr peptides encoded by each of the four Thr codons.  It was used to generate table S3.

3.  fthr_msstats.rmd
This script parses MaqQuant output, and prints a) histograms of estimated FThr incorporation for peptides with at least one FThr-containing PSM across all datasets, conditioned by strain (WT or p0564 KO),  b)  supplementary data tables, which contain measures or which peptides/proteins contain more/less FTHr than expected in the p0564 datasets, and c) assigns FThr-modified and total protein IDs to COG groups.

4.  Data files:
modificationSpecificPeptides.txt
