## This repository contains a series of scripts for extracting information from tabular output from proteomics database searches.  

### 1.  spectrum_counter_MaxQUant.py

This python script accepts as input the MaxQuant output file "msms.txt."  It prints a table with counts for FThr-countaining PSMs, total PSMs, FThr-countaining proteins, and total proteins.  It was used to generate table 1a.  The path to the input file is hard-coded, so you will have to change it if you wish to reproduce.

### 2.  codon_puller_MaxQuant.py

This python script accepts as input the MaxQuant output file and the genbank file with genome sequence and CDSs of *S. cattleya*, and prints a table with the incidence of PSMs corresponding to single-Thr peptides encoded by each of the four Thr codons.  It was used to generate table S3.   The path to the input files are hard-coded, so you will have to change them if you wish to reproduce.

### 3.  fthr_msstats.rmd

This R markdown script parses MaqQuant output file "modificationSpecificPeptides.txt", as well as files for COG assignments and gene name mappings. It prints the following:

-a) histograms of estimated FThr incorporation for peptides with at least one FThr-containing PSM across all datasets, conditioned by strain (WT or p0564 KO)
-b) supplementary data tables, which contain measures or which peptides/proteins contain more/less FTHr than expected in the p0564 datasets
-c) COG group assignments of FThr-modified proteins and all proteins.

Again, all file paths are hard-coded.

### 4.  Data files:

The raw data used as input for the scripts listed above.  The file "msms_truncated.zip" is a truncated version of "msms.txt" that is small enough to upload on github.  It contains all the columns that were used in the above scripts.

modificationSpecificPeptides.txt
-msms_truncated.zip
-cattleya-cchromandplas.gb
-genome_COGs.txt
-plasmid_COGs.txt
-gene-name-comparisons.txt
