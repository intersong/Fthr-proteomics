This repo contains a series of scripts for extracting information from tabular output from proteomics database searches.  Specifically, it is set up to parse output from the Morpheus search tool.
1.  codon_puller.py
Input:
hard-coded genbank files
 XXXX_PSMS.tsv file
 Output:
 a table of codon usage for modified and unmodified positions.
 Notes:
 Potentially returns incorrect results for peptides that occur at multiple locations in the genome.  
 Many key variables are hard coded.

2.  DEseq_cleaner2.py
Input:
XXXX_PSMS.tsv file
Output:
PSMs summed at the protein level for unmodfified and modified peptides that meet specific criteria.
Suitable for input into DESeq, edgeR, or similar.
Not the ideal approach overall, maybe?

3.  number_of_peptides_with_F_vs_total.py
Input:
XXXX_PSMS.tsv file
Output:
Quality PSMs corresponding to unmodified and modified PEPTIDES.

4.  number_of_proteins_with_F_vs_total.py
Input:
XXXX_PSMS.tsv file
Output:
Quality PSMs corresponding to unmodified and modified peptides mapping to specific PROTEINS.

5.  for_getting_peptides_proteins_with_fthr.py
For printing out supplemenental tables with all modified, unmodified peptides.
