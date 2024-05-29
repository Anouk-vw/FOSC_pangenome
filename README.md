# FOSC_pangenome

Repository for building and analysing a pangenome of species containing
accessory genomic regions.

## steps
1. collect and clean genomes
#### genomes/
- genomes/all_meta/clean_data
contain all cleaned genomes
- genomes/all_meta/rename_chroms.py 
script to clean and rename assemblies
Specifically adjusted for FOSC genomes from NCBI adding the formae specialis to the GCA accession number

2. try different settings
#### communities/
- communitites/Snakefile
Pipeline to try different combinations of parameters.
Creates an output directory per combination.
