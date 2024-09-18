# FOSC_pangenome

Repository for building and analysing a pangenome of species containing
accessory genomic regions.

## install
1. clone repository

2. create environment
conda env create -f environment.yml
conda activate mpi

## steps
1. collect and clean genomes
headers as ">genome#chromosome"
2. Concatenate all genomes into one file
3. Add location of all genomes fasta in communities/config.yaml
3. Run snakefile in communities directory

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
