import pandas as pd
import os 
import subprocess

dedup_full = pd.read_csv('FOSC_100-contig_40-N50.csv')

identifier_dict = {}
names_dict = {}
passed = {}
i = 1
for row in dedup_full.iterrows(): #[dedup_full['Assembly Accession'] == 'GCA_001702695.2']['Organism Name']
    GCA_acc = row[1]['Assembly Accession']
    identifier = row[1]['Organism Name'].replace('Fusarium oxysporum', 'Fo').replace('f. sp. ','_').replace(' ', '')
    
    if f'{identifier}_{i}' in names_dict.values():
        passed[identifier] += 1
        i = passed[identifier]
        names_dict[GCA_acc] = f'{identifier}_{i}'
    else:
        i = 1
        names_dict[GCA_acc] = f'{identifier}_{i}'
        passed[f'{identifier}'] = 1

paths = {}
for genome_dir in os.listdir('ncbi_dataset/data/'):
    if genome_dir.startswith('GCA'):
        files = os.listdir(f'ncbi_dataset/data/{genome_dir}')
        f = [x for x in files if x.endswith('_genomic.fna') and x.startswith('GCA')][0]
        paths[genome_dir] = f'ncbi_dataset/data/{genome_dir}/{f}'


for genome, path in paths.items():
	name = f'{genome}_{names_dict[genome]}'
	cmd = f'funannotate sort -i {path} --minlen 500000 -o clean_data/{name}.clean.fa'
	subprocess.run(cmd.split(' '))
	#print(cmd)
