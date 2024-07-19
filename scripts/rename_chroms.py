"""
python rename_chroms.py meta_data.csv genome/dir
"""
import pandas as pd
import os 
import subprocess
from sys import argv

dedup_full = pd.read_csv(argv[1])
data_dir = argv[2]


identifier_dict = {}
names_dict = {}
passed = {}
i = 1
for row in dedup_full.iterrows(): #[dedup_full['Assembly Accession'] == 'GCA_001702695.2']['Organism Name']
    print(row[1]['Organism Name'])
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
#for genome_loc in os.listdir(data_dir):
#	files = [os.listdir(
#	paths[genome_loc] = f'{data_dir}{genome_loc}'

for genome in list(dedup_full['Assembly Accession']):
	genome_loc = [f for f in os.listdir(data_dir) if f.startswith(genome)]
	paths[genome] = f'{data_dir}{genome_loc[0]}'

print(os.listdir(data_dir))

for genome, path in paths.items():
	name = f'{genome}_{names_dict[genome]}'
	cmd = f'funannotate sort -i {path} --minlen 500000 -o clean_data/{name}.clean.fa'
	print(f'running {genome}')
	subprocess.run(cmd.split(' '))
	#print(cmd)
