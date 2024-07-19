"""
Script to split chromosomes falling into multiple communities

python split_chroms.py [community_directory] [genome_directory]

python scripts/split_chroms.py communitites/ 95-2500-50000-0.01/ genomes/new_set/clean_data/ GCA_016801315.1_Fo_capsici_1
"""
import sys
sys.path.append('/home/anouk/anouk2/graph_Fosc2.0/scripts/') 
import pandas as pd
import analyse_communities as a_comm
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sys import argv
import os
import subprocess
import json

def make_bed_paf(indir, scaff_comm_dict, genome):
    """
    From out.paf -> bed file including query and targets
    """
    collection = []
    query = False
    target = False
    with open(f'{indir}/out.paf') as paffile:
        for line in paffile:
            line = line.strip().split('\t')
                
            if genome != 'all':
                #get target and query
                if genome in line[0]:
                    query = True
                    target = False
                elif genome in line[5]:
                    target = True
                    query = False
                else:
                    query = False
                    target = False

                if query:
                    chrom = line[0]
                    start = line[2]
                    stop = line[3]
                    full_len_q = str(int(line[1]) - (int(line[3]) - int(line[2])))
                    full_len_t = str(int(line[6]) - (int(line[8]) - int(line[7])))
                    pid = line[-1].split(':')[-1]
                    try:
                        comm = scaff_comm_dict[line[5]]
                    except KeyError:
                        comm = '-1'
                    info = '-'.join([line[5], line[7], line[8], line[6], full_len_q, full_len_t,pid,comm])
                elif target:
                    chrom = line[5]
                    start = line[7]
                    stop = line[8]
                    full_len_q = str(int(line[6]) - (int(line[8]) - int(line[7])))
                    full_len_t = str(int(line[1]) - (int(line[3]) - int(line[2])))
                    pid = line[-1].split(':')[-1]
                    try:
                        comm = scaff_comm_dict[line[0]]
                    except KeyError:
                        comm = '-1'
                    info = '-'.join([line[0], line[2], line[3], line[1], full_len_q, full_len_t, pid,comm])

                if target or query:
                    collection.append([chrom, start, stop, info])
    return collection

def write_out(to_write, output):
    with open(output, 'w+') as of:
        for line in to_write:
            of.write('\t'.join(line))
            of.write('\n')
    return

def get_windows(genome_dir, genome, window_size, outdir):
    outfile = f'{outdir}/window_bed/w{window_size}_{genome}.bed'
    
    if os.path.exists(outfile):
        print('finished intersect windows')
    else:
        f = open(outfile, "w")
        cmd = ['bedtools', 'makewindows', '-g', \
        f'{genome_dir}/{genome}.clean.fa.fai', '-w', str(window_size)]#, '>', outfile]
        
        subprocess.call(cmd, stdout = f)
        
        f.close()
    
        print(' '.join(cmd))

def intersect_windows(outdir, genome, window_size):
    outfile = f'{outdir}/window_bed/intersect_windows_{genome}.bed'
    
    if os.path.exists(outfile):
        print('finished intersect windows')
    else:
        print('intersect windows')
        
        f = open(outfile, "w")
        cmd = ['bedtools', 'intersect' , '-a', \
        f'{outdir}/window_bed/w{window_size}_{genome}.bed', '-b', \
        f'{outdir}/paf_bed/{genome}_paf.bed', '-wao']# '>', outfile]
        
        subprocess.call(cmd, stdout = f)
        f.close()
    
        print(' '.join(cmd))

def parse_intersect(intersect_bed):
    """
    intersect_bed = path_to_intersect_bed
    
    intersect gives a line per window - mapping. Multiple lines with 
    info on the same window. Merge these in a dict to later get perc. mappings
    
    Key = window
    value = mappings for that window
    """
    paf_dict = {}
    num_chromosomes = set()
    with open(intersect_bed) as infile:
        for line in infile:
            chrom, start, stop, match1, start1, stop1, info1, len_match = line.strip().split('\t')
            num_chromosomes.add(chrom)
            identifier = f'{chrom}-{start}-{stop}'
            try:
                paf_dict[identifier].append(info1)
            except KeyError:
                paf_dict[identifier] = [info1]
    return paf_dict, len(num_chromosomes)

def perc_comm_window(paf_dict, num_chrom):
    """
    Per window get percentage of mappings per community.
    Input: genomic windows and intersect output
    Output: {scaffold1: {window_start1: {comm1 : count_mappings, comm2: },
                        {window_starts2: {comm1 : count_mappings, comm2: }}
    """
    indices = []
    all_dict = {}

    for i in range(1,num_chrom):
        per_comm_dict = {}
        comms = []
        if len(indices) > 0:
            for indi in indices:
                per_comm_dict[indi] = {-1:0}
                    
        for k,v in paf_dict.items():
            if k.startswith(f'{genome}#scaffold{i}-'):
                ind = int(k.split('#')[1].split('-')[1])
                
                if i == 1:
                    indices.append(ind)
                    per_comm_dict[ind] = {}          

                comm_u, counts = np.unique(([x.split('-')[-1] for x in v]), return_counts = True)
                counts = counts/sum(counts)

                for ic, c in enumerate(comm_u):
                    per_comm_dict[ind][c] = counts[ic] 
                
                [comms.extend(x.split('-')[-1]) for x in v]
        comm_u, counts =  np.unique(comms, return_counts = True)
        
        all_dict[f'scaffold{i}'] = per_comm_dict
        
    return all_dict

def plot_perc_windows(all_p_comm_dict):
    fig, ax = plt.subplots(nrows = len(all_p_comm_dict.keys()), figsize = (20, 50), sharex = True)
    i = 0
    for k,v in all_p_comm_dict.items():
        pd.DataFrame(v).T.plot.bar(stacked = True, ax = ax[i])
        ax[i].set_title(f'{genome}#scaffold{i}')
        i += 1
    plt.savefig(f'{genome}_perc_chrom.pdf')

def collect_stats(comm_dict, majority):
    #stats to collect
    majority_expected = 0
    majority_unexpected = 0
    full_support = 0
    inbetween = 0
    low_support = 0
    candidate_comms = {}
    other_windows = {}
    
    first, last = get_first_last(comm_dict)
    #print(str(first))
    #print(str(last))
    
    for window, comms_window in comm_dict.items():
            window_comm = max(comms_window, key=comms_window.get)
            value = max(comms_window.values())
            if window_comm == majority and value >= 0.9:
                majority_expected += 1
                full_support += 1
            elif window_comm == majority and (value >= 0.5 and value < 0.9):
                majority_expected += 1
                inbetween += 1
            elif window_comm == majority and (value < 0.5):
                majority_expected += 1
                low_support += 1    
            elif value == 0:
                'no_mappings'
            #if window is not majority and not unknown, count as unexpected
            elif window_comm != majority and window_comm != -1 and value >= 0.4:
                majority_unexpected += 1
                
                #list alternative comms
                try:
                    candidate_comms[window_comm].append(value)
                except KeyError:
                    candidate_comms[window_comm] = [value]
                
                if str(window) == str(first):
                    print(window)
                    window = f'*{window}'
                if str(window) == str(last):
                    window = f'{window}*'
                other_windows[window] = window_comm
                    
            #current window to previous window before new itteration
            prev_window = window
    
    temp_dir = {'expected': majority_expected, 'unexpected': majority_unexpected,\
                'full support': full_support, 'inbetween' : inbetween, 'low_support': low_support}#, \
                #'candidates':candidate_comms}
    
    return temp_dir, other_windows

def collect_split_candidates(all_dict, scaff_comm_dict, genome, outdir):
    splits = {}
    
    no_majority = False

    scaff_stats_comms = {}

    for scaff, comm_dict in all_dict.items():
        
        switch_window = -1
        switch = False
        
        try:
            majority = scaff_comm_dict[f'{genome}#{scaff}']
            no_majority = False
        except KeyError:
            print(f'{genome}#{scaff}', '?')
            no_majority = True
        
        if not no_majority:
            stats_dict, other_windows = collect_stats(comm_dict, majority)
                    
            scaff_stats_comms[scaff] = stats_dict
            
            if len(other_windows) > 0:
                splits[scaff] = other_windows
        else:
            continue
    write_scaff_splits(f'{outdir}/scaff_split/{genome}.txt', splits)

    return scaff_stats_comms

def write_scaff_splits(scaf_split_out, data):
    with open(scaf_split_out, 'w+') as outfile:
        json.dump(data, outfile)
    
def get_first_last(perc_comm_dict):
    first = min([k for k,v in perc_comm_dict.items() if any(int(str(x).replace('.', '-1')) >= 0 for x in list(v.keys()))])
    last = max([k for k,v in perc_comm_dict.items() if any(int(str(x).replace('.', '-1')) >= 0 for x in list(v.keys()))])
    altf = list(perc_comm_dict.keys())[0]
    altl = list(perc_comm_dict.keys())[-1]
    #print(first, last, altf, altl)
    return first, last


def wrapper(indir, parameters, comm_dir, genome_dir, outdir, genome):
    community_dict = a_comm.get_comm_dict(f'{indir}/flagN/{parameters}')
    community_dict_split = a_comm.get_comm_dict(comm_dir)
    window_size = 100000
    
    #make bed from paf
    paf_lines = make_bed_paf(comm_dir, community_dict, genome)
    write_out(paf_lines, f'{outdir}/paf_bed/{genome}_paf.bed')
    
    get_windows(genome_dir, genome, window_size, f'{outdir}')
    #intersect gives a line per window-mapping pair
    intersect_windows(outdir, genome, window_size)
    
    paf_dict, num_chrom = parse_intersect(f'{outdir}/window_bed/intersect_windows_{genome}.bed')
    pcwd = perc_comm_window(paf_dict, num_chrom)
    
    get_first_last(pcwd['scaffold1'])
    
    plot_perc_windows(pcwd)
    
    scaff_stats = collect_split_candidates(pcwd, community_dict, genome, outdir)
    
    pd.DataFrame(scaff_stats).to_csv(f'{outdir}/scaff_split/{genome}_scaffstats.csv')
    return


if __name__ == '__main__':
    indir = argv[1]
    parameters = argv[2]
    comm_dir = f'{indir}/{parameters}'
    genome_dir = argv[3]
    genome = argv[4]
    outdir = 'chrom_overview_2'
    
    #wrapper(indir, parameters, comm_dir, genome_dir, outdir, genome)
    
    for genome in os.listdir(genome_dir):
        if genome.endswith('.fa'):
            genome = genome.replace('.clean.fa', '')
            wrapper(indir, parameters, comm_dir, genome_dir, outdir, genome)

