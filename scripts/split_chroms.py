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
    """
    run bedtwools makewindows 
    """
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
    """
    Intersect windows and paf bed to get comms per window
    """
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
            #make window identifier
            #append all different info lines (communities) to identifier
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
    
    #per chromosome, per window, determine percentage communities
    for i in range(1,num_chrom):
        per_comm_dict = {}
        comms = []
        
        #if this is NOT the first chromosome (longest), set indices.
        #this is used to plot the chromosomes in relative sizes
        if len(indices) > 0:
            for indi in indices:
                #fill indices with community -1, count 0
                per_comm_dict[indi] = {-1:0}

        #for elements in paf file
        for k,info in paf_dict.items():
            #if chromosome in paf-line
            if k.startswith(f'{genome}#scaffold{i}-'):
                #split k to get starting index of window
                ind = int(k.split('#')[1].split('-')[1])
                
                #if this is the first chromosome (longest, get indices)
                if i == 1:
                    indices.append(ind)
                    per_comm_dict[ind] = {}          

                #split info on '-' and get community (final element)
                comm_u, counts = np.unique(([x.split('-')[-1] for x in info]), return_counts = True)
                #divide counts -> frequency
                counts = counts/sum(counts)
                
                #for the index of the current window (ind)
                #add the unique communites (comm_u) to dictionary
                #{window: {community : frequency}}
                for ic, c in enumerate(comm_u):
                    per_comm_dict[ind][c] = counts[ic] 
                
                [comms.extend(x.split('-')[-1]) for x in info]
                
        comm_u, counts =  np.unique(comms, return_counts = True)
        
        #make a complete dictionary
        #{scaffold: {window: {community : frequency}}}
        all_dict[f'scaffold{i}'] = per_comm_dict
        
    return all_dict

def plot_perc_windows(all_p_comm_dict):
    """
    Plot frequency of community per window per chromosome
    """
    fig, ax = plt.subplots(nrows = len(all_p_comm_dict.keys()), figsize = (20, 50), sharex = True)
    i = 0
    for k,v in all_p_comm_dict.items():
        pd.DataFrame(v).T.plot.bar(stacked = True, ax = ax[i])
        ax[i].set_title(f'{genome}#scaffold{i}')
        i += 1
    plt.savefig(f'{genome}_perc_chrom.pdf')

def collect_stats(comm_dict, majority):
    """
    Collect statistics. How many need splitting?
    """
    #stats to collect
    majority_expected = 0
    majority_unexpected = 0
    full_support = 0
    inbetween = 0
    low_support = 0
    candidate_comms = {}
    other_windows = {}
    
    #get first/last to determine thsi is a flank
    first, last = get_first_last(comm_dict)

    
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
    """
    Get chromosomes that need to be split
    scaff_split/{genome}.txt: dict of split candidates {scaff:{window:comm}}
    """
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
    """
    determine if the window is the firs/last window in chromosome
    Necessary to split only flanks.
    
    Take into account that preceding windows are allowed to be unkown
    """
    first = min([k for k,v in perc_comm_dict.items() if any(int(str(x).replace('.', '-1')) >= 0 for x in list(v.keys()))])
    last = max([k for k,v in perc_comm_dict.items() if any(int(str(x).replace('.', '-1')) >= 0 for x in list(v.keys()))])
    altf = list(perc_comm_dict.keys())[0]
    altl = list(perc_comm_dict.keys())[-1]
    #print(first, last, altf, altl)
    return first, last


def wrapper(indir, parameters, comm_dir, genome_dir, outdir, genome):
    #read communities in directory
    community_dict = a_comm.get_comm_dict(f'{indir}/flagN/{parameters}')
    community_dict_split = a_comm.get_comm_dict(comm_dir)
    window_size = 100000
    
    #make bed from paf
    paf_lines = make_bed_paf(comm_dir, community_dict, genome)
    write_out(paf_lines, f'{outdir}/paf_bed/{genome}_paf.bed')
    
    #get windows
    get_windows(genome_dir, genome, window_size, f'{outdir}')
    #intersect gives a line per window-mapping pair
    intersect_windows(outdir, genome, window_size)
    
    #make dict {scaffold:{window:{comm:frequency}}}
    paf_dict, num_chrom = parse_intersect(f'{outdir}/window_bed/intersect_windows_{genome}.bed')
    pcwd = perc_comm_window(paf_dict, num_chrom)
    
    #make plots
    plot_perc_windows(pcwd)
    
    #gather stats and write dict to out
    scaff_stats = collect_split_candidates(pcwd, community_dict, genome, outdir)
    
    pd.DataFrame(scaff_stats).to_csv(f'{outdir}/scaff_split/{genome}_scaffstats.csv')
    return

"""
filter splits
"""

def find_block(d, start_key, direction, key_mapping):
    """
    find blocks in window dict {window:comm}
    distance should be < 100000 basepairs (or other window value)
    """
    block = []
    comm = []
    start_value = d[start_key]
    #fix '*' indicating start-stop
    start_int_key = int(start_key.strip('*'))
    block.append(start_key)
        
    # Iterate through the sorted keys starting from start_int_key
    current_key = start_int_key
    while True:
        next_key = current_key + (100000 * direction)

        if next_key in key_mapping and d[key_mapping[next_key]] == start_value:
            comm.append(d[key_mapping[next_key]])
            block.append(key_mapping[next_key])
            current_key = next_key
        else:
            break
        
    return block, comm

def find_blocks_from_specific_keys(d):
    """
    for directory {window:comm, window:comm}
    determine blocks beloning to a certain comm
    """
    
    # Define the starting keys
    start_keys = [x for x in d.keys() if '*' in x]
    
    # Convert the keys to integers and create a mapping to original keys
    key_mapping = {int(key.strip('*')): key for key in d.keys()}
    keys = sorted(key_mapping.keys())
    
    # Initialize a list to hold the result
    result = []
    comms = []
    
    # Collect blocks starting from specified keys
    for start_key in start_keys:
        if start_key in d:
            if start_key.startswith('*'):
                blocks,comm = find_block(d, start_key, 1, key_mapping)
                result.append(blocks)  # Forward direction
                comms.append(comm)
            elif start_key.endswith('*'):
                blocks,comm = find_block(d, start_key, -1, key_mapping)
                result.append(blocks) # Backward direction
                comms.append(comm)
                
    return result, comms

def filter_putative_splits(indir):
    """
    '../chrom_overview_2/scaff_split/'
    """
    total_split = []

    for sfile in os.listdir(indir):
        if sfile.endswith('.txt'):
            with open(f'{indir}/{sfile}') as infile:
                ww_dict = json.load(infile)
                    
                for scaff, ww in ww_dict.items():
                    #get consquetive blocks from window dictionary
                    blocks, comms = find_blocks_from_specific_keys(ww)
                    #no blocks found
                    if len(blocks) < 1:
                        splittable = False

                    #blocks found
                    elif len(blocks) >= 1:
                        #if any block more than 10 windows (1 mb in size)
                        splittable = any(len(b) >= 10 for b in blocks)
                        
                        for indblock, b in enumerate(blocks):
                            if len(b) >= 10:
                                blocks = b
                                comms = comms[indblock]
                                comms = comms[0]

                        if splittable:
                            #append into to list 
                            b_temp = [int(x.strip('*')) for x in b]
                            for x in b:
                                if x.startswith('*'):
                                    total_split.append((sfile,scaff, x, str(max(b_temp)), comms))
                                elif x.endswith('*'):
                                    total_split.append((sfile,scaff, str(min(b_temp)), x, comms))
    return total_split

def open_windows(wfile):
    stops = {}
    with open(wfile) as infile:
        for line in infile:
            chrom = line.split('\t')[0].split('#')[1]
            stop = int(line.split('\t')[1])
            stops[chrom] = stop
    return stops

def get_split(paf_in, chrom, split_comm, window):
    """
    determine exact split location from PAF file with splitting enabled
    """
    possible_starts = []
    too_far = []
    prev = 0
    with open(paf_in) as paf:
        for line in paf:
            c,start, stop, info = line.split('\t')
            if c == chrom:
                comm = info.strip().split('-')[-1]
                if comm == split_comm:
                    if abs(int(window) - int(start)) < 100000:
                        possible_starts.append(int(start))
                    else:
                        too_far.append(int(start))
    
    if len(possible_starts) > 0 :
        return min(possible_starts)
    elif len(too_far) > 0:
        #print(min(too_far, key=lambda x:abs(x-window)))
        return min(too_far, key=lambda x:abs(x-window))
    else:
        return 'x'

def splits_to_bed(splits_dir, genome_dir, paf_dir, indir, parameters):
    
    """
    ../genomes/new_set/clean_data
    """
    possible_splits = filter_putative_splits(splits_dir)
    bed_lines = []
    adjustments = {}
    scaff_comm_dict = a_comm.get_comm_dict(f'{indir}/flagN/{parameters}')
    
    for p_split in possible_splits:
        isstart = False
        isstop = False
        
        sfile, chrom, start, stop, new_comm = p_split

        chrom_print = f"{sfile.replace('.txt', '#')}{chrom}"
        paf_file = f"{paf_dir}/{sfile.replace('.txt','')}_paf.bed_sorted"
        stops = open_windows(f"{genome_dir}/{sfile.strip('.txt')}.clean.fa.fai")
        
        current_comm = scaff_comm_dict[chrom_print]
        
        if current_comm not in adjustments.keys():
            adjustments[current_comm] = []
        if new_comm not in adjustments.keys():
            adjustments[new_comm] = []
        
        #determine which side of chrom should be split
        if '*' in start:
            isstart = True
        elif '*' in stop:
            isstop = True

        start = start.strip('*')
        stop = stop.strip('*')
        
        #if 'start' of chromosome should be split
        if isstart and not isstop:
            #determine exact split
            split = get_split(paf_file, chrom_print, str(new_comm), int(stop))
            
            #start must be 0, end is split, write to new_comm
            bed_lines.append([f'{chrom_print}', 0, split, f'{chrom_print}_p1'])
            adjustments[new_comm].append(f'{chrom_print}_p1')
            
            #start is split and end must be stops
            bed_lines.append([chrom_print, split, stops[chrom], f'{chrom_print}_p2'])
            adjustments[current_comm].append(f'{chrom_print}_p2')
           
        elif isstop and not isstart:
            split = get_split(paf_file, chrom_print, str(new_comm), int(start))
            
            #start is split and end must be stops
            bed_lines.append([chrom_print, split, stops[chrom], f'{chrom_print}_p2'])
            adjustments[new_comm].append(f'{chrom_print}_p2')
            
            #start must be 0, end is split
            bed_lines.append([chrom_print, 0, split, f'{chrom_print}_p1'])
            adjustments[current_comm].append(f'{chrom_print}_p1')
            
        elif isstart and isstop:
            print('both')
        
    return bed_lines, adjustments

def write_comm_files(adjustments):
    communities = [x.split('.')[6] for x in os.listdir(f'{indir}/flagN/{parameters}/') if x.startswith('out.paf.edges.weights.txt.community')]
    
    for community in communities:
        try:
            adjusted = adjustments[community]
        except KeyError:
            adjusted = []
        #for community, adjusted in adjustments.items():
        identifiers = [a.split('_p')[0] for a in adjusted]
        
        new_comm_file = f'updated_comms/{community}.txt'
        
        members = []
        with open(f'{indir}/flagN/{parameters}/out.paf.edges.weights.txt.community.{community}.txt') as comm_file:
            for line in comm_file:
                member = line.strip()
                if member in identifiers:
                    index_element = identifiers.index(member)
                    #append adjusted element to list
                    members.append(adjusted[index_element])
                else:
                    members.append(member)
        
        with open(new_comm_file, 'w+') as ncf:
            for m in members:
                ncf.write(m)
                ncf.write('\n')

    return
    

if __name__ == '__main__':
    indir = argv[1]
    parameters = argv[2]
    comm_dir = f'{indir}/{parameters}'
    genome_dir = argv[3]
    genome = argv[4]
    outdir = 'chrom_overview_2'
    
    #for genome in os.listdir(genome_dir):
    #    if genome.endswith('.fa'):
    #        genome = genome.replace('.clean.fa', '')
    #        wrapper(indir, parameters, comm_dir, genome_dir, outdir, genome)
    
    #wrapper(indir, parameters, comm_dir, genome_dir, outdir, genome)
    bl, adjust = splits_to_bed(f'{outdir}/scaff_split', genome_dir, f'{outdir}/paf_bed', indir, parameters)
    print(adjust)
    
    write_comm_files(adjust)
    


    
