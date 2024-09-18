"""
Create a fasta file per community

Include splitting of candidates
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

	for file in os.listdir():
		c_weird_windows = 0
		if file.endswith('.txt'):
			with open(f'../chrom_overview_2/scaff_split/{file}') as infile:
				ww_dict = json.load(infile)
				
				for scaff, ww in ww_dict.items():
				   # if len(ww) > 1:
					   #splittable = False
					blocks, comms = find_blocks_from_specific_keys(ww)

					if len(blocks) < 1:
						splittable = False

					elif len(blocks) >= 1:
						splittable = any(len(b) >= 10 for b in blocks)
						for indblock, b in enumerate(blocks):
							if len(b) >= 10:
								blocks = b
								comms = comms[indblock]
								comms = comms[0]

						if splittable:
							b_temp = [int(x.strip('*')) for x in b]
							for x in b:
								if x.startswith('*'):
									total_split.append((file,scaff, x, str(max(b_temp)), comms))
								elif x.endswith('*'):
									total_split.append((file,scaff, str(min(b_temp)), x, comms))
