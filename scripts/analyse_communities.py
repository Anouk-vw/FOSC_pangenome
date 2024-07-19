import os
from sys import argv

def get_comm_dict(indir):
	"""
	parse communitites to dictionary
	{members:comm}
	"""
	comm_dict = {}
	comms = [x for x in os.listdir(indir) if x.startswith('out.paf.edges.weights.txt.community')]
	for comm in comms:
		com_id = comm.split('.')[6]
		with open(f'{indir}/{comm}') as infile:
			for line in infile:
				member = line.strip()
				comm_dict[member] = com_id
	return comm_dict

def ref_ids(ref, comm_dict):
	"""
	Get all ref chroms
	"""
	ref_ids = []
	ref_dict = {}
	for member, comm in comm_dict.items():
		if member.startswith(ref):
			ref_id = member.split('scaffold')[-1]
			ref_ids.append(ref_id)
			try:
				ref_dict[comm].append(ref_id)
			except KeyError:
				ref_dict[comm] = [ref_id]
	return ref_dict

def get_file(indict, ref):
	"""
	return Node_ID;name;comm;ref_chrom
	"""
	nodes = []
	nodes.append(f'Id;Label;community;ref\n')
	comm_dict = get_comm_dict(indict)
	ref_dict = ref_ids(ref, comm_dict)
	with open(f'{indict}/out.paf.vertices.id2name.txt') as infile:
		for line in infile:
			node, name = line.strip().split(' ')
			community = comm_dict[name]
			try:
				ref = int(ref_dict[community][0])
			except KeyError:
				ref = 0
			nodes.append(f'{node};{name};{community};{ref}\n')
	return nodes

def get_edges(indict):
	edges = []
	weights = []
	with open(f'{indict}/out.paf.edges.list.txt') as infile:
		for line in infile:
			edges.append(line.strip().replace(' ',';'))
			
	with open(f'{indict}/out.paf.edges.weights.txt') as infile:
		for line in infile:
			weights.append(line.strip())

	return edges, weights

def read_paf(in_paf, scaff_comm_dict):
    """
    read paf file
    scaff_comm_dict = {scaffold_id: community}
    """
    lines = []
    with open(in_paf) as ifile:
        for line in ifile:
            query, qlen, qstart, qstop, qstrand, ref, rlen, rstart, rstop, rstraind, matches, blocklen, identity = line.strip().split('\t')
            try:
                community = int(scaff_comm_dict[ref])
                qcomm = int(scaff_comm_dict[query])
            except KeyError:
                community = 'unplaced'
                qcomm = 'unplaced'
            lines.append([query, qlen, qstart, qstop, qstrand, ref, rlen, rstart, rstop, rstraind, matches, blocklen, identity, int(community), int(qcomm)])
    return lines

def write_out(indir):
    """
    write output with comm identities
    """
    cdict = get_comm_dict(indir)
    paf_lines = read_paf(f'{indir}/out.paf',cdict)
    with open(f'{indir}/plot_qascomm.paf', 'w+') as outfile:
        for line in paf_lines:
            line[0] = f'comm-{line[-1]}'
            line = line[0:-2]
            outfile.write('\t'.join(line))
            outfile.write('\n')
            
    paf_lines = read_paf(f'{indir}/out.paf', cdict)
    with open(f'{indir}/plot_tascomm.paf', 'w+') as outfile:
        for line in paf_lines:
            line[5] = f'comm-{line[-2]}'
            line = line[0:-2]
            outfile.write('\t'.join(line))
            outfile.write('\n')
    return

if __name__ == '__main__':
	node_list = get_file(argv[1], 'GCA_000149955.2_Fo_lycopersici4287_1')
	edge_list,weight_list = get_edges(argv[1])
	
	with open(f'{argv[1]}/gephi_edges.csv', 'w+') as outfile_edges:
		outfile_edges.write('Source;Target;Weight\n')
		for edge, weight in zip(edge_list,weight_list):
			outfile_edges.write(f'{edge};{weight}\n')
			
	with open(f'{argv[1]}/gephi_nodes.csv', 'w+') as outfile_nodes:
		for node in node_list:
			outfile_nodes.write(node)

	write_out(argv[1])
