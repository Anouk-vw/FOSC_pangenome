import os
from sys import argv

def get_comm_dict(indict):
	"""
	parse communitites to dictionary
	{members:comm}
	"""
	comm_dict = {}
	comms = [x for x in os.listdir(indict) if x.startswith('out.paf.edges.weights.txt.community')]
	for comm in comms:
		com_id = comm.split('.')[6]
		with open(f'{indict}/{comm}') as infile:
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
	
if __name__ == '__main__':
	node_list = get_file(argv[1], 'GCA_001703175.2_Fo_lycopersici4287_1')
	edge_list,weight_list = get_edges(argv[1])
	
	with open(f'{argv[1]}/gephi_edges.csv', 'w+') as outfile_edges:
		outfile_edges.write('Source;Target;Weight\n')
		for edge, weight in zip(edge_list,weight_list):
			outfile_edges.write(f'{edge};{weight}\n')
			
	with open(f'{argv[1]}/gephi_nodes.csv', 'w+') as outfile_nodes:
		for node in node_list:
			outfile_nodes.write(node)
