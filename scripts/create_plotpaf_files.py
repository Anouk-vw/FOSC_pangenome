from sys import argv
import os

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

def write_out(indir, cdict):
    """
    write output with comm identities
    """
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

if __name__ == '__main__':
    indir = argv[1]
    cdict = get_comm_dict(indir)
    write_out(indir, cdict)
