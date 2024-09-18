"""
Script to match query-target names from paf to community
This can be used to colour chromosome locations based on community
1. For direct plotting in pafR (colors will overlap)
2. For splitting chromosomes, using bedtools intersect -> percentage of
comm per window

#incase chromosomes should be coloured according to flagN communities, indir_comm == flagN indir

Its not a good idea to use flagN as indir to determine coordinates in paf. Because -N makes that
the full chromosome maps, also when only a part should be mapped. Increasing the match size


"""

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
                community = -1
                qcomm = -1
            lines.append([query, qlen, qstart, qstop, qstrand, ref, rlen, rstart, rstop, rstraind, matches, blocklen, identity, int(community), int(qcomm)])
    return lines

def write_out(indir, cdict, oudir):
    """
    write output with comm identities
    """
    paf_lines = read_paf(f'{indir}/out.paf',cdict)
    with open(f'{outdir}/plot_qascomm.paf', 'w+') as outfile:
        for line in paf_lines:
            line[0] = f'comm-{line[-1]}'
            line = line[0:-2]
            outfile.write('\t'.join(line))
            outfile.write('\n')
            
    paf_lines = read_paf(f'{indir}/out.paf', cdict)
    with open(f'{outdir}/plot_tascomm.paf', 'w+') as outfile:
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
    indir_comm = argv[2]
    outdir = argv[3]
    #incase chromosomes should be coloured according to flagN communities
    #but based on full-mappings, indir comm is the flagN comm. ekse indir_comm == indir.
    cdict = get_comm_dict(indir_comm)
    write_out(indir, cdict, outdir)
