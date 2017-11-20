### GLOBALS ####
from __future__ import print_function
import os,json,sys,collections,pysam
import numpy as np

def parse_blacklist(fname):
	with open(fname) as f:
		data = f.read().split('\n')
	blacklist = collections.defaultdict(list)
	for line in data:
		l = line.split('\t')
		if len(l) > 2:
			blacklist[l[0].strip('chr')].append([int(l[1]),int(l[-1])])
	return blacklist



with open(sys.argv[1]) as f:
	data = f.read().split('\n')
	constants = dict()
	for line in data:
		if line and line[0] != '#':
			name,val = line.split('=')
			constants[name] = val


if 'min_mapq' in constants:
	MIN_MAPQ = int(constants['min_mapq'])
else:
	MIN_MAPQ = 40


if 'k' in constants:
	k = int(constants['k'])
else:
	k = 3

if 'bam_file' in constants:
	BAM_FILE=constants['bam_file']
else:
	raise ValueError('Please add input Bam file to config file')


if 'outdir' in constants:
	DIR=constants['outdir']
	if not os.path.exists(DIR):
		os.makedirs(DIR)
else:
	DIR='.'


if 'd' in constants:
	d = int(constants['d'])
else:
	d = 10000

if 'threads' in constants:
	NUM_THREADS = int(constants['threads'])
else:
	NUM_THREADS = 1

if 'sd_mult' in constants:
	SD_MULT = int(constants['sd_mult'])
else:
	SD_MULT = 2


def estimate_lmin_lmax():
	reads = pysam.AlignmentFile(BAM_FILE,"rb")
	mate_pairs = dict()
	ls = []
	length = []
	num = 0
	for i,chrm in enumerate(reads.references):
		for read in reads.fetch(chrm):
			num += 1
			if num > 1000000:
				break
			if read.query_name in mate_pairs:
				if read.reference_start > mate_pairs[read.query_name][0]:
					dist = read.reference_start-mate_pairs[read.query_name][1]
				else:
					try:
						dist = mate_pairs[read.query_name][0]-read.reference_end
					except:
						pass
				if abs(dist) < 2000:
					length.append(read.query_length)
					ls.append(dist)
			else:
				mate_pairs[read.query_name] = (read.reference_start,read.reference_end)
	mean_dist = np.mean(ls)
	std_dist = np.std(ls)
	lmin = max(-int(np.mean(length)),int(mean_dist-std_dist*SD_MULT))
	lmax = int(mean_dist+std_dist*SD_MULT)
	return lmin,lmax

lmin,lmax = estimate_lmin_lmax()


if 'min_sv' in constants:
	min_sv = int(constants['min_sv'])
else:
	min_sv = lmax

if 'min_len' in constants:
	min_len = int(constants['min_len'])
else:
	min_len = 2*lmax

if 'max_len' in constants:
	max_len = int(constants['max_len'])
else:
	max_len = 200000

if 'min_reads' in constants:
	min_reads = int(constants['min_reads'])
else:
	min_reads = 2


if 'min_discs' in constants:
	MIN_DISCS = int(constants['min_discs'])
else:
	MIN_DISCS = 2

if 'candidates' in constants:
	candidates = constants['candidates']
	if not os.path.exists(candidates):
		candidates = ""
else:
	candidates = ""


if 'blacklist' in constants and os.path.exists(constants['blacklist']):	
	if not os.path.exists(constants['blacklist']):
		blacklist = None
	else:
		blacklist = parse_blacklist(constants['blacklist'])
else:
	blacklist = None

R = d

