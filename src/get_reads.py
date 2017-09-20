import os,sys,subprocess,pysam,collections,time,gc,json,random,copy
import multiprocessing as mp
from utils import *
from global_vars import *
from distributions import get_distributions


def inblacklist(cand):
	'''
	Args: candidate novel adjacency
	Returns: True if candidate novel adjacency overlaps with
	interval in blacklist
	'''
	if not blacklist:
		return False
	for s,e in blacklist[cand.chrm]:
		if cand.i > s and cand.i < e:
			return True
	for s,e in blacklist[cand.nextchrm]:
		if cand.j > s and cand.j < e:
			return True	 
	return False


def make_barcodeDict(chrom):
	'''
	Args: chromosome
	Function: Parses postion sorted BAM file and extracts relevant information
	Returns: dict barcoded reads, dict of reads overlapping ref positions, 
	discordant reads, candidate NAs, total coverage
	'''
	print chrom
	cov = 0
	global discs
	global reads_by_barcode
	global barcode_to_num
	global barcode_count
	global discs_by_barcode
	barcode_count = 0
	num_disc = 0
	num_reads = []
	reads = pysam.AlignmentFile(BAM_FILE,"rb")
	count = 0
	r = lmax
	lengths = reads.lengths[list(reads.references).index(chrom)]
	reads_by_LR = collections.defaultdict(list)
	discs_by_barcode = collections.defaultdict(list)
	discs = collections.defaultdict(list)
	interchrom_discs = collections.defaultdict(list)
	barcode_to_num = dict()
	disc_mate_pairs = dict()
	split_reads = collections.defaultdict(list)
	LRs_by_pos = collections.defaultdict(list)
	lengths = sum(reads.lengths)
	iterator = reads.fetch(chrom)
	starttime = time.time()
	prevtime = starttime
	count = 0
	for read in iterator:
		cov += read.query_alignment_length
		## DEBUG
		# if cov > 5000000000:
		# 	break
		if pass_checks(read):
			barcode = get_barcode(read)
			peread = PERead(read)
			if peread.disc:
					peread = PERead(read)
					peread.add_disc(read)
					discs_by_barcode[(peread.chrm,peread.nextchrm,peread.barcode)].append(peread)
					if peread.chrm == peread.nextchrm:
						discs = add_disc(peread,discs)
					else:
						interchrom_discs[(peread.chrm,peread.i/r*r,peread.nextchrm,peread.j/r*r,peread.orient)].append(peread)
			elif read.is_proper_pair and fragment_length(read) > lmin:
					reads_by_LR[(peread.chrm,barcode)].append((peread.start,peread.nextend,peread.hap,peread.mapq))
					if barcode not in LRs_by_pos[(peread.chrm,peread.mid()/R*R)]:
						LRs_by_pos[(peread.chrm,peread.mid()/R*R)].append(barcode)
	return reads_by_LR,LRs_by_pos,discs_by_barcode,discs,interchrom_discs,cov/float(lengths)

def signi(disc):
	if disc.orient[0] == '+':
		return  disc.i+lmin,disc.i+lmax
	else:
		return disc.i-lmax,disc.i-lmin

def signj(disc):
	if disc.orient[1] == '+':
		return  disc.j+lmin,disc.j+lmax
	else:
		return disc.j-lmax,disc.j-lmin


def add_disc(peread,discs):
	r = lmax
	for a in [-r/2,0,r/2]:
		for b in [-r/2,0,r/2]:
			if (a == 0 or (peread.i+a)/r*r != peread.i/r*r) and (b == 0 or (peread.j+b)/r*r != peread.j/r*r):
				discs[(peread.chrm,(peread.i+a)/r*r,peread.nextchrm,(peread.j+b)/r*r,peread.orient)].append(peread)
	return discs

def largest_overlap(items):
	positions_i = collections.defaultdict(int)
	positions_j = collections.defaultdict(int)
	for disc in items:
		start,end = signi(disc)
		for pos in range(start,end):
			positions_i[pos] += 1
		start,end = signj(disc)
		for pos in range(start,end):
			positions_j[pos] += 1
	i = 0
	overlap_i = 0
	for pos,overlap in positions_i.iteritems():
		if overlap > overlap_i:
			overlap_i = overlap
			i = pos
	j = 0
	overlap_j = 0
	for pos,overlap in positions_j.iteritems():
		if overlap > overlap_j:
			overlap_j = overlap
			j = pos
	return i,j,overlap_i,overlap_j

def coordinates(s,e,orient):
	if orient == '+':
		return s
	return e

def disc_intersection(items):	
	positions_i = collections.defaultdict(int)
	positions_j = collections.defaultdict(int)
	intersection = [0,float('Inf'),0,float('Inf')]
	items.sort(key = lambda x: x.i)
	for disc in items:
		starti,endi = signi(disc)
		startj,endj = signj(disc)
		si,ei,sj,ej = intersection
		intersection = [max(si,starti),min(ei,endi),max(sj,startj),min(ej,endj)]
	si,ei,sj,ej = intersection
	if si and sj and si < ei and sj < ej:
		return si,ei,sj,ej
	return 0,0,0,0


def get_candidates(discs,reads_by_LR):
	candidates = []
	r = lmax
	p_len,p_rate,barcode_overlap = get_distributions(reads_by_LR)
	num_cands = 0
	for key,items in discs.iteritems():
		orient = key[4]
		si,ei,sj,ej = disc_intersection(items)
		if si and sj and len(items) >= MIN_DISCS:
			i = coordinates(si,ei,orient[0])
			j = coordinates(sj,ej,orient[1])
			cand = copy.copy(items[0])
			cand.i = i
			cand.j = j
			barcode_overlaps = barcode_overlap[(cand.chrm,cand.i/d*d,cand.nextchrm,cand.j/d*d)]
			if not inblacklist(cand) and ((cand.chrm == cand.nextchrm and cand.j-cand.i < d)\
				or barcode_overlaps >= k):
				already_appended = sum([1 for x in candidates if x.i == cand.i and x.j == cand.j])
				if not already_appended:
					num_cands += 1
					candidates.append(cand)			
	return candidates,p_len,p_rate


def make_barcodeDict_user(candidate):
	'''
	Args: candidate novel adjacency input by user
	Function: Parses postion sorted BAM file and extracts positions indicated by candidate NAs
	Returns: dict barcoded reads, dict of reads overlapping ref positions, 
	discordant reads, total coverage
	'''
	w = max_len-d
	cov = 0
	global discs
	global reads_by_barcode
	global barcode_to_num
	global barcode_count
	global discs_by_barcode
	barcode_count = 0
	num_disc = 0
	num_reads = []
	r = lmax
	reads = pysam.AlignmentFile(BAM_FILE,"rb")
	count = 0
	lengths = 0
	reads_by_LR = collections.defaultdict(list)
	discs_by_barcode = collections.defaultdict(list)
	barcode_to_num = dict()
	disc_mate_pairs = dict()
	split_reads = collections.defaultdict(list)
	LRs_by_pos = collections.defaultdict(list)
	lengths = sum(reads.lengths)
	starttime = time.time()
	prevtime = starttime
	count = 0

	chrm1,break1,chrm2,break2,orientation = candidate
	if break1 > break2 or chrm1 > chrm2:
		chrm1,break1,chrm2,break2 = chrm2,break2,chrm1,break1
	if chrm1 != chrm2 or break1+w < break2-w:
		window = [(chrm1,break1-w,break1+w),(chrm2,break2-w,break2+w)]
	else:
		window = [(chrm1,break1-w,break2+w)]
	for chrm,s,e in window:
		lengths += e-s
		iterator = reads.fetch(chrm,max(0,s),e)
		count = 0
		for read in iterator:
			cov += read.query_alignment_length
			if pass_checks(read):
				barcode = get_barcode(read)
				peread = PERead(read)
				if peread.disc:
						peread = PERead(read)
						peread.add_disc(read)
						discs_by_barcode[(peread.chrm,peread.nextchrm,peread.barcode)].append(peread)
				elif read.is_proper_pair and fragment_length(read) > lmin:
						reads_by_LR[(peread.chrm,barcode)].append((peread.start,peread.nextend,peread.hap,peread.mapq))
						if barcode not in LRs_by_pos[(peread.chrm,peread.mid()/R*R)]:
							LRs_by_pos[(peread.chrm,peread.mid()/R*R)].append(barcode)
	cand = copy.copy(peread)
	cand.chrm = chrm1.strip('chr')
	cand.i = break1
	cand.nextchrm = chrm2.strip('chr')
	cand.j = break2
	cand.orient = orientation
	if (cand.chrm == cand.nextchrm and cand.i > cand.j) or (cand.chrm != cand.nextchrm and cand.chrm > cand.nextchrm):
		cand.chrm = chrm2.strip('chr')
		cand.i = break2
		cand.nextchrm = chrm1.strip('chr')
		cand.j = break1
		cand.orient = orientation[::-1]
	return reads_by_LR,LRs_by_pos,discs_by_barcode,[cand],cov/float(lengths)


def pass_checks(read):
	try:
		good_mapq = read.mapping_quality >= MIN_MAPQ
		not_dup = not read.is_duplicate
		is_chrom = is_proper_chrom(read.reference_name) and is_proper_chrom(read.next_reference_name)
		has_barcode = read.has_tag('BX')
		primary = not read.is_secondary and not read.is_supplementary
		is_first = first_read(read)
		haplotyped = read.has_tag('HP')
	except:
		return False
	return good_mapq and not_dup and is_chrom and has_barcode and primary and is_first



