from global_vars import *
import math,mpmath,copy
import numpy as np
import linecache
import scipy.stats,collections
import multiprocessing as mp

class NA:
	__slots__ = ['score','chrm1','break1','chrm2','break2','orient','haps','disc','pairs','spans']
	def __init__(self,chrm1,chrm2,indi,indj,orient):
		self.score = collections.defaultdict(int)
		self.chrm1 = chrm1
		self.break1 = indi
		self.chrm2 = chrm2
		self.break2 = indj
		self.orient = orient
		self.haps = (0,0)
		self.spans = collections.defaultdict(int)
		self.disc = collections.defaultdict(int)
		self.pairs = collections.defaultdict(int)
	def get_score(self):
		best_score = -float('inf')
		best_haps = (0,0)
		total = 0
		for hap in self.score.keys():
			s = self.score[hap]
			d = self.disc[hap]
			if self.pairs[hap] > self.spans[hap]:
				total += s
			if hap[0] != 0 and self.pairs[(0,hap[1])] >= self.spans[(0,hap[1])]:
				s += max(0,self.score[(0,hap[1])])
				d += self.disc[(0,hap[1])]
			if hap[1] != 0 and self.pairs[(hap[0]),0] >= self.spans[(hap[0]),0]:
				s += max(0,self.score[(hap[0]),0])
				d += self.disc[(hap[0]),0]
			if hap[0] != 0 and hap[1] != 0 and self.pairs[(0,0)] >= self.spans[(0,0)]:
				s += max(0,self.score[(0,0)])
				d += self.disc[(0,0)]
			if s > best_score and self.pairs[hap] > 0 and self.pairs[hap] >= self.spans[hap]:
				best_score = s
				best_haps = hap
		pairs = sum([self.pairs[x] for x in self.pairs.keys()])
		spans = sum([self.spans[x] for x in self.spans.keys()])
		discs = sum([self.disc[x] for x in self.disc.keys()])
		if best_score > -float('inf') and total > best_score and pairs >= spans and\
		((self.score[(1,1)] > 0 and self.score[(2,2)] > 0) or (self.score[(1,2)] > 0 and self.score[(2,1)] > 0)):
			self.haps = (3,3)
			self.score = total
			self.disc = discs
			self.pairs = pairs
			self.spans = spans
		elif best_score > -float('inf'):
			self.haps = best_haps
			score = self.score[best_haps]
			discs = self.disc[best_haps]
			pairs = self.pairs[best_haps]
			spans = self.spans[best_haps]
			if best_haps[0] != 0 and self.score[(0,best_haps[1])] > 0:
				score += self.score[(0,best_haps[1])]
				discs += self.disc[(0,best_haps[1])]
				pairs += self.pairs[(0,best_haps[1])]
				spans += self.spans[(0,best_haps[1])]
			if best_haps[1] != 0 and self.score[(best_haps[0],0)] > 0:
				score += self.score[(best_haps[0],0)]
				discs += self.disc[(best_haps[0],0)]
				pairs += self.pairs[(best_haps[0],0)]
				spans += self.spans[(best_haps[0],0)]
			if best_haps[1] != 0  and best_haps[0] != 0 and self.score[(0,0)] > 0:
				score += self.score[(0,0)]
				discs += self.disc[(0,0)]
				pairs += self.pairs[(0,0)]
				spans += self.spans[(0,0)]
			self.score = score
			self.disc = discs
			self.pairs = pairs
			self.spans = spans
		else:
			self.score = best_score

class LinkedRead:
	__slots__ = ['barcode', 'chrm','start','disc','end','hap','num']
	def __init__(self,barcode):
		self.chrm = 0
		self.start = 0
		self.end = 0
		self.hap = 0
		self.num = 1
		self.disc = []
		self.barcode = barcode
	def add_num(self):
		self.num += 1
	def length(self):
		return float(self.end-self.start)


class PERead:
	__slots__ = ['barcode', 'chrm','orient','start','disc','mapq','end','hap','nextchrm',\
	'nextstart','nextend','nexthap','nextmapq']
	def __init__(self,read):
		read_hap = get_hap(read)
		r = [read.reference_name.strip('chr'),read.reference_start,read.reference_end,\
		read.get_tag('BX'),read.is_reverse,read.mapping_quality]
		m = [read.next_reference_name.strip('chr'),read.next_reference_start,0,0,0,0]
		# if r[1] > m[1] or r[0] > m[0]:
		# 	[chr1,read_start,read_end,barcode1,isrev1,mapq1] = m
		# 	[chr2,mate_start,mate_end,barcode2,isrev2,mapq2] = r
		# else:
		[chr1,read_start,read_end,barcode1,isrev1,mapq1] = r
		[chr2,mate_start,mate_end,barcode2,isrev2,mapq2] = m
		self.barcode = barcode1
		self.chrm = chr1
		self.nextchrm = chr2
		self.start = read_start
		self.mapq = mapq1
		self.orient = get_orient(read)
		self.end = read_end
		self.hap = read_hap
		self.mapq = mapq1
		self.nextstart = mate_start
		self.nextend = mate_start+(read_end-read_start)
		self.disc = self.is_disc() and not read.is_proper_pair

	def is_disc(self):
		if self.orient[0] == '+':
			i = self.end
		else:
			i = self.start
		if self.orient[1] == '+':
			j = self.nextend
		else:
			j = self.nextstart
		return (self.chrm != self.nextchrm or (j-i) > min_sv)

	def mid(self):
		return int((self.end+self.start)/2)
		
	def add_disc(self,read):
		chrom1 = read.reference_name.strip('chr')  
		chrom2 = read.next_reference_name.strip('chr') 
		if read.is_reverse:
			i = read.reference_start
		else:
			i = read.reference_end
		if read.mate_is_reverse:
			j = read.next_reference_start
		else:
			j = read.next_reference_start+read.reference_length
		hap = get_hap(read)
		mapq = read.mapping_quality
		orient = get_orient(read)
		if not first_read(read):
			chrom1,chrom2 = swap(chrom1,chrom2)
			i,j = swap(i,j)
			orient = orient[::-1]
		self.chrm = chrom1
		self.nextchrm = chrom2
		self.i = i
		self.j = j
		self.hap = hap
		self.nexthap = hap
		self.mapq = mapq
		self.nextmapq = mapq
		self.orient = orient


	def add_split(self,read):
		chrm,start,orient,cigar,mapq,hap = read.get_tag('SA').split(',')
		start = int(start)
		mapq = int(mapq)
		hap = int(hap[0])
		if start > self.start or chrm.strip('chr') >= self.chrm:
			self.nextchrm = chrm.strip('chr')
			self.j = start
			self.nexthap = hap
			self.nextmapq = mapq
			self.orient = get_orient(self.orient,orient == '-')
		else:
			self.nextchrm = self.chrm
			self.j = self.i
			self.nexthap = self.hap
			self.nextmapq = self.mapq
			self.chrm = chrm.strip('chr')
			self.i = start
			self.hap = hap
			self.mapq = mapq
			self.orient = get_orient(orient == '-',self.orient)	

def first_read(read):
	cond1 = read.reference_name < read.next_reference_name
	cond2 = read.reference_name == read.next_reference_name and read.reference_start < read.next_reference_start
	return cond1 or cond2


def closest_LR(LR1,LR2,i):
	if not LR1:
		return LR2
	if not LR2:
		return LR1
	LR1_dist = min(abs(LR1[1]-i),abs(LR1[2])-i)
	LR2_dist = min(abs(LR2[1]-i),abs(LR2[2])-i)
	if LR1_dist < LR2_dist:
		return LR1
	return LR2

def swap(a,b):
	c = copy.copy(b)
	d = copy.copy(a)
	a = c
	b = d
	return a,b

def is_proper_chrom(chrom):
	return 'Un' not in chrom and 'random' not in chrom

def get_orient(read):
	a = ''
	if read.is_reverse:
		a += '-'
	else:
		a += '+'
	if read.mate_is_reverse:
		a += '-'
	else:
		a += '+'
	return a

def get_barcode(barcode):
	thefilename = os.path.join(DIR,'LRS.txt')
	line = linecache.getline(thefilename, barcode+1)
	b,line = line.split(':')
	assert b == str(barcode)
	line = line.split('\t')
	line = [[[int(z) for z in y.split(',')] for y in x.split(' ') if len(y.split(',')) == 4] for x in line]
	line = [x for x in line if len(x) > 0]
	return line

def NB(k,m,r):
	return mpmath.gamma(k + r)/(mpmath.gamma(k+1)*mpmath.gamma(r))*np.power(m/(r+m),k)*np.power(r/(r+m),r)

def Pois(k,lam):
	return (np.power(lam,k)*math.exp(-lam))/mpmath.gamma(k+1)
def get_hap(r):
	if r.has_tag('HP'):
		return r.get_tag('HP')
	return 0

def get_hapb(r):
	if r.has_tag('PS'):
		return r.get_tag('PS')
	return 0
	
def get_barcode(r):
	if r.has_tag('BX'):
		return r.get_tag('BX')
	return 0

def to_chr(chrm):
	return chrm
	chrm = chrm.strip('chr')
	# if chrm == '23':
	# 	return 'X'
	# elif chrm == 'Y':
	# 	return '24'
	# return str(chrm)

def plog(x,num):
	ret = 0
	for i in range(num):
		ret += log(x)
	return ret

def log(x):
	try:
		return math.log(x)
	except:
		return math.log(1e-300)

def orient(read):
	if read.is_reverse:
		return '-'
	return '+'



def is_close(index,read_pos,orient):
	if orient == '+':
		return abs(index-read_pos) <= lmax
	else:
		return abs(index-read_pos) <= lmax

def is_convergent(read,mate):
	return not read.is_reverse and mate.is_reverse

def get_orientation(read,mate):
	# if read.is_read1:
	if read.is_reverse:
		if mate.is_reverse:
			return ['--',read.reference_start,mate.reference_start]
		else:
			return ['-+',read.reference_start,mate.reference_end]
	else:
		if read.mate_is_reverse:
			return ['+-',read.reference_end,mate.reference_start]
		else:
			return ['++',read.reference_end,mate.reference_end]

def rnd(num):
	return num/R*R

def flatten(l):
	return [item for sublist in l for item in sublist]


def threshold(cov):
	return 6.943*cov-37.33


def fragment_length(read):
	return max(read.reference_end,read.next_reference_start+read.reference_length)\
	-min(read.reference_start,read.next_reference_start)

def collapse(a,c):
	l = []
	hom = 0
	for linea in a:
		chrom1,s1,chrom2,s2 = linea[0:4]	
		if chrom1 == '24':
			chrom1 = 'Y'
		if chrom1 == '23':
			chrom1 = 'X'
		if chrom2 == '24':
			chrom2 = 'Y'
		if chrom2 == '23':
			chrom2 = 'X'
		split = float(linea[4].split(':')[-1])
		disc = float(linea[5].split(':')[-1])
		orient = linea[6].split(':')[-1]
		haps = linea[7].split(':')[-1]
		score = float(linea[8].split(':')[-1])
		if 'Y' not in chrom1 and 'Y' not in chrom2:
			l.append([chrom1,int(s1),chrom2,int(s2),split,disc,orient,haps,score])
	r = lmax/100*100*5
	l.sort(key = lambda x: x[-1],reverse=True)
	l2 = []
	nas = collections.defaultdict(list)
	for chrom1,s1,chrom2,s2,split,disc,orient,haps,score in l:
		already_appended = False
		for i in [-r,0,r]:
			for j in [-r,0,r]:
				if (chrom1,s1/r*r+i,chrom2,s2/r*r+j) in nas:
					ch1,ds1,ch2,ds2 = nas[(chrom1,s1/r*r+i,chrom2,s2/r*r+j)][0:4]
					if abs(ds1-s1) < r and abs(ds2-s2) < r:
						already_appended = True
		if not already_appended:
			nas[(chrom1,s1/r*r,chrom2,s2/r*r)] = [chrom1,s1,chrom2,s2,split,disc,orient,haps,score]
	for i,elem in nas.iteritems():
		if elem[-1] >= c:
			l2.append(elem+['PASS'])
		else:
			l2.append(elem+['FAIL'])
	l2.sort(key = lambda x: (x[-2]),reverse=True)
	l2 = ['\t'.join([str(y) for y in x])+'\n' for x in l2]
	return l2


def write_scores(scores):
	fname = 'NAIBR_SVs.bedpe'
	print os.path.join(DIR,fname)
	f = open(os.path.join(DIR,fname),'w')
	f.write('\t'.join(['Chr1','Break1','Chr2','Break2','Split molecules','Discordant reads',\
	'Orientation','Haplotype','Score','Pass filter\n']))
	for result in scores:
		f.write(result)
	f.close()
	return

def parallel_execute(function,input_list):
	if NUM_THREADS != 1:
		pool = mp.Pool(NUM_THREADS,maxtasksperchild=1)
		map_fn = pool.map
		data = map_fn(function,input_list)
		pool.close()
		pool.join()
	else:
		data = map(function,input_list)
	return data


