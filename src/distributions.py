from __future__ import print_function,division
from future.utils import iteritems
import os,sys,pysam,collections,time,gc,random,itertools
from utils import *
import numpy as np
from global_vars import *
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pylab import rc
from scipy.stats import norm
from scipy.stats import gamma
import scipy.special as special
import scipy.optimize as optimize
import numpy as np
import mpmath


class negBin(object):
	def __init__(self, p = 0.1, r = 10):
		nbin_mpmath = lambda k, p, r: mpmath.gamma(k + r)/(mpmath.gamma(k+1)\
			*mpmath.gamma(r))*np.power(1-p, r)*np.power(p, k)
		self.nbin = np.frompyfunc(nbin_mpmath, 3, 1)
		self.p = p
		self.r = r
 
	def mleFun(self, par, data, sm):
		p = par[0]
		r = par[1]
		n = len(data)
		f0 = sm/(r+sm)-p
		f1 = np.sum(special.psi(data+r)) - n*special.psi(r) + n*np.log(r/(r+sm))
		return np.array([f0, f1])
 
	def fit(self, data, p = None, r = None):
		if p is None or r is None:
			av = np.average(data)
			va = np.var(data)
			r = (av*av)/(va-av)
			p = (va-av)/(va)
		sm = np.sum(data)/len(data)
		x = optimize.fsolve(self.mleFun, np.array([p, r]), args=(data, sm))
		self.p = x[0]
		self.r = x[1]
 
	def pdf(self, k):
		return self.nbin(k, self.p, self.r).astype('float64')



def plot_distribution(p,distr,xlab,ylab,title):
	fname = '_'.join(title.split(' '))
	nbins = 50
	fig, ax = plt.subplots()
	n, bins, patches = plt.hist(distr,nbins,normed=True, facecolor='blue', alpha=0.70)
	rc('axes', linewidth=1)
	y = [p(b) for b in bins]
	plt.plot(bins,y,color='r',linewidth=5)
	plt.xlabel(xlab)
	plt.ylabel(ylab)
	plt.title(title)
	x0 = max(n)
	xmax = max(bins)
	plt.axis([0, max(bins), 0, max(max(y),max(n))])
	fig.savefig(os.path.join(DIR,fname+'.pdf'),format='pdf')
	plt.close('all')
	return

def linked_reads(reads,chrom):
	reads.sort(key = lambda x: x[0])

	curr_LR = [0,0,0,0]
	LRs = []
	for start,end,hap,mapq in reads:
		if curr_LR[0] == 0 or start-curr_LR[2] > d:
			if curr_LR[0] != 0 and curr_LR[3] >= min_reads and curr_LR[2]-curr_LR[1] >= min_len:
				LRs.append(curr_LR)
			curr_LR = [chrom,start,end,1]
		else:
			curr_LR[2] = max(curr_LR[2],end)
			curr_LR[3] += 1
	if curr_LR[0] != 0 and curr_LR[3] >= min_reads and curr_LR[2]-curr_LR[1] >= min_len:
		LRs.append(curr_LR)
	return LRs


def get_overlap(barcode_LRs):
	global barcode_overlap
	for LR1,LR2 in itertools.combinations(barcode_LRs,2):
		if LR1[0] > LR2[0] or LR1[1] > LR2[1]:
			LR1,LR2 = [LR2,LR1]
		chr1,start1,end1,count1 = LR1
		chr2,start2,end2,count2 = LR2
		index1 = set([int(start1/d)*d,int(end1/d)*d])
		index2 = set([int(start2/d)*d,int(end2/d)*d])
		for id1 in index1:
			for id2 in index2:
				barcode_overlap[(chr1,id1,chr2,id2)] += 1



def get_distributions(reads_by_LR):
	LRs = []
	global barcode_overlap,LRs_by_barcode
	LRs_by_barcode = collections.defaultdict(list)
	barcode_overlap = collections.defaultdict(int)
	for key,reads in iteritems(reads_by_LR):
		chrom,barcode = key
		barcode_LRS = linked_reads(reads,chrom)
		LRs += barcode_LRS
		LRs_by_barcode[barcode] += barcode_LRS
	for barcode,barcode_LRS in iteritems(LRs_by_barcode):
		if len(barcode_LRS) > 1:
			get_overlap(barcode_LRS)
	if len(LRs) < 100:
		return None,None,None
	p_rate = get_rate_distr(LRs)
	p_len = get_length_distr(LRs)	
	return p_len,p_rate,barcode_overlap

def get_length_distr(LRs):
	lengths = [x[2]-x[1] for x in LRs]
	lengths.sort()
	assert len(lengths) >= 100
	assert np.var(lengths) != 0
	b = negBin()
	b.fit(lengths)
	p = b.pdf
	pp = lambda x: max(1e-20,float(p([x])[0]))
	## poisson distribution
	# p = scipy.stats.poisson(np.mean(lengths)).pmf
	# pp = lambda x: max(1e-20,float(p([x])[0]))
	#plot_distribution(pp,lengths,'Linked-read lengths (bp)','Frequency','Linked-read length distribution')
	return pp

def get_rate_distr(LRs):
	rate = [x[3]/float(x[2]-x[1]) for x in LRs]
	rate.sort()
	if len(rate) > 10:
		rate = rate[int(len(rate)/10):int(len(rate)/10*9)]
	alpha,loc,beta = scipy.stats.gamma.fit(rate)
	p = scipy.stats.gamma(alpha,loc,beta).cdf
	pp = lambda x: max(1e-20,float(p([max(x,1e-6)])[0]-p([max(x,1e-6)-1e-6])[0]))
	## normal distribution
	# mu,sig = scipy.stats.norm.fit(rate)
	# p = scipy.stats.norm(mu,sig).cdf
	# pp = lambda x: max(1e-20,float(p([max(x,1e-6)])[0]-p([max(x,1e-6)-1e-6])[0]))
	#plot_distribution(pp,rate,'Per-molecule sequencing rate','Frequency','Per-molecule sequencing rate')
	return pp
