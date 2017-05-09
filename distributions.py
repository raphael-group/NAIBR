import os,sys,pysam,collections,time,gc,random
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
		nbin_mpmath = lambda k, p, r: mpmath.gamma(k + r)/(mpmath.gamma(k+1)*mpmath.gamma(r))*np.power(1-p, r)*np.power(p, k)
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
		return self.nbin(k, self.p, self.r)

def plot_distribution(p,distr,xlab,ylab,title):
	fname = '_'.join(title.split(' '))
	nbins = 50
	fig, ax = plt.subplots()
	n, bins, patches = plt.hist(distr,nbins,normed=1, facecolor='blue', alpha=0.70)
	rc('axes', linewidth=1)
	bins = [int(x) for x in bins]
	y = p(bins)
	plt.plot(bins,y,color='r',linewidth=5)
	plt.xlabel(xlab)
	plt.ylabel(ylab)
	plt.title(title)
	x0 = max(n)
	xmax = max(bins)
	plt.axis([0, max(distr), 0, max(n)])
	fig.savefig(os.path.join(DIR,fname+'.pdf'),format='pdf')
	plt.close('all')
	return

def get_distributions(LENS,RATES):
	p_len = get_length_distr(LENS)
	p_rate = get_rate_distr(RATES)
	return p_len,p_rate

def get_length_distr(lengths):
	lengths = [x for x in lengths if x < 200000]
	b = negBin()
	b.fit(lengths)
	print '------------------------------------------------------------'
	print 'mean of length distr',np.mean(lengths),'sd of length distr',np.std(lengths)
	print 'neg binom fit'
	print 'p',b.p,'r',b.r
	print '------------------------------------------------------------'
	p = b.pdf
	plot_distribution(p,lengths,'Linked-read lengths (bp)','Frequency','Linked-read length distribution')
	return p

def get_rate_distr(rate):
	b = negBin()
	rate = [2.0/max(1.0/100000,x) for x in rate]
	b.fit(rate)
	print '------------------------------------------------------------'
	print 'mean of rate distr',np.mean(rate),'sd of rate distr',np.std(rate)
	print 'gamma fit'
	print 'p',b.p,'r',b.r
	print '------------------------------------------------------------'
	# p = scipy.stats.gamma(alpha,loc,beta).pdf
	# pp = lambda x: p(x*1.0)
	p = b.pdf
	pp = lambda x: p(2.0/max(1.0/100000,x))
	plot_distribution(p,rate,'Per-molecule sequencing rate','Frequency','Per-molecule sequencing rate')
	return pp
