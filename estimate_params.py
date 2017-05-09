import os,sys,subprocess,pysam,collections,time,gc
import multiprocessing as mp
from utils import get_hap,get_orientation,swap,rnd 
import numpy as np
from global_vars import *

def estimate_lmin_lmax():
	reads = pysam.AlignmentFile(BAM_FILE,"rb")
	mate_pairs = dict()
	counter0 = 0
	counter1 = 0
	counter2 = 0
	ls = []
	start = 10000000
	end = 11000000
	for read in reads.fetch('chr10',start,end):
		if get_hap(read) == 0:
			counter0+=read.query_alignment_length
		elif get_hap(read) == 1:
			counter1+=read.query_alignment_length
		elif get_hap(read) == 2:
			counter2+=read.query_alignment_length
		if read.query_name in mate_pairs:
			if abs(read.reference_end-mate_pairs[read.query_name]) < 2000:
				ls.append(read.reference_end-mate_pairs[read.query_name])
		else:
			mate_pairs[read.query_name] = read.reference_start
	mean_dist = max(0,np.median	(ls))
	std_dist = max(0,np.std(ls))
	lmin = max(0,int(mean_dist-std_dist*3))
	lmax = max(0,int(mean_dist+std_dist*3))
	coverage0 = counter0/float(end-start)
	coverage1 = counter1/float(end-start)
	coverage2 = counter2/float(end-start)
	print 'lmin is',lmin
	print 'lmax is',lmax
	print 'coverage of hap0 is',coverage0
	print 'coverage of hap1 is',coverage1
	print 'coverage of hap2 is',coverage2
	return lmin,lmax

def estimate_delta(read_list):
	counter=0
	gaps=[]
	for b,l in read_list:
		counter+=1
		if counter > 10000:
			break
		l.sort(key=lambda x: x[0])
		qname = ''
		prev_end = 0
		for (chrom,start,end,hap,qname2,mapq) in l:
			if qname2 != qname and prev_end > 0 and start-prev_end < max_d:
				gaps.append(start-prev_end)
			qname = qname2
			prev_end = end
	mean = np.mean(gaps)
	gaps.sort()
	d = max(1000,gaps[len(gaps)-int((len(gaps)/100)*5)]/1000*1000)
	print 'd is estimated to be',d
	return d


