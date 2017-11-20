from __future__ import print_function,division
import os,sys,subprocess,pysam,collections,time,gc
import multiprocessing as mp
from utils import *
import numpy as np
from global_vars import *

def estimate_lmin_lmax():
	reads = pysam.AlignmentFile(BAM_FILE,"rb")
	mate_pairs = dict()
	ls = []
	num = 0
	for i,chrm in enumerate(reads.references):
		start = int(reads.lengths[i]/2)
		end = start+1000000
		for read in reads.fetch(chrm,start,end):
			num += 1
			if num > 1000000:
				break
			if read.query_name in mate_pairs:
				if read.reference_start > mate_pairs[read.query_name][0]:
					dist = read.reference_start-mate_pairs[read.query_name][1]
				else:
					dist = mate_pairs[read.query_name][0]-read.reference_end
				if dist < 2000:
					ls.append(dist)
			else:
				mate_pairs[read.query_name] = (read.reference_start,read.reference_end)
	mean_dist = np.median(ls)
	std_dist = np.std(ls)
	lmin = int(mean_dist-std_dist*2)
	lmax = int(mean_dist+std_dist*2)
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
	d = max(1000,int(gaps[len(gaps)-int((len(gaps)/100)*5)]/1000)*1000)
	return d


