import os,sys,subprocess,pysam,collections,time,gc,json
import multiprocessing as mp
from utils import spans,list_from_key,make_key,get_hap,get_orientation,swap,rnd,flatten
import numpy as np
from global_vars import *


def get_candidates(CSMs):
	candidates = []
	r = 10000
	min_overlaps = 5
	for discs in CSMs:
		if len(discs) > 0:
			d_chr1,d_b1,d_hapblock1,d_hap1,d_chr2,d_b2,d_hapblock2,d_hap2 = discs[0]
			#overlaps = len(csmdict[(d_chr1,d_b1/r*r,d_chr2,d_b2/r*r)])
			#print (d_chr1,d_b1/r*r,d_chr2,d_b2/r*r),'discs',discs,'overlaps',overlaps
			if ((d_chr1 != d_chr2) or (abs(d_b2-d_b1) >= min_sv)):# and overlaps >= min_overlaps:
				candidates += [[d_chr1,d_b1,d_hapblock1,d_hap1,d_chr2,d_b2,d_hapblock2,d_hap2]]
	return candidates

def get_discs(lr1,b1,lr2,b2,xx_by_barcode):
	discs = []
	[chr1,hapblock1,hap1,num1,start1,end1,barcode1] = lr1
	if lr2:
		[chr2,hapblock2,hap2,num2,start2,end2,barcode2] = lr2
		for [d_chr1,d_b1,d_hapblock1,d_hap1,d_chr2,d_b2,d_hapblock2,d_hap2] in xx_by_barcode:
			if d_chr1 == chr1 and d_chr2 == chr2 and d_hapblock1 == hapblock1 and d_hapblock2 == hapblock2 and d_hap1 == hap1 and d_hap2 == hap2:
				if abs(b1-d_b1) <= d and abs(b2-d_b2) <= d:
					discs.append([d_chr1,d_b1,d_hapblock1,d_hap1,d_chr2,d_b2,d_hapblock2,d_hap2])
	else:
		for [d_chr1,d_b1,d_hapblock1,d_hap1,d_chr2,d_b2,d_hapblock2,d_hap2] in xx_by_barcode:
			if d_chr1 == chr1 and d_hapblock1 == hapblock1 and d_hap1 == hap1:
				if abs(b1-d_b1) <= d:
					discs.append([d_chr1,d_b1,d_hapblock1,d_hap1,d_chr2,d_b2,d_hapblock2,d_hap2])
	return discs

def get_candidate_split_mols(LRS,pp_by_barcode,pm_by_barcode,mp_by_barcode,mm_by_barcode):
	Mpp_cands = []
	Mpm_cands = []
	Mmp_cands = []
	Mmm_cands = []
	global starts
	global ends
	r = 10000
	s = time.time()
	for i in range(len(LRS)):
		[chr1,hapblock1,hap1,num1,start1,end1,barcode1] = LRS[i]
		lr1 = LRS[i]
		starts[(chr1,start1)].append(lr1)
		ends[(chr1,end1)].append(lr1)
		for j in range(i+1,len(LRS)):
			[chr2,hapblock2,hap2,num2,start2,end2,barcode2] = LRS[j]
			lr2 = LRS[j]
			if chr2 < chr1 or start2 < start1:
				[chr1,hapblock1,hap1,num1,start1,end1,barcode1] = LRS[j]
				[chr2,hapblock2,hap2,num2,start2,end2,barcode2] = LRS[i]	
				lr1 = LRS[j]
				lr2 = LRS[i]
			if lr1[1] != lr2[1] or lr1[2] == lr2[2]: # must be on same hap if on same hap block
				ppd = get_discs(lr1,lr1[-2],lr2,lr2[-2],pp_by_barcode)
				pmd = get_discs(lr1,lr1[-2],lr2,lr2[-3],pm_by_barcode)
				mpd = get_discs(lr1,lr1[-3],lr2,lr2[-2],mp_by_barcode)
				mmd = get_discs(lr1,lr1[-3],lr2,lr2[-3],mm_by_barcode)
				Mpp_cands.append(ppd)
				Mpm_cands.append(pmd)
				Mmp_cands.append(mpd)
				Mmm_cands.append(mmd)
	return Mpp_cands,Mpm_cands,Mmp_cands,Mmm_cands

def get_linked_reads(l):
	s = time.time()
	barcoded_reads,pp_by_barcode,pm_by_barcode,mp_by_barcode,mm_by_barcode,barcode = l
	barcoded_reads.sort(key = lambda x: (x[0],x[1]))
	LRS = []
	prev_end = 0
	prev_chr = 0
	curr_LR = []
	LENS = []
	RATES = []
	for read in barcoded_reads:
		[chr1,start,end,hapblock,hap] = read
		if chr1 == prev_chr and abs(start-prev_end) < d:
			curr_LR[-2] = end
			curr_LR[-4] += 1
		else:
			if curr_LR:
				if curr_LR[-2]-curr_LR[-3] > 1000:
					LENS.append(curr_LR[-2]-curr_LR[-3])
					RATES.append(curr_LR[-4]/float(curr_LR[-2]-curr_LR[-3]))
					LRS.append(curr_LR)
			curr_LR = [chr1,hapblock,hap,1,start,end,barcode]
		prev_end = end
		prev_chr = chr1
	if curr_LR[-2]-curr_LR[-3] >= 1000: #TODO remove constraint?
		LRS.append(curr_LR)
	Mpp,Mpm,Mmp,Mmm = get_candidate_split_mols(LRS,pp_by_barcode,pm_by_barcode,mp_by_barcode,mm_by_barcode)
	return [Mpp,Mpm,Mmp,Mmm,LENS,RATES]

def parallel_get_linked_reads(reads_by_barcode,pp_by_barcode,pm_by_barcode,mp_by_barcode,mm_by_barcode):
	global starts
	global ends
	starts = collections.defaultdict(list)
	ends = collections.defaultdict(list)
	print 'Getting linked-reads...'
	iters = []
	for barcode,item in reads_by_barcode.iteritems():
		iters.append([item,pp_by_barcode[barcode],pm_by_barcode[barcode],mp_by_barcode[barcode],mm_by_barcode[barcode],barcode])
	print len(iters),'barcodes'

	CSMs = map(get_linked_reads,iters)

	mm_cands = []
	mp_cands = []
	pm_cands = []
	pp_cands = []
	LENS = []
	RATES = []
	for [Mpp,Mpm,Mmp,Mmm,lens,rates] in CSMs:
		LENS += lens
		RATES += rates
		mm_cands += get_candidates(Mmm)
		mp_cands += get_candidates(Mmp)
		pm_cands += get_candidates(Mpm)
		pp_cands += get_candidates(Mpp)
	print 'Finished getting linked-reads'
	print 'len pm cands',len(pm_cands)
	return starts,ends,mm_cands,mp_cands,pm_cands,pp_cands,LENS,RATES