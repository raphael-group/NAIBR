import os,sys,subprocess,pysam,collections,time,gc,math,json,copy
import scipy.stats,operator,itertools
import multiprocessing as mp
from utils import *
import numpy as np
from global_vars import *
import linecache






class CandSplitMol:
	__slots__ = ['scoreHA','scoreH0','hap_i','hap_j','discs']
	def __init__(self):
		self.scoreHA = 0
		self.scoreH0 = 0
		self.hap_i = 0
		self.hap_j = 0
		self.discs = 0

	def score_CSM(self,CSM,candScore):
		LRs_i,D,LRs_j = CSM
		#chrm,start,end,num,hap,barcode
		chrm_i,si,ei,mapq_i,hap_i,barcode_i = LRs_i
		chrm_j,sj,ej,mapq_j,hap_j,barcode_j = LRs_j
		num_i,num_j = len(mapq_i),len(mapq_j)
		hap_i = max(hap_i)
		hap_j = max(hap_j)
		len0 = ej-si
		rate0 = (num_j+num_i)/float(len0)
		lenA = (ei-si)+(ej-sj)
		rateA = (num_j+num_i)/float(lenA)
		lenLi = ei-si
		rateLi0 = (num_i)/float(lenLi)
		rateLiA = (num_i)/float(lenLi)
		lenLj = ej-sj
		rateLj0 = (num_j)/float(lenLj)
		rateLjA = (num_j)/float(lenLj)
		mapqiA = self.map_prob(mapq_i)
		mapqi0 = self.map_prob(mapq_i)
		mapqjA = self.map_prob(mapq_j)
		mapqj0 = self.map_prob(mapq_j)
		mapqdiscA = self.map_prob(D)
		mapqdisc0 = self.mismap_prob(D)
		scoreHA = log(plen(lenA)*prate(rateA)*mapqiA*mapqjA*mapqdiscA+\
			plen(lenLi)*prate(rateLiA)*plen(lenLj)*prate(rateLjA)*mapqiA*mapqjA*mapqdiscA)
		if chrm_i == chrm_j:
			scoreH0 = log(plen(len0)*prate(rate0)*mapqi0*mapqj0*mapqdisc0+\
				plen(lenLi)*prate(rateLi0)*plen(lenLj)*prate(rateLj0)*mapqi0*mapqj0*mapqdisc0)
		else:
			scoreH0 = log(0+plen(lenLi)*prate(rateLi0)*plen(lenLj)*prate(rateLj0)*mapqi0*mapqj0*mapqdisc0)
		self.scoreHA = scoreHA
		self.scoreH0 = scoreH0
		self.hap_i = hap_i
		self.hap_j = hap_j
		self.discs = len(D)
		candScore.score[(hap_i,hap_j)] += self.scoreHA-self.scoreH0
		candScore.pairs[(hap_i,hap_j)] += 1
		candScore.disc[(hap_i,hap_j)] += self.discs

	def map_prob(self,mapqs):
		prob = 1
		for mapq in mapqs:
			prob *= 1-math.pow(10,(mapq/-10.0))
		return prob

	def mismap_prob(self,mapqs):
		prob = 1
		for mapq in mapqs:
			prob *= math.pow(10,(mapq/-10.0))
		return prob

def score_pair(i):
	key,CSMs,spans = get_CSMs(i)
	candidate = candidates[i]
	#CSMs,spans = CSMs_by_cand[(candidate.i,candidate.j,candidate.orient)]
	candScore = NA(candidate.chrm,candidate.nextchrm,candidate.i,candidate.j,candidate.orient)
	scoreHA = 0
	scoreH0 = 0
	scores = []
	for CSM in CSMs:
		c = CandSplitMol()
		c.score_CSM(CSM,candScore)
	for hap in spans:
		candScore.spans[hap] += 1
	return candScore,candidate


def near_i(x,candidate):
	chrm = x[0]
	start = x[1]
	end = x[2]
	neari = ((candidate.orient[0] == '+' and end-lmax < candidate.i and end+d > candidate.i)\
	     or (candidate.orient[0] == '-' and start+lmax > candidate.i and start-d < candidate.i)) and\
	 	     candidate.chrm == chrm 
	return neari

def near_j(x,candidate):
	chrm = x[0]
	start = x[1]
	end = x[2]
	nearj = ((candidate.orient[1] == '+' and end-lmax < candidate.j and end+d > candidate.j)\
	     or (candidate.orient[1] == '-' and start+lmax > candidate.j and start-d < candidate.j)) and\
	 	     candidate.nextchrm == chrm 
	return nearj

def spanning(x,candidate):
	start = x[1]
	end = x[2]
	span = start < candidate.i and end > candidate.j	 
	return span

def discs(candidate,barcode):
	ds = discs_by_barcode[(candidate.chrm,candidate.nextchrm,barcode)]
	# for d in ds:
	# 	if d.chrm != d.nextchrm:
	# 		print d.chrm,d.nextchrm,d.i,d.j
	ds = [(read.mapq+read.nextmapq)/2 for read in ds if \
	candidate.chrm == read.chrm and candidate.nextchrm == read.nextchrm\
	and is_close(candidate.i,read.i,read.orient[0]) and \
	    is_close(candidate.j,read.j,read.orient[1]) and candidate.orient == read.orient]
	return ds

def crossing(start,end,i):
	return start < i and end > i
	
def linked_reads(barcode,chrm,candidate):
	span = []
	reads = reads_by_LR[(chrm,barcode)]
	reads.sort(key = lambda x: x[0])
	#chrm,start,end,num,hap,barcode
	curr_LR = [0,0,0,[],[],0]
	LRs = []
	for start,end,hap,mapq in reads:
		if curr_LR[0] == 0 or start-curr_LR[2] > d:
			if curr_LR[0] != 0 and len(curr_LR[3]) >= min_reads and curr_LR[2]-curr_LR[1] >= min_len:
				LRs.append(curr_LR)
			curr_LR = [chrm,start,end,[mapq],[hap],barcode]
		elif candidate.chrm == candidate.nextchrm and curr_LR[0] != 0 and \
		(curr_LR[2] < candidate.i and start > candidate.j):
			LRs.append(curr_LR)
			curr_LR = [chrm,start,end,[mapq],[hap],barcode]
		else:
			curr_LR[2] = max(curr_LR[2],end)
			curr_LR[3].append(mapq)
			curr_LR[4].append(hap)
	if curr_LR[0] != 0 and len(curr_LR[3]) >= min_reads and curr_LR[2]-curr_LR[1] >= min_len:
		LRs.append(curr_LR)
	for i in range(1,len(LRs)):
		if LRs[i][1]-LRs[i-1][2] < d:
			span.append((max(LRs[i-1][-2]),max(LRs[i][-2])))
	return LRs,span

def get_LRs(candidate,barcodes):
	CSMs = []
	spans = []
	for barcode in barcodes:
		span = []
		LRs,s = linked_reads(barcode,candidate.chrm,candidate)
		span += s
		if candidate.chrm != candidate.nextchrm:
			LRs2,s = linked_reads(barcode,candidate.nextchrm,candidate)
			LRs += LRs2
		LRs_i,LRs_j = [[],[]]
		for x in LRs:
			if near_i(x,candidate):
				LRs_i = x
			elif near_j(x,candidate):
				LRs_j = x
			if spanning(x,candidate):
				span += [(x[-2][0],x[-2][-1])]

		if LRs_i and LRs_j and LRs_i[1] < LRs_j[1]:
			D = discs(candidate,barcode)
			CSM = [LRs_i,D,LRs_j]
			CSMs.append(CSM)		
		spans += span
	return CSMs,spans


def get_CSMs(i):
	candidate = candidates[i]
	break1_barcodes = set(LRs_by_pos[(candidate.chrm,candidate.i/R*R)]+\
		LRs_by_pos[(candidate.chrm,candidate.i/R*R-R)]+LRs_by_pos[(candidate.chrm,candidate.i/R*R+R)])
	break2_barcodes = set(list(LRs_by_pos[(candidate.nextchrm,candidate.j/R*R)])+\
		LRs_by_pos[(candidate.nextchrm,candidate.j/R*R-R)]+LRs_by_pos[(candidate.nextchrm,candidate.j/R*R+R)])
	CSMs,spans = get_LRs(candidate,break1_barcodes.intersection(break2_barcodes))
	key = (candidate.i,candidate.j,candidate.orient)
	return key,CSMs,spans

def parse_CSMs(all_CSMs):
	global CSMs_by_cand
	for key,CSMs,spans in all_CSMs:
		CSMs_by_cand[key] = [CSMs,spans]


def get_all_CSMs():
	global CSMs_by_cand
	CSMs_by_cand = dict()
	if not is_interchrom or NUM_THREADS == 1:
		all_CSMs = []
		for i in range(len(candidates)):
			key,CSMs,spans = get_CSMs(i)
			CSMs_by_cand[key] = [CSMs,spans]
	elif is_interchrom and NUM_THREADS != 1:
		pool = mp.Pool(min(5,NUM_THREADS),maxtasksperchild=1)
		map_fn = pool.map
		all_CSMs = map_fn(get_CSMs,range(len(candidates)))
		parse_CSMs(all_CSMs)
		pool.close()
		pool.join()


def get_cand_score(cands):
	global candidates
	candidates = cands

	if not is_interchrom or NUM_THREADS == 1:
		map_fn = map
	elif is_interchrom and NUM_THREADS != 1:
		pool = mp.Pool(min(5,NUM_THREADS),maxtasksperchild=1)
		map_fn = pool.map
	scores = map_fn(score_pair,range(len(candidates)))
	if is_interchrom and NUM_THREADS != 1:
		pool.close()
		pool.join()

	rets = []
	for best,candidate in scores:
		best.get_score()
		if best.pairs > 0 and best.score > -float('inf'):
			ret = [best.chrm1,best.break1,best.chrm2,best.break2]+['split_mols:'+str(best.pairs)]+['discordants:'+str(best.disc)]+\
			['orient:'+candidate.orient,'haps:'+str(best.haps[0])+','+str(best.haps[1]),'score:'+str(best.score)]
			rets.append(ret)
	return rets






def predict_NAs(reads_by_lr,lrs_by_pos,Discs_by_barcode,candidates,p_len,p_rate,cov,interchrom):
	global lmax,plen,prate,discs_by_barcode,LRs_by_pos,chrom_idx,is_interchrom
	global LRs_by_pos,reads_by_LR,discs_by_barcode
	LRs_by_pos,reads_by_LR,discs_by_barcode = lrs_by_pos,reads_by_lr,Discs_by_barcode
	is_interchrom = interchrom
	plen = p_len
	prate = p_rate
	scores = []

	scores += get_cand_score(candidates)
	scores = [x for x in scores if x]
	scores = collapse(scores,threshold(cov))
	return scores



