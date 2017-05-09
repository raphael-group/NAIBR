import os,sys,subprocess,pysam,collections,time,gc,math,json,copy
import scipy.stats
import multiprocessing as multp
from utils import *
import numpy as np
from global_vars import *
import linecache

global pm1
global pm2
def get_key(d_hap1,d_hap2,hap):
	if d_hap1 == d_hap2:
		if hap == 1:
			return 1
		else:
			return 2
	else:
		if hap == 1:
			return 2
		else:
			return 1

def score_pair(pair,t):
	(chr1,hap1,num1,len1,lenover1,b1,discs,chr2,hap2,num2,len2,lenover2,b2) = pair
	num1H0 = num1-discs
	num2H0 = num2-discs
	rate = (num1+num2)/float(len1+len2)
	numover1 = min(5,int(math.ceil(lenover1*rate)))
	numover2 = min(5,int(math.ceil(lenover2*rate)))
	num1 = num1-numover1
	num2 = num2-numover2
	len1 = len1-lenover1
	len2 = len2-lenover2
	scoreHA = 0
	scoreH0 = 0
	if len1 and len2:
		len1H0 = len1+lenover1
		len2H0 = lenover2+len2
		# if t != '+-':
		# 	scoreHA += plog(1e-6,min(1,numover1)+min(1,numover2))
		# else:
		scoreHA += plog(1e-6,numover1+numover2)
		scoreH0 += plog(1e-6,discs)
		
		if chr1 == chr2 and b1-b2 < 50000:
			scoreHA += log(1e-6)
			scoreH0 += log(1e-7)
		else:
			scoreHA += log(1e-6)
			scoreH0 += log(1e-16)
		scoreHA += log(p_len(len1+len2)*p_rate((num1+num2)/float(len1+len2))+p_len(len1)*p_rate(num1/float(len1))*p_len(len2)*p_rate(num2/float(len2)))
		scoreH0 += log(p_len(len1H0+(b2-b1)+len2H0)*p_rate((num1H0+num2H0)/float(len1H0+(b2-b1)+len2H0))+p_len(len1H0)*p_rate(num1H0/float(len1H0))*p_len(len2H0)*p_rate(num2H0/float(len2H0)))
	else:
		return 0,0,0
	return numover1+numover2,discs,scoreHA-scoreH0

def score_single(single,t):
	(chr1,hap1,num1,len1,lenover1,b1) = single
	scoreHA = 0
	scoreH0 = 0
	if lenover1 == 0:
		return 0,0,scoreHA-scoreH0
	rate = (num1)/float(len1)
	numover1 = min(5,int(math.ceil(lenover1*rate)))
	scoreHA += plog(1e-6,numover1)
	scoreH0 += 0
	return numover1,0,scoreHA-scoreH0


def rank_mm(c):
	ret = []
	candidate = c
	[d_chr1,d_b1,d_hapblock1,d_hap1,d_chr2,d_b2,d_hapblock2,d_hap2] = candidate
	ret_score = -float('inf')
	for adj in range(0,lmax/2,100):
		d_b1 = d_b1-adj
		d_b2 = d_b2-adj
		scores = collections.defaultdict(int)
		overs = collections.defaultdict(int)
		discordants = collections.defaultdict(int)
		num_pairs = collections.defaultdict(int)
		num_singles = collections.defaultdict(int)
		CSMs1 = collections.defaultdict(list)
		CSMs2 = collections.defaultdict(list)
		for i in range(d_b1-d/2,min(d_b2,d_b1+d)):
			lrs = starts[(d_chr1,i)]
			for (chr1,hapblock1,hap1,num1,start1,end1,barcode1) in lrs:
				if end1 > d_b1:
					if barcode1 not in CSMs1:
						CSMs1[barcode1] = [chr1,hapblock1,hap1,num1,start1,end1]
					else:
						chrp,hapblockp,happ,nump,startp,endp = CSMs1[barcode1]
						CSMs1[barcode1] = [chr1,hapblock1,hap1,nump+num1,min(startp,start1),max(endp,end1)]
		for i in range(max(d_b1,d_b2-d/2),d_b2+d):
			lrs = starts[(d_chr2,i)]
			for (chr1,hapblock1,hap1,num1,start1,end1,barcode1) in lrs:
				if end1 > d_b2:
					if barcode1 not in CSMs2:
						CSMs2[barcode1] = [chr1,hapblock1,hap1,num1,start1,end1]
					else:
						chrp,hapblockp,happ,nump,startp,endp = CSMs2[barcode1]
						CSMs2[barcode1] = [chr1,hapblock1,hap1,nump+num1,min(startp,start1),max(endp,end1)]
		used_barcodes = set()
		for barcode,(chr1,hapblock1,hap1,num1,start1,end1) in CSMs1.iteritems():
			discs = get_discs(d_chr1,d_b1,d_hapblock1,d_hap1,d_chr2,d_b2,d_hapblock2,d_hap2,mm[barcode])
			used_barcodes.add(barcode)
			lenover1 = max(0,d_b1-start1)
			if barcode in CSMs2:
				(chr2,hapblock2,hap2,num2,start2,end2) = CSMs2[barcode]
				num_pairs[(hap1,hap2)] +=1
				lenover2 = max(0,d_b2-start2)
				over,disc,score = score_pair((chr1,hap1,num1,end1-start1,lenover1,start1,discs,chr2,hap2,num2,end2-start2,lenover2,start2),'--')
				scores[(hap1,hap2)] += score
				overs[(hap1,hap2)] += over
				discordants[(hap1,hap2)] += disc
			else:
				if lenover1:
					num_singles[(hap1,get_key(d_hap1,d_hap2,hap1))] +=1
				over,disc,score = score_single((chr1,hap1,num1,end1-start1,lenover1,start1),'--')
				scores[(hap1,get_key(d_hap1,d_hap2,hap1))] += score
				overs[(hap1,get_key(d_hap1,d_hap2,hap1))]  += over
				discordants[(hap1,get_key(d_hap1,d_hap2,hap1))] += disc	
		for barcode,(chr2,hapblock2,hap2,num2,start2,end2) in CSMs2.iteritems():
			if barcode not in used_barcodes:
				lenover2 = max(0,d_b2-start2)
				if lenover2:
					num_singles[(get_key(d_hap1,d_hap2,hap2),hap2)] +=1
				over,disc,score = score_single((chr2,hap2,num2,end2-start2,lenover2,start2),'--')
				scores[(get_key(d_hap1,d_hap2,hap2),hap2)] += score
				overs[(get_key(d_hap1,d_hap2,hap2),hap2)]  += over
				discordants[(get_key(d_hap1,d_hap2,hap2),hap2)] +=  disc
		best_haps = ()
		best_score = -float('inf')
		hom = 0
		for haps,score in scores.iteritems():
			num_pairs[(3,3)] += num_pairs[haps]
			num_singles[(3,3)] += num_singles[haps]
			overs[(3,3)] += overs[haps]
			discordants[(3,3)] += discordants[haps]
			if score > best_score:
				best_score = score
				best_haps = haps
			hom += score
		if best_score > -float('inf') and hom > best_score:
			best_score = hom
			best_haps = (3,3) 
		if best_score > ret_score and num_pairs[best_haps] > MIN_SPLIT:
			ret_score = best_score
			ret = [d_chr1,d_b1,d_chr2,d_b2]+['split_mols:'+str(num_pairs[best_haps]),'single_mols:'+str(num_singles[best_haps])]+['opposing:'+str(overs[best_haps])]+['discordants:'+str(discordants[best_haps])]+\
			['orient:--','haps:'+str(best_haps[0])+','+str(best_haps[1]),'score:'+str(best_score)]
	return ret

def rank_pp(c):
	ret = []
	candidate = c
	[d_chr1,d_b1,d_hapblock1,d_hap1,d_chr2,d_b2,d_hapblock2,d_hap2] = candidate
	ret_score = -float('inf')
	for adj in range(0,lmax/2,100):
		d_b1 = d_b1+adj
		d_b2 = d_b2+adj
		scores = collections.defaultdict(int)
		overs = collections.defaultdict(int)
		discordants = collections.defaultdict(int)
		num_pairs = collections.defaultdict(int)
		num_singles = collections.defaultdict(int)
		CSMs1 = collections.defaultdict(list)
		CSMs2 = collections.defaultdict(list)
		for i in range(d_b1-d,min(d_b2,d_b1+d/2)):
			lrs = ends[(d_chr1,i)]
			for (chr1,hapblock1,hap1,num1,start1,end1,barcode1) in lrs:
				if start1 < d_b1:
					if barcode1 not in CSMs1:
						CSMs1[barcode1] = [chr1,hapblock1,hap1,num1,start1,end1]
					else:
						chrp,hapblockp,happ,nump,startp,endp = CSMs1[barcode1]
						CSMs1[barcode1] = [chr1,hapblock1,hap1,nump+num1,min(startp,start1),max(endp,end1)]
		for i in range(max(d_b1,d_b2-d),d_b2+d/2):
			lrs = ends[(d_chr2,i)]
			for (chr1,hapblock1,hap1,num1,start1,end1,barcode1) in lrs:
				if start1 < d_b2:
					if barcode1 not in CSMs2:
						CSMs2[barcode1] = [chr1,hapblock1,hap1,num1,start1,end1]
					else:
						chrp,hapblockp,happ,nump,startp,endp = CSMs2[barcode1]
						CSMs2[barcode1] = [chr1,hapblock1,hap1,nump+num1,min(startp,start1),max(endp,end1)]
		used_barcodes = set()
		for barcode,(chr1,hapblock1,hap1,num1,start1,end1) in CSMs1.iteritems():
			discs = get_discs(d_chr1,d_b1,d_hapblock1,d_hap1,d_chr2,d_b2,d_hapblock2,d_hap2,pp[barcode])
			used_barcodes.add(barcode)
			lenover1 = max(0,end1-d_b1)
			if barcode in CSMs2:
				(chr2,hapblock2,hap2,num2,start2,end2) = CSMs2[barcode]
				num_pairs[(hap1,hap2)] +=1
				lenover2 = max(0,end2-d_b2)
				over,disc,score = score_pair((chr1,hap1,num1,end1-start1,lenover1,end1,discs,chr2,hap2,num2,end2-start2,lenover2,end2),'++')
				scores[(hap1,hap2)] += score
				overs[(hap1,hap2)] += over
				discordants[(hap1,hap2)] += disc
			else:
				if lenover1:
					num_singles[(hap1,get_key(d_hap1,d_hap2,hap1))] +=1
				over,disc,score = score_single((chr1,hap1,num1,end1-start1,lenover1,end1),'++')
				scores[(hap1,get_key(d_hap1,d_hap2,hap1))] += score
				overs[(hap1,get_key(d_hap1,d_hap2,hap1))]  += over
				discordants[(hap1,get_key(d_hap1,d_hap2,hap1))] += disc	
		for barcode,(chr2,hapblock2,hap2,num2,start2,end2) in CSMs2.iteritems():
			if barcode not in used_barcodes:
				lenover2 = max(0,end2-d_b2)
				if lenover2:
					num_singles[(get_key(d_hap1,d_hap2,hap2),hap2)] +=1
				over,disc,score = score_single((chr2,hap2,num2,end2-start2,lenover2,end2),'++')
				scores[(get_key(d_hap1,d_hap2,hap2),hap2)] += score
				overs[(get_key(d_hap1,d_hap2,hap2),hap2)]  += over
				discordants[(get_key(d_hap1,d_hap2,hap2),hap2)] +=  disc
		best_haps = ()
		best_score = -float('inf')
		hom = 0
		for haps,score in scores.iteritems():
			num_pairs[(3,3)] += num_pairs[haps]
			num_singles[(3,3)] += num_singles[haps]
			overs[(3,3)] += overs[haps]
			discordants[(3,3)] += discordants[haps]
			if score > best_score:
				best_score = score
				best_haps = haps
			hom += score
		if best_score > -float('inf') and hom > best_score:
			best_score = hom
			best_haps = (3,3) 
		if best_score > ret_score and num_pairs[best_haps] > MIN_SPLIT:
			ret_score = best_score
			ret = [d_chr1,d_b1,d_chr2,d_b2]+['split_mols:'+str(num_pairs[best_haps]),'single_mols:'+str(num_singles[best_haps])]+['opposing:'+str(overs[best_haps])]+['discordants:'+str(discordants[best_haps])]+\
			['orient:++','haps:'+str(best_haps[0])+','+str(best_haps[1]),'score:'+str(best_score)]
	return ret

def rank_mp(c):
	ret = []
	candidate = c
	[d_chr1,d_b1,d_hapblock1,d_hap1,d_chr2,d_b2,d_hapblock2,d_hap2] = candidate
	ret_score = -float('inf')
	for adj in range(0,lmax/2,100):
		d_b1 = d_b1-adj
		d_b2 = d_b2+adj
		scores = collections.defaultdict(int)
		overs = collections.defaultdict(int)
		discordants = collections.defaultdict(int)
		num_pairs = collections.defaultdict(int)
		num_singles = collections.defaultdict(int)
		CSMs1 = collections.defaultdict(list)
		CSMs2 = collections.defaultdict(list)
		for i in range(d_b1-d/2,min(d_b2,d_b1+d)):
			lrs = starts[(d_chr1,i)]
			for (chr1,hapblock1,hap1,num1,start1,end1,barcode1) in lrs:
				if end1 > d_b1:
					if barcode1 not in CSMs1:
						CSMs1[barcode1] = [chr1,hapblock1,hap1,num1,start1,end1]
					else:
						chrp,hapblockp,happ,nump,startp,endp = CSMs1[barcode1]
						CSMs1[barcode1] = [chr1,hapblock1,hap1,nump+num1,min(startp,start1),max(endp,end1)]
		for i in range(max(d_b1,d_b2-d),d_b2+d/2):
			lrs = ends[(d_chr2,i)]
			for (chr1,hapblock1,hap1,num1,start1,end1,barcode1) in lrs:
				if start1 < d_b2:
					if barcode1 not in CSMs2:
						CSMs2[barcode1] = [chr1,hapblock1,hap1,num1,start1,end1]
					else:
						chrp,hapblockp,happ,nump,startp,endp = CSMs2[barcode1]
						CSMs2[barcode1] = [chr1,hapblock1,hap1,nump+num1,min(startp,start1),max(endp,end1)]
		used_barcodes = set()
		for barcode,(chr1,hapblock1,hap1,num1,start1,end1) in CSMs1.iteritems():
			discs = get_discs(d_chr1,d_b1,d_hapblock1,d_hap1,d_chr2,d_b2,d_hapblock2,d_hap2,mp[barcode])
			used_barcodes.add(barcode)
			lenover1 = max(0,d_b1-start1)
			if barcode in CSMs2:
				(chr2,hapblock2,hap2,num2,start2,end2) = CSMs2[barcode]
				num_pairs[(hap1,hap2)] +=1
				lenover2 = max(0,end2-d_b2)
				over,disc,score = score_pair((chr1,hap1,num1,end1-start1,lenover1,start1,discs,chr2,hap2,num2,end2-start2,lenover2,end2),'-+')
				scores[(hap1,hap2)] += score
				overs[(hap1,hap2)] += over
				discordants[(hap1,hap2)] += disc
			else:
				if lenover1:
					num_singles[(hap1,get_key(d_hap1,d_hap2,hap1))] +=1
				over,disc,score = score_single((chr1,hap1,num1,end1-start1,lenover1,start1),'-+')
				scores[(hap1,get_key(d_hap1,d_hap2,hap1))] += score
				overs[(hap1,get_key(d_hap1,d_hap2,hap1))]  += over
				discordants[(hap1,get_key(d_hap1,d_hap2,hap1))] += disc	
		for barcode,(chr2,hapblock2,hap2,num2,start2,end2) in CSMs2.iteritems():
			if barcode not in used_barcodes:
				lenover2 = max(0,d_b2-start2)
				if lenover2:
					num_singles[(get_key(d_hap1,d_hap2,hap2),hap2)] +=1
				over,disc,score = score_single((chr2,hap2,num2,end2-start2,lenover2,end2),'-+')
				scores[(get_key(d_hap1,d_hap2,hap2),hap2)] += score
				overs[(get_key(d_hap1,d_hap2,hap2),hap2)]  += over
				discordants[(get_key(d_hap1,d_hap2,hap2),hap2)] +=  disc
		best_haps = ()
		best_score = -float('inf')
		hom = 0
		for haps,score in scores.iteritems():
			num_pairs[(3,3)] += num_pairs[haps]
			num_singles[(3,3)] += num_singles[haps]
			overs[(3,3)] += overs[haps]
			discordants[(3,3)] += discordants[haps]
			if score > best_score:
				best_score = score
				best_haps = haps
			hom += score
		if best_score > -float('inf') and hom > best_score:
			best_score = hom
			best_haps = (3,3) 
		if best_score > ret_score and num_pairs[best_haps] > MIN_SPLIT:
			ret_score = best_score
			ret = [d_chr1,d_b1,d_chr2,d_b2]+['split_mols:'+str(num_pairs[best_haps]),'single_mols:'+str(num_singles[best_haps])]+['opposing:'+str(overs[best_haps])]+['discordants:'+str(discordants[best_haps])]+\
			['orient:-+','haps:'+str(best_haps[0])+','+str(best_haps[1]),'score:'+str(best_score)]
	return ret


def get_discs(d_chr1,d_b1,d_hapblock1,d_hap1,d_chr2,d_b2,d_hapblock2,d_hap2,discs):
	ret = 0
	for [chr1,b1,hapblock1,hap1,chr2,b2,hapblock2,hap2] in discs:
		if chr1 == d_chr1 and abs(b1-d_b1) < 1000 and d_hapblock1 == hapblock1 and d_hap1 == hap1\
		 and d_chr2 == chr2 and abs(b2-d_b2) < 1000 and d_hapblock2 == hapblock2 and d_hap2 == hap2:
		 ret += 1
	return ret 

def rank_pm(c):
	ret = []
	candidate = c
	[d_chr1,d_b1,d_hapblock1,d_hap1,d_chr2,d_b2,d_hapblock2,d_hap2] = candidate
	ret_score = -float('inf')
	for adj in range(0,lmax/2,100):
		d_b1 = d_b1+adj
		d_b2 = d_b2-adj
		scores = collections.defaultdict(int)
		overs = collections.defaultdict(int)
		discordants = collections.defaultdict(int)
		num_pairs = collections.defaultdict(int)
		num_singles = collections.defaultdict(int)
		CSMs1 = collections.defaultdict(list)
		CSMs2 = collections.defaultdict(list)
		for i in range(d_b1-d,min(d_b2,d_b1+d/2)):
			lrs = ends[(d_chr1,i)]
			for (chr1,hapblock1,hap1,num1,start1,end1,barcode1) in lrs:
				if start1 < d_b1:
					if barcode1 not in CSMs1:
						CSMs1[barcode1] = [chr1,hapblock1,hap1,num1,start1,end1]
					else:
						chrp,hapblockp,happ,nump,startp,endp = CSMs1[barcode1]
						CSMs1[barcode1] = [chr1,hapblock1,hap1,nump+num1,min(startp,start1),max(endp,end1)]
		for i in range(max(d_b1,d_b2-d/2),d_b2+d):
			lrs = starts[(d_chr2,i)]
			for (chr1,hapblock1,hap1,num1,start1,end1,barcode1) in lrs:
				if end1 > d_b2:
					if barcode1 not in CSMs2:
						CSMs2[barcode1] = [chr1,hapblock1,hap1,num1,start1,end1]
					else:
						chrp,hapblockp,happ,nump,startp,endp = CSMs2[barcode1]
						CSMs2[barcode1] = [chr1,hapblock1,hap1,nump+num1,min(startp,start1),max(endp,end1)]
		used_barcodes = set()
		for barcode,(chr1,hapblock1,hap1,num1,start1,end1) in CSMs1.iteritems():
			discs = get_discs(d_chr1,d_b1,d_hapblock1,d_hap1,d_chr2,d_b2,d_hapblock2,d_hap2,pm[barcode])
			used_barcodes.add(barcode)
			lenover1 = max(0,end1-d_b1)
			if barcode in CSMs2:
				(chr2,hapblock2,hap2,num2,start2,end2) = CSMs2[barcode]
				num_pairs[(hap1,hap2)] +=1
				lenover2 = max(0,d_b2-start2)
				over,disc,score = score_pair((chr1,hap1,num1,end1-start1,lenover1,end1,discs,chr2,hap2,num2,end2-start2,lenover2,start2),'+-')
				scores[(hap1,hap2)] += score
				overs[(hap1,hap2)] += over
				discordants[(hap1,hap2)] += disc
			else:
				if lenover1:
					num_singles[(hap1,get_key(d_hap1,d_hap2,hap1))] +=1
				over,disc,score = score_single((chr1,hap1,num1,end1-start1,lenover1,end1),'+-')
				scores[(hap1,get_key(d_hap1,d_hap2,hap1))] += score
				overs[(hap1,get_key(d_hap1,d_hap2,hap1))]  += over
				discordants[(hap1,get_key(d_hap1,d_hap2,hap1))] += disc	
		for barcode,(chr2,hapblock2,hap2,num2,start2,end2) in CSMs2.iteritems():
			if barcode not in used_barcodes:
				lenover2 = max(0,d_b2-start2)
				if lenover2:
					num_singles[(get_key(d_hap1,d_hap2,hap2),hap2)] +=1
				over,disc,score = score_single((chr2,hap2,num2,end2-start2,lenover2,start2),'+-')
				scores[(get_key(d_hap1,d_hap2,hap2),hap2)] += score
				overs[(get_key(d_hap1,d_hap2,hap2),hap2)]  += over
				discordants[(get_key(d_hap1,d_hap2,hap2),hap2)] +=  disc
		best_haps = ()
		best_score = -float('inf')
		hom = 0
		for haps,score in scores.iteritems():
			num_pairs[(3,3)] += num_pairs[haps]
			num_singles[(3,3)] += num_singles[haps]
			overs[(3,3)] += overs[haps]
			discordants[(3,3)] += discordants[haps]
			if score > best_score:
				best_score = score
				best_haps = haps
			hom += score
		if best_score > -float('inf') and hom > best_score:
			best_score = hom
			best_haps = (3,3) 
		if best_score > ret_score and num_pairs[best_haps] > MIN_SPLIT:
			ret_score = best_score
			ret = [d_chr1,d_b1,d_chr2,d_b2]+['split_mols:'+str(num_pairs[best_haps]),'single_mols:'+str(num_singles[best_haps])]+['opposing:'+str(overs[best_haps])]+['discordants:'+str(discordants[best_haps])]+\
			['orient:+-','haps:'+str(best_haps[0])+','+str(best_haps[1]),'score:'+str(best_score)]
	return ret

def parallel_scores(func,iters):
	print len(iters),'discordants'
	if NUM_CORES != 1:
		pool = multp.Pool(NUM_CORES,maxtasksperchild=1)
		map_fn = pool.map
	else:
		map_fn = map
	scores = map_fn(func,iters)
	if NUM_CORES != 1:
		pool.close()
		pool.join()
	return [x for x in scores if x and x[-1] > -float('inf')]

def collapse_scores(scores):
	r = 2000
	print len(scores)
	scores.sort(key = lambda x: (int(x[0]),int(x[2]),int(x[1])/r*r,int(x[3])/r*r,float(x[-1].split(':')[-1])),reverse=True)
	prev_chr1 = 0
	prev_chr2 = 0
	prev_s = 0
	prev_e = 0
	scores2 = []
	for d in scores:
		if prev_chr1 == d[0] and prev_chr2 == d[2] and abs(prev_s-d[1]) <= r and abs(prev_e-d[3]) <= r:
			None
		else:
			scores2.append(d)
		prev_chr1 = d[0]
		prev_chr2 = d[2]
		prev_s = d[1]
		prev_e = d[3]
	scores2.sort(key = lambda x: float(x[-1].split(':')[-1]),reverse=True)
	return scores2

def parallel_rank_cands(pp_by_barcode,pm_by_barcode,mp_by_barcode,mm_by_barcode,s,e,mm_cands,mp_cands,pm_cands,pp_cands,plen,prate):
	global p_len
	global p_rate
	global starts
	global ends
	global pp
	global pm
	global mp
	global mm
	pp = pp_by_barcode
	pm = pm_by_barcode
	mp = mp_by_barcode
	mm = mm_by_barcode
	starts = s
	ends = e
	p_len = plen
	p_rate = prate
	scores = parallel_scores(rank_pm,pm_cands)
	scores += parallel_scores(rank_pp,pp_cands)
	scores += parallel_scores(rank_mp,mp_cands)
	scores += parallel_scores(rank_mm,mm_cands)
	return scores

def predict_NAs(pp_by_barcode,pm_by_barcode,mp_by_barcode,mm_by_barcode,s,e,mm_cands,mp_cands,pm_cands,pp_cands,plen,prate,l):
	global lmax
	lmax = l
	scores = parallel_rank_cands(pp_by_barcode,pm_by_barcode,mp_by_barcode,mm_by_barcode,s,e,mm_cands,mp_cands,pm_cands,pp_cands,plen,prate)
	scores = [x for x in scores if len(x) > 0]
	fname = 'LRSV_SVS.bedpe'
	print os.path.join(DIR,fname)
	f = open(os.path.join(DIR,fname),'w')
	scores = collapse_scores(scores)
	for i,score in enumerate(scores[:100]):
		print i,score
	for result in scores:
		f.write('\t'.join([str(x) for x in result])+'\n')
	f.close()
	return