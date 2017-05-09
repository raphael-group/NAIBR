import os,sys,subprocess,pysam,collections,time,gc,json,random
import multiprocessing as mp
from utils import *
from global_vars import *

def is_discordant(read,mate,lmin,lmax):
	if get_hap(read) == 0 or get_hap(mate) == 0: # TODO: use non-haplotyped reads
		return 0
	global barcode_count
	global reads_by_barcode
	global barcode_to_num
	global pp_by_barcode
	global pm_by_barcode
	global mp_by_barcode
	global mm_by_barcode
	chr1 = to_chr(read.reference_name)
	chr2 = to_chr(mate.reference_name)
	if chr1 != chr2:
		if chr2 < chr1:
			m = mate
			r = read
			mate = r
			read = m
			chr1 = to_chr(read.reference_name)
			chr2 = to_chr(mate.reference_name)

	if read.get_tag('RX') in barcode_to_num:
		barcode = barcode_to_num[read.get_tag('RX')]
	else:
		barcode_count += 1
		barcode = barcode_count
		barcode_to_num[read.get_tag('RX')] = barcode_count
	if not is_convergent(read,mate):
		if not read.is_reverse and not mate.is_reverse and (chr1 != chr2 or mate.reference_end-read.reference_start >= min_sv):
			pp_by_barcode[barcode] += [[chr1,read.reference_end,read.get_tag('PS'),get_hap(read),chr2,mate.reference_end,mate.get_tag('PS'),get_hap(mate)]]
		elif read.is_reverse and mate.is_reverse and (chr1 != chr2 or mate.reference_end-read.reference_start >= min_sv):
			mm_by_barcode[barcode] += [[chr1,read.reference_start,read.get_tag('PS'),get_hap(read),chr2,mate.reference_start,mate.get_tag('PS'),get_hap(mate)]]
		elif (chr1 != chr2 or mate.reference_end-read.reference_start >= min_sv):
			mp_by_barcode[barcode] += [[chr1,read.reference_start,read.get_tag('PS'),get_hap(read),chr2,mate.reference_end,mate.get_tag('PS'),get_hap(mate)]]
	elif mate.reference_end-read.reference_start > lmax and (chr1 != chr2 or mate.reference_end-read.reference_start >= min_sv):
		pm_by_barcode[barcode] += [[chr1,read.reference_end,read.get_tag('PS'),get_hap(read),chr2,mate.reference_start,mate.get_tag('PS'),get_hap(mate)]]
	else:
		reads_by_barcode[barcode] += [[chr1,read.reference_start,max(read.reference_end,mate.reference_end),read.get_tag('PS'),get_hap(read)]]
		return 0
	return 1

def make_barcodeDict(lmin,lmax,CHR):
	global pp_by_barcode
	global pm_by_barcode
	global mp_by_barcode
	global mm_by_barcode
	global reads_by_barcode
	global barcode_to_num
	global barcode_count
	barcode_count = 0
	num_disc = 0
	print 'Processing reads...'
	num_reads = []
	starttime = time.time()
	prevtime = starttime
	reads = pysam.AlignmentFile(BAM_FILE,"rb")
	count = 0
	lasttime = time.time()
	reads_by_barcode = collections.defaultdict(list)
	pp_by_barcode = collections.defaultdict(list)
	pm_by_barcode = collections.defaultdict(list)
	mp_by_barcode = collections.defaultdict(list)
	mm_by_barcode = collections.defaultdict(list)
	barcode_to_num = dict()
	mate_pairs = dict()
	totallen = 0
	if 'chr' in CHR:
		iterator = reads.fetch(CHR)
	else:
		print 'reading until end of file...'
		iterator = reads.fetch(until_eof=True)
	for read in iterator:
		count +=1
		if count%10000000 == 0:
			print count, time.time()-prevtime
			prevtime = time.time()
			if DEBUG:
		 		break
		#TODO: set to hap !=' '0' if require haplotyped
		disc = 0
		if read.mapping_quality >= MIN_MAPQ and not read.is_duplicate and not read.is_secondary and read.is_paired:
			if read.query_name in mate_pairs: # is read 2
				disc = is_discordant(mate_pairs[read.query_name],read,lmin,lmax)
			else:
				mate_pairs[read.query_name] = read # is read 1
			num_disc += disc
	print 'Processed reads in',(time.time()-starttime)/60.0,'minutes'
	print num_disc,'discordant reads identified'
	print 'total num reads',len(mate_pairs)
	mate_pairs = None
	return reads_by_barcode,pp_by_barcode,pm_by_barcode,mp_by_barcode,mm_by_barcode