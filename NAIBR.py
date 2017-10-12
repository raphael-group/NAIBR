import os,sys,time,json,pysam,collections
sys.path.insert(0, 'src')
if len(sys.argv) < 2:
	raise ValueError('No input config file') 
import multiprocessing as mp
from get_reads import *
from global_vars import *
from estimate_params import estimate_lmin_lmax
from utils import *
from rank import *



def run_NAIBR_user(cand):
	'''
	use user input candidate novel adjacencies
	'''
	scores = 0
	reads_by_LR,LRs_by_pos,discs_by_barcode,cands,coverage = make_barcodeDict_user(cand)
	p_len,p_rate,overlap = get_distributions(reads_by_LR)
	if p_len == None:
		return scores
	scores = predict_NAs(reads_by_LR,LRs_by_pos,discs_by_barcode,cands,p_len,p_rate,coverage,False)
	return scores

def run_NAIBR(chrom):
	'''
	automatically identify candidate novel adjacencies
	'''
	reads_by_LR,LRs_by_pos,discs_by_barcode,discs,interchrom_discs,coverage = make_barcodeDict(chrom)
	scores = 0
	if len(reads_by_LR) > 0:
		cands,p_len,p_rate = get_candidates(discs,reads_by_LR)
		if cands == None:
			print 'No candidates from',chrom
			return reads_by_LR,LRs_by_pos,discs_by_barcode,interchrom_discs,coverage,scores
		print 'ranking',len(cands),'candidates from',chrom
		scores = predict_NAs(reads_by_LR,LRs_by_pos,discs_by_barcode,cands,p_len,p_rate,coverage,False)
	else:
		print 'No candidates from',chrom
	return reads_by_LR,LRs_by_pos,discs_by_barcode,interchrom_discs,coverage,scores

def main():
	starttime = time.time()
	if len(candidates) > 0:
		print 'user defined candidates'
		with open(candidates) as f:
			cands = f.read().split('\n')
			cands = [x.split('\t') for x in cands]
			cands = [x for x in cands if x]
			cands = [[x[0],int(x[1]),x[3],int(x[4]),x[-1]] for x in cands if len(x) >= 4]
		scores = flatten(parallel_execute(run_NAIBR_user,cands))
		write_scores(scores)

	else:
		reads = pysam.AlignmentFile(BAM_FILE,"rb")
		chroms = reads.references
		chroms = [x for x in chroms if is_proper_chrom(x)]
		data = parallel_execute(run_NAIBR,chroms)
		reads_by_LR = collections.defaultdict(list)
		LRs_by_pos = collections.defaultdict(list)
		discs_by_barcode = collections.defaultdict(list)
		discs = collections.defaultdict(list)
		coverage = []
		scores = []
		for reads_by_LR_chrom,LRs_by_pos_chrom,discs_by_barcode_chrom,discs_chrom,cov_chrom,scores_chrom in data:
			if scores_chrom:
				reads_by_LR.update(reads_by_LR_chrom)
				LRs_by_pos.update(LRs_by_pos_chrom)
				discs_by_barcode.update(discs_by_barcode_chrom)
				discs.update(discs_chrom)
				coverage.append(cov_chrom)
				scores += scores_chrom
		cands,p_len,p_rate = get_candidates(discs,reads_by_LR)
		if cands != None:
			print 'ranking',len(cands),'interchromosomal candidates'
			scores += predict_NAs(reads_by_LR,LRs_by_pos,discs_by_barcode,cands,p_len,p_rate,np.mean(coverage),True)
		else:
			print 'No interchromosomal candidates'
		write_scores(scores)



	print 'Finished in',(time.time()-starttime)/60.0,'minutes'
	return

if __name__ == "__main__":
    main()
