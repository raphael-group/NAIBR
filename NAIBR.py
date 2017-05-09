import os,sys,time,json
from get_reads import make_barcodeDict
from global_vars import *
from estimate_params import estimate_lmin_lmax
from linked_reads import parallel_get_linked_reads
from distributions import get_distributions
from utils import *
from rank import predict_NAs



def main():
	starttime = time.time()
	lmin,lmax = estimate_lmin_lmax()
	reads_by_barcode,pp_by_barcode,pm_by_barcode,mp_by_barcode,mm_by_barcode = make_barcodeDict(lmin,lmax,'')
	starts,ends,mm_cands,mp_cands,pm_cands,pp_cands,LENS,RATES = parallel_get_linked_reads(reads_by_barcode,pp_by_barcode,pm_by_barcode,mp_by_barcode,mm_by_barcode)
	reads_by_barcode = None
	p_len,p_rate = get_distributions(LENS,RATES)
	predict_NAs(pp_by_barcode,pm_by_barcode,mp_by_barcode,mm_by_barcode,starts,ends,mm_cands,mp_cands,pm_cands,pp_cands,p_len,p_rate,lmax)
	print 'Finished in',(time.time()-starttime)/60.0,'minutes'

if __name__ == "__main__":
    main()
