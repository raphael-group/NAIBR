from global_vars import *
import math,mpmath
import numpy as np
import linecache
import scipy.stats

def get_barcode(barcode):
	thefilename = os.path.join(DIR,'LRS.txt')
	line = linecache.getline(thefilename, barcode+1)
	b,line = line.split(':')
	assert b == str(barcode)
	line = line.split('\t')
	line = [[[int(z) for z in y.split(',')] for y in x.split(' ') if len(y.split(',')) == 4] for x in line]
	line = [x for x in line if len(x) > 0]
	return line

def NB(k,m,r):
	return mpmath.gamma(k + r)/(mpmath.gamma(k+1)*mpmath.gamma(r))*np.power(m/(r+m),k)*np.power(r/(r+m),r)

def Pois(k,lam):
	return (np.power(lam,k)*math.exp(-lam))/mpmath.gamma(k+1)
def get_hap(r):
	if r.has_tag('HP'):
		return r.get_tag('HP')
	return 0

def to_chr(chrm):
	chrm = chrm.strip('chr')
	if chrm == 'X':
		return '23'
	elif chrm == 'Y':
		return '24'
	return str(chrm)

def plog(x,num):
	ret = 0
	for i in range(num):
		ret += log(x)
	return ret

def log(x):
	try:
		return math.log(x)
	except:
		return math.log(1e-300)


def pow(x,y):
	return max(1e-300,math.pow(x,y))


def orient(read):
	if read.is_reverse:
		return '-'
	return '+'

def write_file(fname,read,mate,barcode):
	if mate.reference_end-read.reference_start >= 4000:
		fname.write('\t'.join([str(barcode),to_chr(read.reference_name),str(read.reference_start),str(read.reference_end),orient(read),\
			to_chr(read.reference_name),str(mate.reference_start),str(mate.reference_end),orient(mate)+'\n']))
	return

def is_convergent(read,mate):
	return not read.is_reverse and mate.is_reverse

def get_orientation(read,mate):
	# if read.is_read1:
	if read.is_reverse:
		if mate.is_reverse:
			return ['--',read.reference_start,mate.reference_start]
		else:
			return ['-+',read.reference_start,mate.reference_end]
	else:
		if read.mate_is_reverse:
			return ['+-',read.reference_end,mate.reference_start]
		else:
			return ['++',read.reference_end,mate.reference_end]
	# elif mate.is_read1:
	# 	if mate.is_reverse:
	# 		if read.is_reverse:
	# 			return ['--',mate.reference_start,read.reference_start]
	# 		else:
	# 			return ['-+',mate.reference_start,read.reference_end]
	# 	else:
	# 		if read.mate_is_reverse:
	# 			return ['+-',mate.reference_end,read.reference_start]
	# 		else:
	# 			return ['++',mate.reference_end,read.reference_end]
	# else:
	# 	raise ValueError('Nothing is read1')
def swap(read,mate,break1,break2):
	b1 = break2
	b2 = break1
	read1 = mate
	mate1 = read
	return(read1,mate1,b1,b2)

def make_key(name_list):
	return '.'.join([str(x) for x in name_list])

def list_from_key(key):
	return key.split('.')

def rnd(num):
	return num/R*R

def flatten(l):
	return [item for sublist in l for item in sublist]

def spans(Window,point):
	(start,end) = Window
	if point >= start-LMAX and point <= end+LMAX:
		return True
	return False