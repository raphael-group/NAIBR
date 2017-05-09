### GLOBALS ####
import os,json
with open('config.json') as f:
	constants = json.load(f)
MIN_MAPQ = constants['MIN_MAPQ']
BAM_FILE=constants['BAM_FILE']
DIR=constants['DIR']
if not os.path.exists(DIR):
    os.makedirs(DIR)
d = constants['d']
DEBUG = False
if DEBUG:
	NUM_CORES = 1
else:
	NUM_CORES = constants['NUM_CORES']
MIN_DISCS = constants['MIN_DISCS']
MIN_SPLIT = constants['MIN_SPLIT']
min_sv = constants['min_sv']