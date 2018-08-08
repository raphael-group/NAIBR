## Overview

NAIBR (Novel Adjacency Identification with Barcoded Reads) identifies novel adjacencies created by structural variation events such as deletions, duplications, inversions, and complex rearrangements using linked-read whole-genome sequencing data produced by 10X Genomics. Please refer to the [publication](https://doi.org/10.1093/bioinformatics/btx712) for details about the method.


NAIBR takes as in put a BAM file produced by 10X Genomic's Long Ranger pipeline and outputs a BEDPE file containing predicted novel adjacencies and a likelihood score for each adjacency.

## Installing NAIBR
```
git clone https://github.com/raphael-group/NAIBR.git
```

NAIBR is written in python 2.7 and requires the following dependencies: pysam, numpy, scipy, and matplotlib

In order to install all dependencies, please execute the following command via `pip`

```
pip install -r requirements.txt
```

## Running NAIBR

NAIBR can be run using the following command:

```
python NAIBR.py <configfile>
```

A template config file can be found in example/example.config. The following parameters can be set in the config file:

* bam_file: Input BAM file < required >
* min_mapq: Minimum mapping quality for a read to be included in analysis (default: 40)
* outdir: Output directory (default: . )
* d: The maximum distance between reads in a linked-read
* blacklist: tap separated list of regions to be excluded from analysis (default: None)
* candidates: List in BEDPE format of novel adjacencies to be scored by NAIBR. This will override automatic detection of candidate novel adjacencies. 
* threads: Number of threads (default: 1)
* min_sv: Minimum size of a structural variant to be detected (default: lmax, the 95th percentile of the paired-end read insert size distribution)
* k: minimum number of barcode overlaps supporting a candidate NA (default = 3)

## Output

NAIBR outputs a BEDPE file containing all novel scored novel adjacencies. Predicted novel adjacencies with scores greater than the threshold c are labelled 'PASS' and others are labelled 'FAIL'. 

## Example
Example files are provided in the 'example' directory. Running

```
python NAIBR.py example/example.config
```

will produce the file 'example/NAIBR_SVs.bedpe'.

### Citing NAIBR
Elyanow, Rebecca, Hsin-Ta Wu, and Benjamin J. Raphael. "Identifying structural variants using linked-read sequencing data." Bioinformatics (2017).
```
@article{elyanow2017identifying,
  title={Identifying structural variants using linked-read sequencing data},
  author={Elyanow, Rebecca and Wu, Hsin-Ta and Raphael, Benjamin J},
  journal={Bioinformatics},
  year={2017}
}
```


