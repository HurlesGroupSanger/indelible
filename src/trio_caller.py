"""
    Author: Alejandro Sifrim
    Affiliation: Wellcome Trust Sanger Institute
    
    Script to call inheritance mode of an event

    Parameters
    ----------
    1) child score file
	2) mother bam
	3) father bam

    Returns
    -------
	1) File with putative denovo calls

"""

import os
import sys
import pysam
import csv
from indelible_lib import *
import scipy.stats as sp
import math
import pprint

SCORE_THRESHOLD = 0.6
SR_THRESHOLD = 1
COV_THRESHOLD = 9
def compute_stats(bam_file,chrom,position,relationship):
	res = {}
	res['coverage'] = coverage_at_position(bam_file,chrom,position)
	res['indel_context'] = reads_with_indels_in_neighbourhood(bam_file,chrom,position,20)
	res['sr_context'] = split_reads_in_neighbourhood(bam_file,chrom,position,5)
	# res['relationship'] = relationship
	return res

db = read_database("/nfs/users/nfs_a/as33/Projects/MEISRFinder/y2_parents_database_sorted.tsv")

mum_bam = pysam.Samfile(sys.argv[2],'rb')
dad_bam = pysam.Samfile(sys.argv[3],'rb')
scored_file = csv.DictReader(open(sys.argv[1],'r'), delimiter="\t")
new_fieldnames = scored_file.fieldnames
new_fieldnames.extend(("mum_sr","dad_sr",'mum_indel_context','dad_indel_context',"mum_cov","dad_cov"))
output_file = csv.DictWriter(open(sys.argv[1]+".denovo",'w'),fieldnames=new_fieldnames,delimiter="\t")
output_file.writeheader()
for v in csv.DictReader(open(sys.argv[1],'r'), delimiter="\t"):
	if float(v['prob_Y']) >= SCORE_THRESHOLD:

	 	mum_stats = compute_stats(mum_bam,v['chrom'],int(v['position']),'mum')
		dad_stats = compute_stats(dad_bam,v['chrom'],int(v['position']),'dad')
		
		v['mum_sr'] =  mum_stats['sr_context']
		v['dad_sr'] =  dad_stats['sr_context']
		v['mum_indel_context'] =  mum_stats['indel_context']['deletions'] +  mum_stats['indel_context']['insertions']
		v['dad_indel_context'] =  dad_stats['indel_context']['deletions'] + dad_stats['indel_context']['insertions']
		v['mum_cov'] = mum_stats['coverage']
		v['dad_cov'] = dad_stats['coverage']

		if v['mum_sr'] <= SR_THRESHOLD and v['dad_sr'] <= SR_THRESHOLD and v['mum_cov'] >= COV_THRESHOLD and v['dad_cov'] >= COV_THRESHOLD:
			output_file.writerow(v)
