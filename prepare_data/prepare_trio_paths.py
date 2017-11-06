import pandas as pd
import csv
import os
from signal import signal, SIGPIPE, SIG_DFL

signal(SIGPIPE,SIG_DFL) 
"""
    Author: Alejandro Sifrim
    Affiliation: Wellcome Trust Sanger Institute
    
    Generates a file with the paths to the BAMs across the 4k DDD trios

    Parameters
    ----------
    
    Returns
    -------
    prints tab delimited paths to bam files for trios (child,mother,father) to data/DDD4k_trio_bam_paths.txt
    
"""
    

trios = csv.DictReader(open("/nfs/users/nfs_a/as33/Projects/Indelible/prepare_data/trios_main_ids.txt",'r'),fieldnames=("proband","father","mother"), delimiter="\t")
output_file = open("/nfs/users/nfs_a/as33/Projects/Indelible/data/DDD4k_trio_bam_paths.txt",'w')
bam_file_paths = "/nfs/users/nfs_a/as33/Projects/Indelible/prepare_data/all_ddd_bams.fofn"
bam_paths = {}

for bam_path in open(bam_file_paths):
	bam_path = bam_path.rstrip()
	bam_paths[os.path.basename(bam_path).strip(".bam")] = bam_path

for trio in trios:
	paths = {}
	for person in ("proband","father","mother"):
		main_id = trio[person]
		for m_id in main_id.split(","):
			if m_id in bam_paths:
				paths[person] = bam_paths[m_id]
			# elif os.path.lexists(bam_file_path_2+m_id+'.bam'):
			# 	paths[person] = bam_file_path_2+m_id+'.bam'
	if "proband" in paths and "father" in paths and "mother" in paths:			
		output_file.write("%(proband)s\t%(mother)s\t%(father)s\n" % paths)
