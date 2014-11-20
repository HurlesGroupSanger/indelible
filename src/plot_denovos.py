import sys
import csv
import string
import random
import igv
import os
import time
import subprocess
import numpy as np

"""
    Author: Alejandro Sifrim
    Affiliation: Wellcome Trust Sanger Institute
    
    Annotate unlabeled examples in a trio with label after manual verification by IGV

    Parameters
    ----------
    
	1) Putatitve denovo events file
	2) Training file

    Returns
    -------
	Nothing, outputs a file with the same name but suffix .annotated

"""
    

def random_string(length):
	characters = string.ascii_uppercase + string.digits
	random_characters = [random.choice(characters) for _ in range(length)]
	return ''.join(random_characters)

def read_relationship_file(relationship_path,rotated=False):
	relationships = {}
	for line in open(relationship_path,'r'):
		data = line.split("\t")
		child = os.path.basename(data[0]).strip(".bam")
		if rotated:
			child = os.path.basename(data[1]).strip(".bam")
		else:
			child = os.path.basename(data[0]).strip(".bam")
		relationships[child] = {"child_bam": data[0],"mum_bam": data[1],"dad_bam": data[2]}
	return relationships


def select_denovo(denovo):
	if denovo['ddg2p'] and int(float(denovo['mum_sr'])) <= 1 and int(float(denovo['dad_sr'])) <= 1:
		if denovo['maf']:
			if float(denovo['maf']) < 0.01:
				return True
			else:
				return False
		else:
			return True
	else:
		return False


def decipher_translation_hash(file):
	thash = {}
	for line in open(file,'r'):
		data = line.rstrip().split("\t")
		if len(data) != 3:
			continue
		else:
			thash[data[1]] = data[2]
			thash[data[2]] = data[1]
	return thash	


denovo_file = open("/nfs/users/nfs_a/as33/Projects/Indelible/results/ddd4k_denovos.tsv")
relationships = read_relationship_file("/nfs/users/nfs_a/as33/Projects/Indelible/data/DDD4k_trio_bam_paths.txt",rotated=False)
id_thash = decipher_translation_hash("/nfs/users/nfs_a/as33/Projects/Indelible/data/personid_decipher_id_sangerid.txt")
screenshot_dir = "/nfs/users/nfs_a/as33/Projects/Indelible/results/igv_plots"
igv = igv.IGV()

for row in csv.DictReader(denovo_file,delimiter="\t"):
	if not select_denovo(row):
		continue
	sid = row["sample_id"]

	if sid in relationships:
	 	igv.send('new')
		igv.load(relationships[sid]["child_bam"])
		igv.load(relationships[sid]["mum_bam"])
		igv.load(relationships[sid]["dad_bam"])
		igv.send('collapse')
		pos = int(float(row["position"]))
		igv.go("%s:%s-%s"%(row["chrom"],pos-100,pos+100))
		screenshot_path = "%s/%s_%s_-_%s_%s_%s-%s.png"%(screenshot_dir,sid,id_thash[sid],row["ddg2p"],row["chrom"],pos-100,pos+100)
		screenshot = igv.save(screenshot_path)
		while not os.path.isfile(screenshot_path):
			time.sleep(2)
		
