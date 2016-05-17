"""
Author: Alejandro Sifrim
Affiliation: Wellcome Trust Sanger Institute

Library to blast and interpret results of soft-clipped segments: 
- Output segments to fasta file
- Blast them
- Read result

Parameters
----------

Returns
-------

"""
import sys
import os
import csv
import time
import pyblast
import yaml

today = time.strftime('%Y%m%d') 


"""
PARAMETERS
"""

MINIMUM_LENGTH = 20
BLASTdb = "/lustre/scratch113/projects/ddd/users/as33/resources/resources_striped/hs37d5.fa"
WINDOWMASKERdb = "/lustre/scratch113/projects/ddd/users/as33/resources/resources_striped/wmasker.obinary"
REPEATdb = "/lustre/scratch113/projects/ddd/users/as33/Indelible/repbase_fixed.fasta"
BLAST_FIELDS = [
    'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
    'sstart', 'send', 'sstrand','evalue', 'bitscore']
# 
# for i in range(1,100):
# 	testfa = open("/nfs/users/nfs_a/as33/Projects/Python_libs/test_data/test.fa")
# 	print i
# 	for r in pyblast.blastn(testfa, db=blastdb ):
# 		msg = 'query {} has {} hits'.format(r.id, len(r.hits))
# 		if r.hits:
# 		    for h in r.hits:
# 		    	print h.sseqid
# 		    	print h.sstart
# 		    	print h.send
# 		    	print h.length
# 		    	print h.pident
# 	testfa.close()
# 	    # Use Hit.evalue as a float for comparison
# 		# min_evalue = sorted([h.evalue for h in r.hits])[0]

# 		# msg += '; minimum evalue {:f}'.format(min_evalue)
# 	# print msg

def generate_fasta(scored_file):
	output_fasta = open(scored_file+".fasta",'w')
	for row in csv.DictReader(open(scored_file,'r'),delimiter="\t"):
		if len(row['seq_longest']) >= MINIMUM_LENGTH:
			output_fasta.write(">%s_%s_%s\n" % (row['chrom'],row['position'],len(row['seq_longest'])))
			output_fasta.write("%s\n" % row['seq_longest'])
	output_fasta.close()
	return scored_file+".fasta"

def blast_fasta(fasta_file,output_path):
	fa = open(fasta_file)
	fieldnames = ["chrom","pos","query_length","target_chrom","target_start","target_end","target_identity","target_strand","evalue"]
	writer = csv.DictWriter(open(output_path,'w'), fieldnames=fieldnames,delimiter="\t",extrasaction='ignore')
	writer.writeheader()
	for r in pyblast.blastn(fa, window_masker_db=WINDOWMASKERdb,word_size=15,max_target_seqs=100,penalty=-3,evalue=0.001,reward=1,db=BLASTdb,pb_fields=BLAST_FIELDS):
		msg = 'query {} has {} hits'.format(r.id, len(r.hits))
		if len(r.hits) < 10: 
			for hit in r.hits:
				res = {}
				res["chrom"] = r.id.split("_")[0]
				res["pos"] = r.id.split("_")[1]
				res["query_length"] = r.id.split("_")[2]
				res["target_chrom"] = hit.sseqid
				res["target_start"] = hit.sstart
				res["target_end"] = hit.send
				res["target_identity"] = hit.pident
				res["target_strand"] = hit.sstrand
				res["evalue"] = hit.evalue
				writer.writerow(res)
		# print msg

def blast_repeats(fasta_file,output_path):
	fa = open(fasta_file)
	fieldnames = ["chrom","pos","query_length","target_id","target_start","target_end","target_identity","target_strand","evalue"]
	writer = csv.DictWriter(open(output_path,'w'), fieldnames=fieldnames,delimiter="\t",extrasaction='ignore')
	writer.writeheader()
	for r in pyblast.blastn(fa,word_size=15,max_target_seqs=100,penalty=-3,evalue=0.001,reward=1,db=REPEATdb,pb_fields=BLAST_FIELDS):
		msg = 'query {} has {} hits'.format(r.id, len(r.hits))
		if len(r.hits) < 10: 
			for hit in r.hits:
				res = {}
				res["chrom"] = r.id.split("_")[0]
				res["pos"] = r.id.split("_")[1]
				res["query_length"] = r.id.split("_")[2]
				res["target_id"] = hit.sseqid
				res["target_start"] = hit.sstart
				res["target_end"] = hit.send
				res["target_identity"] = hit.pident
				res["target_strand"] = hit.sstrand
				res["evalue"] = hit.evalue
				writer.writerow(res)
		# print msg

def blast(input_path, config):
	fasta_path = generate_fasta(input_path)
	blast_fasta(fasta_path,fasta_path+".hits_nonrepeats")
	blast_repeats(fasta_path,fasta_path+".hits_repeats")
