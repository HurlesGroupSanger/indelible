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
import csv
import time
from Bio.Blast.Applications import NcbiblastnCommandline
import StringIO

today = time.strftime('%Y%m%d')


"""
PARAMETERS
"""

MINIMUM_LENGTH = 20
MINIMUM_SCORE = 0.5

BLAST_FIELDS = [
	'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
	'sstart', 'send', 'sstrand','evalue', 'bitscore']


def generate_fasta(scored_file):
	output_fasta = open(scored_file+".fasta",'w')
	for row in csv.DictReader(open(scored_file,'r'),delimiter="\t"):
		if len(row['seq_longest']) >= MINIMUM_LENGTH and float(row['prob_Y']) > MINIMUM_SCORE:
			output_fasta.write(">%s_%s_%s\n" % (row['chrom'],row['position'],len(row['seq_longest'])))
			output_fasta.write("%s\n" % row['seq_longest'])
	output_fasta.close()
	return scored_file+".fasta"


def run_blast(fasta_file, db, WINDOWMASKERdb=None):

	blastn_clin = NcbiblastnCommandline(query=fasta_file, db=db, word_size=15,
										max_target_seqs=100, penalty=-3, evalue=0.001, reward=1,
										outfmt="\'6 " + " ".join(BLAST_FIELDS) + "\'")
	if not WINDOWMASKERdb is None:
		blastn_clin.set_parameter("window_masker_db",WINDOWMASKERdb)

	stdout, stderr = blastn_clin()
	input = StringIO.StringIO(stdout)

	hits = {}

	for r in input:
		data = r.rstrip().split("\t")
		if data[0] in hits:
			hits[data[0]].append(data)
		else:
			hits[data[0]] = [data]

	return hits


def blast_fasta(fasta_file,output_path,db,WINDOWMASKERdb = None):

	fieldnames = ["chrom","pos","query_length","target_chrom","target_start","target_end","target_identity","target_strand","evalue"]
	writer = csv.DictWriter(open(output_path,'w'), fieldnames=fieldnames,delimiter="\t",extrasaction='ignore')
	writer.writeheader()

	hits = run_blast(fasta_file,db,WINDOWMASKERdb)

	for seq in hits:
		msg = 'query {} has {} hits'.format(seq, len(hits[seq]))
		#print msg
		if len(hits[seq]) < 10:
			for result in hits[seq]:
				res = {}
				res["chrom"] = result[0].split("_")[0]
				res["pos"] = result[0].split("_")[1]
				res["query_length"] = result[0].split("_")[2]
				res["target_chrom"] = result[1]
				res["target_start"] = result[8]
				res["target_end"] = result[9]
				res["target_identity"] = result[2]
				res["target_strand"] = result[10]
				res["evalue"] = result[11]
				writer.writerow(res)

def blast(input_path, config):
	BLASTdb = config['blastdb']
	WINDOWMASKERdb = config['windowmaskerdb']
	REPEATdb = config['repeatdb']
	fasta_path = generate_fasta(input_path)
	blast_fasta(fasta_path, fasta_path + ".hits_nonrepeats", BLASTdb, WINDOWMASKERdb)
	blast_fasta(fasta_path, fasta_path + ".hits_repeats", REPEATdb)
