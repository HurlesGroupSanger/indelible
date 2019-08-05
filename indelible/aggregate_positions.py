import sys
import pysam
import csv
import math
from collections import Counter
from pyfaidx import Fasta
import swalign
import numpy
import re
from indelible.indelible_lib import *
"""
Author: Alejandro Sifrim
Affiliation: Wellcome Trust Sanger Institute

Takes the output of fetch_reads.py and computes stats per
position such as sr_coverage, coverage overall, shannon entropy
of surrounding regions, sequence similarity between split segments
and other good stuff...

Parameters
----------
1) fetch_reads.py outputfile
2) original BAM file
3) reference fasta
Returns
-------
File with stats for each position

"""


#CONSTANTS
#=========

def dedup(sr_reads = []):
	tmp = {"3":{},"5":{}}
	for read in sr_reads:
		tmp[read["prime"]][read["split_length"]] = read
	res = []
	res.extend(tmp["3"].values())
	res.extend(tmp["5"].values())
	return res


def sr_coverage(sr_reads = [],cutoff=10):
	total = 0
	total_long = 0
	total_short = 0
	short_3 = 0
	long_3 = 0
	short_5 = 0
	long_5 = 0
	double_reads = 0

	deduped_reads = dedup(sr_reads)
	for read in deduped_reads:
		total += 1

		if read["is_double"] == "True":
			double_reads += 1

		if int(read["split_length"]) < cutoff:
			total_short += 1
			if int(read["prime"]) == 3:
				short_3 += 1
			else:
				short_5 += 1
		else:
			total_long += 1
			if int(read["prime"]) == 3:
				long_3 += 1
			else:
				long_5 += 1

	return (total,total_long,total_short,long_5,short_5,long_3,short_3,double_reads)


def entropy_longest_sr(sr_reads=[]):
	longest = sorted(sr_reads,key=lambda x:int(x["split_length"]),reverse=True)[0]
	return entropy(longest["seq"])


def seq_longest(sr_reads=[]):
	longest = sorted(sr_reads,key=lambda x:int(x["split_length"]),reverse=True)[0]
	return longest["seq"]


def sequence_similarity_score(sr_reads=[]):
	match = 2
	mismatch = -1
	scoring = swalign.NucleotideScoringMatrix(match, mismatch)
	sw = swalign.LocalAlignment(scoring)
	sequences = [x["seq"] for x in sr_reads]
	sequences = sorted(sequences,key=len,reverse=True)
	aln_scores = []
	for seq1_idx in range(len(sequences)):
		for seq2_idx in range(len(sequences)):
			if seq1_idx != seq2_idx:
				seq1 = sequences[seq1_idx]
				seq2 = sequences[seq2_idx]
				min_length = min(len(seq1),len(seq2))
				matches = sw.align(seq1,seq2).matches
				aln_scores.append(float(matches)/min_length)
	return numpy.mean(aln_scores)


def avg_mapq(sr_reads=[]):
	mapqs = [int(x["mapq"]) for x in sr_reads]
	return numpy.mean(mapqs)


def avg_avg_sr_qual(sr_reads=[]):
	avg_sr_quals = list([float(x["avg_sr_qual"]) for x in sr_reads])
	return numpy.mean(avg_sr_quals)


def aggregate_positions(input_path, input_bam, output_path, reference_path, config):

	#ARGUMENTS ARE 1) sr_reads file 2) BAM file 3) reference fasta file

	chr_dict = {}

	splitfile = open(input_path,'r')

	bam_file = bam_open(input_bam)

	reference = Fasta(reference_path,as_raw=True)
	splitreader = csv.DictReader(splitfile,delimiter="\t",quoting=csv.QUOTE_NONE)
	header = ("chrom","position","coverage","insertion_context","deletion_context",
				  "sr_total","sr_total_long","sr_total_short",
				  "sr_long_5","sr_short_5","sr_long_3",
				  "sr_short_3","sr_entropy","context_entropy",
				  "entropy_upstream","entropy_downstream","sr_sw_similarity",
				  "avg_avg_sr_qual","avg_mapq","seq_longest","pct_double_split")
	outputfile = open(output_path,'w')
	splitwriter = csv.DictWriter(outputfile,fieldnames=header,delimiter="\t",lineterminator="\n")
	splitwriter.writeheader()

	for row in splitreader:
		if row['chr'] is "hs37d5":
			continue
		if re.search("N",row['seq']):
			continue
		if not row['chr'] in chr_dict:
			chr_dict[row['chr']] = {}
			chr_dict[row['chr']][row['split_position']] = []
		if not row['split_position'] in chr_dict[row['chr']]:
			chr_dict[row['chr']][row['split_position']] = []

        chr_dict[row['chr']][row['split_position']].append(row)

	countint = 0
	count = 0

	for chrom in chr_dict:
		for position in chr_dict[chrom]:
			countint += 1
			pos = int(position)
			print len(chr_dict[chrom][position])
			if len(chr_dict[chrom][position]) >= config["MINIMUM_SR_COVERAGE"]:
				sr_reads = chr_dict[chrom][position]

				count += 1

				# res = {}
				# res["chrom"] = chrom
				# res["position"] = pos
				#
				# # covPU = coverage_at_position_pileup(bam_file, chrom, pos)
				# covFE = coverage_at_position_pileup(bam_file, chrom, pos)
				#
				# res["coverage"] = covFE
				#
				# indel_counts = reads_with_indels_in_neighbourhood(bam_file,chrom,pos,config)
				# res["insertion_context"] = indel_counts["insertions"]
				# res["deletion_context"] = indel_counts["deletions"]
				# sr_cov = sr_coverage(sr_reads,config["SHORT_SR_CUTOFF"])
				# res["sr_total"] = sr_cov[0]
				# res["sr_total_long"] = sr_cov[1]
				# res["sr_total_short"] = sr_cov[2]
				# res["sr_long_5"] = sr_cov[3]
				# res["sr_short_5"] = sr_cov[4]
				# res["sr_long_3"] = sr_cov[5]
				# res["sr_short_3"] = sr_cov[6]
				# res["pct_double_split"] = float(sr_cov[7]) / float(res["sr_total"])
				# res["sr_entropy"] = entropy_longest_sr(sr_reads)
				# seq_context = reference[str(chrom)][pos-20:pos+20]
				# res["context_entropy"] = entropy(seq_context)
				# res["entropy_upstream"] = entropy(seq_context[0:20])
				# res["entropy_downstream"] = entropy(seq_context[20:])
				# res["sr_sw_similarity"] = sequence_similarity_score(sr_reads)
				# res["avg_mapq"] = avg_mapq(sr_reads)
				# res["avg_avg_sr_qual"] = avg_avg_sr_qual(sr_reads)
				# res["seq_longest"] = seq_longest(sr_reads)
				# splitwriter.writerow(res)
				# outputfile.flush()
	print "Total Init Reads " + str(countint)
	print "Total Reads " + str(count)