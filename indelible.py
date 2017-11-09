#!/usr/bin/env python

"""
    Author: Alejandro Sifrim
    Affiliation: Wellcome Trust Sanger Institute

    Indelible main program

"""
import time
import yaml
import os
import argparse
import indelible
import shutil

today = time.strftime('%Y%m%d')

def timestamp():
	return time.strftime('%d/%m/%y - %H:%M:%S')

usage_text = """

Indelible
=========
Alejandro Sifrim - Wellcome Trust Sanger Institute
as33@sanger.ac.uk

Usage:
./indelible command

Where command can be:

fetch - fetches reads from BAM file
aggregate - aggregate information per position
score - score positions using Random Forest model
annotate - annotate positions with additional information
denovo - searches for de novo events

"""
parser = argparse.ArgumentParser(prog='indelible')
subparsers = parser.add_subparsers(help='One of the following commands:',dest="command",metavar="<command>")
#Fetch parser
subparser = subparsers.add_parser('fetch', help='fetch reads from BAM file')
subparser.add_argument('--i', help='path to input BAM file', metavar="<input_path>", required=True, dest="input_path")
subparser.add_argument('--o', help='path to output file', metavar="<output_path>", required=True, dest="output_path")

subparser = subparsers.add_parser('aggregate', help='aggregate information per position')
subparser.add_argument('--i', help='path to input file (output of fetch command)', metavar="<input_path>", required=True, dest="input_path")
subparser.add_argument('--b', help='path to input bam file which was used in the fetch command', metavar="<input_bam>", required=True, dest="input_bam")
subparser.add_argument('--o', help='path to output file', metavar="<output_path>", required=True, dest="output_path")
subparser.add_argument('--r', help='path to reference genome', metavar="<reference_path>", required=True, dest="reference_path")

subparser = subparsers.add_parser('score', help='score positions using Random Forest model')
subparser.add_argument('--i', help='input file (output of aggregate command)', metavar="<input_path>", required=True, dest="input_path")
subparser.add_argument('--o', help='path to output file', metavar="<output_path>", required=True, dest="output_path")

subparser = subparsers.add_parser('blast', help='blast clipped sequences')
subparser.add_argument('--i', help='input file (output of score command)', metavar="<input_path>", required=True, dest="input_path")

subparser = subparsers.add_parser('annotate', help='annotate positions with additional information')
subparser.add_argument('--i', help='input file (output of score command)', metavar="<input_path>", required=True, dest="input_path")
subparser.add_argument('--o', help='path to output file', metavar="<output_path>", required=True, dest="output_path")

subparser = subparsers.add_parser('denovo', help='searches for de novo events')
subparser.add_argument('--c', help='path to scored/annotated calls in the child', metavar="<child_indelible_path>", required=True, dest="child_input")
subparser.add_argument('--m', help='path to maternal BAM file', metavar="<mother_bam_path>", required=True, dest="mother_bam")
subparser.add_argument('--p', help='path to paternal BAM file', metavar="<father_bam_path>", required=True, dest="father_bam")
subparser.add_argument('--o', help='path to output file', metavar="<output_path>", required=True, dest="output_path")

subparser = subparsers.add_parser('train', help='trains the Random Forest model on a bunch of examples')
subparser.add_argument('--i', help='Input file with labeled examples', metavar="<input_path>", required=True, dest="input_path")
subparser.add_argument('--o', help='Output path for RF model', metavar="<output_path>", required=True, dest="output_path")

subparser = subparsers.add_parser('complete', help='Performs the complete Indelible analysis')
subparser.add_argument('--i', help='path to input BAM file', metavar="<input_path>", required=True, dest="input_path")
subparser.add_argument('--o', help='path to output directory', metavar="<output_path>", required=True, dest="output_path")
subparser.add_argument('--r', help='path to reference genome', metavar="<reference_path>", required=True, dest="reference_path")
subparser.add_argument('--keeptmp', action='store_const', const=True,  dest="keep_tmp")
args = parser.parse_args()


"""
Read config file
"""
config_path = os.path.join(os.path.dirname(__file__), 'config.yml')
config = yaml.load(open(config_path))

required = ["MINIMUM_SR_COVERAGE","SHORT_SR_CUTOFF",
"MINIMUM_LENGTH_SPLIT_READ","MINIMUM_MAPQ",
"MININUM_AVERAGE_BASE_QUALITY_SR","HC_THRESHOLD",'random_forest_model']

for r in required:
	if r not in config:
		print("ERROR: %s not specified in config file!")
		exit(1)

"""
FETCH command
"""
if args.command == "fetch":
	indelible.fetch_reads(args.input_path, args.output_path, config)

"""
AGGREGATE command
"""

if args.command == "aggregate":
	if not os.path.isfile(args.input_path):
		print("ERROR: Input file does not exist or cannot be accessed!")
		exit(1)
	if not os.path.isfile(args.input_bam):
		print("ERROR: Input BAM does not exist or cannot be accessed!")
		exit(1)
	if not os.path.isfile(args.reference_path):
		print("ERROR: Reference path does not exist or cannot be accessed!")
		exit(1)
	indelible.aggregate_positions(args.input_path, args.input_bam, args.output_path, args.reference_path, config)

"""
SCORE command
"""
if args.command == "score":
	if not os.path.isfile(config['random_forest_model']):
		print("ERROR: The path specified for random_forest_model in config does not exist!")
		exit(1)
	if not os.path.isfile(args.input_path):
		print("ERROR: Input file does not exist!")
		exit(1)
	indelible.score_positions(args.input_path, args.output_path, config)

"""
BLAST
"""
if args.command == "blast":
	if not os.path.isfile(args.input_path):
		print("ERROR: Input file does not exist!")
		exit(1)
	indelible.blast(args.input_path,config)

"""
ANNOTATE command
"""
if args.command == "annotate":
	if not os.path.isfile(args.input_path):
		print("ERROR: Input file does not exist!")
		exit(1)
	indelible.annotate(args.input_path, args.output_path, config)

"""
DENOVO command
"""
if args.command == "denovo":
	if not os.path.isfile(args.child_input):
		print("ERROR: Input file does not exist!")
		exit(1)
	if not os.path.isfile(args.mother_bam):
		print("ERROR: Maternal BAM file does not exist!")
		exit(1)
	if not os.path.isfile(args.father_bam):
		print("ERROR: Paternal BAM file does not exist!")
		exit(1)
	indelible.denovo_caller(args.child_input, args.mother_bam, args.father_bam, args.output_path, config)

"""
TRAIN command
"""
if args.command == "train":
	if not os.path.isfile(args.input_path):
		print("ERROR: Input file does not exist!")
		exit(1)
	indelible.train(args.input_path, args.output_path)

"""
COMPLETE command
"""
if args.command == "complete":
	if not os.path.isfile(args.input_path):
		print("ERROR: Input file does not exist!")
		exit(1)
	if not os.path.isfile(args.reference_path):
		print("ERROR: Reference path does not exist or cannot be accessed!")
		exit(1)
	if not os.path.isdir(args.output_path):
		print("ERROR: Output directory does not exist!")
		exit(1)

	reads_path = args.output_path+"/"+os.path.basename(args.input_path)+".sc_reads"
	counts_path = reads_path+".aggregated"
	scored_path = counts_path+".scored"
	annotated_path = counts_path+".annotated"
	final_path = args.output_path+"/"+os.path.basename(args.input_path)+".indelible.tsv"
	print("%s: Fetching reads..." % timestamp())
	indelible.fetch_reads(args.input_path, reads_path, config)
	print("%s: Aggregating across positions..." % timestamp())
	indelible.aggregate_positions(reads_path, args.input_path, counts_path, args.reference_path, config)
	print("%s: Scoring positions..." % timestamp())
	indelible.score_positions(counts_path, scored_path, config)
	print("%s: Blasting soft-clipped segments..." % timestamp())
	indelible.blast(scored_path, config)
	print("%s: Annotating positions..." % timestamp())
	indelible.annotate(scored_path, annotated_path, config)
	shutil.copy(annotated_path,final_path)
	if args.keep_tmp != True:
		print("%s: Removing temporary files..." % timestamp())
		os.remove(reads_path)
		os.remove(counts_path)
		os.remove(scored_path)
		os.remove(annotated_path)
		os.remove(scored_path+".fasta.hits_nonrepeats")
		os.remove(scored_path+".fasta.hits_repeats")
		os.remove(scored_path+".fasta")
