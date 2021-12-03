#!/usr/bin/env python

"""
    Author: Alejandro Sifrim
    Author: Eugene Gardner
    Affiliation: Wellcome Sanger Institute

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

version = "1.1.3"

parser = argparse.ArgumentParser(prog='indelible',description="InDelible v" + version + " -- Structural variant discovery with split reads.\n\n")
subparsers = parser.add_subparsers(help='One of the following commands:',dest="command",metavar="<command>")

subparser = subparsers.add_parser('fetch', help='fetch reads from BAM file')
subparser.add_argument('--config', help='path to config file', metavar="<config_path>", required=True, dest="config_path", default = None)
subparser.add_argument('--i', help='path to input BAM file', metavar="<input_path>", required=True, dest="input_path")
subparser.add_argument('--o', help='path to output file', metavar="<output_path>", required=True, dest="output_path")

subparser = subparsers.add_parser('aggregate', help='aggregate information per position')
subparser.add_argument('--config', help='path to config file', metavar="<config_path>", required=True, dest="config_path", default = None)
subparser.add_argument('--i', help='path to input file (output of fetch command)', metavar="<input_path>", required=True, dest="input_path")
subparser.add_argument('--b', help='path to input bam file which was used in the fetch command', metavar="<input_bam>", required=True, dest="input_bam")
subparser.add_argument('--o', help='path to output file', metavar="<output_path>", required=True, dest="output_path")
subparser.add_argument('--r', help='path to reference genome', metavar="<reference_path>", required=True, dest="reference_path")

subparser = subparsers.add_parser('score', help='score positions using Random Forest model')
subparser.add_argument('--config', help='path to config file', metavar="<config_path>", required=True, dest="config_path", default = None)
subparser.add_argument('--i', help='input file (output of aggregate command)', metavar="<input_path>", required=True, dest="input_path")
subparser.add_argument('--o', help='path to output file', metavar="<output_path>", required=True, dest="output_path")

subparser = subparsers.add_parser('database', help='build SR allele frequency database')
subparser.add_argument('--config', help='path to config file', metavar="<config_path>", required=True, dest="config_path", default = None)
subparser.add_argument('--f', help='File of files from the score command representing a complete dataset', metavar="<fof>", required=True, dest="fof")
subparser.add_argument('--r', help='path to reference genome', metavar="<reference_path>", required=True, dest="reference_path")
subparser.add_argument('--o', help='path to output file', metavar="<output_path>", required=True, dest="output_path")
subparser.add_argument('--priors', help='path to priors MAF database [optional]', metavar="<priors>", required=False, dest="priors", default=None)
subparser.add_argument('--old-maf', help='Use MAF from the priors file for rediscovered loci?', metavar="<old_maf>", dest="old_maf", action='store_const', const = True)
subparser.add_argument('--tb', help='number of threads to use for bwa alignment. [1]', metavar="<bwa_thread>", required=False, dest="bwa_thread",default=1)

subparser = subparsers.add_parser('annotate', help='annotate positions with additional information')
subparser.add_argument('--config', help='path to config file', metavar="<config_path>", required=True, dest="config_path", default = None)
subparser.add_argument('--i', help='input file (output of score command)', metavar="<input_path>", required=True, dest="input_path")
subparser.add_argument('--o', help='path to output file', metavar="<output_path>", required=True, dest="output_path")
subparser.add_argument('--d', help='path to indelible frequency database', metavar="<database>", required=True, dest="database")

subparser = subparsers.add_parser('denovo', help='searches for de novo events')
subparser.add_argument('--config', help='path to config file', metavar="<config_path>", required=True, dest="config_path", default = None)
subparser.add_argument('--c', help='path to scored/annotated calls in the child', metavar="<child_indelible_path>", required=True, dest="child_input")
subparser.add_argument('--m', help='path to maternal BAM file', metavar="<mother_bam_path>", required=False, dest="mother_bam")
subparser.add_argument('--p', help='path to paternal BAM file', metavar="<father_bam_path>", required=False, dest="father_bam")
subparser.add_argument('--o', help='path to output file', metavar="<output_path>", required=True, dest="output_path")

subparser = subparsers.add_parser('complete', help='Performs the complete Indelible analysis')
subparser.add_argument('--config', help='path to config file', metavar="<config_path>", required=True, dest="config_path", default=None)
subparser.add_argument('--i', help='path to input BAM file', metavar="<input_path>", required=True, dest="input_path")
subparser.add_argument('--o', help='path to output directory', metavar="<output_path>", required=True, dest="output_path")
subparser.add_argument('--r', help='path to reference genome', metavar="<reference_path>", required=True, dest="reference_path")
subparser.add_argument('--m', help='path to maternal BAM file', metavar="<mother_bam_path>", required=False, dest="mother_bam")
subparser.add_argument('--p', help='path to paternal BAM file', metavar="<father_bam_path>", required=False, dest="father_bam")
subparser.add_argument('--priors', help='path to priors MAF database [optional]', metavar="<priors>", required=False, dest="priors", default=None)
subparser.add_argument('--old-maf', help='Use MAF from the priors file for rediscovered loci?', metavar="<old_maf>", dest="old_maf", action='store_const', const = True)
subparser.add_argument('--tb', help='number of threads to use for bwa alignment. [1]', metavar="<bwa_thread>", required=False, dest="bwa_thread",default=1)
subparser.add_argument('--keeptmp', action='store_const', const=True,  dest="keep_tmp")

subparser = subparsers.add_parser('train', help='trains the Random Forest model on a bunch of examples')
subparser.add_argument('--config', help='path to config file', metavar="<config_path>", required=True, dest="config_path", default = None)
subparser.add_argument('--i', help='Input file with labeled examples', metavar="<input_path>", required=True, dest="input_path")
subparser.add_argument('--o', help='Output path for RF model', metavar="<output_path>", required=True, dest="output_path")
subparser.add_argument('--k', help='value of k hyperparameter for training the Random Forest [75].', metavar="<k>", required=False, default=75, type=int, dest="k")
subparser.add_argument('--s', help='value of the stop parameter for training the Random Forest [0.01].', metavar="<stop_parameter>", required=False, default=0.01, type=float, dest="stop_parameter")

args = parser.parse_args()

"""
Read config file
"""

config_path = args.config_path
config = yaml.load(open(config_path), Loader=yaml.SafeLoader)

required = ["MINIMUM_SR_COVERAGE","SHORT_SR_CUTOFF",
"MINIMUM_LENGTH_SPLIT_READ","MINIMUM_MAPQ",
"MININUM_AVERAGE_BASE_QUALITY_SR","HC_THRESHOLD",'random_forest_model']

for r in required:
    if r not in config:
        print("ERROR: %s not specified in config file!")
        exit(1)

if args.command is None:
    print("\nMust specify one of 'fetch', 'aggregate', 'score', 'database', 'annotate', 'denovo', 'complete', 'train'\n")
    parser.print_help()
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
DATABASE command
"""
if args.command == "database":
    indelible.build_database(args.fof, args.output_path, args.reference_path, config, args.priors, args.old_maf, args.bwa_thread)

"""
ANNOTATE command
"""
if args.command == "annotate":
    if not os.path.isfile(args.input_path):
        print("ERROR: Input file does not exist!")
        exit(1)
    indelible.annotate(args.input_path, args.output_path, args.database, config)

"""
DENOVO command
"""
if args.command == "denovo":

    if not os.path.isfile(args.child_input):
        print("ERROR: Input file does not exist!")
        exit(1)

    indelible.denovo_caller_trio(args.child_input, args.mother_bam, args.father_bam, args.output_path, config)

"""
TRAIN command
"""
if args.command == "train":
    if not os.path.isfile(args.input_path):
        print("ERROR: Input file does not exist!")
        exit(1)
    indelible.train(args.input_path, args.output_path, args.k, args.stop_parameter)

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
    score_list_path = args.output_path + "/" + "score.lst"
    score_writer = open(score_list_path,'w')
    score_writer.write(scored_path + "\n")
    score_writer.close()
    db_path = args.output_path+"/"+"maf_db.txt"
    final_path = args.output_path+"/"+os.path.basename(args.input_path)+".indelible.tsv"
    denovo_path = args.output_path+"/"+os.path.basename(args.input_path)+".indelible.denovo.tsv"

    print(("%s: Fetching reads..." % timestamp()))
    indelible.fetch_reads(args.input_path, reads_path, config)
    print(("%s: Aggregating across positions..." % timestamp()))
    indelible.aggregate_positions(reads_path, args.input_path, counts_path, args.reference_path, config)
    print(("%s: Scoring positions..." % timestamp()))
    indelible.score_positions(counts_path, scored_path, config)
    print(("%s: Building InDelible database and adding priors..." % timestamp()))
    indelible.build_database(score_list_path, db_path, args.reference_path, config, args.priors, args.old_maf, args.bwa_thread)
    print(("%s: Annotating positions..." % timestamp()))
    indelible.annotate(scored_path, annotated_path, db_path, config)
    shutil.copy(annotated_path,final_path)
    print(("%s: Calling de novo variants..." % timestamp()))
    indelible.denovo_caller_trio(final_path, args.mother_bam, args.father_bam, denovo_path, config)
    if args.keep_tmp != True:
        print(("%s: Removing temporary files..." % timestamp()))
        os.remove(reads_path)
        os.remove(counts_path)
        os.remove(scored_path)
        os.remove(annotated_path)
        os.remove(scored_path+".fasta.hits_nonrepeats")
        os.remove(scored_path+".fasta.hits_repeats")
        os.remove(scored_path+".fasta")