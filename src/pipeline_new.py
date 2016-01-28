from ruffus import *
import os
import sys
import time
from bsub import bsub

"""
 Author: Alejandro Sifrim
 Affiliation: Wellcome Trust Sanger Institute



 Parameters
 ----------

 Returns
 -------

"""

OUTPUT_DIR = "/lustre/scratch113/projects/ddd/users/as33/Indelible/DDD4k/new/"
SAMTOOLS_BIN = "/software/hgi/pkglocal/samtools-0.1.19/bin/samtools"
SRC_DIR = "/nfs/users/nfs_a/as33/Projects/Indelible/src"
REF = "/lustre/scratch113/projects/ddd/resources/v1.2/hs37d5.fasta"
RF_MODEL = "/nfs/users/nfs_a/as33/Projects/Indelible/data/random_forest.pkl"
#Input arguments
trio = {}
child_base = os.path.basename(sys.argv[1]).strip(".bam")
trio[os.path.basename(sys.argv[1]).strip(".bam")] = sys.argv[1]
trio[os.path.basename(sys.argv[2]).strip(".bam")] = sys.argv[2]
trio[os.path.basename(sys.argv[3]).strip(".bam")] = sys.argv[3]

#Create the outputdir and change to it
trio_dir = OUTPUT_DIR+child_base
if not os.path.exists(trio_dir):
	os.makedirs(trio_dir)
os.chdir(trio_dir)
starting_files = sys.argv[1:]

@transform(starting_files,regex(r".+/(.*).bam"),r"%s/\1.sr_reads"%trio_dir)
def fetch_reads(input_file,output_file):
	command = "%s view -b %s" % (SAMTOOLS_BIN,input_file)
	command += " | python %s/fetch_reads.py > %s" % (SRC_DIR,output_file)
	print command
	sub = bsub("Indelible_-_Fetch_Reads",R=" select[mem>1000] rusage[mem=1000]",M=1000)
	sub(command)
	time.sleep(10)
	bsub.poll(sub.job_id)

@transform(fetch_reads,suffix(".sr_reads"),".sr_reads.counts")
def aggregate_positions(input_file,output_file):
	base = os.path.basename(input_file).strip(".sr_reads")
	bam_path = trio[base]
	command = "python %s/aggregate_positions.py %s %s %s" % (SRC_DIR,input_file,bam_path,REF)
	print command
	sub = bsub("Indelible_-_Aggregate_positions",R="select[mem>1000] rusage[mem=1000]",M=1000)
	sub(command)
	time.sleep(10)
	bsub.poll(sub.job_id)

@transform(aggregate_positions,suffix(".sr_reads.counts"),".sr_reads.counts.scored")
def score_positions(input_file,output_file):
	command = "python %s/score_positions.py %s %s" % (SRC_DIR,input_file,RF_MODEL)
	print command
	sub = bsub("Indelible_-_Score_positions",R="select[mem>1000] rusage[mem=1000]",M=1000)
	sub(command)
	time.sleep(10)
	bsub.poll(sub.job_id)

@merge(score_positions,child_base+".sr_reads.counts.scored.denovo")
def call_denovos(input_files,output_file):
	mum_bam = sys.argv[2]
	dad_bam = sys.argv[3]
	child_calls = child_base+".sr_reads.counts.scored"
	command = "python %s/trio_caller.py %s %s %s" % (SRC_DIR,child_calls,mum_bam,dad_bam)
	print command
	sub = bsub("Indelible_-_DeNovo",R="select[mem>1000] rusage[mem=1000]",M=1000)
	sub(command)
	time.sleep(10)
	bsub.poll(sub.job_id)

# @transform(call_denovos,suffix(".sr_reads.counts.scored.denovo"),".sr_reads.counts.scored.denovo.annotated")
# def annotate_positions(input_file,output_file):
# 	#First blast the clipped segments:
# 	command = "python %s/blast.py %s" % (SRC_DIR,input_file)
# 	print command
# 	sub = bsub("Indelible_-_Blast_positions",R="select[mem>1000] rusage[mem=1000]",M=1000)
# 	sub(command)
# 	time.sleep(10)
# 	bsub.poll(sub.job_id)

# 	#Then annotate:
# 	command = "python %s/annotate.py %s" % (SRC_DIR,input_file)
# 	print command
# 	sub = bsub("Indelible_-_Annotate_positions",R="select[mem>1000] rusage[mem=1000]",M=1000)
# 	sub(command)
# 	time.sleep(10)
# 	bsub.poll(sub.job_id)


# @merge(annotate_positions,child_base+".sr_reads.counts.scored.annotated.denovo")
# def call_denovos(input_files,output_file):
# 	mum_bam = sys.argv[2]
# 	dad_bam = sys.argv[3]
# 	child_calls = child_base+".sr_reads.counts.scored.annotated"
# 	command = "python %s/trio_caller.py %s %s %s" % (SRC_DIR,child_calls,mum_bam,dad_bam)
# 	print command
# 	sub = bsub("Indelible_-_DeNovo",R="select[mem>1000] rusage[mem=1000]",M=1000)
# 	sub(command)
# 	time.sleep(10)
# 	bsub.poll(sub.job_id)

# @follows(call_denovos)
# def clean():
# 	command ="%s/run_clean.sh %s" % (SRC_DIR,trio_dir)
# 	os.popen(command).read()



pipeline_run(forcedtorun_tasks=[],multiprocess=3)
