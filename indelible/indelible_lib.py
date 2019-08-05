import pysam
from collections import Counter
import sys
import csv
import math
from collections import defaultdict
from pyfaidx import Fasta
import swalign
import numpy
import re

BASE_QUALITIES = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"

def bam_open(bam_file):

    if bam_file is None:
        bam_reader = None
    elif 'bam' in bam_file:
        bam_reader = pysam.Samfile(bam_file, 'rb')
    elif 'cram' in bam_file:
        bam_reader = pysam.Samfile(bam_file, 'rc')

    return bam_reader

def reads_with_indels_in_neighbourhood(bam_file,chrom,pos,config):
    counts = {"insertions": 0, "deletions": 0}
    window_size = config['WINDOW_SIZE']
    fetch_start = pos-(window_size/2)
    fetch_end = pos+(window_size/2)
    if fetch_start <= 0:
        fetch_start = 0

    for alignedread in bam_file.fetch(chrom,fetch_start,fetch_end):
        cigar = alignedread.cigar
        cigar_types = [c[0] for c in cigar]
        if 1 in cigar_types: counts["insertions"] += 1
        if 2 in cigar_types: counts["deletions"] += 1
    return counts

def hard_clip(seq,qual,threshold=10):
    res_seq = ""
    res_qual = ""

    offset = 0
    q = BASE_QUALITIES.index(qual[offset])
    while q < threshold and offset < len(seq):
        offset += 1
        if offset < len(seq):
            q = BASE_QUALITIES.index(qual[offset])

    return [seq[offset:],qual[offset:]]

def average_quality(qual):
    if len(qual) != 0:
        s = 0
        for q in qual:
            s += BASE_QUALITIES.index(q)
        return float(s)/len(qual)
    else:
        return 0

def entropy(s):
    p, lns = Counter(s), float(len(s))
    return -sum( count/lns * math.log(count/lns, 2) for count in p.values())


def read_database(path):
    db = defaultdict(dict)
    for v in csv.DictReader(open(path,'r'),fieldnames=("chrom","position","maf","count","total"),delimiter="\t"):
        db[v["chrom"]][int(v["position"])] = float(v["maf"])
    return db

def split_reads_in_neighbourhood(bam_file,chrom,pos,config):
    count = 0
    window_size = config["WINDOW_SIZE"]

    for s in bam_file.fetch(chrom,pos-(window_size/2),pos+(window_size/2)):

        cigar = s.cigar
        sr = {}
        if len(cigar) == 2:
            #5' Split Reads
            if cigar[0][0] == 4 and cigar[1][0] == 0: # 4 = SOFT CLIP/ 0 = MATCH

                sr["chr"] = bam_file.getrname(s.tid)
                sr["split_position"] = s.pos
                sr["prime"] = 5
                sr["seq"],sr["qual"] = hard_clip(s.seq[0:cigar[0][1]],s.qual[0:cigar[0][1]],config['HC_THRESHOLD'])
                sr["length"] = len(sr["seq"])
                sr["mapq"] = s.mapq
                sr["avg_sr_qual"] = average_quality(sr["qual"])

            #3' Split Reads
            if cigar[0][0] == 0 and cigar[1][0] == 4:

                sr["chr"]  = bam_file.getrname(s.tid)
                sr["split_position"] = s.pos + cigar[0][1]
                sr["prime"] = 3
                seq = s.seq[-cigar[1][1]:]
                qual = s.qual[-cigar[1][1]:]
                seq,qual = hard_clip(seq[::-1],qual[::-1],config['HC_THRESHOLD'])
                sr["seq"] = seq[::-1]
                sr["qual"] = qual[::-1]
                sr["length"] = len(sr["seq"])
                sr["mapq"] = s.mapq
                sr["avg_sr_qual"] = average_quality(sr["qual"])

        elif len(cigar) == 3:

            # These are reads with Soft-clips on both sides, likely due to dropped quality at the end
            if cigar[0][0] == 4 and cigar[1][0] == 0 and cigar[2][0] == 4:
                #1st split-segment
                if cigar[0][1] >= cigar[2][1]:
                    sr["chr"] = bam_file.getrname(s.tid)
                    sr["split_position"] = s.pos
                    sr["prime"] = 5
                    sr["seq"],sr["qual"] = hard_clip(s.seq[0:cigar[0][1]],s.qual[0:cigar[0][1]],config['HC_THRESHOLD'])
                    sr["length"] = len(sr["seq"])
                    sr["mapq"] = s.mapq
                    sr["avg_sr_qual"] = average_quality(sr["qual"])
                #2nd split segment
                else:
                    sr["chr"]  = bam_file.getrname(s.tid)
                    sr["split_position"] = s.pos + cigar[1][1] #alignment_start + length of matching segment = start of 3' split segment
                    sr["prime"] = 3
                    seq = s.seq[-cigar[2][1]:]
                    qual = s.qual[-cigar[2][1]:]
                    seq,qual = hard_clip(seq[::-1],qual[::-1],config['HC_THRESHOLD'])
                    sr["seq"] = seq[::-1]
                    sr["qual"] = qual[::-1]
                    sr["length"] = len(sr["seq"])
                    sr["mapq"] = s.mapq
                    sr["avg_sr_qual"] = average_quality(sr["qual"])
        if sr:
            if sr['length'] > config['MINIMUM_LENGTH_SPLIT_READ']:
                count += 1
    return count

