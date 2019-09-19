import pysam
from collections import Counter
import sys
import csv
import math
from collections import defaultdict
from pyfaidx import Fasta
import numpy
import re
import subprocess

BASE_QUALITIES = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"


def bam_open(bam_file):

    is_cram = False

    if bam_file is None:
        bam_reader = None
    elif 'bam' in bam_file:
        bam_reader = pysam.Samfile(bam_file, 'rb')
    elif 'cram' in bam_file:
        bam_reader = pysam.Samfile(bam_file, 'rc')
        is_cram = True

    return {"reader":bam_reader, "is_cram": is_cram}


def hard_clip(seq, qual, threshold=10):
    res_seq = ""
    res_qual = ""

    offset = 0
    q = BASE_QUALITIES.index(qual[offset])
    while q < threshold and offset < len(seq):
        offset += 1
        if offset < len(seq):
            q = BASE_QUALITIES.index(qual[offset])

    return [seq[offset:], qual[offset:]]


def average_quality(qual):
    if len(qual) != 0:
        s = 0
        for q in qual:
            s += BASE_QUALITIES.index(q)
        return float(s) / len(qual)
    else:
        return 0


def entropy(s):
    p, lns = Counter(s), float(len(s))
    return -sum(count / lns * math.log(count / lns, 2) for count in list(p.values()))


def read_database(path):
    db = defaultdict(dict)
    for v in csv.DictReader(open(path, 'r'), fieldnames=("chrom", "position", "maf", "count", "total"), delimiter="\t"):
        db[v["chrom"]][int(v["position"])] = float(v["maf"])
    return db

def bgzip_and_tabix(path):

    # bgzip
    sigint = subprocess.call(["bgzip", "-f", path])
    if sigint != 0:
        raise Exception("bgzip on the file " + path + " did not run properly... Exiting!")
    # tabix index
    sigint = subprocess.call(["tabix", "-f", "-p", "bed", path + ".gz"])
    if sigint != 0:
        raise Exception("tabix on the file " + path + ".gz did not run properly... Exiting!")