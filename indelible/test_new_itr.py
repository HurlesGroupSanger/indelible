"""
    Author: Eugene Gardner
    Affiliation: Wellcome Sanger Institute

    Object which handles various coverage calculations for aggregate_positions.py and trio_caller.py

    ----------
    Parameters
    1.) chr_dictionary of split reads from fetch
    2.) input bam path
    3.) config file of parameters for all of Indelible


"""

import pybedtools as bedtools
import subprocess
from indelible.indelible_lib import *
import pysam
import sys

bam_file = sys.argv[1]

bam_reader = pysam.Samfile(bam_file, 'rb')

for read in bam_reader.fetch():
    cigar = read.cigartuples
    cigar_types = [c[0] for c in cigar]

    insertion = 0
    deletion = 0

    if 1 in cigar_types: insertion += 1
    if 2 in cigar_types: deletion += 1

    print()