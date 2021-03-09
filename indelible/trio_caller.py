"""
    Author: Alejandro Sifrim & Eugene Gardner
    Affiliation: Wellcome Sanger Institute

    Script to call inheritance mode of an event

    Parameters
    ----------
    1) child score file
    2) mother bam
    3) father bam

    Returns
    -------
    1) File with putative denovo calls

"""

import os
import sys
import pysam
import csv
from indelible.indelible_lib import *
import scipy.stats as sp
import math
from .coverage_calculator import CoverageCalculator


def compute_stats(cov_calc, chrom, position):

    res = {}
    if cov_calc.input_bam is None:
        res = {'coverage': "NA", 'sr_context': "NA", 'indel_context': {'deletions': "NA", 'insertions': ""}}
    else:
        res['coverage'] = cov_calc.calculate_coverage(chrom, position)
        res['indel_context'] = cov_calc.reads_with_indels_in_neighbourhood(chrom, position)
        res['sr_context'] = cov_calc.split_reads_in_neighbourhood(chrom, position)

    return res


def denovo_caller_trio(child_input, mother_bam, father_bam, output_path, config):

    csv.field_size_limit(sys.maxsize)

    mum_cov = CoverageCalculator({}, mother_bam, "", config)
    dad_cov = CoverageCalculator({}, father_bam, "", config)

    scored_file = csv.DictReader(open(child_input, 'r'), delimiter="\t")
    new_fieldnames = scored_file.fieldnames
    new_fieldnames.extend(("mum_sr", "dad_sr", 'mum_indel_context', 'dad_indel_context', "mum_cov", "dad_cov"))
    output_file = csv.DictWriter(open(output_path, 'w'), fieldnames=new_fieldnames, delimiter="\t", lineterminator="\n")
    output_file.writeheader()

    for v in csv.DictReader(open(child_input, 'r'), delimiter="\t"):
        mum_stats = compute_stats(mum_cov, v['chrom'], int(v['position']))
        dad_stats = compute_stats(dad_cov, v['chrom'], int(v['position']))

        v['mum_sr'] = mum_stats['sr_context']
        v['dad_sr'] = dad_stats['sr_context']
        v['mum_indel_context'] = mum_stats['indel_context']['deletions'] + mum_stats['indel_context']['insertions']
        v['dad_indel_context'] = dad_stats['indel_context']['deletions'] + dad_stats['indel_context']['insertions']
        v['mum_cov'] = mum_stats['coverage']
        v['dad_cov'] = dad_stats['coverage']

        # Only filter here if we have BOTH mom and dad bams
        if mother_bam is None or father_bam is None:
            output_file.writerow(v)
        else:
            if v['mum_sr'] <= config['SR_THRESHOLD'] and v['dad_sr'] <= config['SR_THRESHOLD'] and v['mum_cov'] >= \
                    config['COV_THRESHOLD'] and v['dad_cov'] >= config['COV_THRESHOLD']:
                output_file.writerow(v)
