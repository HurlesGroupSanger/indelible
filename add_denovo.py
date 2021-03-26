#!/usr/bin/env python

import csv
import sys
import os
import re

def build_denovo(denovo):

    csv.field_size_limit(sys.maxsize)
    denovo_dict = csv.DictReader(open(denovo), delimiter="\t")

    dnm = {}

    for val in denovo_dict:

        position_key = val['chrom'] + "_" + val['position']
        dnm[position_key] = {'mum_sr': val['mum_sr'],
                             'dad_sr': val['dad_sr'],
                             'mum_indel_context': val['mum_indel_context'],
                             'dad_indel_context': val['dad_indel_context'],
                             'mum_cov': val['mum_cov'],
                             'dad_cov': val['dad_cov']}

    return dnm


file_num = int(os.getenv('LSB_JOBINDEX'))
files = sys.argv[1]
file_list = open(files).readlines()

annotated_file = file_list[file_num - 1].rstrip()
annotated_file = annotated_file + ".new_annotated"

match = re.search(r"(\S*\.cram)\.counts\.scored.new_annotated", annotated_file)
file_prefix = match.group(1)
denovo_file = file_prefix + ".indelible.denovo.tsv"
denovo_out = file_prefix + ".indelible.denovo.revision.tsv"

dnm_db = build_denovo(denovo_file)
annotated_dict = csv.DictReader(open(annotated_file), delimiter="\t")

new_fieldnames = annotated_dict.fieldnames
new_fieldnames.extend(('mum_sr','dad_sr','mum_indel_context','dad_indel_context',
                      'mum_cov','dad_cov'))

denovo_writer = csv.DictWriter(open(denovo_out, 'w'), fieldnames=new_fieldnames, delimiter = "\t", lineterminator = "\n")
denovo_writer.writeheader()

for line in annotated_dict:

    key = line['chrom'] + "_" + line['position']
    if key in dnm_db:
        print(key)
        to_add = dnm_db[key]
        for k,v in to_add.items():
            line[k] = v
        denovo_writer.writerow(line)
