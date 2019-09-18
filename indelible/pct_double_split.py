import csv
import pandas as pd
import re
from indelible.aggregate_positions import sr_coverage

def split_positions(input_path_sr, input_agg, output_path, config):

    chr_dict = {}
    splitfile = open(input_path_sr, 'r')

    splitreader = csv.DictReader(splitfile, delimiter="\t", quoting=csv.QUOTE_NONE)
    header = ("chrom", "position", "coverage", "insertion_context", "deletion_context",
              "sr_total", "sr_total_long", "sr_total_short",
              "sr_long_5", "sr_short_5", "sr_long_3",
              "sr_short_3", "sr_entropy", "context_entropy",
              "entropy_upstream", "entropy_downstream", "sr_sw_similarity",
              "avg_avg_sr_qual", "avg_mapq", "seq_longest", "pct_double_split")

    outputfile = open(output_path, 'w')
    splitwriter = csv.DictWriter(outputfile, fieldnames=header, delimiter="\t", lineterminator="\n")
    splitwriter.writeheader()

    for row in splitreader:
        if row['chr'] == "hs37d5":
            continue
        if re.search("N", row['seq']):
            continue
        if not row['chr'] in chr_dict:
            chr_dict[row['chr']] = {}
            chr_dict[row['chr']][row['split_position']] = []
        if not row['split_position'] in chr_dict[row['chr']]:
            chr_dict[row['chr']][row['split_position']] = []

        chr_dict[row['chr']][row['split_position']].append(row)

    df = pd.read_csv(input_agg, sep="\t")

    for row in csv.DictReader(open(input_agg, 'r'), delimiter="\t"):

            sr_reads = chr_dict[row["chrom"]][row["position"]]

            sr_cov = sr_coverage(sr_reads, config["SHORT_SR_CUTOFF"])

            row["pct_double_split"] = float(sr_cov[7]) / float(row["sr_total"])

            splitwriter.writerow(row)