"""
    Author: Eugene Gardner
    Affiliation: Wellcome Sanger Institute

    Script to build to AF database necessary for annotation:

    Parameters
    ----------
    1) input fof of all *.scored files from dataset
    2) output coordinate file
    3) score threshold from the config.yml file

    Returns
    -------
    1) Coordinate file with allele frequencies

"""

import pandas
import os
import pysam
from indelible.bwa_runner import BWARunner
from indelible.blast_repeats import BlastRepeats
import csv


def add_alignment_information(curr_bp, decisions, repeat_info):
    name = "%s_%s" % (curr_bp["chrom"], curr_bp["pos"])
    if name in repeat_info:
        curr_bp["otherside"] = "repeats_hit"
        curr_bp["mode"] = "BLAST_REPEAT"
        curr_bp["svtype"] = "INS_" + repeat_info[name]["target_chrom"]
        curr_bp["size"] = "NA"
        curr_bp["aln_length"] = repeat_info[name]["query_length"]
        curr_bp["otherside_found"] = "NA"
        curr_bp["is_primary"] = "NA"
        curr_bp["variant_coord"] = "NA"
    elif name in decisions:
        curr_bp["otherside"] = decisions[name]["otherside"]
        curr_bp["mode"] = decisions[name]["mode"]
        curr_bp["svtype"] = decisions[name]["svtype"]
        curr_bp["size"] = decisions[name]["size"]
        curr_bp["aln_length"] = decisions[name]["aln_length"]
        curr_bp["otherside_found"] = "NA"
        curr_bp["is_primary"] = "NA"
        curr_bp["variant_coord"] = "NA"
        ## Search for otherside in actual hits
        if curr_bp["otherside"] != "NA":
            if curr_bp["otherside"] in decisions:
                other_coord = curr_bp["otherside"].split("_")
                if int(other_coord[1]) < curr_bp["pos"]:
                    curr_bp["is_primary"] = "false"
                    curr_bp["variant_coord"] = "%s:%s-%s" % (curr_bp["chrom"],other_coord[1], curr_bp["pos"])
                else:
                    curr_bp["is_primary"] = "true"
                    curr_bp["variant_coord"] = "%s:%s-%s" % (curr_bp["chrom"], curr_bp["pos"], other_coord[1])
    else:
        #"otherside", "mode", "svtype", "size", "aln_length"
        curr_bp["otherside"] = "NA"
        curr_bp["mode"] = "FAIL_ALIGNMENT"
        curr_bp["svtype"] = "UNK"
        curr_bp["size"] = "NA"
        curr_bp["aln_length"] = "NA"
        curr_bp["otherside_found"] = "NA"
        curr_bp["is_primary"] = "NA"
        curr_bp["variant_coord"] = "NA"
    return curr_bp


def find_bwa():
    paths = os.environ.get("PATH").split(":")
    for path in paths:
        search_path = path + "/bwa"
        if os.path.isfile(search_path):
            return(search_path)


def decide_direction(left, right):
    if (left * 0.5) > right:
        dir = "left"
    elif (right * 0.5) > left:
        dir = "right"
    else:
        dir = "uncer"
    return dir


def build_database(score_files, output_path, fasta, config, bwa_threads):

    # Pull stuff out of arguments/config that we need
    fasta = pysam.FastaFile(fasta)
    score_threshold = config['SCORE_THRESHOLD']
    REPEATdb = config['repeatdb']

    # Ensure bwa is path and find it:
    bwa_loc = find_bwa()
    if bwa_loc is None:
        raise Exception("bwa could not be found in path")

    # Do basic mashing together to get "allele frequencies"
    data = []
    allele_count = float(0)

    # Open all scored files and mash them together while filtering low quality variants
    for file in open(score_files, 'r'):
        allele_count += 1
        file = file.rstrip()
        frame = pandas.read_csv(
            file,
            sep="\t",
            header=0)
        is_pos = frame["prob_Y"] >= score_threshold
        frame = frame[is_pos][["chrom", "position", "seq_longest", "sr_long_5", "sr_short_5", "sr_long_3", "sr_short_3"]]
        frame["left"] = frame["sr_long_5"] + frame["sr_short_5"]
        frame["right"] = frame["sr_long_3"] + frame["sr_short_3"]
        frame.drop(["sr_long_5", "sr_short_5", "sr_long_3", "sr_short_3"], axis=1)
        data.append(frame)

    data_joined = pandas.concat(data)
    data_joined["coord"] = data_joined["chrom"].astype(str) + "_" + data_joined["position"].astype(str)

    final_frame = data_joined.groupby('coord').agg(chrom = ('chrom','first'),
                                                   pos = ('position','first'),
                                                   counts=('coord', len),
                                                   tot_left = ('left','sum'),
                                                   tot_right = ('right','sum'),
                                                   longest = ('seq_longest','max'))

    # Set direction value:
    final_frame['dir'] = final_frame.apply(lambda x: decide_direction(x.tot_left, x.tot_right), axis=1)

    # Generate "Allele Frequencies"
    final_frame["pct"] = final_frame["counts"] / allele_count
    final_frame["tot"] = allele_count

    # Blast repeat database to mask:
    repeat_blast = BlastRepeats(output_path, final_frame, REPEATdb)
    repeat_info = repeat_blast.get_repeat_blast_info()

    # Generate bwa files and perform alignment:
    bwa_parser = BWARunner(final_frame, output_path, fasta, bwa_loc, bwa_threads)
    decisions = bwa_parser.get_decisions()

    final_frame["chrom"] = pandas.Categorical(final_frame["chrom"], fasta.references)
    final_frame = final_frame.sort_values(by=["chrom", "pos"])

    # Write final database file:
    header = ["chrom", "pos", "pct", "counts", "tot", "otherside", "mode", "svtype", "size", "aln_length",
              "otherside_found", "is_primary", "variant_coord"]
    output_file = csv.DictWriter(open(output_path, 'w'), fieldnames=header, delimiter="\t",
                                 lineterminator="\n")

    for index, row in final_frame.iterrows():

        row_dict = row.to_dict()
        v = {}
        v["chrom"] = row_dict["chrom"]
        v["pos"] = row_dict["pos"]
        v["pct"] = row_dict["pct"]
        v["counts"] = row_dict["counts"]
        v["tot"] = row_dict["tot"]

        v = add_alignment_information(v, decisions, repeat_info)

        output_file.writerow(v)

    fasta.close()
