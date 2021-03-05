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
from indelible.bwa_runner import BWARunner
from indelible.blast_repeats import BlastRepeats
import csv

def decide_split_information(curr_data):

    left_total = sum(curr_data["left"])
    right_total = sum(curr_data["right"])

    # Decide direction:
    if (left_total * 0.5) > right_total:
        dir = "left"
    elif (right_total * 0.5) > left_total:
        dir = "right"
    else:
        dir = "uncer"

    # Decide longest sequence:
    fin_seq = ""
    for seq in curr_data["seq_longest"]:
        if len(seq) > len(fin_seq):
            fin_seq = seq

    returnable = {'lefts': left_total, 'rights': right_total, 'dir': dir, 'longest_seq': fin_seq}
    return returnable


def add_alignment_information(curr_bp, decisions, repeat_info):

    coord = curr_bp["coord"]
    if coord in repeat_info:

        curr_bp["otherside"] = "repeats_hit"
        curr_bp["mode"] = "BLAST_REPEAT"
        curr_bp["svtype"] = "INS_" + repeat_info[coord]["target_chrom"]
        curr_bp["size"] = "NA"
        curr_bp["aln_length"] = repeat_info[coord]["query_length"]
        curr_bp["otherside_found"] = "NA"
        curr_bp["is_primary"] = "NA"
        curr_bp["variant_coord"] = "NA"

    elif coord in decisions:

        curr_bp["otherside"] = decisions[coord]["otherside"]
        curr_bp["mode"] = decisions[coord]["mode"]
        curr_bp["svtype"] = decisions[coord]["svtype"]
        curr_bp["size"] = decisions[coord]["size"]
        curr_bp["aln_length"] = decisions[coord]["aln_length"]
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


def build_database(score_files, output_path, fasta, config, bwa_threads):

    # Pull stuff out of config that we need
    score_threshold = config['SCORE_THRESHOLD']
    REPEATdb = config['repeatdb']

    # Ensure bwa is path and find it:
    bwa_loc = find_bwa()
    if bwa_loc is None:
        raise Exception("bwa could not be found in path")

    # Do basic mashing together to get "allele frequencies"
    data = []
    allele_count = float(0)

    for file in open(score_files, 'r'):
        allele_count += 1
        file = file.rstrip()

        sample_id = os.path.basename(file).split(".")[0]

        frame = pandas.read_csv(
            file,
            sep="\t",
            header=0)

        is_pos = frame["prob_Y"] >= score_threshold
        frame = frame[is_pos][["chrom", "position","seq_longest","sr_long_5","sr_short_5","sr_long_3","sr_short_3"]]
        frame["left"] = frame["sr_long_5"] + frame["sr_short_5"]
        frame["right"] = frame["sr_long_3"] + frame["sr_short_3"]
        frame.drop(["sr_long_5", "sr_short_5", "sr_long_3", "sr_short_3"], axis=1)

        data.append(frame)

    data_joined = pandas.concat(data)
    data_joined["coord"] = data_joined["chrom"].astype(str) + "_" + data_joined["position"].astype(str)

    counts = data_joined["coord"].value_counts()

    final_frame = pandas.DataFrame()

    final_frame["coord"] = counts.index.values
    final_frame["counts"] = counts.values

    # Assess each coordinate for properties we need to examine for SV type:
    lefts = []
    rights = []
    dirs = []
    longest = []

    for index, row in final_frame.iterrows():

        curr_data = data_joined.loc[data_joined["coord"] == row["coord"]]
        curr_info = decide_split_information(curr_data)

        lefts.append(curr_info["lefts"])
        rights.append(curr_info["rights"])
        dirs.append(curr_info["dir"])
        longest.append(curr_info["longest_seq"])

    # Add all of this information to the final_frame
    final_frame["tot_left"] = lefts
    final_frame["tot_right"] = rights
    final_frame["dir"] = dirs
    final_frame["longest"] = longest

    # Convert the coordinate to positions:
    split = final_frame["coord"].str.split("_", n=1, expand=True)
    final_frame["chrom"] = split[0]
    final_frame["pos"] = split[1]
    final_frame["pos"] = final_frame["pos"].astype(int)
    final_frame = final_frame.sort_values(by=["chrom", "pos"])

    # Generate "Allele Frequencies"
    final_frame["pct"] = final_frame["counts"] / allele_count
    final_frame["tot"] = allele_count

    # Blast repeat database to mask:
    repeat_blast = BlastRepeats(output_path, final_frame, REPEATdb)
    repeat_info = repeat_blast.get_repeat_blast_info()

    # Generate bwa files and perform alignment:
    bwa_parser = BWARunner(final_frame, output_path, fasta, bwa_loc, bwa_threads)
    decisions = bwa_parser.get_decisions()

    # Write final database file:
    header = ["chrom", "pos", "pct", "counts", "tot", "otherside", "mode", "svtype", "size", "aln_length",
              "otherside_found", "is_primary", "variant_coord"]
    # header = ["coord","chrom", "pos", "otherside", "mode", "svtype", "size", "aln_length"]
    output_file = csv.DictWriter(open(output_path, 'w'), fieldnames=header, delimiter="\t",
                                 lineterminator="\n")

    for index, row in final_frame.iterrows():

        row_dict = row.to_dict()
        v = {}
        v["coord"] = row_dict["coord"]
        v["chrom"] = row_dict["chrom"]
        v["pos"] = row_dict["pos"]
        v["pct"] = row_dict["pct"]
        v["counts"] = row_dict["counts"]
        v["tot"] = row_dict["tot"]

        v = add_alignment_information(v, decisions, repeat_info)
        v.pop('coord', None)

        output_file.writerow(v)
