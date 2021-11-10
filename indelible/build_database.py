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
import numpy as np
import os
import pysam
from difflib import SequenceMatcher
from indelible.bwa_runner import BWARunner
from indelible.blast_repeats import BlastRepeats
import csv
import statistics
import random


def build_priors(prior_file, final_frame, fasta):
    priors_frame = pandas.read_csv(prior_file,
                                   sep="\t",
                                   header=None,
                                   index_col=False,
                                   low_memory=False,
                                   names = ("chrom", "pos", "pct", "counts", "tot", "coverage", "otherside", "mode", "svtype",
                                            "size", "aln_length", "otherside_found", "is_primary", "variant_coord"),
                                   dtype = {'chrom': pandas.CategoricalDtype(categories = fasta.references, ordered=True), 'pos': np.int, 'pct': np.float, 'counts': np.int,
                                            'tot': np.int, 'coverage': np.float, 'otherside': np.str, 'mode': np.str, 'svtype': np.str,
                                            'size': np.float, 'aln_length': np.float, 'otherside_found': "boolean", 'is_primary': np.str,
                                            'variant_coord': np.str})
    priors_frame["coord"] = priors_frame["chrom"].astype(str) + ":" + priors_frame["pos"].astype(str)
    priors_frame = priors_frame.set_index("coord")

    # Flag sites that have already been identified
    already_found = final_frame.index.to_list()
    priors_frame["already_found"] = priors_frame.index.isin(already_found)
    # priors_frame = priors_frame.drop(index=already_found)

    return(priors_frame)


def add_alignment_information(curr_bp, decisions, repeat_info):
    name = "%s:%s" % (curr_bp["chrom"], curr_bp["pos"])
    if name in repeat_info:
        curr_bp["otherside"] = "repeats_hit"
        curr_bp["mode"] = "BLAST_REPEAT"
        curr_bp["svtype"] = "INS_" + repeat_info[name]["target_chrom"]
        curr_bp["size"] = "na"
        curr_bp["aln_length"] = repeat_info[name]["query_length"]
        curr_bp["otherside_found"] = False
        curr_bp["is_primary"] = "nan"
        curr_bp["variant_coord"] = "nan"
    elif name in decisions:
        curr_bp["otherside"] = decisions[name]["otherside"]
        curr_bp["mode"] = decisions[name]["mode"]
        curr_bp["svtype"] = decisions[name]["svtype"]
        curr_bp["size"] = decisions[name]["size"]
        curr_bp["aln_length"] = decisions[name]["aln_length"]
        curr_bp["otherside_found"] = False
        curr_bp["is_primary"] = "nan"
        curr_bp["variant_coord"] = "nan"
        ## Search for otherside in actual hits
        if curr_bp["otherside"] != "NA":
            if curr_bp["otherside"] in decisions:
                other_coord = curr_bp["otherside"].split(":")
                curr_bp["otherside_found"] = True
                if int(other_coord[1]) < curr_bp["pos"]:
                    curr_bp["is_primary"] = False
                    curr_bp["variant_coord"] = "%s:%s-%s" % (curr_bp["chrom"], other_coord[1], curr_bp["pos"])
                else:
                    curr_bp["is_primary"] = False
                    curr_bp["variant_coord"] = "%s:%s-%s" % (curr_bp["chrom"], curr_bp["pos"], other_coord[1])
    else:
        # "otherside", "mode", "svtype", "size", "aln_length"
        curr_bp["otherside"] = "NA"
        curr_bp["mode"] = "FAIL_ALIGNMENT"
        curr_bp["svtype"] = "UNK"
        curr_bp["size"] = "nan"
        curr_bp["aln_length"] = "nan"
        curr_bp["otherside_found"] = False
        curr_bp["is_primary"] = "nan"
        curr_bp["variant_coord"] = "nan"
    return curr_bp


def find_bwa():
    paths = os.environ.get("PATH").split(":")
    for path in paths:
        search_path = path + "/bwa"
        if os.path.isfile(search_path):
            return(search_path)


def seq_similarity(seq1, seq2):
    seq = SequenceMatcher(a=seq1, b=seq2)
    return seq.ratio()


def decide_longest(seqs):

    if len(seqs) == 1:
        return seqs.to_list()[0]
    else:

        returnable = None
        for s in sorted(seqs, key = len, reverse=True):
            if len(seqs) >= 20:
                # Set seed for random to ensure identical results when variant calling:
                random.seed(1234)
                seq_frame = pandas.DataFrame(data={'orig': np.repeat(s, 20), 'old': random.sample(seqs.to_list(), 20)})
            else:
                seq_frame = pandas.DataFrame(data={'orig': np.repeat(s, len(seqs)), 'old': seqs})

            seq_frame['similarity'] = seq_frame.apply(lambda x: seq_similarity(x[0], x[1]), axis=1)
            seq_frame[seq_frame['similarity'] != 1] # I don't know what this does, but am too much of a coward to delete it
            if statistics.mean(seq_frame['similarity']) >= 0.6:
                returnable = s
                break

        if returnable is None:
            return max(seqs)
        else:
            return returnable


def decide_direction(left, right):
    if (left * 0.5) > right:
        dir = "left"
    elif (right * 0.5) > left:
        dir = "right"
    else:
        dir = "uncer"
    return dir


def build_database(score_files, output_path, fasta, config, priors, old_maf, bwa_threads):

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
            header=0,
        dtype = {'coverage': np.float})
        # Calculating mean coverage and then adding a fairly lenient threshold. This prevents areas with a ton of reads
        # from "poisoning the well" when deciding which read to align later in this algorithm
        # The 'A' base is a placeholder that will automatically fail alignment if no other sites are found
        coverage_cutoff = statistics.mean(frame["coverage"]) * 10
        frame.loc[frame['coverage'] >= coverage_cutoff, 'seq_longest'] = 'A'

        is_pos = frame["prob_Y"] >= score_threshold
        frame = frame[is_pos][["chrom", "position", "coverage", "seq_longest", "sr_long_5", "sr_short_5", "sr_long_3", "sr_short_3"]]
        frame["left"] = frame["sr_long_5"] + frame["sr_short_5"]
        frame["right"] = frame["sr_long_3"] + frame["sr_short_3"]
        frame.drop(["sr_long_5", "sr_short_5", "sr_long_3", "sr_short_3"], axis=1)
        data.append(frame)

    data_joined = pandas.concat(data)
    data_joined["coord"] = data_joined["chrom"].astype(str) + ":" + data_joined["position"].astype(str)

    print("Total number of sites to iterate through: %s" % len(data_joined.groupby('coord').agg(counts=('coord', len))))

    final_frame = data_joined.groupby('coord').agg(chrom = ('chrom','first'),
                                                   pos = ('position','first'),
                                                   counts=('coord', len),
                                                   tot_left = ('left','sum'),
                                                   tot_right = ('right','sum'),
                                                   coverage = ('coverage','mean'),
                                                   longest = ('seq_longest',decide_longest))

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

    final_frame["chrom"] = pandas.Categorical(final_frame["chrom"].astype(str), categories = fasta.references, ordered=True)
    final_frame = final_frame.sort_values(by=["chrom", "pos"])

    # Write final database file:
    header = ["chrom", "pos", "pct", "counts", "tot", "coverage", "otherside", "mode", "svtype", "size", "aln_length",
              "otherside_found", "is_primary", "variant_coord"]
    output_maf_db = open(output_path, 'w')
    output_file = csv.DictWriter(output_maf_db, fieldnames=header, delimiter="\t", lineterminator="\n")

    # Build priors frame if required:
    if priors is not None:
        priors_frame = build_priors(priors, final_frame, fasta)

    for index, row in final_frame.iterrows():

        row_dict = row.to_dict()
        v = {}
        v["chrom"] = row_dict["chrom"]
        v["pos"] = row_dict["pos"]
        v["coverage"] = row_dict["coverage"]

        # Update MAF information from priors frame if requested
        # Otherwise just use precalculated value from current data
        if priors is not None and old_maf is True:
            if index in priors_frame.index:
                prior_record = priors_frame.loc[index]
                v["pct"] = prior_record["pct"]
                v["counts"] = prior_record["counts"]
                v["tot"] = prior_record["tot"]
            else:
                v["pct"] = row_dict["pct"]
                v["counts"] = row_dict["counts"]
                v["tot"] = row_dict["tot"]
        else:
            v["pct"] = row_dict["pct"]
            v["counts"] = row_dict["counts"]
            v["tot"] = row_dict["tot"]

        v = add_alignment_information(v, decisions, repeat_info)

        output_file.writerow(v)
    output_maf_db.flush()

    ## Check priors and append to the bottom:
    if priors_frame is not None:
        priors_frame = priors_frame[priors_frame["already_found"] == False]
        priors_frame = priors_frame.drop("already_found", axis = 1)
        for index, row in priors_frame.iterrows():
            row_dict = row.to_dict()
            output_file.writerow(row_dict)

    output_maf_db.close()
    fasta.close()
