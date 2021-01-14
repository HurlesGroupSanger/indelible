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


def build_database(score_files, output_file, score_threshold):

    output_writer = open(output_file, "w")

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
        frame = frame[is_pos][["chrom", "position"]]
        data.append(frame)

    data_joined = pandas.concat(data)
    data_joined["coord"] = data_joined["chrom"].astype(str) + ":" + data_joined["position"].astype(str)

    counts = data_joined["coord"].value_counts()

    final_frame = pandas.DataFrame()

    final_frame["coord"] = counts.index.values
    final_frame["counts"] = counts.values

    split = final_frame["coord"].str.split(":", n=1, expand=True)
    final_frame["chrom"] = split[0]
    final_frame["pos"] = split[1]

    final_frame = final_frame.drop(["coord"], axis=1)

    final_frame["pos"] = final_frame["pos"].astype(int)
    final_frame = final_frame.sort_values(by=["chrom", "pos"])

    final_frame["pct"] = final_frame["counts"] / allele_count

    final_frame["tot"] = allele_count

    header = ["chrom", "pos", "pct", "counts", "tot"]
    csv_obj = final_frame.to_csv(sep="\t", index=False, header=False, columns=header)

    output_writer.write(csv_obj)
