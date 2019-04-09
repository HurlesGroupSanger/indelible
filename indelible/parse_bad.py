"""
Author: Eugene Gardner
Affiliation: Wellcome Trust Sanger Institute

Parses aggregated positions and checks whether split reads are

Parameters
----------

Returns
-------

"""

from indelible.indelible_lib import *


def check_sr(input_path, output_path, config):

    denovo_file = csv.DictReader(open(input_path, 'r'), delimiter="\t")
    new_fieldnames = denovo_file.fieldnames
    new_fieldnames.extend(("total_one", "total_two", "total_reads", "percent"))
    output_file = csv.DictWriter(open(output_path, 'w'), fieldnames=new_fieldnames, delimiter="\t")
    output_file.writeheader()

    for v in csv.DictReader(open(input_path,'r'), delimiter="\t"):

        infile = bam_open(v['cram'])

        start_coord = int(v['position'])-config['WINDOW_SIZE']
        end_coord = int(v['position']) + config['WINDOW_SIZE']
        iter = infile.fetch(v['chrom'],start_coord,end_coord)

        # outfile = open(output_path, 'w')
        # outfile.write("chr\tsplit_position\tprime\tsplit_length\tseq\tqual\tmapq\tavg_sr_qual\treverse_strand\n")

        total_one = 0.0
        total_two = 0.0
        total = 0.0

        for s in iter:
            cigar = s.cigartuples
            if cigar is not None:
                refname = infile.getrname(s.tid)
                if refname is not "hs37d5":
                    if len(cigar) >= 2:
                        total += 1
                        if len(cigar) == 2:
                            # 5' Split Reads
                            total_one += 1
                        elif len(cigar) == 3:
                            total_two += 1
                            # sr = {}
                            # sr["chr"] = infile.getrname(s.tid)
                            # sr["split_position"] = s.pos
                            # sr["prime"] = 5
                            # sr["seq"], sr["qual"] = hard_clip(s.seq[0:cigar[0][1]], s.qual[0:cigar[0][1]],
                            #                                   config["HC_THRESHOLD"])
                            # sr["length"] = len(sr["seq"])
                            # sr["mapq"] = s.mapq
                            # sr["avg_sr_qual"] = average_quality(sr["qual"])
                            # sr["strand"] = s.is_reverse
                            # These are reads with Soft-clips on both sides, likely due to dropped quality at the end
                            # if cigar[0][0] == 4 and cigar[1][0] == 0 and cigar[2][0] == 4:
        v['total_one'] = total_one
        v['total_two'] = total_two
        v['total_reads'] = total
        v['percent'] = total_two / total
        output_file.writerow(v)