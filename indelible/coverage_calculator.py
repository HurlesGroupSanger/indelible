"""
    Author: Eugene Gardner
    Affiliation: Wellcome Sanger Institute

    Object which handles various coverage calculations for aggregate_positions.py and trio_caller.py

    ----------
    Parameters
    1.) chr_dictionary of split reads from fetch
    2.) input bam path
    3.) input path
    4.) config file of parameters for all of Indelible


"""

import subprocess
from indelible.indelible_lib import *
import pysam
from Bio import SeqIO


class CoverageCalculator:

    def __init__(self, chr_dict, input_bam, input_path, config):

        self.chr_dict = chr_dict
        self.input_bam = input_bam
        self.__input_path = input_path
        opener = bam_open(self.input_bam)
        self.__bam_file = opener["reader"]
        self.__is_cram = opener["is_cram"]
        self.__tabix_file = self.__input_path + ".bg.gz"
        self.__config = config
        self.__use_bam = True

        self.__decide_coverage_method()

    def __decide_coverage_method(self):

        count = 0

        for chrom in self.chr_dict:
            for position in self.chr_dict[chrom]:
                pos = int(position)
                if len(self.chr_dict[chrom][position]) >= self.__config["MINIMUM_SR_COVERAGE"]:
                    count += 1

        # It is much faster to calculate bam-wide coverage all at once:
        # (__calculate_coverage_bam + __coverage_at_position_tabix)
        # than doing the bam iterator method (__coverage_at_position_pileup) individual for each coordinate if the
        # number of reads is excessive. I haven't exactly pinpointed when the number of coordinates to examine reaches
        # this threshold, but I don't think it does too much damage anyway as the time to calculate WES coverage is
        # relatively low

        if count <= 30000:
            self.__use_bam = True
        else:
            self.__use_bam = False
            output_file = self.__input_path + ".bg"
            self.__calculate_coverage_bam(output_file)

    def __calculate_coverage_bam(self, output_file):

        cmd = "bedtools genomecov -bg -ibam " + self.input_bam
        proc = subprocess.Popen(cmd, shell=True, stdout=open(output_file, "w"), stderr=subprocess.PIPE)
        stderr = proc.communicate()

        if proc.returncode != 0:
            print("The following bedtools command failed:")
            print(cmd)
            print("STDERROR follows\n")
            print(stderr.decode('utf-8'))
            raise

        bgzip_and_tabix(output_file)

    def calculate_coverage(self, chr, pos):

        if self.__use_bam is True:
            return self.__coverage_at_position_pileup(chr, pos)
        else:
            return self.__coverage_at_position_tabix(chr, pos)

    def __coverage_at_position_pileup(self, chr, pos):

        if pos <= 0:
            pileup = self.__bam_file.pileup(chr, 0, pos + 1, truncate=True, multiple_iterators=False)
        else:
            pileup = self.__bam_file.pileup(chr, pos - 1, pos + 1, truncate=True, multiple_iterators=False)

        cov = 0

        for pileupcolumn in pileup:
            if pileupcolumn.pos == pos:
                cov = pileupcolumn.n
                break

        return cov

    def __coverage_at_position_tabix(self, chr, pos):

        cov = 0
        tbx = pysam.TabixFile(self.__tabix_file)

        if pos <= 1:
            itr = tbx.fetch(chr, 1, 2, parser=pysam.asTuple())
        else:
            itr = tbx.fetch(chr, pos - 1, pos, parser=pysam.asTuple())

        for row in itr:
            cov = row[3]
        return cov

    def reads_with_indels_in_neighbourhood(self, chrom, pos):
        counts = {"insertions": 0, "deletions": 0}
        window_size = self.__config['WINDOW_SIZE']
        fetch_start = pos - (window_size / 2)
        fetch_end = pos + (window_size / 2)
        if fetch_start <= 0:
            fetch_start = 0

        for alignedread in self.__bam_file.fetch(chrom, fetch_start, fetch_end):
            cigar = alignedread.cigar
            cigar_types = [c[0] for c in cigar]
            if 1 in cigar_types: counts["insertions"] += 1
            if 2 in cigar_types: counts["deletions"] += 1
        return counts

    def split_reads_in_neighbourhood(self, chrom, pos):
        count = 0
        window_size = self.__config["WINDOW_SIZE"]

        start = pos - (window_size / 2)
        end = pos + (window_size / 2)

        if start <= 0:
            start = 1

        for s in self.__bam_file.fetch(chrom, start, end):

            cigar = s.cigar
            sr = {}
            if len(cigar) == 2:
                # 5' Split Reads
                if cigar[0][0] == 4 and cigar[1][0] == 0:  # 4 = SOFT CLIP/ 0 = MATCH

                    sr["chr"] = self.__bam_file.getrname(s.tid)
                    sr["split_position"] = s.pos
                    sr["prime"] = 5
                    sr["seq"], sr["qual"] = hard_clip(s.seq[0:cigar[0][1]], s.qual[0:cigar[0][1]],
                                                      self.__config['HC_THRESHOLD'])
                    sr["length"] = len(sr["seq"])
                    sr["mapq"] = s.mapq
                    sr["avg_sr_qual"] = average_quality(sr["qual"])

                # 3' Split Reads
                if cigar[0][0] == 0 and cigar[1][0] == 4:
                    sr["chr"] = self.__bam_file.getrname(s.tid)
                    sr["split_position"] = s.pos + cigar[0][1]
                    sr["prime"] = 3
                    seq = s.seq[-cigar[1][1]:]
                    qual = s.qual[-cigar[1][1]:]
                    seq, qual = hard_clip(seq[::-1], qual[::-1], self.__config['HC_THRESHOLD'])
                    sr["seq"] = seq[::-1]
                    sr["qual"] = qual[::-1]
                    sr["length"] = len(sr["seq"])
                    sr["mapq"] = s.mapq
                    sr["avg_sr_qual"] = average_quality(sr["qual"])

            elif len(cigar) == 3:

                # These are reads with Soft-clips on both sides, likely due to dropped quality at the end
                if cigar[0][0] == 4 and cigar[1][0] == 0 and cigar[2][0] == 4:
                    # 1st split-segment
                    if cigar[0][1] >= cigar[2][1]:
                        sr["chr"] = self.__bam_file.getrname(s.tid)
                        sr["split_position"] = s.pos
                        sr["prime"] = 5
                        sr["seq"], sr["qual"] = hard_clip(s.seq[0:cigar[0][1]], s.qual[0:cigar[0][1]],
                                                          self.__config['HC_THRESHOLD'])
                        sr["length"] = len(sr["seq"])
                        sr["mapq"] = s.mapq
                        sr["avg_sr_qual"] = average_quality(sr["qual"])
                    # 2nd split segment
                    else:
                        sr["chr"] = self.__bam_file.getrname(s.tid)
                        sr["split_position"] = s.pos + cigar[1][
                            1]  # alignment_start + length of matching segment = start of 3' split segment
                        sr["prime"] = 3
                        seq = s.seq[-cigar[2][1]:]
                        qual = s.qual[-cigar[2][1]:]
                        seq, qual = hard_clip(seq[::-1], qual[::-1], self.__config['HC_THRESHOLD'])
                        sr["seq"] = seq[::-1]
                        sr["qual"] = qual[::-1]
                        sr["length"] = len(sr["seq"])
                        sr["mapq"] = s.mapq
                        sr["avg_sr_qual"] = average_quality(sr["qual"])
            if sr:
                if sr['length'] > self.__config['MINIMUM_LENGTH_SPLIT_READ']:
                    count += 1
        return count
