"""
    Author: Eugene Gardner
    Affiliation: Wellcome Sanger Institute

    Object which handles various coverage calculations for aggregate_positions.py

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


class CoverageCalculator:

    def __init__ (self, chr_dict, input_bam, config):

        self.chr_dict = chr_dict
        self.input_bam = input_bam
        self.__bam_file = bam_open(self.input_bam)
        self.__tabix_file = self.input_bam + ".bg.gz"
        self.config = config
        self.__use_bam = True

        self.__decide_coverage_method()

        print self.__use_bam

    def __decide_coverage_method(self):

        count = 0

        for chrom in self.chr_dict:
            for position in self.chr_dict[chrom]:
                pos = int(position)
                if len(self.chr_dict[chrom][position]) >= self.config["MINIMUM_SR_COVERAGE"]:
                    count += 1

        # It is much faster to calculate bam-wide coverage all at once:
        # (__calculate_coverage_bam + __coverage_at_position_tabix)
        # than doing the bam iterator method (__coverage_at_position_pileup) individual for each coordinant if the
        # number of reads is excessive. I haven't exactly pinpointed when the number of coordinates to examine reaches
        # this threshold, but I don't think it does too much damage anyway as the time to calculate WES coverage is
        # relatively low

        if count <= 30000:
            print "Low-ish number of reads (" + str(count) + ")... Using pysam to extract coverage..."
            self.__use_bam = True
        else:
            print "Excessive numbers of reads (" + str(count) + ")... Using bedtools + tabix method to calculate per-site coverage..."
            self.__use_bam = False
            output_file = self.input_bam + ".bg"
            self.__calculate_coverage_bam(output_file)

    def __calculate_coverage_bam(self, output_file):

        bt = bedtools.BedTool(self.input_bam)
        bg_file = bt.genome_coverage(bg=True)
        bg_file.moveto(output_file)

        # bgzip
        sigint = subprocess.call(["bgzip", output_file])
        if sigint != 0:
            raise Exception("bgzip on the file " + output_file + " did not run properly... Exiting!")
        # tabix index
        sigint = subprocess.call(["tabix", "-p", "bed", output_file + ".gz"])
        if sigint != 0:
            raise Exception("tabix on the file " + output_file + ".gz did not run properly... Exiting!")

    def calculate_coverage(self, chr, pos):

        if self.__use_bam is True:
            return self.__coverage_at_position_pileup(chr, pos)
        else:
            return self.__coverage_at_position_tabix(chr, pos)

    def __coverage_at_position_pileup(self, chr, pos):

        if pos <= 0:
            pileup = self.__bam_file.pileup(chr, 0, pos + 1, truncate=True)
        else:
            pileup = self.__bam_file.pileup(chr, pos - 1, pos + 1, truncate=True)

        for pileupcolumn in pileup:
            if pileupcolumn.pos == pos:
                return pileupcolumn.n
        else:
            return 0

    def __coverage_at_position_tabix(self, chr, pos):

        cov = 0

        tbx = pysam.TabixFile(self.__tabix_file)
        for row in tbx.fetch(chr,pos,pos+1):
            cov = row[0]

        return cov

    def reads_with_indels_in_neighbourhood (self, chrom, pos):
        counts = {"insertions": 0, "deletions": 0}
        window_size = self.config['WINDOW_SIZE']
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