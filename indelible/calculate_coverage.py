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

import pybedtools as bedtools
import subprocess
from indelible.indelible_lib import *
import pysam

class CoverageCalculator:

    def __init__ (self, chr_dict, bam_file, config):

        self.chr_dict = chr_dict
        self.bam_file = bam_file
        self.__tabix_file = bam_file + "bg.gz"
        self.config = config
        self.use_bam

        self.__decide_coverage_method()


    def calculate_coverage(self, chr, pos):

        if self.use_bam is True:
            return self.__coverage_at_position_pileup(chr, pos)
        else:
            return self.__coverage_at_position_tabix(chr, pos)

    def __decide_coverage_method(self):

        count = 0

        for chrom in self.chr_dict:
            for position in self.chr_dict[chrom]:
                pos = int(position)
                if len(self.chr_dict[chrom][position]) >= self.config["MINIMUM_SR_COVERAGE"]:
                    count += 1

        if count <= 10:
            self.use_bam = True
        else:
            print "Using tabix method..."
            self.use_bam = False
            output_file = self.bam_file + ".bg"
            self.__calculate_coverage_bam(output_file)


    def __calculate_coverage_bam (self, output_file):

        bt = bedtools.BedTool(self.bam_file)
        bg_file = bt.genome_coverage(bg=True)
        bg_file.saveas(output_file)

        # bgzip
        sigint = subprocess.call(["bgzip",output_file])
        if sigint != 0:
            raise Exception("bgzip on the file " + output_file + " did not run properly... Exiting!")
        # tabix index
        sigint = subprocess.call(["tabix","-p","bed",output_file + ".gz"])
        if sigint != 0:
            raise Exception("tabix on the file " + output_file + ".gz did not run properly... Exiting!")


    def __coverage_at_position_pileup(self, chr, pos):

        for pileupcolumn in CoverageCalculator.bam_file.pileup(chr,pos-1,pos+1,max_depth=100,truncate=True):
            if pileupcolumn.pos == pos:
                return pileupcolumn.n
        else:
            return 0


    def __coverage_at_position_tabix(tabix_file, chr, pos):

        cov = 0

        tbx = pysam.TabixFile(tabix_file)
        for row in tbx.fetch(chr,pos,pos):
            cov = row[0]

        return cov