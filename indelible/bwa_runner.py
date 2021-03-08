"""
    Author: Eugene Gardner
    Affiliation: Wellcome Sanger Institute

    Object which handles various coverage calculations for aggregate_positions.py and trio_caller.py

    ----------
    Parameters
    1.)


"""

import pysam
import subprocess

class BWARunner:

    def __init__(self, final_frame, output_file, fasta, bwa_loc, bwa_threads):

        self.__fasta = fasta
        self.__bwa_loc = bwa_loc
        self.__bwa_threads = bwa_threads
        self.final_frame = final_frame

        # dict of the structure:
        # readname:
        #   otherside:
        #   PASS/FAIL state (enums possible with Python...?)
        #   SV Type
        #   SV size
        #   Length of split read alignment
        self.__final_decision = {}

        # First – create fastq files with read information:
        split_fq = output_file + ".1.fastq"
        ref_fq = output_file + ".2.fastq"
        split_pairs = open(split_fq, "w")
        ref_pairs = open(ref_fq, "w")

        for index, row in final_frame.iterrows():

            ref_seq = self.__get_ref_string(row["dir"], row["chrom"], row["pos"])
            self.__write_fastq_entry(split_pairs, row["chrom"], row["pos"], row["longest"])
            self.__write_fastq_entry(ref_pairs, row["chrom"], row["pos"], ref_seq)

        split_pairs.close()
        ref_pairs.close()

        # Now run bwa mem:
        self.__bwa_engine(split_fq, ref_fq, output_file + ".sam")

        # And process alignments...
        self.__process_alignment(output_file + ".sam")

    def get_decisions(self):

        return self.__final_decision

    def __write_fastq_entry(self, writer, chrom, position, sequence):

        qual = 'Z' * len(sequence)
        name = "%s_%s" % (chrom, position)

        writer.write('@' + name + "\n")
        writer.write(sequence + "\n")
        writer.write('+' + "\n")
        writer.write(qual + "\n")

    def __get_ref_string(self, dir, chrom, position):

        chrom = str(chrom)
        if dir == "left":

            left = position
            right = position + 100

        elif dir == "right":

            left = position - 100
            right = position

        else:

            left = position - 50
            right = position + 50

        ref_seq = self.__fasta.fetch(reference=chrom, start=left, end=right)

        return(ref_seq)

    def __bwa_engine(self, split_fq, ref_fq, outsam):

        cmd = self.__bwa_loc + " mem -t " + str(self.__bwa_threads) + " -T 10 -k 10 -o " + outsam + " " + self.__fasta.filename + " " + split_fq + " " + ref_fq
        p = subprocess.Popen(cmd, shell=True)
        p.wait()
        print("bwa ran with code:", p.returncode)

    def __process_alignment(self, insam):

        samfile = pysam.AlignmentFile(insam, "r")

        current_id = None
        read_one = None
        read_two = None

        for read in samfile:

            name = read.qname

            # We have hit the next set of IDs, pause, check, and add to dict if necessary before continuing
            if name != current_id:

                self.__process_current_group(current_id, read_one, read_two)

                # reset loop parameters:
                current_id = None
                read_one = None
                read_two = None

            if read.is_supplementary or read.is_duplicate or read.is_qcfail or read.is_secondary or read.is_unmapped:
                # Toss in the bin (have this nothing statement because it's more understandable than negating everything above...)
                stuff = 1
            else:
                current_id = name
                if read.is_read1:
                    read_one = read
                else:
                    read_two = read

        # Process the last read if possible:
        self.__process_current_group(current_id, read_one, read_two)

    def __process_current_group(self, current_id, read_one, read_two):

        ## Only fill the dictionary if we get valid information, will default to None for empties
        if read_one is not None and read_two is not None:

            curr_data = self.final_frame.loc[current_id].to_dict()
            if curr_data["dir"] == "uncer":

                self.__fill_dict(current_id, "NA", "FAIL_MULTISPLIT", "UNK", "NA", read_one.query_alignment_length)

            elif read_one.mapq == 0:

                self.__fill_dict(current_id, "NA", "FAIL_LOWMAPQ", "UNK", "NA", read_one.query_alignment_length)

            elif read_one.reference_name != read_two.reference_name:

                self.__fill_dict(current_id, read_one.reference_name + "_" + str(read_one.reference_start), "REALN_CHR",
                                 "TRANS_SEGDUP", "NA", read_one.query_alignment_length)

            else:

                size = None
                sv_type = None
                bp = None
                aln_len = read_one.query_alignment_length

                if curr_data["dir"] == "left":

                    size = abs((read_one.reference_end) - curr_data["pos"])
                    bp = "%s_%s" % (read_one.reference_name,read_one.reference_end)

                    if curr_data["pos"] < read_one.reference_end:
                        sv_type = "DUP"
                    else:
                        sv_type = "DEL"

                else:

                    size = abs((read_one.reference_start) - curr_data["pos"])
                    bp = "%s_%s" % (read_one.reference_name, read_one.reference_start)

                    if curr_data["pos"] > read_one.reference_start:
                        sv_type = "DUP"
                    else:
                        sv_type = "DEL"

                if read_one.is_proper_pair and read_two.is_proper_pair:
                    self.__fill_dict(current_id, bp, "REALN", sv_type, size, aln_len)
                else:
                    self.__fill_dict(current_id, bp, "REALN_XL", sv_type, size, aln_len)

    def __fill_dict(self, breakpoint_id, otherside, mode, svtype, size, aln_length):

        self.__final_decision[breakpoint_id] = {'otherside': otherside,
                                                'mode': mode,
                                                'svtype': svtype,
                                                'size': size,
                                                'aln_length': aln_length}


