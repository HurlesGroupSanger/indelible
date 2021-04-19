"""
Author: Alejandro Sifrim & Eugene Gardner
Affiliation: Wellcome Trust Sanger Institute

Library to blast and interpret results of soft-clipped segments:
- Output segments to fasta file
- Blast them
- Read result

Parameters
----------

Returns
-------

"""
import csv
import time
from typing import Dict, Any

from Bio.Blast.Applications import NcbiblastnCommandline
import io

today = time.strftime('%Y%m%d')

"""
PARAMETERS
"""

MINIMUM_LENGTH = 20

BLAST_FIELDS = [
    'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
    'sstart', 'send', 'sstrand', 'evalue', 'bitscore']


class BlastRepeats:

    def __init__(self, output_path, final_frame, REPEATdb):

        self.__repeat_blast_info = {}
        fasta_path = self.__generate_fasta(output_path, final_frame)
        self.__blast_fasta(fasta_path, REPEATdb)


    def get_repeat_blast_info(self):

        return self.__repeat_blast_info

    def __generate_fasta(self, output_path, final_frame):
        output_fasta = open(output_path + ".fasta", 'w')
        for index, row in final_frame.iterrows():
            if len(row['longest']) >= MINIMUM_LENGTH:
                output_fasta.write(">%s:%s:%s\n" % (row['chrom'], row['pos'], len(row['longest'])))
                output_fasta.write("%s\n" % row['longest'])
        output_fasta.close()
        return output_path + ".fasta"

    def __run_blast(self, fasta_file, db):
        blastn_clin = NcbiblastnCommandline(query=fasta_file, db=db, word_size=15,
                                            max_target_seqs=100, penalty=-3, evalue=0.001, reward=1,
                                            outfmt='6 ' + " ".join(BLAST_FIELDS))

        stdout, stderr = blastn_clin()
        input = io.StringIO(stdout)

        hits = {}

        for r in input:
            data = r.rstrip().split("\t")
            if data[0] in hits:
                hits[data[0]].append(data)
            else:
                hits[data[0]] = [data]

        return hits

    def __blast_fasta(self, fasta_file, db):

        hits = self.__run_blast(fasta_file, db)

        for seq in hits:
            if len(hits[seq]) < 10:
                for result in hits[seq]:

                    id_split = result[0].split(":")
                    current_id = id_split[0] + ":" + id_split[1]
                    res = {}
                    res["target_chrom"] = result[1]
                    res["evalue"] = float(result[11])
                    res["query_length"] = res["query_length"] = id_split[2]

                    # This will effectively take the first hit if all of the values are the same,
                    # which should be roughly the same repeat type anyway...
                    # BLAST sorts by e.value, so this code is actually kind of just wasting time...
                    # Need to warn that family information for MEIs in unlikely to be accurate
                    if current_id in self.__repeat_blast_info:
                        past_evalue = self.__repeat_blast_info[current_id]["evalue"]
                        if res["evalue"] < past_evalue:
                            self.__repeat_blast_info[current_id] = res
                            print(res["evalue"], past_evalue)
                    else:
                        self.__repeat_blast_info[current_id] = res
