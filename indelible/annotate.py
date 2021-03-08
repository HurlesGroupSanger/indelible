import csv
import pybedtools as bedtools
import re
import sys
from indelible.indelible_lib import *
import re
from intervaltree import Interval, IntervalTree


def create_exon_intervaltree(exon_file):

    exons = {}

    header = ['chr','start','stop','transcript','score','strand','coding_start','coding_stop','score','num_exons','exon_lengths','exon_starts']
    exon_intervals = csv.DictReader(open(exon_file), delimiter="\t", fieldnames=header)

    for transcript in exon_intervals:

        lengths_array = transcript['exon_lengths'].split(",")
        lengths_array.pop()
        starts_array = transcript['exon_starts'].split(",")
        starts_array.pop()

        for i in range(0, len(lengths_array), 1):
            current_start = (int(transcript['start']) + int(starts_array[i])) + 1 ## +1 is to correct for 0-based bed format
            current_stop = current_start + (int(lengths_array[i]))

            current_chr = normalize_chr(transcript['chr'])

            if current_chr in exons:

                exons.get(current_chr).addi(current_start, current_stop, transcript['transcript'])

            else:
                current_tree = IntervalTree()
                current_tree.addi(current_start, current_stop, i)
                exons[current_chr] = current_tree

    return exons


def read_database(path):
    db = defaultdict(dict)

    header = ["chrom", "pos", "pct", "counts", "tot", "otherside", "mode", "svtype", "size", "aln_length",
     "otherside_found", "is_primary", "variant_coord"]

    for v in csv.DictReader(open(path, 'r'), fieldnames=header, delimiter="\t"):
        db[normalize_chr(v["chrom"]) + "_" + int(v["position"])] = {'maf': float(v["pct"]),
                                                                    'otherside': v["otherside"],
                                                                    'mode': v["mode"],
                                                                    'svtype': v["svtype"],
                                                                    'size': int(v["size"]),
                                                                    'aln_length': int(v["aln_length"]),
                                                                    'otherside_found': bool(v["otherside_found"]),
                                                                    'is_primary': bool(v["is_primary"]),
                                                                    'variant_coord': v["variant_coord"]}
    return db


def create_exac_constraint_hash(config):
    synonym_hash = create_gene_synonym_hash(config['hgnc_synonyms'])
    constraint_hash = {}
    exac_constraints_file = config['exac_constraint']

    for row in csv.DictReader(open(exac_constraints_file, 'r'), delimiter="\t"):
        if row['gene'] in synonym_hash:
            row['gene'] = synonym_hash[row['gene']]
        constraint_hash[row['gene']] = float(row["pLI"])

    return constraint_hash


def read_ddg2p(ddg2p_bed):
    ddg2p_db = {}

    for d in open(ddg2p_bed, 'r'):
        data = d.rstrip().split("\t")

        current_chr = normalize_chr(data[0])
        current_start = int(data[1])
        current_end = int(data[2])
        current_gene = data[3]

        if current_chr in ddg2p_db:
            ddg2p_db[current_chr].addi(current_start, current_end, current_gene)
        else:
            ddg2p_db[current_chr] = IntervalTree()
            ddg2p_db[current_chr].addi(current_start, current_end, current_gene)

    return ddg2p_db


def read_hgnc_genes(hgnc_bed):
    hgnc_db = {}
    for row in csv.DictReader(open(hgnc_bed, 'r'), delimiter="\t"):
        hgnc_db[row['HGNC']] = {"chrom": normalize_chr(row['CHROM']), "start": int(row["START"]), "end": int(row["END"])}
    return hgnc_db


def create_gene_synonym_hash(hgnc_synonyms):
    syn_hash = {}
    syn_file = hgnc_synonyms
    for row in csv.DictReader(open(syn_file, 'r'), delimiter="\t"):
        if row['Status'] != "Approved":
            continue
        synonyms = []
        synonyms.extend(row['Previous Symbols'].split(", "))
        # synonyms.extend(row['Synonyms'].split(", "))
        for synonym in synonyms:
            syn_hash[synonym] = row['Approved Symbol']
    return syn_hash


def attach_db(v, db):

        key = v["chrom"] + "_" + v["pos"]
        # While this is _slightly_ dangerous there should be a 0% chance that the key is not contained within this db.
        # (so long as the user didn't change the score cutoff during runtime...)
        db_entry = db[key]
        for key,value in db_entry.iteritems():
            v[key] = value


def find_hgnc_genes(chrom, start, end, hgnc_db):
    res = []
    for (gene, coords) in list(hgnc_db.items()):
        if chrom == normalize_chr(coords["chrom"]):
            if interval_overlap(start, end, coords["start"], coords["end"]):
                res.append(gene)

    if res != []:
        return res
    else:
        return None


def find_ddg2p_gene(chrom, start, end, ddg2p_db):
    res = []
    if chrom in ddg2p_db:
        overlaps = ddg2p_db.get(chrom).overlap(start, end)
        for d in overlaps:
            res.append(d.data)
        if res != []:
            return res
        else:
            return None
    else:
        return None


def hgnc_constrained_subset(genes, constraint_hash):
    constrained_hgnc = []
    for hg in genes:
        if hg in constraint_hash:
            if constraint_hash[hg] > 0.9:
                constrained_hgnc.append(hg)

    if constrained_hgnc != []:
        return constrained_hgnc
    else:
        return None


def find_protein_coding_ensembl_exon(chrom, start, end, ensembl_exons):

    res_exons = set()
    if chrom in ensembl_exons:
        for v in ensembl_exons[chrom].overlap(start - 10, end + 10):
            res_exons.add(v[2])

    if len(res_exons) > 0:
        return res_exons
    else:
        return None


def interval_overlap(start1, end1, start2, end2):
    s1 = min(int(start1), int(end1))
    e1 = max(int(start1), int(end1))
    s2 = min(int(start2), int(end2))
    e2 = max(int(start2), int(end2))

    overlap = max(0, min(e1, e2) - max(s1, s2))
    if overlap > 0:
        return True
    else:
        return False


def annotate(input_path, output_path, database, config):

    scored_file = csv.DictReader(open(input_path), delimiter="\t")

    # Prepare outputfile
    new_fieldnames = scored_file.fieldnames
    new_fieldnames.extend(("ddg2p", 'hgnc', 'hgnc_constrained', "exonic", "transcripts", "maf", "otherside","sv_type"))
    output_file = csv.DictWriter(open(output_path, 'w'), fieldnames=scored_file.fieldnames, delimiter="\t", lineterminator="\n")
    output_file.writeheader()

    # Prepare searchable hashes
    ensembl_exons = create_exon_intervaltree(config['ensembl_exons'])
    db = read_database(database)
    constraint_hash = create_exac_constraint_hash(config)
    ddg2p_db = read_ddg2p(config['ddg2p_bed'])
    hgnc_db = read_hgnc_genes(config['hgnc_file'])

    for v in scored_file:

        if v["prob_y"] >= config["SCORE_THRESHOLD"]:

            # Keys that I have removed:
            # hit["blast_hit"] = "no_hit"
            # hit["blast_dist"] = "NA"
            # hit["blast_identity"] = "NA"
            # hit["blast_strand"] = "NA"
            # hit["blast_hgnc"] = "NA"

            chrom = normalize_chr(v["chrom"])
            pos = int(v["position"])

            v = attach_db(v, db)

            # Decide if we can use additional info from the db for following steps:
            left_search = pos
            right_search = pos + 1
            if v["otherside"] != "NA":
                other_coord = v["otherside"].split("_")
                qchr = normalize_chr(other_coord[0])
                qpos = int(other_coord[1])

                if qchr == chrom:
                    if qpos < pos:
                        left_search = qpos
                        right_search = pos
                    else:
                        left_search = pos
                        right_search = qpos

            hgnc_genes = find_hgnc_genes(chrom, left_search, right_search, hgnc_db)
            ddg2p_genes = find_ddg2p_gene(chrom, left_search, right_search, ddg2p_db)

            if hgnc_genes != None:
                v["hgnc"] = ";".join(hgnc_genes)
                hgnc_constrained = hgnc_constrained_subset(hgnc_genes, constraint_hash)
                if hgnc_constrained != None:
                    v["hgnc_constrained"] = ";".join(hgnc_constrained)
                else:
                    v["hgnc_constrained"] = "NA"
            else:
                v["hgnc"] = "NA"

            if ddg2p_genes != None:
                v["ddg2p"] = ";".join(ddg2p_genes)
            else:
                v["ddg2p"] = "NA"

            exons = find_protein_coding_ensembl_exon(chrom, left_search, right_search, ensembl_exons)
            if exons == None:
                v["exonic"] = False
                v["transcripts"] = "NA"
            elif len(exons) > 10:
                v["exonic"] = True
                v["transcripts"] = "multiple_transcripts"
            else:
                v["exonic"] = True
                v["transcripts"] = ";".join(str(x) for x in exons)

            if pos in db[chrom]:
                v["maf"] = db[chrom][pos]
            else:
                v["maf"] = "NA"

            output_file.writerow(v)
