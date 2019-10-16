import csv
import pybedtools as bedtools
import re
import sys
from indelible.indelible_lib import *
import re
from intervaltree import Interval, IntervalTree

CHROMOSOMES = [str(x) for x in range(1, 23)] + ["X", "Y"]


def create_exon_intervaltree(exon_file):

    exons = {}

    header = ['chr','start','stop','transcript','score','strand','coding_start','coding_stop','score','num_exons','exon_lengths','exon_starts']
    exon_intervals = csv.DictReader(open(exon_file), delimiter="\t", fieldnames=header)

    for transcript in exon_intervals:

        lengths_array = transcript['exon_lengths'].split(",")
        lengths_array.pop()
        starts_array = transcript['exon_starts'].split(",")
        starts_array.pop()

        exon_starts = []
        exon_stops = []

        for i in range(0, len(lengths_array), 1):
            current_start = (int(transcript['start']) + int(starts_array[i])) + 1 ## +1 is to correct for 0-based bed format
            current_stop = current_start + (int(lengths_array[i]))

            if transcript['chr'] in exons:

                exons.get(transcript['chr']).addi(current_start, current_stop, transcript['transcript'])

            else:
                current_tree = IntervalTree()
                current_tree.addi(current_start, current_stop, i)
                exons[transcript['chr']] = current_tree

    return exons


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


def create_exac_constraint_hash(config):
    synonym_hash = create_gene_synonym_hash(config['hgnc_synonyms'])
    constraint_hash = {}
    exac_constraints_file = config['exac_constraint']

    for row in csv.DictReader(open(exac_constraints_file, 'r'), delimiter="\t"):
        if row['gene'] in synonym_hash:
            row['gene'] = synonym_hash[row['gene']]
        constraint_hash[row['gene']] = float(row["pLI"])

    return constraint_hash


def find_protein_coding_ensembl_exon(chrom, pos, blast_hit, ensembl_exons):

    # This is simply because of a weird scenario where SRs can be found at the beginning of a contig/chromosome...
    if pos == 0:
        pos = 1
    else:
        pos = pos

    p = re.compile("(\S+):(\d+)\-(\d+)")
    m = p.match(blast_hit)

    if m:
        if m.group(1) is chrom:
            if int(m.group(2)) < pos:
                start_coord = int(m.group(2))
                end_coord = pos + 10
            else:
                start_coord = pos - 10
                end_coord = int(m.group(2))
        else:
            start_coord = pos - 10
            end_coord = pos + 10
    else:
        start_coord = pos - 10
        end_coord = pos + 10

    res_exons = set()
    if chrom in ensembl_exons:
        for v in ensembl_exons[chrom].overlap(start_coord, end_coord):
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


def read_ddg2p(ddg2p_bed):
    ddg2p_db = {}
    for c in CHROMOSOMES:
        ddg2p_db[c] = []
    for d in open(ddg2p_bed, 'r'):
        data = d.rstrip().split("\t")
        if data[0] in CHROMOSOMES:
            ddg2p_db[data[0]].append({"start": int(data[1]), "end": int(data[2]), "gene": data[3]})
    return ddg2p_db


def find_ddg2p_gene(chrom, start, end, ddg2p_db):
    res = []
    if chrom in ddg2p_db:
        for d in ddg2p_db[chrom]:
            if interval_overlap(start, end, d["start"], d["end"]):
                res.append(d["gene"])
        if res != []:
            return res
        else:
            return None
    else:
        return None


def read_hgnc_genes(hgnc_bed):
    hgnc_db = {}
    for row in csv.DictReader(open(hgnc_bed, 'r'), delimiter="\t"):
        hgnc_db[row['HGNC']] = {"chrom": row['CHROM'], "start": int(row["START"]), "end": int(row["END"])}
    return hgnc_db


def find_hgnc_genes(chrom, start, end, hgnc_db):
    res = []
    for (gene, coords) in list(hgnc_db.items()):
        if chrom == coords["chrom"]:
            if interval_overlap(start, end, coords["start"], coords["end"]):
                res.append(gene)
    if res != []:
        return res
    else:
        return None


def hgnc_constrained_subset(genes, constraint_hash):
    constrained_hgnc = []
    for hg in genes:
        if hg in constraint_hash:
            if constraint_hash[hg] > 0.9:
                constrained_hgnc.append(hg)
    return constrained_hgnc


def create_blast_hash(scored_file):
    blast_nonrepeats_path = scored_file + ".fasta.hits_nonrepeats"
    blast_repeats_path = scored_file + ".fasta.hits_repeats"

    h = {}

    for row in csv.DictReader(open(blast_nonrepeats_path, 'r'), delimiter="\t"):
        key = row['chrom'] + "_" + row["pos"] + "_" + row["query_length"]
        if key not in h:
            h[key] = {}
            h[key]["nonrepeats"] = []
            h[key]["repeats"] = []
        h[key]["nonrepeats"].append(row)

    for row in csv.DictReader(open(blast_repeats_path, 'r'), delimiter="\t"):
        key = row['chrom'] + "_" + row["pos"] + "_" + row["query_length"]
        if key not in h:
            h[key] = {}
            h[key]["nonrepeats"] = []
            h[key]["repeats"] = []

        h[key]["repeats"].append(row)
    return h


# THIS IS POSSIBLY THE UGLIEST CODE I'VE EVER WRITTEN
def annotate_blast(hit, blast_hash, ddg2p_db, constraint_hash, hgnc_db):
    key = hit["chrom"] + "_" + hit["position"] + "_" + str(len(hit["seq_longest"]))
    if key in blast_hash:
        blast_hit = blast_hash[key]
        if blast_hit["repeats"] == []:
            if len(blast_hit["nonrepeats"]) == 1:
                blast_hit = blast_hit["nonrepeats"][0]
                hit["blast_hit"] = "%s:%s-%s" % (
                blast_hit['target_chrom'], blast_hit['target_start'], blast_hit['target_end'])
                hit["blast_strand"] = blast_hit["target_strand"]

                if blast_hit["target_chrom"] == hit["chrom"]:
                    hit["blast_dist"] = min(
                        abs(int(hit["position"]) - int(blast_hit["target_start"])),
                        abs(int(hit["position"]) - int(blast_hit["target_end"]))
                    )
                    # Choose smallest interval
                    if abs(int(hit["position"]) - int(blast_hit["target_start"])) > abs(
                            int(hit["position"]) - int(blast_hit["target_end"])):
                        hit["blast_hgnc"] = find_hgnc_genes(hit["chrom"], hit["position"], blast_hit["target_end"],
                                                            hgnc_db)
                    else:
                        hit["blast_hgnc"] = find_hgnc_genes(hit["chrom"], hit["position"], blast_hit["target_start"],
                                                            hgnc_db)
                    if hit["blast_hgnc"] != None:
                        hit["blast_hgnc"] = ";".join(hit["blast_hgnc"])
                else:
                    hit["blast_dist"] = "other_chrom"
                    hit["blast_hgnc"] = "NA"
                hit["blast_identity"] = blast_hit["target_identity"]
            else:
                hit["blast_hit"] = "multi_hit"
                hit["blast_dist"] = "NA"
                hit["blast_identity"] = "NA"
                hit["blast_strand"] = "NA"
                hit["blast_hgnc"] = "NA"
        else:
            hit["blast_hit"] = "repeats_hit"
            hit["blast_dist"] = "NA"
            hit["blast_identity"] = "NA"
            hit["blast_strand"] = "NA"
            hit["blast_hgnc"] = "NA"
    else:
        hit["blast_hit"] = "no_hit"
        hit["blast_dist"] = "NA"
        hit["blast_identity"] = "NA"
        hit["blast_strand"] = "NA"
        hit["blast_hgnc"] = "NA"
    return hit


def annotate(input_path, output_path, database, config):

    scored_file = csv.DictReader(open(input_path), delimiter="\t")

    # Prepare outputfile
    # new_fieldnames = scored_file.fieldnames
    # new_fieldnames.extend(("ddg2p", 'hgnc', 'hgnc_constrained', "exonic", "transcripts", "maf",
    #                        "blast_hit", "blast_strand", "blast_identity", "blast_dist", "blast_hgnc"))
    # output_file = csv.DictWriter(open(output_path, 'w'), fieldnames=new_fieldnames, delimiter="\t", lineterminator="\n")
    output_file = csv.DictWriter(open(output_path, 'w'), fieldnames=scored_file.fieldnames, delimiter="\t", lineterminator="\n")
    output_file.writeheader()
    # Prepare searchable hashes
    # bhash = create_blast_hash(input_path)
    ensembl_exons = create_exon_intervaltree(config['ensembl_exons'])
    db = read_database(database)
    constraint_hash = create_exac_constraint_hash(config)
    ddg2p_db = read_ddg2p(config['ddg2p_bed'])
    hgnc_db = read_hgnc_genes(config['hgnc_file'])

    for v in scored_file:

        # v = annotate_blast(v, bhash, ddg2p_db, constraint_hash, hgnc_db)
        chrom = v["chrom"]
        pos = int(v["position"])

        hgnc_genes = find_hgnc_genes(chrom, pos, pos + 1, hgnc_db)

        ddg2p_genes = find_ddg2p_gene(chrom, pos, pos + 1, ddg2p_db)

        if hgnc_genes != None:
            v["hgnc"] = ";".join(hgnc_genes)
            hgnc_constrained = hgnc_constrained_subset(hgnc_genes, constraint_hash)
            v["hgnc_constrained"] = ";".join(hgnc_constrained)
        else:
            v["hgnc"] = "NA"

        if ddg2p_genes != None:
            v["ddg2p"] = ";".join(ddg2p_genes)
        else:
            v["ddg2p"] = "NA"

        exons = find_protein_coding_ensembl_exon(chrom, pos, v["blast_hit"], ensembl_exons)
        if exons == None:
            v["exonic"] = False
            v["transcripts"] = "NA"
        else:
            v["exonic"] = True
            v["transcripts"] = ";".join(str(x) for x in exons)

        if pos in db[chrom]:
            v["maf"] = db[chrom][pos]
        else:
            v["maf"] = "NA"

        output_file.writerow(v)
