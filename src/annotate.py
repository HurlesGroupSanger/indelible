import csv
import pybedtools as bedtools
import re
import sys
from indelible_lib import *
ddg2p_bed = "/nfs/users/nfs_a/as33/Projects/Indelible/data/DDG2P_18112014.bed"
ensembl_exons = bedtools.IntervalFile("/nfs/users/nfs_a/as33/Projects/Indelible/data/protein_coding_ensembl_exons_gencode17.gff")
db = read_database("/nfs/users/nfs_a/as33/Projects/Indelible/data/ddd4k_database.tsv")
scored_file = sys.argv[1]
CHROMOSOMES = [str(x) for  x in range(1,23)]+["X","Y"]

def find_protein_coding_ensembl_exon(chrom,pos):
	query = bedtools.Interval(chrom,int(pos-1),int(pos+1))
	res_exons = []
	for v in ensembl_exons.all_hits(query):
		res_exons.append(v.attrs)
	if res_exons != []:
		return [map(lambda x:x["transcript_id"].strip("\""),res_exons),map(lambda x:x["exon_number"].strip("\""),res_exons)]
		# return res_exons
	else:
		return None

def read_ddg2p(ddg2p_bed):
	ddg2p_db = {}
	for c in CHROMOSOMES:
		ddg2p_db[c] = []
	for d in open(ddg2p_bed,'r'):
		data = d.rstrip().split("\t")
		if data[0] in CHROMOSOMES:
			ddg2p_db[data[0]].append({"start":int(data[1]),"end":int(data[2]),"gene":data[3]})
	return ddg2p_db

def find_ddg2p_gene(ddg2p_db,chrom,pos):
	res = []
	for d in ddg2p_db[chrom]:
		if d["start"] <= pos and d["end"] >= pos:
			res.append(d["gene"])
	if res != []:
		return res
	else:
		return None

# def find_ddg2p_gene(chrom,pos):
# 	query = bedtools.Interval(chrom,int(pos-1),int(pos+1))
# 	res_genes = []
# 	for v in ddg2p.all_hits(query):
# 		print v
# 		res_genes.append(v.name)
# 	if res_genes != []:
# 		return res_genes
# 	else:
# 		return None

scored_file = csv.DictReader(open(scored_file),delimiter="\t")

new_fieldnames = scored_file.fieldnames
new_fieldnames.extend(("ddg2p","exonic","transcripts","exon_numbers","maf"))
output_file = csv.DictWriter(open(sys.argv[1]+".annotated",'w'),fieldnames=new_fieldnames,delimiter="\t")
output_file.writeheader()
ddg2p_db = read_ddg2p(ddg2p_bed)

for v in scored_file:
	chrom = v["chrom"]
	pos = int(v["position"])
	ddg2p_genes = find_ddg2p_gene(ddg2p_db,chrom,pos)
	if ddg2p_genes != None:
		v["ddg2p"] = ";".join(ddg2p_genes)
	else:
		v["ddg2p"] = ""

	exons = find_protein_coding_ensembl_exon(chrom,pos)
	if exons == None:
		v["exonic"] = False
		v["transcripts"] = ""
		v["exon_numbers"] = ""
	else:
		v["exonic"] = True
		v["transcripts"] = ";".join(exons[0])
		v["exon_numbers"] = ";".join(exons[1])
	if pos in db[chrom]:
		v["maf"] = db[chrom][pos]
	else:
		v["maf"] = ""
	output_file.writerow(v)


