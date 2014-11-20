import sys
import csv
import json
from collections import defaultdict
import os.path

"""
    Author: Alejandro Sifrim
    Affiliation: Wellcome Trust Sanger Institute
    
    Aggregates over .score files and prints a tab delimited file with the frequency of the position being called and the samples it's being called in.

    Parameters
    ----------
    1) Score files (can use wildcards)
    
    Returns
    -------
    2) prints to stdout
    
"""


SCORE_THRESHOLD = 0.6
CHROMOSOMES = map(str,range(1,23))+["X","Y"]

db = defaultdict(list)
count = 0
for f in open(sys.argv[1],'r'):
	f = f.rstrip()
	if os.path.isfile(f):
		sample = os.path.basename(f).replace(".sr_reads.counts.scored","")
		count += 1
		for v in csv.DictReader(open(f,'r'), delimiter="\t"):
			if float(v["prob_Y"]) >= SCORE_THRESHOLD:
				db[v["chrom"]+"_"+v["position"]].append(sample)
				#print json.dumps(v,sort_keys=True,indent=4, separators=(',', ': '))

for k,v in db.iteritems():
	s = ""
	s += "\t".join(k.split("_"))
	s += "\t%.5f" % (float(len(v))/count)
	s += "\t%s" % len(v)
	s += "\t%s" % count
	# s += "\t"+",".join(v)
	print s