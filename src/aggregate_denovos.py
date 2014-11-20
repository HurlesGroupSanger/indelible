import sys
import os
import pandas as pd

datas = []
for f in open(sys.argv[1]):
	sample_id = f
	data = pd.read_table(f.rstrip())
	data["sample_id"] = os.path.basename(sample_id).strip(".sr_reads.counts.scored.annotated.denovo\n")
	datas.append(data)
all_data = pd.concat(datas)
print all_data
all_data["position"] = all_data["position"].astype(int)
all_data.to_csv("/nfs/users/nfs_a/as33/Projects/Indelible/results/ddd4k_denovos.tsv",sep="\t",index=False)
