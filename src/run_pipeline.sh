#!/bin/bash

TRIO=$(sed "$LSB_JOBINDEX q;d" /nfs/users/nfs_a/as33/Projects/Indelible/data/DDD4k_trio_bam_paths.txt)
CHILD=$(echo $TRIO | cut -f 1 -d " ")
MUM=$(echo $TRIO | cut -f 2 -d " ")
DAD=$(echo $TRIO | cut -f 3 -d " ")

python /nfs/users/nfs_a/as33/Projects/Indelible/src/pipeline.py $CHILD $MUM $DAD