#!/bin/bash

cd $1
mkdir logs
mv *err *out logs
rm *.sr_reads
rm *.counts