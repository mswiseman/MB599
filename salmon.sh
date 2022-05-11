#!/bin/bash

for i in `ls -1 *.fastq | cut -c -17 |uniq`
do
        echo "working with sample $i"
        salmon quant -l A -i ../../reference/cascade_index/ -r ${i}_cleaned.fastq -p 4 -o ../trx/${i} --numAuxModelSamples 100000 --numPreAuxModelSamples 100000 --writeUnmappedNames  ../trx/${i}_unmappedReads.txt
done
