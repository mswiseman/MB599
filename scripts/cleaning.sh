#!/bin/bash

for i in $(ls -1 *fastq | cut -c -44 | uniq )
do
echo "cleaning sample $i"
cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o ${i}_cleaned.fastq -p ${i}.fastq >> ~/project/cleaned/cutadapt.log
done
