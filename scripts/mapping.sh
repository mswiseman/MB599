#!/bin/bash

for i in `ls -1 *.fastq| cut -c -17 |uniq`
do
  	echo "Input: sample $i"
        echo "Output: sample $i.gtf"
echo " Running hisat"
touch results.txt
hisat2 --threads 4 -x ../../reference/cascade_index -U "$i"_cleaned.fastq -S ../aligned/$i.sam >> results.txt

echo " Running samtools"
samtools sort -@ 4 -o ../aligned/$i.sam.bam ../aligned/$i.sam

echo " Running stringtie"
stringtie -p 4 -e -G ../../reference/combinedGeneModels.tenScaffolds.repeatFiltered.gtf -o ../aligned/$i.gtf ../aligned/$i.sam.bam
done

SGE_Batch -c "bash mapping.sh" >>results.txt
