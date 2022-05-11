---
Title: "RNAseq Pipeline"
Authors: "Michele Wiseman, Jenny Hesser, Riely White, Sophie Baur"
Date: "May 2022"
---
# Acknowledgements
A starting point for much of the code included here was provided by Renee Erikssen and Andrew Black (as part of the [CQLS RNAseq course](https://cqls.oregonstate.edu/training/workshops)). Thank you both for your help.

# Programs used:

* [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [multiqc](https://multiqc.info/)
* [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
* [hisat2](https://ccb.jhu.edu/software/hisat2/manual.shtml)
* [samtools](http://www.htslib.org/doc/samtools-1.2.html)
* [stringtie](http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual)
* [salmon](https://salmon.readthedocs.io/en/latest/salmon.html)

# General set-up
The entirety of this pipeline will be done using the [CQLS Cybercomputing infastructure](https://cqls.oregonstate.edu/). For sake of security, I have excluded specific commands of logging in and/or submitting jobs to our cluster; however, hopefully you will find this code to be useful and largely reproducible. 

## QC

```shell
#in <project> directory, have subdirectories:

mkdir aligned
mkdir cleaned
mkdir rawReads
mkdir reference
mkdir trx
mkdir denovo

#in ~/<project>/rawReads, download fastq files
gunzip *.gz

#checking quality of FASTQ files using fastqc
#in ~/<project>/rawReads 
mkdir fastqc_results

#this code deposits quality check output results (.html and .zip files) into ~/<project>/raw/fastqc_results
fastqc *.fastq -o ./fastqc_results

#in ~/<project>/raw/fastqc_results/
mkdir ./multQC 
```

## Removal of adaptors

```shell
#code for script can be found in scripts folder
chmod +x cleaning.sh
bash cleaning.sh
```
## Download reference genome and annotation files

```shell
#note: some of the necessary files werent on hopbase as of 5/11/22, so I requested them from Dr. Henning and used CyberDuck to put them in our reference folder. 

wget "http://hopbase.cqls.oregonstate.edu/content/cascadeDovetail/assemblyData/dovetailCascade10ScaffoldsMasked.fasta.gz" -O dovetailCascade10ScaffoldsMasked.fasta.gz #masked version
wget "http://hopbase.cqls.oregonstate.edu/content/cascadeDovetail/assemblyData/dovetailCascade10ScaffoldsUnmasked.fasta.gz" -O dovetailCascade10ScaffoldsUnmasked.fasta.gz #unmasked
wget "http://hopbase.cqls.oregonstate.edu/content/cascadeDovetail/geneData/combinedGeneModels/combinedGeneModels.tenScaffolds.pep.fasta.gz" -O combinedGeneModels.tenScaffolds.pep.fasta.gz
wget "http://hopbase.cqls.oregonstate.edu/content/cascadeDovetail/geneData/combinedGeneModels/combinedGeneModels.txt" -O combinedGeneModels.txt
wget "http://hopbase.cqls.oregonstate.edu/content/cascadeDovetail/repeatData/compiledRepeats.gff.gz" -O compiledRepeats.gff.gz
wget "http://hopbase.cqls.oregonstate.edu/content/cascadeDovetail/geneData/maker/makerGenes.gff3.gz" -O makerGenes.gff3.gz
wget "http://hopbase.cqls.oregonstate.edu/content/cascadeDovetail/geneData/transdecoder/completeTransdecoderGeneModels/transdecoder.tenScaffolds.cds.fasta.gz" -O transdecoder.tenScaffolds.cds.fasta.gz
wget "http://hopbase.cqls.oregonstate.edu/content/cascadeDovetail/geneData/transdecoder/completeTransdecoderGeneModels/transdecoder.tenScaffolds.pep.fasta.gz" -O transdecoder.tenScaffolds.pep.fasta.gz
wget "http://hopbase.cqls.oregonstate.edu/content/cascadeDovetail/geneData/transdecoder/transdecoderOutput/transcripts.fasta.transdecoder.genomeCentric.gff3.gz" -O transcripts.fasta.transdecoder.genomeCentric.gff3.gz
wget "http://hopbase.cqls.oregonstate.edu/content/cascadeDovetail/geneData/maker/all.maker.proteins.fasta.gz" -O all.maker.proteins.fasta.gz
wget "http://hopbase.cqls.oregonstate.edu/content/cascadeDovetail/geneData/maker/all.maker.transcripts.fasta.gz" -O all.maker.transcripts.fasta.gz
wget "http://hopbase.cqls.oregonstate.edu/content/cascadeDovetail/geneData/maker/all.maker.cds.fasta.gz" -O all.maker.cds.fasta.gz

#decompress files
gunzip *.gz

#determine length of scaffolds
cat dovetailCascade10ScaffoldsMasked.fasta | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > polishedScaffoldLengths.txt

#For some reason gffread was picky about the extra features, so I removed them
grep -v "UTR" combinedGeneModels.tenScaffolds.repeatFiltered.gff > combinedGeneModels.tenScaffolds.repeatFiltered.gff
grep -v "gene" combinedGeneModels.tenScaffolds.repeatFiltered.gff > combinedGeneModels.tenScaffolds.repeatFiltered.gff

#convert to .gtf for hisat2
gffread combinedGeneModels.tenScaffolds.repeatFiltered.gff -T -o combinedGeneModels.tenScaffolds.repeatFiltered.gtf
```
## Mapping

```shell
#code for script can be found in scripts folder
chmod +x mapping.sh
bash mapping.sh
```

## Visualize results
```shell 
#to visualize results::
grep 'overall alignment rate' results.txt | cut -c -5 > rate.tmp
ls -1 *fastq | cut -c -17 | uniq > names.tmp
paste names.tmp rate.tmp > mapping_genome.txt
cat mapping_genome.txt

#can export and visualize in R
```

## to merge transcripts
```shell
#this takes a list of gtf files [gtf_list] and constructs a gtf file of merged reads into transcripts
stringtie --merge -G ../../reference/combinedGeneModels.tenScaffolds.repeatFiltered.gtf -o stringtie_merged.gtf gtf_list.txt
```
## Prep for DESeq
```
#make a file with the sample names
ls -1 *.gtf
ls -1 *.gtf > samplelist.txt
cat samplelist.txt

#prepend sampleX to the files
awk '{print "sample"NR, "\t", $0}' < samplelist.txt > Samples.txt
rm samplelist.txt
nano Samples.txt
#make sure that sample1 reads sample01, sample2 reads sample02
```
## Running DE in String tie
```shell
#make sure prepDE.py is installed on the cluster
#if not: https://ccb.jhu.edu/software/stringtie/dl/prepDE.py

chmod +x prepDE.py
```
## Preview DE output
```shell
less -S gene_count_matrix.csv
less -S transcript_count_matrix.csv
```
## Running salmon
```shell
#script can be found in scripts folder
chmod +x salmon.sh
bash salmon.sh
```

## Summarize mapping rate
```shell
#in trx directory

#first pull out the sample names
grep "Mapping rate" * -r | cut -c -21 | uniq > names
cat names
#remove /log at end if it's there
cut -b -17 names > names1

#next pull out the mapping rate
grep "Mapping rate" * -r | cut -c 100- | uniq > rate
cat rate

#create a column of factors to assign group
echo Symphony Nugget Symphony Nugget Symphony Nugget Symphony Nugget Nugget Nugget Symphony Symphony Symphony Symphony Nugget Nugget Nugget Nugget Symphony  Symphony Symphony Symphony Nugget Nugget Nugget Nugget Nugget Symphony Symphony Symphony Symphony Symphony Symphony Symphony Symphony Nugget Nugget Nugget Nugget| tr " " "\n" > group

#paste columns into new text document
paste names1 rate group > mapping_trx.txt
cat mapping_trx.txt

#tidy up space
rm -f names rate group
```






