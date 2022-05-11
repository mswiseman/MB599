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
```
