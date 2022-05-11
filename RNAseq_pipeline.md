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

