---
Title: "RNAseq Pipeline"
Author: "Michele Wiseman"
Date: "May 2022"
---

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

```shell
#in <project> directory, have subdirectories:

mkdir aligned
mkdir cleaned
mkdir rawReads
mkdir reference
mkdir trx
mkdir denovo
```
