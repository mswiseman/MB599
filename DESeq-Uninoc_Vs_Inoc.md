Continued from [Nugget vs. Symphony](DESeq-Nugget_vs_Sym.md). 

# Unioculated plants vs. Inoculated plants. 

![UninocVsInoc.png](images/UninocVsInoc.png | width = 300)

Since the time points seem to cluster closely, it might be more informative to look just at uninoculated/inoculated. So, I made a new dds object to examine just that. 

```r load packages, message=FALSE, warning=FALSE, include=FALSE

#Remember you can install packages from cran using install.packages("packagename")
#For Bioc packages, it's a little different
#BiocManager::install('apeglm')

library(PCAtools)
library(genefilter)
library(pheatmap)
library(RColorBrewer)
library(DESeq2)
library(tximportData)
library(tximport)
library(dplyr)
library(ggplot2)
library(hexbin)
library(readr)
library(adegenet)
library(readxl)
library(janitor)
library(gplots)
library(factoextra)
library(cowplot)
library(apeglm)
library(ggpubr)
```
# Set up DESeq changed metadata. 

```r
#Change design to inoc/no inoc
#Need to set up meta data for DESeq2. Changed from time to condition. 

condition = as.factor(c(rep("Uninoculated", 7),rep("Inoculated", 30))) 

genotype=as.factor(c(rep("Nugget", 4), rep("Symphony",3), rep("Nugget", 3), rep("Symphony", 4), rep("Nugget", 4), rep("Symphony", 4), rep("Nugget", 4), rep("Symphony", 4), rep("Nugget", 4), rep("Symphony", 3)))

colData = data.frame(col.names = c(" 01t00p1Nug"," 02t00p2Nug"," 03t00p3Nug"," 04t00p4Nug", " 05t00p1Sym", " 06t00p2Sym","07t00p3Sym","09t12p1Nug","10t12p2Nug","11t12p3Nug","12t12p1Sym","13t12p2Sym","14t12p3Sym","15t12p4Sym","16t24p1Nug","17t24p2Nug","18t24p3Nug","19t24p4Nug","20t24p1Sym","21t24p2Sym","22t24p3Sym","23t24p4Sym","24t48p1Nug","25t48p2Nug","26t48p3Nug","27t48p4Nug","28t48p1Sym","29t48p2Sym","30t48p3Sym","31t48p4Sym","32t72p1Nug","33t72p2Nug","34t72p3Nug","35t72p4Nug","36t72p1Sym","37t72p2Sym","39t72p4Sym"), genotype, condition) #condition

#create DESEQ2 object genotype and condition plus their interaction
dds <- DESeqDataSetFromMatrix(countData = countData2_NugSym_mtx, colData = colData, ~ genotype + condition + genotype:condition)

#delete rows with less than 10 reads
dds <-dds[ rowSums(counts(dds)) > 10, ]
```

# Load in annotations for later

```r
goTerms_biological <- read_tsv('/Users/michelewiseman/Desktop/goTerms_biologicalProcesses.tsv', col_names=FALSE)
goTerms_cellular <- read_tsv('/Users/michelewiseman/Desktop/goTerms_cellularComponents.tsv', col_names=FALSE)
goTerms_molecular <- read_tsv('/Users/michelewiseman/Desktop/goTerms_molecularFunction.tsv', col_names=FALSE)
kegg_pathway <- read_csv('/Users/michelewiseman/Desktop/kegg_pathway_map.csv', col_names=TRUE)

names(goTerms_biological) <- c('geneID', 'UNIProt_ID', 'GO ID', 'Description', 'Key terms')
names(goTerms_cellular) <- c('geneID', 'UNIProt_ID', 'GO ID', 'Description', 'Key terms')
names(goTerms_molecular) <- c('geneID', 'UNIProt_ID', 'GO ID', 'Description', 'Key terms')
```
