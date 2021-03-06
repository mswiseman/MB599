---
title: "Phyloseq work"
authors: "Riely White and Sophie Baur"
date: June 2022
---


```r
#create samples_by_metadata for phyloseq
library(readxl)
countData_NugSym <- read_excel("C:/Users/sophi/Downloads/NEWVERSION_dovetailAssemblyFullNonRepeatAssociatedGeneList.xlsx", 
                               col_types = c("text", "skip", "skip", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))

#store GeneID info because the matrix doesnt like it. 
GeneID <- countData_NugSym$GeneID

#convert counts to matrix and subtract geneID column. 
countData2_NugSym <- as.matrix(countData_NugSym[ , -1])

#check work if you want
#view(countData2_NugSym)

#add gene IDs back. 
row.names(countData2_NugSym) <- GeneID   

#will need this later
Sample <- c(" 01t00p1Nug"," 02t00p2Nug"," 03t00p3Nug"," 04t00p4Nug", " 05t00p1Sym", " 06t00p2Sym","07t00p3Sym","08t00p4Sym","09t12p1Nug","10t12p2Nug","11t12p3Nug","12t12p1Sym","13t12p2Sym","14t12p3Sym","15t12p4Sym","16t24p1Nug","17t24p2Nug","18t24p3Nug","19t24p4Nug","20t24p1Sym","21t24p2Sym","22t24p3Sym","23t24p4Sym","24t48p1Nug","25t48p2Nug","26t48p3Nug","27t48p4Nug","28t48p1Sym","29t48p2Sym","30t48p3Sym","31t48p4Sym","32t72p1Nug","33t72p2Nug","34t72p3Nug","35t72p4Nug","36t72p1Sym","37t72p2Sym","38t72p3Sym","39t72p4Sym") 

#can check everything by uncommenting and running commands below. 
#View(countData2_NugSym)
#str(countData2_NugSym)

#need to set up meta data for DESeq2
Condition = as.factor(c(rep("Uninoculated",8),rep("Inoculated",31))) 

Cultivar=as.factor(c(rep("Nugget", 4), rep("Symphony",4), rep("Nugget", 3), rep("Symphony", 4), rep("Nugget", 4), rep("Symphony", 4), rep("Nugget", 4), rep("Symphony", 4), rep("Nugget", 4), rep("Symphony", 4)))

Time=as.factor(c(rep("0 hr", 8), rep("12 hr", 7), rep("24 hr", 8), rep("48 hr", 8), rep("72 hr", 8)))

#colData = data.frame(col.names = c(" 01t00p1Nug"," 02t00p2Nug"," 03t00p3Nug"," 04t00p4Nug", " 05t00p1Sym", " 06t00p2Sym","07t00p3Sym","08t00p4Sym","09t12p1Nug","10t12p2Nug","11t12p3Nug","12t12p1Sym","13t12p2Sym","14t12p3Sym","15t12p4Sym","16t24p1Nug","17t24p2Nug","18t24p3Nug","19t24p4Nug","20t24p1Sym","21t24p2Sym","22t24p3Sym","23t24p4Sym","24t48p1Nug","25t48p2Nug","26t48p3Nug","27t48p4Nug","28t48p1Sym","29t48p2Sym","30t48p3Sym","31t48p4Sym","32t72p1Nug","33t72p2Nug","34t72p3Nug","35t72p4Nug","36t72p1Sym","37t72p2Sym","38t72p3Sym","39t72p4Sym"), condition, genotype, time)

samples_by_metadata <- data.frame(Sample, Condition, Cultivar, Time)

```


```r
library(DESeq2)
## create countData object for deseq analysis and table of sig genes
countData_NugSym <- read_excel("C:/Users/sophi/Downloads/NEWVERSION_dovetailAssemblyFullNonRepeatAssociatedGeneList.xlsx", 
                               col_types = c("text", "skip", "skip", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))

View(countData_NugSym)
#Store GeneID info because the matrix doesnt like it. 
GeneID <- countData_NugSym$GeneID

#convert counts to matrix and subtract geneID column. 
countData2_NugSym <- as.matrix(countData_NugSym[ , -1])

#check work if you want
#view(countData2_NugSym)

#add gene IDs back... I knowm it's ridiculous. 
row.names(countData2_NugSym) <- GeneID   

#will need this later
Samples <- c(" 01t00p1Nug"," 02t00p2Nug"," 03t00p3Nug"," 04t00p4Nug", " 05t00p1Sym", " 06t00p2Sym","07t00p3Sym","08t00p4Sym","09t12p1Nug","10t12p2Nug","11t12p3Nug","12t12p1Sym","13t12p2Sym","14t12p3Sym","15t12p4Sym","16t24p1Nug","17t24p2Nug","18t24p3Nug","19t24p4Nug","20t24p1Sym","21t24p2Sym","22t24p3Sym","23t24p4Sym","24t48p1Nug","25t48p2Nug","26t48p3Nug","27t48p4Nug","28t48p1Sym","29t48p2Sym","30t48p3Sym","31t48p4Sym","32t72p1Nug","33t72p2Nug","34t72p3Nug","35t72p4Nug","36t72p1Sym","37t72p2Sym","38t72p3Sym","39t72p4Sym") 

#can check everything by uncommenting and running commands below. 
#View(countData2_NugSym)
#str(countData2_NugSym)

#need to set up meta data for DESeq2
condition = as.factor(c(rep("Uninoculated",8),rep("Inoculated",31))) 

genotype=as.factor(c(rep("Nugget", 4), rep("Symphony",4), rep("Nugget", 3), rep("Symphony", 4), rep("Nugget", 4), rep("Symphony", 4), rep("Nugget", 4), rep("Symphony", 4), rep("Nugget", 4), rep("Symphony", 4)))

time=as.factor(c(rep("0 hr", 8), rep("12 hr", 7), rep("24 hr", 8), rep("48 hr", 8), rep("72 hr", 8)))

colData = data.frame(col.names = c(" 01t00p1Nug"," 02t00p2Nug"," 03t00p3Nug"," 04t00p4Nug", " 05t00p1Sym", " 06t00p2Sym","07t00p3Sym","08t00p4Sym","09t12p1Nug","10t12p2Nug","11t12p3Nug","12t12p1Sym","13t12p2Sym","14t12p3Sym","15t12p4Sym","16t24p1Nug","17t24p2Nug","18t24p3Nug","19t24p4Nug","20t24p1Sym","21t24p2Sym","22t24p3Sym","23t24p4Sym","24t48p1Nug","25t48p2Nug","26t48p3Nug","27t48p4Nug","28t48p1Sym","29t48p2Sym","30t48p3Sym","31t48p4Sym","32t72p1Nug","33t72p2Nug","34t72p3Nug","35t72p4Nug","36t72p1Sym","37t72p2Sym","38t72p3Sym","39t72p4Sym"), condition, genotype, time)

#head(colData)

#create DESEQ2 object
dds_nug_sym <- DESeqDataSetFromMatrix(countData = countData2_NugSym, colData = colData, ~ genotype + condition + genotype:condition)
nrow(dds_nug_sym)

dds_nug_sym2 <- DESeqDataSetFromMatrix(countData= countData2_NugSym, colData = colData, ~ genotype + genotype:condition)
nrow(dds_nug_sym2)

```

```r
dds <- estimateSizeFactors(dds_nug_sym)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

deseq2Data <- DESeq(dds_nug_sym)

deseq2Results <- results(deseq2Data)
deseq2Results
```



```r
gene_counts_by_sample <- read_excel("C:/Users/sophi/Downloads/NEWVERSION_dovetailAssemblyFullNonRepeatAssociatedGeneList.xlsx", 
                               col_types = c("text", "skip", "skip", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))

gene_counts_by_sample
```

```r
library(phyloseq)
library(ggplot2)
theme_set(theme_bw())
###########otu table###########
###############################
GeneID <- gene_counts_by_sample$GeneID
countData2_NugSym <- as.matrix(gene_counts_by_sample[ , -1])
row.names(countData2_NugSym) <- GeneID
OTU1<-otu_table(countData2_NugSym, taxa_are_rows = TRUE)
head(OTU1)
```
```r
tax <- read.csv("C:/Users/sophi/Downloads/doveTail Gene Count With GO Terms.csv", header = TRUE)
#head(tax)

GeneID <- tax$geneID
tax2 <- as.matrix(tax[,-1])
row.names(tax2) <- GeneID
head(tax2)
#view(tax2)
tax3<-tax2[,c(-11,-12,-13,-14,-15,-16,-17,-18,-19,-20,-21,-22,-23,-24,-25,-26,-27,-28,-29,-30,-31,-32,-33,-34,-35,-36,-37,-38,-39,-40,-41,-42,-43,-44,-45,-46,-47,-48,-49)]
#view(tax3)
#head(tax3) # matrix

taxa_df<-data.frame(tax3) #taxa df
#taxa_df
tax4 <- as.matrix(taxa_df)
#head(tax4)
```


```r
TAX = tax_table(tax4)
```

```r
head(TAX)
```



```r
#TAX
class(TAX)
```

```r
physeq = phyloseq(OTU1, TAX)
```

```r
physeq
```



```r
Sample <- samples_by_metadata$Sample
samples <- as.matrix(samples_by_metadata[ , -1])
row.names(samples) <- Sample
```

```r
samples
```



```r
#starting after youve done DESeq analysis.

# only select significant values
res.05 <- subset(deseq2Results, padj<.05)

#convert to dataframe
res.05 <- as.data.frame(res.05)
#View(res.05)

#make row namess a column
res.05 <- tibble::rownames_to_column(res.05, "geneID")
```

```r
head(countData_NugSym)
```
```r
#head(res.05)

sig_genes <- res.05[1]
```

```r
sig_genes
```
```r
names(sig_genes)[1] <- "GeneID"
```



```r
library(dplyr)
res.05_gene_counts <- left_join(sig_genes, countData_NugSym, by = 'GeneID') 
```

```r
tail(res.05_gene_counts)
```

[HopBase annotation data](http://hopbase.cqls.oregonstate.edu/dovetailDownloads/dovetailCascadeMasked.php)

```r
#import molecular function df - you can download this off HopBase if you dont have it
library(readr)

goTerms_molecular <- read_tsv('C:/Users/sophi/Downloads/goTerms_molecularFunction.tsv.gz', col_names=FALSE)
names(goTerms_molecular) <- c('geneID', 'UNIProt_ID', 'GO ID', 'Description', 'Key terms')

#merge go terms to results df
res.05<-merge(res.05, goTerms_molecular, by.a="geneID", by.b="geneID", all.x=FALSE, all.y=FALSE)

#use dplyr to pull out the rows with uniq geneID... this drops the other molecular key terms though... can't figure out how to combine them. 
res.05_unique <- distinct(res.05, geneID, .keep_all = TRUE)
```

```r
library(tidyr)
#drop na values
res.05_unique <- res.05_unique %>% 
  drop_na()

#subset only important stuff
res.05_subset<- res.05_unique[c("geneID","Key terms","log2FoldChange","padj")]
```

```r
res.05_unique
```
```r
library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
plot(random_tree)
```

```r
#format sample data table

sampledata = sample_data(samples_by_metadata)

row.names(sampledata) <- Sample

sample_data_ = sampledata[,-1]

sampledata = sample_data_
```


```r
#physeq with all 3 tables, no tree
physeq1 = merge_phyloseq(physeq, sampledata)
physeq1
```

```r
#plot a quick heatmap to visualize the data
plot_heatmap(physeq1)
```

```r
library(DESeq2)
dds <- estimateSizeFactors(dds_nug_sym)
#sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
#head(dds)
```

```r
#write normalized_counts for use in the future
write.csv(normalized_counts,"C:/Users/sophi/Downloads/normalized_counts1.csv", row.names = TRUE)
```

## Make new matrix with sig genes and remake phyloseq object


```r
ord_pcoa = ordinate(physeq1, "PCoA", "bray")
```

```r
plot_ordination(physeq1, ord_pcoa, type="taxa")
```

```r
head(TAX)
```


```r
plot_ordination(physeq1, ord_pcoa, type="taxa", color="Key.terms.cellular", title = "Cellular Key Terms (All Genes)") + 
  theme(legend.position="none")
```

```r
plot_ordination(physeq1, ord_pcoa, color="Cultivar", title = "Samples (Raw Data)")
```
```r
plot_ordination(physeq1, ord_pcoa, color="Time", title = "Samples (Raw Data)")
```

```r
# take a look at this function?
phyloseq_to_deseq2()
```


```r
norm_otu = otu_table(normalized_counts, taxa_are_rows = TRUE)
```

```r
head(norm_otu)
```

```r
physeq_norm = phyloseq(norm_otu, TAX, sampledata)
```


```r
physeq_norm
```

```r
ord_pcoa_norm = ordinate(physeq_norm, "PCoA", "bray")
```

```r
plot_ordination(physeq_norm, ord_pcoa_norm, type="taxa", color="Key.terms.cellular", title = "Taxa (Normalized Data)") + 
  theme(legend.position="none")
```
```r
plot_ordination(physeq_norm, ord_pcoa_norm, color="Time", title = "Samples (Normalized Data")
```

```r
plot_ordination(physeq_norm, ord_pcoa_norm, color="Cultivar", title = "Samples (Normalized Data)")
```

```r
wh0 = genefilter_sample(physeq1, filterfun_sample(function(x) x > 5), A=0.5*nsamples(physeq1))
physeq_pruned = prune_taxa(wh0, physeq1)
```

```r
physeq_pruned
```


```r
physeq_pruned = transform_sample_counts(physeq_pruned, function(x) 1E6 * x/sum(x))
```

```r
cellular.sum = tapply(taxa_sums(physeq_pruned), tax_table(physeq_pruned)[, "Key.terms.cellular"], sum, na.rm=TRUE)
top5cellular = names(sort(cellular.sum, TRUE))[1:5]
physeq_pruned = prune_taxa((tax_table(physeq_pruned)[, "Key.terms.cellular"] %in% top5cellular), physeq_pruned)
```

```r
physeq_pruned.ord <- ordinate(physeq_pruned, "PCoA", "bray")
p1 = plot_ordination(physeq_pruned, physeq_pruned.ord, type="taxa", color="Key.terms.cellular", title="Cellular Key Terms (Pruned and Significant Genes Only)")
print(p1)
```
```r
plot_ordination(physeq_pruned, physeq_pruned.ord, color="Cultivar", title = "Samples")
```

```r
p1 + facet_wrap(~Key.terms.cellular, 3)
```
```r
res.05_gene_counts
```


```r
# make sig genes table into OTU table 

Gene_ID <- res.05_gene_counts$GeneID
res.05_gene_counts_ma <- as.matrix(res.05_gene_counts[ , -1])
row.names(res.05_gene_counts_ma) <- Gene_ID
OTU_sig_genes<-otu_table(res.05_gene_counts_ma, taxa_are_rows = TRUE)
head(OTU_sig_genes)
```


```r
#new phyloseq with significant genes only
physeq_sig_genes = phyloseq(OTU_sig_genes, TAX, sampledata)
```

```r
physeq_sig_genes
```
```r
wh1 = genefilter_sample(physeq_sig_genes, filterfun_sample(function(x) x > 5), A=0.5*nsamples(physeq_sig_genes))
physeq_pruned1 = prune_taxa(wh1, physeq_sig_genes)
```

```r
# pruned and significant genes only
physeq_pruned1
```


```r
cellular.sum1 = tapply(taxa_sums(physeq_pruned1), tax_table(physeq_pruned1)[, "Key.terms.cellular"], sum, na.rm=TRUE)
top5cellular1 = names(sort(cellular.sum1, TRUE))[1:5]
physeq_pruned1 = prune_taxa((tax_table(physeq_pruned1)[, "Key.terms.cellular"] %in% top5cellular1), physeq_pruned1)
```

```r
physeq_pruned1.ord <- ordinate(physeq_pruned1, "PCoA", "bray")
p2 = plot_ordination(physeq_pruned1, physeq_pruned1.ord, type="taxa", color="Key.terms.cellular", title="Cellular Key Terms (Raw Data, Significant Genes Only)")
print(p2)
```
```r
p2 + facet_wrap(~Key.terms.cellular, 3)
```

```r
plot_ordination(physeq_pruned1, physeq_pruned1.ord, type="taxa", color="Key.terms.biological", title = "Biological Key Terms (Significant Genes Only)") + 
  theme(legend.position="none")
```

```r
bio.sum1 = tapply(taxa_sums(physeq_pruned1), tax_table(physeq_pruned1)[, "Key.terms.biological"], sum, na.rm=TRUE)
top5bio1 = names(sort(bio.sum1, TRUE))[1:5]
physeq_pruned2 = prune_taxa((tax_table(physeq_pruned1)[, "Key.terms.biological"] %in% top5bio1), physeq_pruned1)
```

```r
physeq_pruned2.ord <- ordinate(physeq_pruned2, "PCoA", "bray")
p3 = plot_ordination(physeq_pruned1, physeq_pruned2.ord, type="taxa", color="Key.terms.biological", title="Biological Key Terms (Significant Genes Only)")
print(p3)
```

```r
library(ggplot2)
```

```r
data("midwest", package = "ggplot2")
options(scipen=999)
```
```r
dim(midwest)
```

```r
ggplot(midwest, aes(x=area, y=poptotal)) + 
  geom_point() +
  geom_smooth(method="lm")
```

```r
ggplot(midwest, aes(x=area, y=poptotal)) + 
  geom_point(col="steelblue", size=2) +   # Set static color and size for points
  geom_smooth(method="lm", col="firebrick") +  # change the color of line
  coord_cartesian(xlim=c(0, 0.1), ylim=c(0, 1000000)) + 
  labs(title="Area Vs Population", subtitle="From midwest dataset", y="Population", x="Area", caption="Midwest Demographics")
```

```r
ggplot(midwest, aes(x=area, y=poptotal)) + 
  geom_point(aes(col=state), size=3) +  # Set color to vary based on state categories.
  geom_smooth(method="lm", col="firebrick", size=2) + 
  coord_cartesian(xlim=c(0, 0.1), ylim=c(0, 1000000)) + 
  labs(title="Area Vs Population", subtitle="From midwest dataset", y="Population", x="Area", caption="Midwest Demographics")
```
```r
ggplot(data=countData_NugSym, aes(x = GeneID, y=01t00p1Nug)) + 
  geom_point()
```

