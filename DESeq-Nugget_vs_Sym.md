---
title: "Nugget vs Sym"
author: "Michele Wiseman"
date: "May 11th, 2022"
---

# Nugget vs Sym

```{r load packages, message=FALSE, warning=FALSE, include=FALSE}
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



```{r loading data}

#load count data into R. We have to skipa annotation columns for now since this needs to be a matrix.
countData_NugSym <- read_excel("Desktop/NEWVERSION_dovetailAssemblyFullNonRepeatAssociatedGeneList.xlsx", 
                        col_types = c("text", "skip", "skip", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "skip", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "skip", "numeric"))

#View(countData_NugSym)
#Store GeneID info because the matrix doesnt like it. 
GeneID <- countData_NugSym$GeneID

#convert counts to matrix and subtract geneID column. 
countData2_NugSym_mtx <- as.matrix(countData_NugSym[ , -1])

#check work if you want
#View(countData2_NugSym)

#add gene IDs back... I knowm it's ridiculous. 
row.names(countData2_NugSym_mtx) <- GeneID  

```

# Some data visualization 
## Check raw reads for trends

```{r pre-processing visuals}
row.names(countData_NugSym) <- GeneID   
countData_NugSym = data.frame(countData_NugSym)
countData_NugSym <- countData_NugSym[,-1]
#View(countData_NugSym)

#Quick view of raw reads
par(mar=c(8,4,4,1)+0.1)
barplot(colSums(countData_NugSym)/1e6, las =3)

#clearly need to transform
hist(countData_NugSym$X01t00p1Nug, br=100)

#log transformation seems to normalize it pretty well 
logcountdata = log2(1+countData_NugSym)
hist(logcountdata$X01t00p1Nug, br=100)

#pretty strongly correlated within nugget uninoculated treatments
plot(logcountdata[,1], logcountdata[,2])

#A lot more differences between treatments (less correlation)
plot(logcountdata[,1], logcountdata[,12])

```

# Back to DESeq prep

```{r deseq setup}

#will need this later
Samples <- c(" 01t00p1Nug"," 02t00p2Nug"," 03t00p3Nug"," 04t00p4Nug", " 05t00p1Sym", " 06t00p2Sym","07t00p3Sym","09t12p1Nug","10t12p2Nug","11t12p3Nug","12t12p1Sym","13t12p2Sym","14t12p3Sym","15t12p4Sym","16t24p1Nug","17t24p2Nug","18t24p3Nug","19t24p4Nug","20t24p1Sym","21t24p2Sym","22t24p3Sym","23t24p4Sym","24t48p1Nug","25t48p2Nug","26t48p3Nug","27t48p4Nug","28t48p1Sym","29t48p2Sym","30t48p3Sym","31t48p4Sym","32t72p1Nug","33t72p2Nug","34t72p3Nug","35t72p4Nug","36t72p1Sym","37t72p2Sym","39t72p4Sym") 

#set up experiment metadata
genotype=as.factor(c(rep("Nugget", 4), rep("Symphony",3), rep("Nugget", 3), rep("Symphony", 4), rep("Nugget", 4), rep("Symphony", 4), rep("Nugget", 4), rep("Symphony", 4), rep("Nugget", 4), rep("Symphony", 3)))
time=as.factor(c(rep("0hr", 7), rep("12hr", 7), rep("24hr", 8), rep("48hr", 8), rep("72hr", 7)))
colData = data.frame(col.names = c(" 01t00p1Nug"," 02t00p2Nug"," 03t00p3Nug"," 04t00p4Nug", " 05t00p1Sym", " 06t00p2Sym","07t00p3Sym","09t12p1Nug","10t12p2Nug","11t12p3Nug","12t12p1Sym","13t12p2Sym","14t12p3Sym","15t12p4Sym","16t24p1Nug","17t24p2Nug","18t24p3Nug","19t24p4Nug","20t24p1Sym","21t24p2Sym","22t24p3Sym","23t24p4Sym","24t48p1Nug","25t48p2Nug","26t48p3Nug","27t48p4Nug","28t48p1Sym","29t48p2Sym","30t48p3Sym","31t48p4Sym","32t72p1Nug","33t72p2Nug","34t72p3Nug","35t72p4Nug","36t72p1Sym","37t72p2Sym","39t72p4Sym"), genotype, time)

#create DESEQ2 object genotype and condition plus their interaction
dds <- DESeqDataSetFromMatrix(countData = countData2_NugSym_mtx, colData = colData, ~ genotype + time + genotype:time)

#delete rows with less than 10 reads
dds <-dds[ rowSums(counts(dds)) > 10, ]

```

```{r running deseq}

#run DESeq... this performs the median of ratios normalization method
dds <- DESeq(dds)

#check out the normalization
sizeFactors(dds)

## Total number of raw counts per sample
colSums(counts(dds))

## Total number of normalized counts per sample
colSums(counts(dds, normalized=T))

#extract normalized counts
normalized_counts <- counts(dds, normalized=TRUE)
#View(normalized_counts)

# write to table
write.table(normalized_counts, file="normalized_counts.txt", sep="\t", quote=F, col.names=NA)

#NOTE: DESeq2 doesn’t actually use normalized counts, rather it uses the raw counts and models the normalization inside the Generalized Linear Model (GLM). These normalized counts will be useful for downstream visualization of results, but cannot be used as input to DESeq2 or any other tools that perform differential expression analysis which use the negative binomial model.

## Plot dispersion estimates
plotDispEsts(dds)

```

```{r printing deseq results, still as dds object}

#lets just compare nugget and symphony
contrast_genotypes <- c("genotype", "Nugget", "Symphony")
res005_genotypes <- results(dds, contrast=contrast_genotypes, alpha=0.0025)

#lets look at time differences (this unfortunately clumps in nugget and symphony)
contrast_time <- c("time", "0hr", "72hr")
res005_time <- results(dds, contrast=contrast_time, alpha=0.0025)
#View(res005_time)


#visulize sig points
plotMA(res005_genotypes, ylim=c(-2,2)) 
plotMA(res005_time, ylim=c(-2,2))
 
```

```{r convert to df}
#order by pvalue, genotype
resOrdered <- as.data.frame(res005_genotypes)

#for time
resTime <- as.data.frame(res005_time)

#import molecular function df
goTerms_molecular <- read_tsv('/Users/michelewiseman/Desktop/goTerms_molecularFunction(1).tsv', col_names=FALSE)
names(goTerms_molecular) <- c('geneID', 'UNIProt_ID', 'GO ID', 'Description', 'Key terms')

#import kegg pathway mapp
kegg_pathway <- read_csv('/Users/michelewiseman/Desktop/kegg_pathway_map.csv', col_names=TRUE)

#merge go terms to results df - genotype
resOrdered <- tibble::rownames_to_column(resOrdered, "geneID")
resOrdered<-merge(resOrdered, goTerms_molecular, by.a="geneID", by.b="geneID", all.x=TRUE, all.y=FALSE)
resOrdered<-merge(resOrdered, kegg_pathway, by.a="UNIProt_ID", by.b="UNIProt_ID", all.x=TRUE, all.y=FALSE)

#merge go terms to results df - time
resTime <- tibble::rownames_to_column(resTime, "geneID")
resTime<-merge(resTime, goTerms_molecular, by.a="geneID", by.b="geneID", all.x=TRUE, all.y=FALSE)
resTime<-merge(resTime, kegg_pathway, by.a="UNIProt_ID", by.b="UNIProt_ID", all.x=TRUE, all.y=FALSE)

#sort by p-value and set cut-off to 0.05
resOrdered <- resOrdered[order(resOrdered$pvalue),]
resOrdered_pval_cutoff <- resOrdered %>%
  filter(padj < 0.05)

#subsample only results with less than 0.01 p-value and greater than 1-fold change or less than -1 fold change
resSig <- resOrdered_pval_cutoff[ resOrdered_pval_cutoff$padj < 0.01 & (resOrdered_pval_cutoff$log2FoldChange >1| resOrdered_pval_cutoff$log2FoldChange < -1), ]

#remove any rows that lack UNIProt_ID or Kegg_accession
resSig<-resSig %>%
    filter_at(vars(UNIProt_ID, Kegg_accession), all_vars(!is.na(.)))

#write to csv
write.csv(resOrdered_pval_cutoff, "Nugget_vs_Sym_results_with_Kegg_path_and_goTermsMolecular.csv")

```


```{r transformation of counts}
#transform data
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)

#rename row names with descriptive names
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(genotype, time, sep=" : "))

```

# PCA plots

```{r PCA Plots}
rv <- rowVars(assay(rld))
select <- order(rv, decreasing=T)[seq_len(min(500,length(rv)))]
pc <- prcomp(t(assay(vsd)[select,]))

condition <- genotype
scores <- data.frame(pc$x, condition)

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(mat, trace = "none", col = rev(hmcol), margin = c(13,13))

#ggplot pcaplot
(pcaplot <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(condition))))
  + geom_point(size = 5)
  + ggtitle("Principal Components")
  + scale_colour_brewer(name = " ", palette = "Set1")
  + theme(
    plot.title = element_text(face = 'bold'),
    legend.position = c(.9,.2),
    legend.key = element_rect(fill = 'NA'),
    legend.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(colour = "Black"),
    axis.text.x = element_text(colour = "Black"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = 'bold'),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(color = 'black',fill = NA)
  ))

#transpose and then do comphrensive PCA.
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))

#Make a new df that binds metadata and PCA values
df <- cbind(colData, pca$x)

#vst version
vst_mat<-assay(vsd)
pca <- prcomp(t(vst_mat))
vst_cor<-cor(vst_mat)
pheatmap(vst_cor, border_color=NA, fontsize = 10, 
  		fontsize_row = 6, height=20)

#fancy PCA by individuals
fviz_pca_ind(pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )
#fancy PCA by variables
fviz_pca_var(pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )

fviz_pca_biplot(pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
                )

groups <- as.factor(time)

fviz_pca_ind(pca,
             col.ind = groups, # color by groups
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
             )
```

```{r heat maps}
#correlation matrix for heat map
rld_cor<-cor(rld_mat)

#heat map
pheatmap(rld_cor, border_color=NA, fontsize = 10, 
  		fontsize_row = 6, height=20)

#scree plot
fviz_eig(pca) 

```



# Another form of heatmap

```{r}
topGene <- resOrdered$UNIProt_ID[which.min(resOrdered$padj)]
View(topGene)
res

plotMA(resOrdered, ylim=c(-5,5)) 

with(resOrdered[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

plotDispEsts(dds)

#showing p-value across all genes
hist(res005_genotypes$padj, breaks=20, col="grey50", border="white")

#showing p-value after subsampling
hist(resSig$padj, breaks=20, col="grey50", border="white")

#looking at p-values were mean count is greater than 10
hist(res005_genotypes$padj[res005_genotypes$baseMean > 10], breaks=20, col="grey50", border="white")
hist(resSig$padj[resSig$baseMean > 10], breaks=20, col="grey50", border="white")

colors <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255)
sidecols <- c("grey","dodgerblue")[ rld$col.names ]
topVarGenes <- head(order(-rowVars(assay(rld))),35)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
colnames(mat) <- paste0(rld$genotype,"-",rld$time)
heatmap.2(mat, trace="none", col=colors, ColSideColors=sidecols,
          labRow=FALSE, mar=c(10,2), scale="row")


```

# Reduced model (trying to remove effect of genotype). 

```{r}
ddsLRT <- DESeq(dds, test="LRT", reduced = ~ genotype + time)
resLRT <- results(ddsLRT)
resLRT$symbol <- mcols(ddsLRT)$symbol
head(resLRT[order(resLRT$pvalue),],4)

#show most differentially expressed gene over time
data <- plotCounts(ddsLRT, which.min(resLRT$pvalue), 
                   intgroup=c("time","genotype"), returnData=TRUE)

View(data)
ggplot(data, aes(x=time, y=count, color=genotype, group=genotype)) + 
  geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10()

```
