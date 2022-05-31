Continued from [Nugget vs. Symphony](DESeq-Nugget_vs_Sym.md). 

# Unioculated plants vs. Inoculated plants. 

<img src="images/UninocVsInoc.png" width="400">

Since the time points seem to cluster closely, it might be more informative to look just at uninoculated/inoculated (it also makes for much simpler comparisons as opposed to our previous [comparisons](images/allway.png)). So, we made a new dds object to examine just that. 

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

#transform
#log transformation
rld <- rlog(dds) #log and vst transformations look pretty similiar. 
rld_mat <- assay(rld)

#vst transformation
vst <- vst(dds)
vst_mat <- assay(vst)

```

# Just look at symphony vs. nugget

```{r}
#difference between nugget and symphony
base_differences <- results(dds, contrast=c("genotype","Nugget","Symphony"))
base_differences<-data.frame(base_differences)

#separate out up-regulated genes (ie upregulated in nugget, down in sym)
base_differences_up <- base_differences %>%
  filter(padj < 0.05) %>%
  filter(log2FoldChange > 0)

#separate out down-regulated genes (ie upregulated in symphony, down in nugget)
base_differences_down <- base_differences %>%
  filter(padj < 0.05) %>%
  filter(log2FoldChange < 0)

#prepping annotation
colnames(base_differences_down)[1] <- "geneID"
colnames(base_differences_up)[1] <- "geneID" 
base_differences_down$geneID <- rownames(base_differences_down)
base_differences_up$geneID <- rownames(base_differences_up)

#add annotation terms
base_differences_down<-merge(base_differences_down, goTerms_biological, by.a="geneID", by.b="geneID", all.x=TRUE, all.y=FALSE) 
base_differences_down<-merge(base_differences_down, kegg_pathway, by.a="UNIProt_ID", by.b="UNIProt_ID", all.x=TRUE, all.y=FALSE)
base_differences_up<-merge(base_differences_up, goTerms_biological, by.a="geneID", by.b="geneID", all.x=TRUE, all.y=FALSE)
base_differences_up<-merge(base_differences_up, kegg_pathway, by.a="UNIProt_ID", by.b="UNIProt_ID", all.x=TRUE, all.y=FALSE)

base_differences_up <- base_differences_up[order(base_differences_up$padj),]
base_differences_down <- base_differences_down[order(base_differences_down$padj),]

View(base_differences_up) #upregulated in Nugget over symphony https://biit.cs.ut.ee/gplink/l/GW4RMMCBTo

#top gene nug over symphony
plotCounts(dds, gene=which.min(base_differences$padj), intgroup="genotype")

#quick heatmap comparing nug and sym
select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","genotype")])
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)

```
![heatmap](https://github.com/mswiseman/MB599/blob/main/images/heatmap3.png)

# Looking at the effect of the inoculation on both cultivars

```r
both_cult_treatment_diff = results(dds, contrast=c("condition","Uninoculated","Inoculated")) #number reflects uninoc as baseline
ix = which.min(both_cult_treatment_diff$padj) # most significant
both_cult_treatment_diff <- both_cult_treatment_diff[order(both_cult_treatment_diff$padj),] # sort
kable(both_cult_treatment_diff[1:100,-(3:4)])

#convert to df
both_cult_treatment_diff<-data.frame(both_cult_treatment_diff)

#separate out down-regulated genes (ie downreg in uninoc, up in inoc)
both_cult_treatment_diff_down <- both_cult_treatment_diff %>%
  filter(padj < 0.05) %>%
  filter(log2FoldChange < 0)

#separate out up-regulated genes (ie upregulated in uninoc, down in inoc)
both_cult_treatment_diff_up <- both_cult_treatment_diff %>%
  filter(padj < 0.05) %>%
  filter(log2FoldChange > 0)

#prepping annotation
colnames(both_cult_treatment_diff_down)[1] <- "geneID"
colnames(both_cult_treatment_diff_up)[1] <- "geneID" 
both_cult_treatment_diff_down$geneID <- rownames(both_cult_treatment_diff_down)
both_cult_treatment_diff_up$geneID <- rownames(both_cult_treatment_diff_up)

#add annotation terms
both_cult_treatment_diff_down<-merge(both_cult_treatment_diff_down, goTerms_molecular, by.a="geneID", by.b="geneID", all.x=TRUE, all.y=FALSE)
both_cult_treatment_diff_down<-merge(both_cult_treatment_diff_down, kegg_pathway, by.a="UNIProt_ID", by.b="UNIProt_ID", all.x=TRUE, all.y=FALSE)
both_cult_treatment_diff_up<-merge(both_cult_treatment_diff_up, goTerms_molecular, by.a="geneID", by.b="geneID", all.x=TRUE, all.y=FALSE)
both_cult_treatment_diff_up<-merge(both_cult_treatment_diff_up, kegg_pathway, by.a="UNIProt_ID", by.b="UNIProt_ID", all.x=TRUE, all.y=FALSE)

both_cult_treatment_diff_up <- both_cult_treatment_diff_up[order(both_cult_treatment_diff_up$padj),]
both_cult_treatment_diff_down <- both_cult_treatment_diff_down[order(both_cult_treatment_diff_down$padj),]

#print top 100 in table form
kable(both_cult_treatment_diff_up[1:100,-(3:4)]) #(ie upregulated in uninoc, down in inoc)
kable(both_cult_treatment_diff_down[1:100,-(3:4)]) #(ie downreg in uninoc, up in inoc)

```

# Looking at the effect of the inoculation on sym

```{r}
Sym_treat_effects <- results(dds, list( c("condition_Uninoculated_vs_Inoculated","genotypeSymphony.conditionUninoculated") ))
ix = which.min(Sym_treat_effects$padj) # most significant
Sym_treat_effects <- Sym_treat_effects[order(Sym_treat_effects$padj),] # sort

#convert to df
Sym_treat_effects <- data.frame(Sym_treat_effects)

#separate out down-regulated genes (ie downreg in uninoc, up in inoc)
Sym_treat_effects_down <- Sym_treat_effects %>%
  filter(padj < 0.05) %>%
  filter(log2FoldChange < 0)

#separate out up-regulated genes (ie upregulated in uninoc, down in inoc)
Sym_treat_effects_up <- Sym_treat_effects %>%
  filter(padj < 0.05) %>%
  filter(log2FoldChange > 0)

#prepping annotation
colnames(Sym_treat_effects_down)[1] <- "geneID"
colnames(Sym_treat_effects_up)[1] <- "geneID" 
Sym_treat_effects_down$geneID <- rownames(Sym_treat_effects_down)
Sym_treat_effects_up$geneID <- rownames(Sym_treat_effects_up)

#add annotation terms; downreg = downreg in uninoc
Sym_treat_effects_down<-merge(Sym_treat_effects_down, goTerms_molecular, by.a="geneID", by.b="geneID", all.x=TRUE, all.y=FALSE) 
Sym_treat_effects_down<-merge(Sym_treat_effects_down, kegg_pathway, by.a="UNIProt_ID", by.b="UNIProt_ID", all.x=TRUE, all.y=FALSE)
Sym_treat_effects_up<-merge(Sym_treat_effects_up, goTerms_molecular, by.a="geneID", by.b="geneID", all.x=TRUE, all.y=FALSE)
Sym_treat_effects_up<-merge(Sym_treat_effects_up, kegg_pathway, by.a="UNIProt_ID", by.b="UNIProt_ID", all.x=TRUE, all.y=FALSE)

#View(Sym_treat_effects_down)
```

#Volcano plot

```r

#Symphony

vol_data1 <- data.frame(geneID=row.names(Sym_treat_effects), pval=-log10(Sym_treat_effects$padj), lfc=Sym_treat_effects$log2FoldChange)
vol_data1 <- na.omit(vol_data1)
View(vol_data1)

vol_data1 <- mutate(vol_data1, color=case_when(
    vol_data1$lfc > 0 & vol_data1$pval > 1.3 ~ "Increased",
    vol_data1$lfc < 0 & vol_data1$pval > 1.3 ~ "Decreased",
    vol_data1$pval < 1.3 ~ "Nonsignificant"))
vol1 <- ggplot(vol_data1, aes(x=lfc, y=pval, color=color))

vol_data2 <- data.frame(geneID=Sym_treat_effects_up$geneID, pval=-log10(Sym_treat_effects_up$padj), lfc=Sym_treat_effects_up$log2FoldChange)
# remove na
vol_data2 <- na.omit(vol_data2)
# set upper and lower threshold
vol_data2 <- mutate(vol_data2, color=case_when(
    vol_data2$lfc > 0 & vol_data2$pval > 1.3 ~ "Increased",
    vol_data2$lfc < 0 & vol_data2$pval > 1.3 ~ "Decreased",
    vol_data2$pval < 1.3 ~ "nonsignificant"))
vol2 <- ggplot(vol_data2, aes(x=lfc, y=pval, color=color))

plot_grid(vol1 +  ggtitle(label="Volcano Plot") +
            geom_point(size=2.5, alpha=0.8, na.rm=T) +
            scale_color_manual(name="Directionality",
                               values=c(Increased="#008B00", Decreased="#CD4F39", nonsignificant="darkgray")) +
            theme_bw(base_size=14) +
            theme(legend.position="right") +
            xlab(expression(log[2]("Symphony Uninoc vs. Inoc"))) +
            ylab(expression(-log[10]("adjusted p-value"))) +
            geom_hline(yintercept=1.3, colour="darkgrey") +
            scale_y_continuous(trans="log1p"), vol2 + geom_point(size=2.5, alpha=0.8, na.rm=T) +
            scale_color_manual(name="Directionality",
                     values=c(Increased="#008B00", Decreased="#CD4F39",
                              nonsignificant="darkgray")) +
            theme_bw(base_size=14) +
            theme(legend.position="right") +
            xlab(expression(log[2]("Nugget Uninoc vs. Inoc"))) +
            ylab(expression(-log[10]("adjusted p-value"))) +
            geom_hline(yintercept=1.3, colour="darkgrey") +
            scale_y_continuous(trans="log1p"), byrow = TRUE, nrow = 2)
dev.off()
```

![volcanoplot](images/volcanoplot.png)

```r
#order by pvalue, genotype
resOrdered <- as.data.frame(res05_genotypes)

#merge go terms to results df - genotype
resOrdered <- tibble::rownames_to_column(resOrdered, "geneID")
resOrdered<-merge(resOrdered, goTerms_molecular, by.a="geneID", by.b="geneID", all.x=TRUE, all.y=FALSE)
resOrdered<-merge(resOrdered, kegg_pathway, by.a="UNIProt_ID", by.b="UNIProt_ID", all.x=TRUE, all.y=FALSE)

#sort by p-value and set cut-off to 0.05 for genotype comparison
resOrdered <- resOrdered[order(resOrdered$pvalue),]
resOrdered_pval_cutoff <- resOrdered %>%
  filter(padj < 0.05)

#subsample only results with less than 0.05 p-value and greater than 1-fold change or less than -1 fold change
resSig <- resOrdered_pval_cutoff[ resOrdered_pval_cutoff$padj < 0.05 & (resOrdered_pval_cutoff$log2FoldChange >1| resOrdered_pval_cutoff$log2FoldChange < -1), ]

#remove any rows that lack UNIProt_ID or Kegg_accession
resSig<-resSig %>%
    filter_at(vars(UNIProt_ID, Kegg_accession), all_vars(!is.na(.)))

#genotype comparison
write.csv(resSig, "resSig.csv")
```



# Venn diagram

```r
venn_prep_inoc <- results.condition %>%
  filter(padj < 0.05) %>%
  mutate(Change = if_else(log2FoldChange < 0, "Downreg_Uninoc", "Upreg_Uninoc"))

venn_prep_geno <- results.interaction.symvnug %>%
  filter(padj < 0.05) %>%
  mutate(Change = if_else(log2FoldChange < 0, "Downreg_Genotype", "Upreg_Genotype"))

View(venn_prep_geno)

x <- list(up_uninoc=rownames(venn_prep_inoc)[venn_prep_inoc$Change=="Upreg_Uninoc"],
          down_uninoc=rownames(venn_prep_inoc)[venn_prep_inoc$Change=="Downreg_Uninoc"],
          up_gen=rownames(venn_prep_geno)[venn_prep_geno$Change=="Upreg_Genotype"],
          down_gen=rownames(venn_prep_geno)[venn_prep_geno$Change=="Downreg_Genotype"])

ggVennDiagram(x, label_alpha = 0, category.names = c("Upreg Uninoc","Downreg Uninoc","Upreg Nugget", "Downreg Nugget"))
```

[Genes upregulated in Nugget as compared to Symphony](https://biit.cs.ut.ee/gplink/l/aNUtJlb0Ts)

[Genes downregulated in Nugget as compared to Symphony](https://biit.cs.ut.ee/gplink/l/uPNRUSYIQq)

```r

#pathogenesis related genes from Bhardwaj 2011
pathogenesis_proteins <- read_excel("Downloads/pathogenesis_proteins.xlsx")

# negative = downregulated in uninoculated; positive = upregulated in uninoculated
Condition_effects1_path <- merge(Condition_effects1, pathogenesis_proteins, by.a="UNIProt_ID", by.b="geneID", all.x=FALSE, all.y=FALSE)
View(Condition_effects1_path)

#subset pathology-related genes
#with UNIProt ID
genes_UniProt_Path <- unique(Condition_effects1_path$UNIProt_ID)

#with cascade dovetail ID
genes_Path <- unique(Condition_effects1_path$geneID)
View(genes_Path)

#top inoc upregulated gene
NPR1 <- plotCounts(dds, gene="HUMLU_CAS0027044.t1.p1", intgroup=c("condition", "genotype"), returnData=TRUE)

#top differentially expressed (pre/post inoculation) gene that's pathogenesis-related
RPM1 <- plotCounts(dds, gene="HUMLU_CAS0061636.t1.p1", intgroup=c("condition", "genotype"), returnData=TRUE)

#define our color pallette
colourPallette <- c("#7145cd","#bbcfc4","#90de4a","#cd46c1","#77dd8e","#592b79","#d7c847","#6378c9","#619a3c","#d44473","#63cfb6","#dd5d36","#5db2ce","#8d3b28","#b1a4cb","#af8439","#c679c0","#4e703f","#753148","#cac88e","#352b48","#cd8d88","#463d25","#556f73")

#plot the normalized gene counts of RPM1. 
ggplot(RPM1, aes(x=genotype, y=count, fill=time)) + geom_col(position = "dodge2") +  theme_bw() +  theme(axis.text.x=element_text(angle=15, hjust=1)) + scale_colour_manual(values=colourPallette) + guides(colour=guide_legend(ncol=3)) + ggtitle("Disease resistance protein RPM1")

#plot the normalized gene counts of NPR1. NPR1 is a key player in the SA pathway, so this makes sense. 
ggplot(NPR1, aes(x=genotype, y=count, fill=time)) + geom_col(position = "dodge2") +  theme_classic() +  theme(axis.text.x=element_text(angle=15, hjust=1)) + scale_colour_manual(values=colourPallette) + guides(colour=guide_legend(ncol=3)) + ggtitle("ARABIDOPSIS NONEXPRESSER OF PR GENES 1")

RPM1

#SA pathway: https://www.genome.jp/dbget-bin/www_bget?ath:AT1G64280
```
