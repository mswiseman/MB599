library(PCAtools)
install.packages('PCAtools')
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


#load count data into R. We have to skipa annotation columns for now since this needs to be a matrix.
countData_NugSym <- read_excel("/users/rielywhiter/Desktop/capstone/NEWVERSION_dovetailAssemblyFullNonRepeatAssociatedGeneList.xlsx", 
                               col_types = c("text", "skip", "skip", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))


#Store GeneID info because the matrix doesnt like it. 
GeneID <- countData_NugSym$GeneID

#convert counts to matrix and subtract geneID column. 
countData2_NugSym <- as.matrix(countData_NugSym[ , -1])

#add gene IDs back... I knowm it's ridiculous. 
row.names(countData2_NugSym) <- GeneID   

#will need this later
Samples <- c(" 01t00p1Nug"," 02t00p2Nug"," 03t00p3Nug"," 04t00p4Nug", " 05t00p1Sym", " 06t00p2Sym","07t00p3Sym","08t00p4Sym","09t12p1Nug","10t12p2Nug","11t12p3Nug","12t12p1Sym","13t12p2Sym","14t12p3Sym","15t12p4Sym","16t24p1Nug","17t24p2Nug","18t24p3Nug","19t24p4Nug","20t24p1Sym","21t24p2Sym","22t24p3Sym","23t24p4Sym","24t48p1Nug","25t48p2Nug","26t48p3Nug","27t48p4Nug","28t48p1Sym","29t48p2Sym","30t48p3Sym","31t48p4Sym","32t72p1Nug","33t72p2Nug","34t72p3Nug","35t72p4Nug","36t72p1Sym","37t72p2Sym","38t72p3Sym","39t72p4Sym") 


#need to set up meta data for DESeq2
condition = as.factor(c(rep("Uninoculated",8),rep("Inoculated",31))) 

genotype=as.factor(c(rep("Nugget", 4), rep("Symphony",4), rep("Nugget", 3), rep("Symphony", 4), rep("Nugget", 4), rep("Symphony", 4), rep("Nugget", 4), rep("Symphony", 4), rep("Nugget", 4), rep("Symphony", 4)))

time=as.factor(c(rep("0 hr", 8), rep("12 hr", 7), rep("24 hr", 8), rep("48 hr", 8), rep("72 hr", 8)))

colData = data.frame(col.names = c(" 01t00p1Nug"," 02t00p2Nug"," 03t00p3Nug"," 04t00p4Nug", " 05t00p1Sym", " 06t00p2Sym","07t00p3Sym","08t00p4Sym","09t12p1Nug","10t12p2Nug","11t12p3Nug","12t12p1Sym","13t12p2Sym","14t12p3Sym","15t12p4Sym","16t24p1Nug","17t24p2Nug","18t24p3Nug","19t24p4Nug","20t24p1Sym","21t24p2Sym","22t24p3Sym","23t24p4Sym","24t48p1Nug","25t48p2Nug","26t48p3Nug","27t48p4Nug","28t48p1Sym","29t48p2Sym","30t48p3Sym","31t48p4Sym","32t72p1Nug","33t72p2Nug","34t72p3Nug","35t72p4Nug","36t72p1Sym","37t72p2Sym","38t72p3Sym","39t72p4Sym"), condition, genotype, time)


#create DESEQ2 object
dds_nug_sym <- DESeqDataSetFromMatrix(countData = countData2_NugSym, colData = colData, ~ genotype + condition + genotype:condition)
nrow(dds_nug_sym)

dds_nug_sym2 <- DESeqDataSetFromMatrix(countData= countData2_NugSym, colData = colData, ~ genotype + genotype:condition)
nrow(dds_nug_sym2)



dds <- estimateSizeFactors(dds_nug_sym)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
deseq2Data <- DESeq(dds_nug_sym)
deseq2Results <- results(deseq2Data)

#dds_nug_sym <- DESeqDataSetFromMatrix(countData = countData2_NugSym, colData = colData, ~ genotype + condition + genotype:condition)

############# OTU ##################
norm_otu = otu_table(normalized_counts, taxa_are_rows = TRUE)

########### sample ###########
condition = as.factor(c(rep("Uninoculated",8),rep("Inoculated",31))) 

genotype=as.factor(c(rep("Nugget", 4), rep("Symphony",4), rep("Nugget", 3), rep("Symphony", 4), rep("Nugget", 4), rep("Symphony", 4), rep("Nugget", 4), rep("Symphony", 4), rep("Nugget", 4), rep("Symphony", 4)))

time=as.factor(c(rep("0 hr", 8), rep("12 hr", 7), rep("24 hr", 8), rep("48 hr", 8), rep("72 hr", 8)))

rowData = data.frame(row.names = c(" 01t00p1Nug"," 02t00p2Nug"," 03t00p3Nug"," 04t00p4Nug", " 05t00p1Sym", " 06t00p2Sym","07t00p3Sym","08t00p4Sym","09t12p1Nug","10t12p2Nug","11t12p3Nug","12t12p1Sym","13t12p2Sym","14t12p3Sym","15t12p4Sym","16t24p1Nug","17t24p2Nug","18t24p3Nug","19t24p4Nug","20t24p1Sym","21t24p2Sym","22t24p3Sym","23t24p4Sym","24t48p1Nug","25t48p2Nug","26t48p3Nug","27t48p4Nug","28t48p1Sym","29t48p2Sym","30t48p3Sym","31t48p4Sym","32t72p1Nug","33t72p2Nug","34t72p3Nug","35t72p4Nug","36t72p1Sym","37t72p2Sym","38t72p3Sym","39t72p4Sym"), condition, genotype, time)


SAMPLE<-sample_data(rowData, errorIfNULL = TRUE)


######### taxa #############

goTerms_biological <- read_excel('/Users/rielywhiter/Desktop/capstone/biologicalprocesses.xlsx', col_names=FALSE)
goTerms_cellular <- read_excel('/Users/rielywhiter/Desktop/capstone/cellular.xlsx', col_names=FALSE)
goTerms_molecular <- read_excel('/Users/rielywhiter/Desktop/capstone/molecular.xlsx', col_names=FALSE)

names(goTerms_biological) <- c('geneID', 'UNIProt_ID', 'GO ID', 'Description', 'Key terms')
names(goTerms_cellular) <- c('geneID', 'UNIProt_ID', 'GO ID', 'Description', 'Key terms')
names(goTerms_molecular) <- c('geneID', 'UNIProt_ID', 'GO ID', 'Description', 'Key terms')


df.desc<-merge(goTerms_biological,goTerms_cellular, by.x="geneID", by.y="geneID", all.x=FALSE, all.y=TRUE)
df.desc2<-merge(df.desc, goTerms_molecular, by.a="geneID", by.b="geneID", all.x=FALSE, all.y=TRUE)

#removing extra Uniprot columns
df.desc2 <- df.desc2 %>% select(-UNIProt_ID.y)  
df.desc2 <- df.desc2  %>% select(-UNIProt_ID)

#renaming columns...
names(df.desc2)[names(df.desc2) == 'UNIProt_ID.x'] <- 'UNIProt_ID'
names(df.desc2)[names(df.desc2) == 'GO ID.x'] <- 'GO ID.biological'
names(df.desc2)[names(df.desc2) == 'Description.x'] <- 'Description.biological'
names(df.desc2)[names(df.desc2) == 'Key terms.x'] <- 'Key terms.biological'
names(df.desc2)[names(df.desc2) == 'GO ID.y'] <- 'GO ID.cellular'
names(df.desc2)[names(df.desc2) == 'Descriptio
                n.y'] <- 'Description.cellular'
names(df.desc2)[names(df.desc2) == 'Key terms.y'] <- 'Key terms.cellular'
names(df.desc2)[names(df.desc2) == 'GO ID'] <- 'GO ID.molecular'
names(df.desc2)[names(df.desc2) == 'Description'] <- 'Description.molecular'
names(df.desc2)[names(df.desc2) == 'Key terms'] <- 'Key terms.molecular'

#importing count matrix with Uniprot IDs (note, I didn't skip any columns)
countData_NugSym_with_annotations <- read_excel("/Users/rielywhiter/Desktop/capstone/NEWVERSION_dovetailAssemblyFullNonRepeatAssociatedGeneList.xlsx")

#match geneID case format
names(countData_NugSym_with_annotations)[names(countData_NugSym_with_annotations) == 'GeneID'] <- 'geneID'

#merge the GO terms file with the count file
df.desc_all<-merge(df.desc2, countData_NugSym_with_annotations, by.a="geneID", by.b="geneID", all.a=FALSE, all.b=TRUE)

#remove a few more names
df.desc_all <- df.desc_all  %>% select(-`UniProt Identifier`) 
df.desc_all <- df.desc_all  %>% select(-`Description`)

#change a few column names 
names(df.desc_all)[names(df.desc_all) == 'UniProt_ID.x'] <- 'UniProt_ID'
names(df.desc_all)[names(df.desc_all) == 'GO ID.x'] <- 'GO ID.biological'
names(df.desc_all)[names(df.desc_all) == 'Description.x'] <- 'Description.biological'
names(df.desc_all)[names(df.desc_all) == 'Key terms'] <- 'Key terms.biological'
names(df.desc_all)[names(df.desc_all) == 'GO ID.y'] <- 'GO ID.cellular'
names(df.desc_all)[names(df.desc_all) == 'Description.y'] <- 'Description.cellular'
names(df.desc_all)[names(df.desc_all) == 'Key terms.y'] <- 'Key terms.cellular'

#write to CSV
write.csv(df.desc_all,"/Users/rielywhiter/Desktop/capstone/Gene_Count_With_GO_Terms3.csv", row.names = FALSE)

tax<- read.csv("/Users/rielywhiter/Desktop/capstone/Gene_Count_With_GO_Terms3.csv", header = TRUE)


GeneID <- tax$geneID
tax2 <- as.matrix(tax[,-1])
row.names(tax2) <- GeneID

tax3<-tax2[,c(-11,-12,-13,-14,-15,-16,-17,-18,-19,-20,-21,-22,-23,-24,-25,-26,-27,-28,-29,-30,-31,-32,-33,-34,-35,-36,-37,-38,-39,-40,-41,-42,-43,-44,-45,-46,-47,-48,-49)]


tax5<-data.frame(tax3)
tax4 <- as.matrix(tax5)
TAX = tax_table(tax4)


example(enterotype, ask=FALSE)

physeq3 = phyloseq(norm_otu, TAX,SAMPLE ) 

#starting after youve done DESeq analysis.
# only select significant values
res.05 <- subset(deseq2Results, padj<.05)
#convert to dataframe
res.05 <- as.data.frame(res.05)

#make row namess a column
res.05 <- tibble::rownames_to_column(res.05, "geneID")

head(countData_NugSym)
sig_genes <- res.05[1]
sig_genes

names(sig_genes)[1] <- "GeneID"
library(dplyr)
res.05_gene_counts <- left_join(sig_genes, countData_NugSym, by = 'GeneID')
tail(res.05_gene_counts)
Gene_ID <- res.05_gene_counts$GeneID
res.05_gene_counts_ma <- as.matrix(res.05_gene_counts[ , -1])
row.names(res.05_gene_counts_ma) <- Gene_ID
OTU_sig_genes<-otu_table(res.05_gene_counts_ma, taxa_are_rows = TRUE)
physeq_sig_genes = phyloseq(OTU_sig_genes, TAX, SAMPLE)


#################### pcoa for non-normalized data ############
#compare normalized vs non normalized

GP.ord<-ordinate(physeq_sig_genes, "PCoA","bray") #another distance + ordination 

ord_pcoa=plot_ordination(physeq_sig_genes, GP.ord,title = "PCOA Plot of Samples Cultivar ",color="genotype") +
  theme(legend.position="right")
ord_pcoa


################### PERMANOVA #######################
install.packages("vegan")
library(vegan)

my_phyloseq_bray<-phyloseq::distance (physeq_sig_genes, method = "bray")
sampledf<-data.frame(sample_data(physeq_sig_genes))
table<-adonis2(my_phyloseq_bray~genotype, data=sampledf)


####### principle of random forest #############
install.packages('e1071')
install.packages('randomForest')
install.packages('caret')
library(ape)
install.packages('glmnet')
library(glmnet)
library(caret)
library(phyloseq)
library(randomForest)
library(e1071)

########## randomforest #####
ps<-physeq_sig_genes

set.seed(1)
index_train <- createDataPartition(ps@sam_data$genotype, p = 0.7)[[1]]
x_train <- ps@otu_table[, index_train]
x_test <- ps@otu_table[,-index_train]

#split the phyloseq objects into training and testing to make our lives easier later on
ps_train <- phyloseq(otu_table(x_train, taxa_are_rows = TRUE), ps@sam_data[index_train, ])
ps_test <- phyloseq(otu_table(x_test, taxa_are_rows = TRUE), ps@sam_data[-index_train, ])


#set.seed(1)
#ttest<-data.frame(t(physeq3@otu_table))
x_train<-data.frame(t(ps_train@otu_table))
g_role = ps_train@sam_data$genotype
control <- trainControl(method='repeatedcv', 
                        number=3, 
                        repeats=3,
                        allowParallel = F)


tunegrid <- expand.grid(.mtry=c(1:20)) #mtry is the depth of each decision tree. We'll be trying out models where each tree is 3 to 20 splits deep
rf <- train(x_train,
            g_role,
            method='rf', 
            metric='Accuracy', 
            tuneGrid=tunegrid, 
            trControl=control)
print(rf)


mtry_best = as.numeric(rf$bestTune)
model = randomForest(x_train, y = as.factor(ps_train@sam_data$genotype), mtry = mtry_best)

#Performance on test set
preds = predict(model, t(x_test)) ###### stuck here ######
print(paste("Accuracy: ", sum(preds == t(as.factor(ps_test@sam_data$genotype))) / t(nsamples(ps_test))))




#Visualize on test dataset

ord <- ordinate(ps_test, method = "PCoA", distance = 'bray')
ps_test@sam_data$rf_predictions = predict(model,t(ps_test@otu_table)) ##stuck here

plot_ordination(ps_test, ord, 'samples', color = 'genotype', shape = 'rf_predictions') + geom_point(size = 4)


model = randomForest(t(ps@otu_table), y = as.factor(ps@sam_data$genotype), mtry = mtry_best) #####something wrong
varImpPlot(model, type = 2)

imp_list <- list()
for(i in 1:50){
  model = randomForest(t(ps@otu_table), y = as.factor(ps@sam_data$genotype), mtry = mtry_best)
  imp_list[i] <- varImp(model)
}

imp_df <- do.call(rbind.data.frame, imp_list)
colnames(imp_df) <- colnames(x_train)
view(colMeans(imp_df))

barplot(cex.names=0.28,cex.axis = 0.5,sort(colMeans(imp_df)), horiz = T, las =1 , xlab = "Mean variable importance")
