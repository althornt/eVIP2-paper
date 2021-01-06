#https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#kallisto

library(tximportData)
library(tximport)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(DESeqDataSetFromTximport)
library(DESeq2)


dir <- "/Users/alexis/Desktop/eVIP2/tutorial_files/RNF43_kallisto_outputs"

#make conversion data frame
edb <- EnsDb.Hsapiens.v86
tx2gene = transcripts(edb , columns=c("tx_id", "gene_name"),return.type="DataFrame")
head(tx2gene)


################
# GFP vs 659   #
################

#set up 
sample_names = c("GFP_1","GFP_2","GFP_3","GFP_4","RNF43_G659fs_1","RNF43_G659fs_2","RNF43_G659fs_3","RNF43_G659fs_4")
samples <- data.frame(samples=sample_names)
samples
condition <- c("GFP", "GFP","GFP","GFP","659","659","659","659")
sampleTable <- data.frame(samples=samples, condition=condition)
sampleTable

files <- file.path(dir, samples$sample, "abundance.tsv")
names(files) <- samples$sample

files
all(file.exists(files))

#import 
txi.kallisto <- tximport(files, type = "kallisto",tx2gene = tx2gene, txOut = FALSE,ignoreTxVersion=TRUE)
head(txi.kallisto$counts)

rownames(sampleTable) <- colnames(txi.kallisto$counts)
sampleTable

#DESeq
dds <- DESeqDataSetFromTximport(txi.kallisto, sampleTable, ~condition)

colData(dds)$condition<-factor(colData(dds)$condition, levels=c('GFP','659'))
GFP_659_dds<-DESeq(dds)
GFP_659_res<-results(GFP_659_dds)
GFP_659_res<-GFP_659_res[order(GFP_659_res$padj),]
head(GFP_659_res)

## Merge with normalized count data
GFP_659_resdata <- merge(as.data.frame(GFP_659_res), as.data.frame(counts(GFP_659_dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(GFP_659_resdata)[1] <- "Gene"
head(GFP_659_resdata)
mcols(GFP_659_res,use.names=TRUE)
write.csv(as.data.frame(GFP_659_res),file='deseq2_outputs/GFP_vs_659_results_deseq2.csv')
#write.csv(as.data.frame(GFP_659_resdata),file='GFP_vs_659_results_deseq2_counts.csv')

###############
# GFP vs WT  #
#############

#set up 
sample_names = c("GFP_1","GFP_2","GFP_3","GFP_4","RNF43_WT_1","RNF43_WT_2","RNF43_WT_3","RNF43_WT_4")
samples <- data.frame(samples=sample_names)
samples
condition <- c("GFP", "GFP","GFP","GFP","WT","WT","WT","WT")
sampleTable <- data.frame(samples=samples, condition=condition)
sampleTable

files <- file.path(dir, samples$sample, "abundance.h5")
names(files) <- samples$sample

files
all(file.exists(files))

#import 
txi.kallisto <- tximport(files, type = "kallisto",tx2gene = tx2gene, txOut = FALSE,ignoreTxVersion=TRUE)
head(txi.kallisto$counts)

rownames(sampleTable) <- colnames(txi.kallisto$counts)
sampleTable

#DESeq
dds <- DESeqDataSetFromTximport(txi.kallisto, sampleTable, ~condition)

colData(dds)$condition<-factor(colData(dds)$condition, levels=c('GFP','WT'))
GFP_WT_dds<-DESeq(dds)
GFP_WT_res<-results(GFP_WT_dds)
GFP_WT_res<-GFP_WT_res[order(GFP_WT_res$padj),]
head(GFP_WT_res)

## Merge with normalized count data
GFP_WT_resdata <- merge(as.data.frame(GFP_WT_res), as.data.frame(counts(GFP_WT_dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(GFP_WT_resdata)[1] <- "Gene"
head(GFP_WT_resdata)
mcols(GFP_WT_res,use.names=TRUE)
write.csv(as.data.frame(GFP_WT_res),file='deseq2_outputs/GFP_vs_WT_results_deseq2.csv')
#write.csv(as.data.frame(GFP_WT_resdata),file='GFP_vs_WT_results_deseq2_counts.csv')


###############
# WT vs 659  #
#############

#set up 
sample_names = c("RNF43_WT_1","RNF43_WT_2","RNF43_WT_3","RNF43_WT_4","RNF43_659fs_1","RNF43_659fs_2","RNF43_659fs_3","RNF43_659fs_4")
samples <- data.frame(samples=sample_names)
samples
condition <- c("WT","WT","WT","WT","659","659","659","659")
sampleTable <- data.frame(samples=samples, condition=condition)
sampleTable

files <- file.path(dir, samples$sample, "abundance.h5")
names(files) <- samples$sample

files
all(file.exists(files))

#import 
txi.kallisto <- tximport(files, type = "kallisto",tx2gene = tx2gene, txOut = FALSE,ignoreTxVersion=TRUE)
head(txi.kallisto$counts)

rownames(sampleTable) <- colnames(txi.kallisto$counts)
sampleTable

#DESeq
dds <- DESeqDataSetFromTximport(txi.kallisto, sampleTable, ~condition)

colData(dds)$condition<-factor(colData(dds)$condition, levels=c('WT','659'))
WT_659_dds<-DESeq(dds)
WT_659_res<-results(WT_659_dds)
WT_659_res<-WT_659_res[order(WT_659_res$padj),]
head(WT_659_res)

## Merge with normalized count data
WT_659_resdata <- merge(as.data.frame(WT_659_res), as.data.frame(counts(WT_659_dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(WT_659_resdata)[1] <- "Gene"
head(WT_659_resdata)
mcols(WT_659_res,use.names=TRUE)
write.csv(as.data.frame(WT_659_res),file='deseq2_outputs/WT_vs_659_results_deseq2.csv')
#write.csv(as.data.frame(WT_659_resdata),file='WT_vs_659_results_deseq2_counts.csv')

################
# GFP vs 117   #
################

#set up 
sample_names = c("GFP_1","GFP_2","GFP_3","GFP_4","RNF43_R117fs_1","RNF43_R117fs_2","RNF43_r117fs_3","RNF43_r117fs_4")
samples <- data.frame(samples=sample_names)
samples
condition <- c("GFP", "GFP","GFP","GFP","117","117","117","117")
sampleTable <- data.frame(samples=samples, condition=condition)
sampleTable

files <- file.path(dir, samples$sample, "abundance.h5")
names(files) <- samples$sample

files
all(file.exists(files))

#import 
txi.kallisto <- tximport(files, type = "kallisto",tx2gene = tx2gene, txOut = FALSE,ignoreTxVersion=TRUE)
head(txi.kallisto$counts)

rownames(sampleTable) <- colnames(txi.kallisto$counts)
sampleTable

#DESeq
dds <- DESeqDataSetFromTximport(txi.kallisto, sampleTable, ~condition)

colData(dds)$condition<-factor(colData(dds)$condition, levels=c('GFP','117'))
GFP_117_dds<-DESeq(dds)
GFP_117_res<-results(GFP_117_dds)
GFP_117_res<-GFP_117_res[order(GFP_117_res$padj),]
head(GFP_117_res)

## Merge with normalized count data
GFP_117_resdata <- merge(as.data.frame(GFP_117_res), as.data.frame(counts(GFP_117_dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(GFP_117_resdata)[1] <- "Gene"
head(GFP_117_resdata)
mcols(GFP_117_res,use.names=TRUE)
write.csv(as.data.frame(GFP_117_res),file='deseq2_outputs/GFP_vs_117_results_deseq2.csv')
#write.csv(as.data.frame(GFP_117_resdata),file='GFP_vs_117_results_deseq2_counts.csv')


###############
# WT vs 117  #
#############

#set up 
sample_names = c("RNF43_WT_1","RNF43_WT_2","RNF43_WT_3","RNF43_WT_4","RNF43_117fs_1","RNF43_117fs_2","RNF43_117fs_3","RNF43_117fs_4")
samples <- data.frame(samples=sample_names)
samples
condition <- c("WT","WT","WT","WT","117","117","117","117")
sampleTable <- data.frame(samples=samples, condition=condition)
sampleTable

files <- file.path(dir, samples$sample, "abundance.h5")
names(files) <- samples$sample

files
all(file.exists(files))

#import 
txi.kallisto <- tximport(files, type = "kallisto",tx2gene = tx2gene, txOut = FALSE,ignoreTxVersion=TRUE)
head(txi.kallisto$counts)

rownames(sampleTable) <- colnames(txi.kallisto$counts)
sampleTable

#DESeq
dds <- DESeqDataSetFromTximport(txi.kallisto, sampleTable, ~condition)

colData(dds)$condition<-factor(colData(dds)$condition, levels=c('WT','117'))
WT_117_dds<-DESeq(dds)
WT_117_res<-results(WT_117_dds)
WT_117_res<-WT_117_res[order(WT_117_res$padj),]
head(WT_117_res)

## Merge with normalized count data
WT_117_resdata <- merge(as.data.frame(WT_117_res), as.data.frame(counts(WT_117_dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(WT_117_resdata)[1] <- "Gene"
head(WT_117_resdata)
mcols(WT_117_res,use.names=TRUE)
write.csv(as.data.frame(WT_117_res),file='deseq2_outputs/WT_vs_117_results_deseq2.csv')
#write.csv(as.data.frame(WT_117_resdata),file='WT_vs_117_results_deseq2_counts.csv')