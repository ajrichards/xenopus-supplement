#!/usr/bin/Rscript
## install notes
# > source("http://bioconductor.org/biocLite.R")
# > biocLite("DESeq2")

## samples notes
#Group 1: A,B,E,H
#Group 2: C
#Group 3: D,F,G

library("DESeq2")

args <-commandArgs(TRUE)

if (length(args)==0){
    print("ERROR: Did not specify counts file e.g 'Rscript runDESeq.R gg-genes-counts.csv'")    
    q()
}

countsFile = args[1]
outFile = args[2]

if (!file.exists(countsFile)){
    print("ERROR: invalid counts file")    
    q()
}

#### from counts file (count matrix)
#countData <- read.csv(countsFile,header=TRUE,row.names=1)
countData <- read.table(countsFile,header=TRUE,row.names=1)
countData <- round(countData)

colData <- data.frame(
    row.names=colnames(countData),
    condition=c("group1", "group1", "group2", "group3", "group1", "group3","group3","group1"),
    libType=c("paired", "paired", "paired","paired","paired","paired","paired","paired")
)

## create a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData=countData,colData=colData,design=~condition)

## set 'non-endurant' as the control
dds$condition <- factor(dds$condition,levels=c("group1","group2","group3"))

dds <- DESeq(dds)
rld <- rlog(dds, blind=TRUE)
res <- results(dds)
resOrdered <-res[order(res$padj),]
print(head(resOrdered))

## export results to csv
write.csv(as.data.frame(resOrdered),file=outFile,quote=FALSE)
write.csv(as.data.frame(assay(rld)),file=gsub("\\.csv","_samples.csv",outFile),quote=FALSE)
