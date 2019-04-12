#!/usr/bin/Rscript
## install notes
# > source("http://bioconductor.org/biocLite.R")
# > biocLite("DESeq2")

## samples notes
#    A = endurant
#    B = non-endurant
#    C = non-endurant
#    D = endurant
#    E = endurant
#    F = non-endurant
#    G = endurant
#    H = non-endurant

library("DESeq2")
library("pasilla")
library("Biobase")

#### from custom counts file
countsFile <- file.path(path.expand("~"),'documents','hts-integrate-docs','xenopus','data','transcript-count-summary.csv')
deFile <- file.path(path.expand("~"),'documents','hts-integrate-docs','xenopus','data','deseq.csv')
countData <- read.csv(countsFile,header=T,row.names=1)
colData <- data.frame(
    row.names=colnames(countData),
    condition=c("endurant", "non-endurant", "non-endurant", "endurant", "endurant", "non-endurant", "endurant","non-endurant"),
    libType=c("paired", "paired", "paired", "paired","paired","paired","paired","paired")
    )
## create a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData=countData,colData=colData,design=~condition)

## set 'non-endurant' as the control
dds$condition <- factor(dds$condition,levels=c("non-endurant","endurant"))

dds <- DESeq(dds)
res <- results(dds)
resOrdered <-res[order(res$padj),]
print(head(resOrdered))

## plot the values
pdf("deseq2.pdf")
plotMA(res,main="DESeq2",ylim=c(-2.5,2.5))
dvo <- dev.off()

## export results to csv
write.csv(as.data.frame(resOrdered),file=deFile)
