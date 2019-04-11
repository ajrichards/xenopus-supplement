#!/usr/bin/Rscript
## install notes
# > source("http://bioconductor.org/biocLite.R")
# > biocLite("edgeR")

## samples notes
#    A = endurant
#    B = non-endurant
#    C = non-endurant
#    D = endurant
#    E = endurant
#    F = non-endurant
#    G = endurant
#    H = non-endurant
## this script was modified from Trinity 

library(edgeR)

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

data <- read.csv(countsFile,header=TRUE,row.names=1,com='')
#data <- read.table(countsFile,header=TRUE,row.names=1,com='')
data <- round(data)

#data = read.table("/home/adam/sequencing/xenopus/dn-trinity/features/Trinity_genes.counts.matrix", header=T, row.names=1, com='')
col_ordering <- c(1,4,5,7,2,3,6,8)
rnaseqMatrix <- data[,col_ordering]
rnaseqMatrix <- round(rnaseqMatrix)
rnaseqMatrix <- rnaseqMatrix[rowSums(rnaseqMatrix)>=2,]
conditions <- factor(c(rep("endurant", 4), rep("non-endurant", 4)))

exp_study <- DGEList(counts=rnaseqMatrix,group=conditions)
exp_study <- calcNormFactors(exp_study)
exp_study <- estimateCommonDisp(exp_study)
exp_study <- estimateTagwiseDisp(exp_study)

et = exactTest(exp_study)
tTags = topTags(et,n=NULL)

#write.csv(as.data.frame(),file=outFile)
write.csv(tTags,file=outFile,quote=FALSE,row.names=TRUE)
print('complete')
