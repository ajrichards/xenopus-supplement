#!/usr/bin/Rscript
## install notes
#   > install.packages("GSA")

## example from documentation
##### two class unpaired comparison
#   > y must take values 1,2
#   > set.seed(100)
#   > x<-matrix(rnorm(1000*20),ncol=20)
#   > dd<-sample(1:1000,size=100)
#   > u<-matrix(2*rnorm(100),ncol=10,nrow=100)
#   > x[dd,11:20]<-x[dd,11:20]+u
#   > y<-c(rep(1,10),rep(2,10))
#   > genenames=paste("g",1:1000,sep="")
##### create some random gene sets
#   > genesets=vector("list",50)
#   > for(i in 1:50){
#    > genesets[[i]]=paste("g",sample(1:1000,size=30),sep="")
# }
#geneset.names=paste("set",as.character(1:50),sep="")
#GSA.obj<-GSA(x,y, genenames=genenames, genesets=genesets, resp.type="Two class unpaired", nperms=100)
#GSA.listsets(GSA.obj, geneset.names=geneset.names,FDRcut=.5)

# samples notes
#    A = endurant
#    B = non-endurant
#    C = non-endurant
#    D = endurant
#    E = endurant
#    F = non-endurant
#    G = endurant
#    H = non-endurant

library("GSA")

## import the counts
countsFile <- file.path(path.expand("~"),'documents','hts-integrate-docs','xenopus','data','transcript-count-summary.csv')
countsData <- read.csv(countsFile,header=T)
geneNames <- countsData$transcript
x <- data.matrix(countsData[2:9])
sx <- apply(x,2,function(y) y - mean(y))
y <- c(2,1,1,2,2,1,2,1)

## import the gene sets
gmtFile <- file.path(path.expand("~"),'documents','hts-integrate-docs','xenopus','data','xeno-noiea.gmt')
geneset.obj<- GSA.read.gmt(gmtFile)

#### import the gene sets
#genesets=vector("list",50)
#for(i in 1:50){
#    genesets[[i]]=paste("g",sample(1:1000,size=30),sep="")
#}
#geneset.names=paste("set",as.character(1:50),sep="")
#
#print(typeof(genesets))

q()

## run GSA
GSA.obj<-GSA(x,y, genenames=genenames, genesets=genesets, resp.type="Two class unpaired", nperms=100)
GSA.listsets(GSA.obj, geneset.names=geneset.names,FDRcut=.5)


#to use "real" gene set collection, we read it in from a gmt file:
#
# geneset.obj<- GSA.read.gmt("file.gmt")
#
# where file.gmt is a gene set collection from GSEA collection or
# or the website http://www-stat.stanford.edu/~tibs/GSA, or one
# that you have created yourself. Then
# GSA.obj<-GSA(x,y, genenames=genenames, genesets=geneset.obj$genesets, resp.type="Two class unpaired", nperms=100)
#


q()


#colData <- data.frame(
#    row.names=colnames(countData),
#    condition=c("endurant", "non-endurant", "non-endurant", "endurant", "endurant", "non-endurant", "endurant","non-endurant"),
#    libType=c("paired", "paired", "paired", "paired","paired","paired","paired","paired")
#    )

#dds <- DESeq(dds)
#res <- results(dds)
#resOrdered <-res[order(res$padj),]
#print(head(resOrdered))

## plot the values
#pdf("deseq2.pdf")
#plotMA(res,main="DESeq2",ylim=c(-2.5,2.5))
#dvo <- dev.off()

## export results to csv
#write.csv(as.data.frame(resOrdered),file=deFile)
