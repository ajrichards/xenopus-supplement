#!/usr/bin/Rscript
## install notes
#   > install.packages("GSA")

library("GSA")
library("stats")

# samples notes
#    A = endurant
#    B = non-endurant
#    C = non-endurant
#    D = endurant
#    E = endurant
#    F = non-endurant
#    G = endurant
#    H = non-endurant

args <-commandArgs(TRUE)

if (length(args)==0){
    print("ERROR: Did not specify name-aspect e.g 'Rscript runGSA.R x-noiea'")
    q()
}

nameBase <- args[1]

run_gsa <- function(assembly,aspect){

    geneSetsFile <- file.path(path.expand("~"),'sequencing','xenopus','gsa',nameBase,
                              paste(assembly,"-",aspect,".gmt",sep=""))

    if(assembly == 'ref'){
        countsFile <- file.path(path.expand("~"),'sequencing','xenopus','reference',
                                "features","deseq-samples.csv")
    }
    else{
        countsFile <- file.path(path.expand("~"),'sequencing','xenopus',paste(assembly,'-trinity',sep=""),
                                "features",paste("deseq-samples.csv",sep=""))
    }
    if (!file.exists(geneSetsFile)){
        print("ERROR: invalid gene set file")
        print(geneSetsFile)
        q()
    }

    if (!file.exists(countsFile)){
        print("ERROR: invalid counts file")
        print(countsFile)
        q()
    }

    ## load count matrix and specify phenotype labels
    countsData <- read.csv(countsFile,header=TRUE)
    names(countsData) <- c("id","A","B","C","D","E","F","G","H")
    print(names(countsData))
    geneNames <- countsData$id
    print(paste('gene name:', geneNames[1]))
    x <- data.matrix(countsData[2:9])
    sx <- apply(x,2,function(y) y - mean(y))
    y <- c(2,1,1,2,2,1,2,1)
    fdrCut <- 0.8
    
    ## read in gene sets and test
    geneset.obj <- GSA.read.gmt(geneSetsFile)
    GSA.obj<-GSA(sx,y, genenames=geneNames, genesets=geneset.obj$genesets,resp.type="Two class unpaired",
                 nperms=400,restand=TRUE,,minsize=10) ## data|catalog
    print(paste('assembly', assembly))
    print(GSA.listsets(GSA.obj,geneset.names=geneset.obj$geneset.names,FDRcut=fdrCut))
    print('done')

    return(GSA.listsets(GSA.obj,geneset.names=geneset.obj$geneset.names,FDRcut=fdrCut))
}

## create a function to play with the columns
process_results <- function(x, outFile,assembly,direction,aspect) {
    names(x) <- lapply(names(x),function(z) gsub("\\-","_",tolower(z)))
    geneSet <- x[[2]]
    score <- x[[3]]
    pValue <-  x[[4]]
    fdr <- x[[5]]
    rowToOutput <- paste(assembly,direction,geneSet,aspect,score,pValue,fdr,sep=",")
    #print(paste(assembly,transcript,direction,pValue, sep=","))
    cat(rowToOutput,file=outFile,append=T,fill=T)
}

## prepare the outfile
outFile <- file.path(path.expand("~"),'sequencing','xenopus','gsa',paste(nameBase,"-","geneset-results",".csv",sep=""))
newHeader <- paste("assembly","direction","gene_set","aspect","score","p-value","FDR",sep=",")
cat(newHeader,file=outFile,append=F,fill=T)

aspect = "bp"
assembly = "dn"

results = run_gsa(assembly,aspect)
apply(results$negative,1,process_results,outFile=outFile,assembly=assembly,direction='neg',aspect=aspect)
apply(results$positive,1,process_results,outFile=outFile,assembly=assembly,direction='pos',aspect=aspect)

#results = run_gsa("gg")
#apply(results$negative,1,process_results,outFile=outFile,assembly="gg",transcript="isoforms",direction='neg')
#apply(results$positive,1,process_results,outFile=outFile,assembly="gg",transcript="isoforms",direction='pos')

#results = run_gsa("dn")
#apply(results$negative,1,process_results,outFile=outFile,assembly="dn",transcript="isoforms",direction='neg')
#apply(results$positive,1,process_results,outFile=outFile,assembly="dn",transcript="isoforms",direction='pos')

print("done")
