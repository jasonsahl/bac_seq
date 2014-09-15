library("DESeq")
options(echo=TRUE)
Args <- commandArgs(TRUE);
countsTable <- read.table(Args[1], header=TRUE)
rownames(countsTable) <- countsTable$gene
condsTable <- read.table(Args[2])
conds <- as.numeric(condsTable[1,])
cds<-newCountDataSet(countsTable,conds)
cds<-estimateSizeFactors(cds)
cds = estimateDispersions(cds)
res<-nbinomTest(cds, Args[3], Args[4])
write.table(res,file="all_results.txt",sep="\t")


