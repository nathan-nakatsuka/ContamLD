#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(scipen=999)

files=read.table(paste(args[1],"/",args[2],sep=""),header=F)
Table=data.frame(1:nrow(files))
for(k in 1:nrow(files)){
args1=as.character(files[k,1])
args3=as.character(files[k,2])
setwd(paste(args[1],"/directories/",args1,"_jackknife",sep=""))
score=read.table(paste(args1,"_",args3,"panel_scores_reads_ExtCorr.jout",sep=""),header=F)
Table[k,1]=args1
Table[k,2]=score[1,1]-as.numeric(args[3])
Table[k,3]=score[1,2]
Table[k,4]=args3
}
colnames(Table) <- c("Sample_ID","ContamLD_ExtCorr_Estimate","ContamLD_ExtCorr_Estimate_se","Panel")
write.table(Table,file=paste(args[1],"/FinalContamScores_ExtCorr_",args[2],sep=""),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
