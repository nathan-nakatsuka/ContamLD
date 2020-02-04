#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(scipen=999)


files=read.table(paste(args[1],"/",args[2],sep=""),header=F)
for(j in 1:nrow(files)){
args1=files[j,1]
args3=files[j,2]
setwd(paste(args[1],"/directories/",args1,"_jackknife",sep=""))
jinfile=read.table(paste(args1,"_",args3,"_",args[3],"panel_fullreadpairs.jin",sep=""),header=F)
### Loop through chromosomes; get values
for(i in 1:22){
score=read.table(paste(args1,"_",args3,"_",args[3],"panel_1ind_noChr",i,"_topscore_fullreads.txt",sep=""),header=F)
jinfile[i,3]=score[1,1]
}
num=mean(jinfile[,3])
writejin <- file(paste(args1,"_",args3,"_",args[3],"panel_scores_reads_ExtCorr.jin",sep=""),"w")
write.table(jinfile,writejin,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
flush(writejin)
close(writejin)
writenum <- file(paste(args1,"_",args3,"_",args[3],"panel_scores_reads_ExtCorr_meanscore.txt",sep=""),"w")
write.table(num,writenum,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
flush(writenum)
close(writenum)
}
