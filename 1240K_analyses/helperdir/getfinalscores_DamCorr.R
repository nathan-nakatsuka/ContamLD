#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(scipen=999)

files=read.table(paste(args[1],"/",args[2],sep=""),header=F)
Table=data.frame(1:nrow(files))
for(k in 1:nrow(files)){
args1=as.character(files[k,1])
args3=as.character(files[k,2])
setwd(paste(args[1],"/directories/",args1,"_jackknife",sep=""))
score=read.table(paste(args1,"_",args3,"panel_scores_reads_DamCorr.jout",sep=""),header=F)
score2=read.table(paste(args1,"_",args3,"panel_kscores_reads_DamCorr.jout",sep=""),header=F)
Uncorr_score=read.table(paste(args1,"_",args3,"panel_scores_reads_ExtCorr.jout",sep=""),header=F)
damageratio=read.table(paste(args1,"_damageratio.txt",sep=""),header=F)
UUscore = read.table(paste(args1,"_UU_",args3,"panel_1ind_noChr0_topscore_fullreads.txt",sep=""),header=F)
DUscore = read.table(paste(args1,"_DU_",args3,"panel_1ind_noChr0_topscore_fullreads.txt",sep=""),header=F)
Table[k,1]=args1
Table[k,2]=score[1,1]
Table[k,3]=score[1,2]
Table[k,4]=score2[1,1]
Table[k,5]=score2[1,2]
Table[k,6]=1-damageratio[1,2]
Table[k,7]=args3
Table[k,8]="None"
if(DUscore[1,1]>(UUscore[1,1]+0.02)){
Table[k,8]="Model_Misspecified"
}
if(Uncorr_score[1,1]>0.15){
Table[k,8]="Very_High_Contamination"
}}
colnames(Table) <- c("Sample_ID","ContamLD_DamCorr_Estimate","ContamLD_DamCorr_Estimate_se","K_score","K_score_se","Damage_Ratio","Panel","Warnings")
write.table(Table,file=paste(args[1],"/FinalContamScoresK_damageratio_",args[2],sep=""),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
