#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(scipen=999)

### Make a table for storing values (first pass)
TableA=data.frame(0:600)
TableA[,1]=seq(-0.1,0.5,by=0.001)

alpha=rep(TableA[1:601,1],601)
k23=rep(TableA[1:601,1],each=601)
BigTable=data.frame(alpha,k23)
names(BigTable)=c("alpha","k")

## Extra table for storing more specific values.
number=401*401
BigTable2=data.frame(1:number)

files=read.table(paste(args[1],"/",args[2],sep=""),header=F)
for(s in 1:nrow(files)){
args1=as.character(files[s,1])
args3=as.character(files[s,2])
setwd(paste(args[1],"/directories/",args1,"_jackknife",sep=""))

jinfile=read.table(paste(args1,"_",args3,"_",args[3],"panel_fullreadpairs.jin",sep=""),header=F)
jinfile[,3]=0
jinfile2=jinfile

damageratio=read.table(paste(args1,"_",args3,"_",args[3],"panel_damageratio.txt",sep=""),header=F)

j=1
d=damageratio[(j+1),2]
## All possible alpha+k values
BigTable[,3]=round(BigTable[,1]+BigTable[,2],digits=4)
## All possible alpha/2+k values
BigTable[,4]=round(BigTable[,1]/2+BigTable[,2],digits=4)

TableA = read.table(paste(args1,"_UU_",args3,"_",args[3],"panel_1ind_noChr",j,"_alphatable_fullreads.txt",sep=""),header=F)
TableB = read.table(paste(args1,"_DU_",args3,"_",args[3],"panel_1ind_noChr",j,"_alphatable_fullreads.txt",sep=""),header=F)
TableC = read.table(paste(args1,"_DD_",args3,"_",args[3],"panel_1ind_noChr",j,"_alphatable_fullreads.txt",sep=""),header=F)

## Get log likelihoods of A(alpha+k) table (UU pairs)
BigTable[,5]=TableA[match(as.numeric(BigTable[,3]),as.numeric(TableA[,1])),2]
## Get log likelihoods of B(alpha/2+k) table (DU pairs)
BigTable[,6]=TableB[match(as.numeric(BigTable[,4]),as.numeric(TableB[,1])),2]
## Get log likelihoods of C(k) table (DD pairs)
BigTable[,7]=TableC[match(as.numeric(BigTable[,2]),as.numeric(TableC[,1])),2]
BigTable[,8]=BigTable[,5]+BigTable[,6]+BigTable[,7]

### Get new, more detailed Table based on broad estimates from above.
fun1=seq((BigTable[which.max(BigTable[,8]),1]-0.02),(BigTable[which.max(BigTable[,8]),1]+0.02),by=0.0001)
fun3=seq((BigTable[which.max(BigTable[,8]),2]-0.02),(BigTable[which.max(BigTable[,8]),2]+0.02),by=0.0001)
alpha_2 = rep(fun1,401)
k23_2 =rep(fun3,each=401)

BigTable2[,1]=alpha_2
BigTable2[,2]=k23_2

## All possible alpha+k values
BigTable2[,3]=round(BigTable2[,1]+BigTable2[,2],digits=4)
## All possible alpha/2+k values
BigTable2[,4]=round(BigTable2[,1]/2+BigTable2[,2],digits=4)
## Get log likelihoods of A(alpha+k) table (UU pairs)
BigTable2[,5]=TableA[match(as.numeric(BigTable2[,3]),as.numeric(TableA[,1])),2]
## Get log likelihoods of B(alpha/2+k) table (DU pairs)
BigTable2[,6]=TableB[match(as.numeric(BigTable2[,4]),as.numeric(TableB[,1])),2]
## Get log likelihoods of C(k) table (DD pairs)
BigTable2[,7]=TableC[match(as.numeric(BigTable2[,2]),as.numeric(TableC[,1])),2]
BigTable2[,8]=BigTable2[,5]+BigTable2[,6]+BigTable2[,7]


#### Store the values in jinfile.
jinfile[j,3]=BigTable2[which.max(BigTable2[,8]),1]*d
jinfile2[j,3]=BigTable2[which.max(BigTable2[,8]),2]

for(j in 2:22){
TableA = read.table(paste(args1,"_UU_",args3,"_",args[3],"panel_1ind_noChr",j,"_alphatable_fullreads.txt",sep=""),header=F)
TableB = read.table(paste(args1,"_DU_",args3,"_",args[3],"panel_1ind_noChr",j,"_alphatable_fullreads.txt",sep=""),header=F)
TableC = read.table(paste(args1,"_DD_",args3,"_",args[3],"panel_1ind_noChr",j,"_alphatable_fullreads.txt",sep=""),header=F)

## Get log likelihoods of A(alpha+k) table (UU pairs)
BigTable2[,5]=TableA[match(as.numeric(BigTable2[,3]),as.numeric(TableA[,1])),2]
## Get log likelihoods of B(alpha/2+k) table (DU pairs)
BigTable2[,6]=TableB[match(as.numeric(BigTable2[,4]),as.numeric(TableB[,1])),2]
## Get log likelihoods of C(k) table (DD pairs)
BigTable2[,7]=TableC[match(as.numeric(BigTable2[,2]),as.numeric(TableC[,1])),2]
BigTable2[,8]=BigTable2[,5]+BigTable2[,6]+BigTable2[,7]

#### Store the values in jinfile.
jinfile[j,3]=BigTable2[which.max(BigTable2[,8]),1]*d
jinfile2[j,3]=BigTable2[which.max(BigTable2[,8]),2]
}
## make jin files for the particular contamination level.
writejin <- file(paste(args1,"_",args3,"_",args[3],"panel_scores_reads_DamCorr.jin",sep=""),"w")
write.table(jinfile,writejin,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
flush(writejin)
close(writejin)
writejin2 <- file(paste(args1,"_",args3,"_",args[3],"panel_kscores_reads_DamCorr.jin",sep=""),"w")
write.table(jinfile2,writejin2,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
flush(writejin2)
close(writejin2)
num=mean(jinfile[,3])
writenum <- file(paste(args1,"_",args3,"_",args[3],"panel_scores_reads_DamCorr_meanscore2.txt",sep=""),"w")
write.table(num,writenum,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
flush(writenum)
close(writenum)
num=mean(jinfile2[,3])
writenum2 <- file(paste(args1,"_",args3,"_",args[3],"panel_kscores_reads_DamCorr_meanscore2.txt",sep=""),"w")
write.table(num,writenum2,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
flush(writenum2)
close(writenum2)

### Get meanscores
j=0
TableA = read.table(paste(args1,"_UU_",args3,"_",args[3],"panel_1ind_noChr",j,"_alphatable_fullreads.txt",sep=""),header=F)
TableB = read.table(paste(args1,"_DU_",args3,"_",args[3],"panel_1ind_noChr",j,"_alphatable_fullreads.txt",sep=""),header=F)
TableC = read.table(paste(args1,"_DD_",args3,"_",args[3],"panel_1ind_noChr",j,"_alphatable_fullreads.txt",sep=""),header=F)

## Get log likelihoods of A(alpha+k) table (UU pairs)
BigTable2[,5]=TableA[match(as.numeric(BigTable2[,3]),as.numeric(TableA[,1])),2]
## Get log likelihoods of B(alpha/2+k) table (DU pairs)
BigTable2[,6]=TableB[match(as.numeric(BigTable2[,4]),as.numeric(TableB[,1])),2]
## Get log likelihoods of C(k) table (DD pairs)
BigTable2[,7]=TableC[match(as.numeric(BigTable2[,2]),as.numeric(TableC[,1])),2]
BigTable2[,8]=BigTable2[,5]+BigTable2[,6]+BigTable2[,7]

alphascore=BigTable2[which.max(BigTable2[,8]),1]
kscore=BigTable2[which.max(BigTable2[,8]),2]

writealpha <- file(paste(args1,"_",args3,"_",args[3],"panel_scores_reads_DamCorr_meanscore.txt",sep=""),"w")
write.table(alphascore,writealpha,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
flush(writealpha)
close(writealpha)
writek <- file(paste(args1,"_",args3,"_",args[3],"panel_kscores_reads_DamCorr_meanscore.txt",sep=""),"w")
write.table(kscore,writek,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
flush(writek)
close(writek)
}
