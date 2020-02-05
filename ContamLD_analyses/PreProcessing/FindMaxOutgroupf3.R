#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(scipen=999)

setwd(args[1])
test = read.table(args[2],header=F)
test2 = read.table(args[3],header=F)
Table=data.frame(1:nrow(test2))
Table[,2]=0
for(i in 0:(nrow(test2)-1)){
temp=test[(1+i*26):((i+1)*26),]
Table[(i+1),1]=as.character(temp[which.max(temp[,5]),3])
Table[(i+1),2]=as.character(temp[which.max(temp[,5]),2])
Table[(i+1),3]=as.character(temp[which.max(temp[,5]),5])
}
Table2=Table
Table2[,1]=test2[,1]
Table2[,3]=NULL
Table2$V2 <- substr(Table2$V2, 0, 3)
write.table(Table2,file=args[4],sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
