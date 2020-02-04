#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(scipen=999)

setwd(paste(args[1],sep=""))
LDfile = read.table(paste(args[2],".ld",sep=""),header=T)
snpfile=read.table(paste(args[2],".snp",sep=""),header=F)
fileName <- paste(args[2],".eigenstratgeno",sep="")
odd=readChar(fileName, file.info(fileName)$size)
odd2=substr(odd[1],1,nchar(odd[1])-1)
a2 <- unlist(strsplit(odd2, "\n"))
mm <- do.call(rbind, strsplit(a2, ""))
dd<-data.frame(mm)


countfreq <- function(dataframerow){
temp=as.data.frame(table(t(dataframerow)))
## Count number of 0s
if(nrow(temp[temp[,1]==0,]) !=0){
	count0=temp[temp[,1]==0,][2]
	count0=count0[1,1]}
if(nrow(temp[temp[,1]==0,]) ==0){
	count0=0}
## Count number of 1s
if(nrow(temp[temp[,1]==1,]) !=0){
	count1=temp[temp[,1]==1,][2]
	count1=count1[1,1]}
if(nrow(temp[temp[,1]==1,]) ==0){
	count1=0}
## Count number of 2s
if(nrow(temp[temp[,1]==2,]) !=0){
	count2=temp[temp[,1]==2,][2]
	count2=count2[1,1]}
if(nrow(temp[temp[,1]==2,]) ==0){
	count2=0}

return(c(count0,count1,count2))
}

test=apply(dd,1,countfreq)
SNPfreq=data.frame(t(test))

odd=rowSums(SNPfreq)
SNPfreq$total=odd

SNPfreq$freq0 = SNPfreq[,1]/SNPfreq$total
SNPfreq$freq1 = SNPfreq[,2]/SNPfreq$total
SNPfreq$freq2 = SNPfreq[,2]/SNPfreq$total


### match on rsid, then plug in value from other dataframe (dd); sum across 
### 

LDfile[,1]=match(LDfile[,3],snpfile[,1])
LDfile[,4]=match(LDfile[,6],snpfile[,1])


## Given a dataframerow of locations in genotype file (positions in the snp file),
## this function returns a list of haplotypes ("00","02", etc.)
gethaplotypes <- function(dataframerow,dd){
test=unlist(dd[dataframerow[1],])
test2=unlist(dd[dataframerow[2],])
return(paste0(test,test2))
}


## Get all the SNP file locations of the SNPs in high LD
new=data.frame(LDfile[,1],LDfile[,4])
### Get all haplotypes 
test=apply(new,1,gethaplotypes,dd=dd)
Haplotypes=data.frame(t(test))


#### Find frequency of Haplotypes in snp file

### Function to count the haplotypes in a file after attaining them.
counthaplotypes <- function(dataframerow){
temp=as.data.frame(table(t(dataframerow)))
## Count number of 00s
if(nrow(temp[temp[,1]=="00",]) !=0){
	count00=temp[temp[,1]=="00",][2]
	count00=count00[1,1]}
if(nrow(temp[temp[,1]=="00",]) ==0){
	count00=0}
## Count number of 01s
if(nrow(temp[temp[,1]=="01",]) !=0){
	count01=temp[temp[,1]=="01",][2]
	count01=count01[1,1]}
if(nrow(temp[temp[,1]=="01",]) ==0){
	count01=0}
## Count number of 02s
if(nrow(temp[temp[,1]=="02",]) !=0){
	count02=temp[temp[,1]=="02",][2]
	count02=count02[1,1]}
if(nrow(temp[temp[,1]=="02",]) ==0){
	count02=0}
## Count number of 10s
if(nrow(temp[temp[,1]=="10",]) !=0){
	count10=temp[temp[,1]=="10",][2]
	count10=count10[1,1]}
if(nrow(temp[temp[,1]=="10",]) ==0){
	count10=0}
## Count number of 11s
if(nrow(temp[temp[,1]=="11",]) !=0){
	count11=temp[temp[,1]=="11",][2]
	count11=count11[1,1]}
if(nrow(temp[temp[,1]=="11",]) ==0){
	count11=0}
## Count number of 12s
if(nrow(temp[temp[,1]=="12",]) !=0){
	count12=temp[temp[,1]=="12",][2]
	count12=count12[1,1]}
if(nrow(temp[temp[,1]=="12",]) ==0){
	count12=0}
## Count number of 20s
if(nrow(temp[temp[,1]=="20",]) !=0){
	count20=temp[temp[,1]=="20",][2]
	count20=count20[1,1]}
if(nrow(temp[temp[,1]=="20",]) ==0){
	count20=0}
## Count number of 21s
if(nrow(temp[temp[,1]=="21",]) !=0){
	count21=temp[temp[,1]=="21",][2]
	count21=count21[1,1]}
if(nrow(temp[temp[,1]=="21",]) ==0){
	count21=0}
## Count number of 22s
if(nrow(temp[temp[,1]=="22",]) !=0){
	count22=temp[temp[,1]=="22",][2]
	count22=count22[1,1]}
if(nrow(temp[temp[,1]=="22",]) ==0){
	count22=0}
return(c(count00,count01,count02,count10,count11,count12,count20,count21,count22))
}

test2=apply(Haplotypes,1,counthaplotypes)
Haplotypefreq=data.frame(t(test2))

odd=rowSums(Haplotypefreq)
Haplotypefreq$total=odd

Haplotypefreq$freq00 = Haplotypefreq[,1]/Haplotypefreq$total
Haplotypefreq$freq01 = Haplotypefreq[,2]/Haplotypefreq$total
Haplotypefreq$freq02 = Haplotypefreq[,3]/Haplotypefreq$total
Haplotypefreq$freq10 = Haplotypefreq[,4]/Haplotypefreq$total
Haplotypefreq$freq11 = Haplotypefreq[,5]/Haplotypefreq$total
Haplotypefreq$freq12 = Haplotypefreq[,6]/Haplotypefreq$total
Haplotypefreq$freq20 = Haplotypefreq[,7]/Haplotypefreq$total
Haplotypefreq$freq21 = Haplotypefreq[,8]/Haplotypefreq$total
Haplotypefreq$freq22 = Haplotypefreq[,9]/Haplotypefreq$total

### Calculate expected frequencies
SNPfreq$expected00=SNPfreq$freq0*SNPfreq$freq0
SNPfreq$expected01=SNPfreq$freq0*SNPfreq$freq1
SNPfreq$expected02=SNPfreq$freq0*SNPfreq$freq2
SNPfreq$expected10=SNPfreq$freq1*SNPfreq$freq0
SNPfreq$expected11=SNPfreq$freq1*SNPfreq$freq1
SNPfreq$expected12=SNPfreq$freq1*SNPfreq$freq2
SNPfreq$expected20=SNPfreq$freq2*SNPfreq$freq0
SNPfreq$expected21=SNPfreq$freq2*SNPfreq$freq1
SNPfreq$expected22=SNPfreq$freq2*SNPfreq$freq2

chrs=c(1:22)
chrs=as.character(chrs)
LDfile$CHR_A = as.character(LDfile$CHR_A)

write.table(SNPfreq,file=paste(args[2],"_",args[3]"_snpfreqs_heldout.txt",sep=""),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
write.table(Haplotypefreq,file=paste(args[2],"_",args[3]"_haplotypefreqs_heldout.txt",sep=""),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)

## Remove anything except chromosomes 1-22 from .geno, snpfreq, and .snp files
SNPfreq_no23=SNPfreq[(snpfile$V2 %in% chrs),]
write.table(SNPfreq_no23,file=paste(args[2],"_",args[3]"_snpfreqs_heldout_noX.txt",sep=""),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)

## Remove Chromosome 23 from LD and Haplotypefreq files
LDfile_no23 = LDfile[LDfile$CHR_A %in% chrs,]
write.table(LDfile_no23,file=paste(args[2],"_",args[3]"_eig_heldout_noX.ld",sep=""),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)

Haplotypefreq_no23=Haplotypefreq[(LDfile$CHR_A %in% chrs),]
write.table(Haplotypefreq_no23,file=paste(args[2],"_",args[3]"_haplotypefreqs_heldout_noX.txt",sep=""),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)

snpfile_no23 = snpfile[(snpfile$V2 %in% chrs),]
write.table(snpfile_no23,file=paste(args[2],"_",args[3]"_eig_heldout_noX.snp",sep=""),sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
