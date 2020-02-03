#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(scipen=999)

setwd(paste(args[2],"/directories/",args[3],"_jackknife",sep=""))

### snpfile must be in Chr, pos, readdepth#1, readdepth#2 format. (All Hg19)
snpfile=read.table(paste(args[3],"_reads_107k.txt",sep=""),header=F)

SNPfreq_panel=read.table(gzfile(paste(args[1],"/panels/",args[4],"_snpfreqs_heldout_noX.txt.gz",sep="")),header=T,stringsAsFactors=F)
### These are the positions used to make the SNPfreq_panel and HaplotypeFreq_panel files.
SNPfreq_panel_snpfile=read.table(gzfile(paste(args[1],"/panels/",args[4],"_eig_heldout_noX_ChromPos.snp.gz",sep="")),header=F,stringsAsFactors=F)

#### Intersect the SNPfreq_snpfile and snpfile to only get the places where they match. (if the input snpfile is different from the file used to make the SNPfreq_panel file).
### Concatenate chromosome and position for intersection.
snpfile$V5=paste(snpfile$V1,"-",snpfile$V2,sep="")
### Narrow down the SNPfreq panel and SNPfreq panel snpfile to the places where there exists the same snp in the snpfile.
SNPfreq_panel=SNPfreq_panel[(SNPfreq_panel_snpfile[,1] %in% snpfile[,5]),]
SNPfreq_panel_snpfile=SNPfreq_panel_snpfile[(SNPfreq_panel_snpfile[,1] %in% snpfile[,5]),]
SNPfreq_panel_snpfile=data.frame(SNPfreq_panel_snpfile)

### Shrink the snpfile down to the places where there exists a SNP in the SNPfreq file.
snpfile=snpfile[(snpfile[,5] %in% SNPfreq_panel_snpfile[,1]),]

### Find the matches of SNP pairs of snpfile in LDfile (will be used to search in HaplotypeFreq_panel also)
LDfile = read.table(gzfile(paste(args[1],"/panels/",args[4],"_eig_heldout_noX.ld.gz",sep="")),header=T)
### Keep the chromosome number for later (jackknife resampling)
LDfile[,8]=LDfile[,1]
### Concatenate chromosome and position for intersection.
LDfile[,9]=paste(LDfile[,1],"-",LDfile[,2],sep="") 
LDfile[,10]=paste(LDfile[,4],"-",LDfile[,5],sep="") 

HaplotypeFreq_panel=read.table(gzfile(paste(args[1],"/panels/",args[4],"_haplotypefreqs_heldout_noX.txt.gz",sep="")),header=T,stringsAsFactors=F,colClasses="character")

### Goal:  Narrow down LDfile and HaplotypeFreq_panel to the places where it matches snpfile
HaplotypeFreq_panel=HaplotypeFreq_panel[((LDfile[,9] %in% snpfile[,5]) & (LDfile[,10] %in% snpfile[,5])),]
LDfile=LDfile[((LDfile[,9] %in% snpfile[,5]) & (LDfile[,10] %in% snpfile[,5])),]

LDfile[,1]=match(LDfile[,9],snpfile[,5])
LDfile[,4]=match(LDfile[,10],snpfile[,5])

for(i in 1:ncol(HaplotypeFreq_panel)){
HaplotypeFreq_panel[,i]=as.numeric(HaplotypeFreq_panel[,i])
}
for(i in 1:ncol(SNPfreq_panel)){
SNPfreq_panel[,i]=as.numeric(SNPfreq_panel[,i])
}

### Add a little to make sure no 0s.
SNPfreq_panel[,1]=SNPfreq_panel[,1]+0.000000001
SNPfreq_panel[,2]=SNPfreq_panel[,2]+0.000000001
SNPfreq_panel[,3]=SNPfreq_panel[,3]+0.000000001
SNPfreq_panel[,4]=SNPfreq_panel[,4]+0.000000003
SNPfreq_panel[,5]=SNPfreq_panel[,1]/SNPfreq_panel[,4]
SNPfreq_panel[,6]=SNPfreq_panel[,2]/SNPfreq_panel[,4]
SNPfreq_panel[,7]=SNPfreq_panel[,3]/SNPfreq_panel[,4]

HaplotypeFreq_panel[,1]=HaplotypeFreq_panel[,1]+0.000000001
HaplotypeFreq_panel[,2]=HaplotypeFreq_panel[,2]+0.000000001
HaplotypeFreq_panel[,3]=HaplotypeFreq_panel[,3]+0.000000001
HaplotypeFreq_panel[,4]=HaplotypeFreq_panel[,4]+0.000000001
HaplotypeFreq_panel[,5]=HaplotypeFreq_panel[,5]+0.000000001
HaplotypeFreq_panel[,6]=HaplotypeFreq_panel[,6]+0.000000001
HaplotypeFreq_panel[,7]=HaplotypeFreq_panel[,7]+0.000000001
HaplotypeFreq_panel[,8]=HaplotypeFreq_panel[,8]+0.000000001
HaplotypeFreq_panel[,9]=HaplotypeFreq_panel[,9]+0.000000001
HaplotypeFreq_panel[,10]=HaplotypeFreq_panel[,10]+0.000000009

HaplotypeFreq_panel$freq00 = HaplotypeFreq_panel[,1]/HaplotypeFreq_panel$total
HaplotypeFreq_panel$freq01 = HaplotypeFreq_panel[,2]/HaplotypeFreq_panel$total
HaplotypeFreq_panel$freq02 = HaplotypeFreq_panel[,3]/HaplotypeFreq_panel$total
HaplotypeFreq_panel$freq10 = HaplotypeFreq_panel[,4]/HaplotypeFreq_panel$total
HaplotypeFreq_panel$freq11 = HaplotypeFreq_panel[,5]/HaplotypeFreq_panel$total
HaplotypeFreq_panel$freq12 = HaplotypeFreq_panel[,6]/HaplotypeFreq_panel$total
HaplotypeFreq_panel$freq20 = HaplotypeFreq_panel[,7]/HaplotypeFreq_panel$total
HaplotypeFreq_panel$freq21 = HaplotypeFreq_panel[,8]/HaplotypeFreq_panel$total
HaplotypeFreq_panel$freq22 = HaplotypeFreq_panel[,9]/HaplotypeFreq_panel$total

### Get just frequency of SNPs (0,1,2)
SNPfreq_summarized=data.frame(1:nrow(SNPfreq_panel))
SNPfreq_summarized[,1]=SNPfreq_panel$freq0
SNPfreq_summarized[,2]=SNPfreq_panel$freq1
SNPfreq_summarized[,3]=SNPfreq_panel$freq2

## Get all the SNP file locations of the SNPs in high LD
new=data.frame(LDfile[,1],LDfile[,4])

### Goal:  Get frequencies of all SNPs at positions where haplotype frequencies are calculated (i.e. SNP LD is high)
freq0col1=SNPfreq_summarized[new[,1],]
freq0col2=SNPfreq_summarized[new[,2],]
HaplotypeFreq_panel$exp00=freq0col1[,1]*freq0col2[,1]
HaplotypeFreq_panel$exp01=freq0col1[,1]*freq0col2[,2]
HaplotypeFreq_panel$exp02=freq0col1[,1]*freq0col2[,3]
HaplotypeFreq_panel$exp10=freq0col1[,2]*freq0col2[,1]
HaplotypeFreq_panel$exp11=freq0col1[,2]*freq0col2[,2]
HaplotypeFreq_panel$exp12=freq0col1[,2]*freq0col2[,3]
HaplotypeFreq_panel$exp20=freq0col1[,3]*freq0col2[,1]
HaplotypeFreq_panel$exp21=freq0col1[,3]*freq0col2[,2]
HaplotypeFreq_panel$exp22=freq0col1[,3]*freq0col2[,3]

## Get counts of all pair-wise combinations of possible haplotypes at each SNP pair (based on number of reads).
odd=snpfile[new[,1],]
odd2=snpfile[new[,2],]
Haplotypes_testindiv=data.frame(1:nrow(new))
## 1st column = # of 00 haplotypes
Haplotypes_testindiv[,1]=odd[,4]*odd2[,4]
## 2nd column = # of 02 haplotypes
Haplotypes_testindiv[,2]=odd[,4]*odd2[,3]
## 3rd column = # of 20 haplotypes
Haplotypes_testindiv[,3]=odd[,3]*odd2[,4]
## 4th column = # of 22 haplotypes
Haplotypes_testindiv[,4]=odd[,3]*odd2[,3]

HaplotypeFreq_panel_testindiv=HaplotypeFreq_panel

## Remove all places where there are no haplotypes.
HaplotypeFreq_panel_testindiv=HaplotypeFreq_panel_testindiv[!(Haplotypes_testindiv[,1]==0 & Haplotypes_testindiv[,2]==0 & Haplotypes_testindiv[,3]==0 & Haplotypes_testindiv[,4]==0),]
LDfile=LDfile[!(Haplotypes_testindiv[,1]==0 & Haplotypes_testindiv[,2]==0 & Haplotypes_testindiv[,3]==0 & Haplotypes_testindiv[,4]==0),]
Haplotypes_testindiv=Haplotypes_testindiv[!(Haplotypes_testindiv[,1]==0 & Haplotypes_testindiv[,2]==0 & Haplotypes_testindiv[,3]==0 & Haplotypes_testindiv[,4]==0),]


########## EM algorithm
### Set initial haplotype frequency estimates (00, 01, 10, 11)
HaploEstimate=data.frame(1:nrow(Haplotypes_testindiv))
HaploEstimate[,1]=0.25
HaploEstimate[,2]=0.25
HaploEstimate[,3]=0.25
HaploEstimate[,4]=0.25

### Calculate P(i,j|h). (this is the probability distribution of diploid counts; i.e. probability of getting a particular count 00, 02, etc.)
HaploEstimate$freq00=HaploEstimate[,1]*HaploEstimate[,1]
HaploEstimate$freq01=HaploEstimate[,1]*HaploEstimate[,2]+HaploEstimate[,2]*HaploEstimate[,1]
HaploEstimate$freq02=HaploEstimate[,2]*HaploEstimate[,2]
HaploEstimate$freq10=HaploEstimate[,1]*HaploEstimate[,3]+HaploEstimate[,3]*HaploEstimate[,1]
HaploEstimate$freq11=HaploEstimate[,2]*HaploEstimate[,3]+HaploEstimate[,3]*HaploEstimate[,2]+HaploEstimate[,1]*HaploEstimate[,4]+HaploEstimate[,4]*HaploEstimate[,1]
HaploEstimate$freq12=HaploEstimate[,2]*HaploEstimate[,4]+HaploEstimate[,4]*HaploEstimate[,2]
HaploEstimate$freq20=HaploEstimate[,3]*HaploEstimate[,3]
HaploEstimate$freq21=HaploEstimate[,3]*HaploEstimate[,4]+HaploEstimate[,4]*HaploEstimate[,3]
HaploEstimate$freq22=HaploEstimate[,4]*HaploEstimate[,4]

## Calculate n, the expected number of times the (a,b) configuration from one parent is contributing to the diploid count.
n=data.frame(1:nrow(HaploEstimate))
n[,1]=HaplotypeFreq_panel_testindiv$freq00*HaploEstimate[,1]*HaploEstimate[,1]/HaploEstimate$freq00 ## (0,0) is contributing to (0,0)
n[,2]=HaplotypeFreq_panel_testindiv$freq01*HaploEstimate[,1]*HaploEstimate[,2]/HaploEstimate$freq01 ## (0,0) is contributing to (0,1)
n[,3]=0
n[,4]=HaplotypeFreq_panel_testindiv$freq10*HaploEstimate[,1]*HaploEstimate[,3]/HaploEstimate$freq10 ## (0,0) is contributing to (1,0)
n[,5]=HaplotypeFreq_panel_testindiv$freq11*HaploEstimate[,1]*HaploEstimate[,4]/HaploEstimate$freq11 ## (0,0) is contributing to (1,1)
n[,6]=0
n[,7]=0
n[,8]=0
n[,9]=0
n[,10]=0
n[,11]=HaplotypeFreq_panel_testindiv$freq01*HaploEstimate[,2]*HaploEstimate[,1]/HaploEstimate$freq01 ## (0,1) is contributing to (0,1)
n[,12]=HaplotypeFreq_panel_testindiv$freq02*HaploEstimate[,2]*HaploEstimate[,2]/HaploEstimate$freq02 ## (0,1) is contributing to (0,2)
n[,13]=0
n[,14]=HaplotypeFreq_panel_testindiv$freq11*HaploEstimate[,2]*HaploEstimate[,3]/HaploEstimate$freq11 ## (0,1) is contributing to (1,1)
n[,15]=HaplotypeFreq_panel_testindiv$freq12*HaploEstimate[,2]*HaploEstimate[,4]/HaploEstimate$freq12 ## (0,1) is contributing to (1,2)
n[,16]=0
n[,17]=0
n[,18]=0
n[,19]=0
n[,20]=0
n[,21]=0
n[,22]=HaplotypeFreq_panel_testindiv$freq10*HaploEstimate[,3]*HaploEstimate[,1]/HaploEstimate$freq10 ## (1,0) is contributing to (1,0)
n[,23]=HaplotypeFreq_panel_testindiv$freq11*HaploEstimate[,3]*HaploEstimate[,2]/HaploEstimate$freq11 ## (1,0) is contributing to (1,1)
n[,24]=0
n[,25]=HaplotypeFreq_panel_testindiv$freq20*HaploEstimate[,3]*HaploEstimate[,3]/HaploEstimate$freq20 ## (1,0) is contributing to (2,0)
n[,26]=HaplotypeFreq_panel_testindiv$freq21*HaploEstimate[,3]*HaploEstimate[,4]/HaploEstimate$freq21 ## (1,0) is contributing to (2,1)
n[,27]=0
n[,28]=0
n[,29]=0
n[,30]=0
n[,31]=0
n[,32]=HaplotypeFreq_panel_testindiv$freq11*HaploEstimate[,4]*HaploEstimate[,1]/HaploEstimate$freq11 ## (1,1) is contributing to (1,1)
n[,33]=HaplotypeFreq_panel_testindiv$freq12*HaploEstimate[,4]*HaploEstimate[,2]/HaploEstimate$freq12 ## (1,1) is contributing to (1,2)
n[,34]=0
n[,35]=HaplotypeFreq_panel_testindiv$freq21*HaploEstimate[,4]*HaploEstimate[,3]/HaploEstimate$freq21 ## (1,1) is contributing to (2,1)
n[,36]=HaplotypeFreq_panel_testindiv$freq22*HaploEstimate[,4]*HaploEstimate[,4]/HaploEstimate$freq22 ## (1,1) is contributing to (2,2)

### Maximize the expected value of haplotype.  For each haplotype configuration, sum up the probability of that configuration over all diploid configurations.
d=data.frame(1:nrow(HaploEstimate))
d[,1]=(n[,1]+n[,1])*HaplotypeFreq_panel_testindiv$freq00+(n[,2]+n[,2])*HaplotypeFreq_panel_testindiv$freq01+(n[,4]+n[,4])*HaplotypeFreq_panel_testindiv$freq10+(n[,5]+n[,5])*HaplotypeFreq_panel_testindiv$freq11
d[,2]=(n[,11]+n[,11])*HaplotypeFreq_panel_testindiv$freq01+(n[,12]+n[,12])*HaplotypeFreq_panel_testindiv$freq02+(n[,14]+n[,14])*HaplotypeFreq_panel_testindiv$freq11+(n[,15]+n[,15])*HaplotypeFreq_panel_testindiv$freq12
d[,3]=(n[,22]+n[,22])*HaplotypeFreq_panel_testindiv$freq10+(n[,23]+n[,23])*HaplotypeFreq_panel_testindiv$freq11+(n[,25]+n[,25])*HaplotypeFreq_panel_testindiv$freq20+(n[,26]+n[,26])*HaplotypeFreq_panel_testindiv$freq21
d[,4]=(n[,32]+n[,32])*HaplotypeFreq_panel_testindiv$freq11+(n[,33]+n[,33])*HaplotypeFreq_panel_testindiv$freq12+(n[,35]+n[,35])*HaplotypeFreq_panel_testindiv$freq21+(n[,36]+n[,36])*HaplotypeFreq_panel_testindiv$freq22

### Get the haplotype distributions.
h=data.frame(1:nrow(HaploEstimate))
h[,1]=d[,1]/(d[,1]+d[,2]+d[,3]+d[,4])
h[,2]=d[,2]/(d[,1]+d[,2]+d[,3]+d[,4])
h[,3]=d[,3]/(d[,1]+d[,2]+d[,3]+d[,4])
h[,4]=d[,4]/(d[,1]+d[,2]+d[,3]+d[,4])

total1=h[,1]-HaploEstimate[,1]
total2=h[,2]-HaploEstimate[,2]
total3=h[,3]-HaploEstimate[,3]
total4=h[,4]-HaploEstimate[,4]
total5=(sum(total1))^2+(sum(total2))^2+(sum(total3))^2+(sum(total4))^2

##### Repeat the process until convergence.
while(total5>.001){
HaploEstimate[,1]=h[,1]
HaploEstimate[,2]=h[,2]
HaploEstimate[,3]=h[,3]
HaploEstimate[,4]=h[,4]

### Calculate P(i,j|h). (this is the distribution of diploid counts)
HaploEstimate$freq00=HaploEstimate[,1]*HaploEstimate[,1]
HaploEstimate$freq01=HaploEstimate[,1]*HaploEstimate[,2]+HaploEstimate[,2]*HaploEstimate[,1]
HaploEstimate$freq02=HaploEstimate[,2]*HaploEstimate[,2]
HaploEstimate$freq10=HaploEstimate[,1]*HaploEstimate[,3]+HaploEstimate[,3]*HaploEstimate[,1]
HaploEstimate$freq11=HaploEstimate[,2]*HaploEstimate[,3]+HaploEstimate[,3]*HaploEstimate[,2]+HaploEstimate[,1]*HaploEstimate[,4]+HaploEstimate[,4]*HaploEstimate[,1]
HaploEstimate$freq12=HaploEstimate[,2]*HaploEstimate[,4]+HaploEstimate[,4]*HaploEstimate[,2]
HaploEstimate$freq20=HaploEstimate[,3]*HaploEstimate[,3]
HaploEstimate$freq21=HaploEstimate[,3]*HaploEstimate[,4]+HaploEstimate[,4]*HaploEstimate[,3]
HaploEstimate$freq22=HaploEstimate[,4]*HaploEstimate[,4]

## Calculate n, the expected number of times the (a,b) configuration is contributing to the diploid count.
n=data.frame(1:nrow(HaploEstimate))
n[,1]=HaplotypeFreq_panel_testindiv$freq00*HaploEstimate[,1]*HaploEstimate[,1]/HaploEstimate$freq00 ## (0,0) is contributing to (0,0)
n[,2]=HaplotypeFreq_panel_testindiv$freq01*HaploEstimate[,1]*HaploEstimate[,2]/HaploEstimate$freq01 ## (0,0) is contributing to (0,1)
n[,3]=0
n[,4]=HaplotypeFreq_panel_testindiv$freq10*HaploEstimate[,1]*HaploEstimate[,3]/HaploEstimate$freq10 ## (0,0) is contributing to (1,0)
n[,5]=HaplotypeFreq_panel_testindiv$freq11*HaploEstimate[,1]*HaploEstimate[,4]/HaploEstimate$freq11 ## (0,0) is contributing to (1,1)
n[,6]=0
n[,7]=0
n[,8]=0
n[,9]=0
n[,10]=0
n[,11]=HaplotypeFreq_panel_testindiv$freq01*HaploEstimate[,2]*HaploEstimate[,1]/HaploEstimate$freq01 ## (0,1) is contributing to (0,1)
n[,12]=HaplotypeFreq_panel_testindiv$freq02*HaploEstimate[,2]*HaploEstimate[,2]/HaploEstimate$freq02 ## (0,1) is contributing to (0,2)
n[,13]=0
n[,14]=HaplotypeFreq_panel_testindiv$freq11*HaploEstimate[,2]*HaploEstimate[,3]/HaploEstimate$freq11 ## (0,1) is contributing to (1,1)
n[,15]=HaplotypeFreq_panel_testindiv$freq12*HaploEstimate[,2]*HaploEstimate[,4]/HaploEstimate$freq12 ## (0,1) is contributing to (1,2)
n[,16]=0
n[,17]=0
n[,18]=0
n[,19]=0
n[,20]=0
n[,21]=0
n[,22]=HaplotypeFreq_panel_testindiv$freq10*HaploEstimate[,3]*HaploEstimate[,1]/HaploEstimate$freq10 ## (1,0) is contributing to (1,0)
n[,23]=HaplotypeFreq_panel_testindiv$freq11*HaploEstimate[,3]*HaploEstimate[,2]/HaploEstimate$freq11 ## (1,0) is contributing to (1,1)
n[,24]=0
n[,25]=HaplotypeFreq_panel_testindiv$freq20*HaploEstimate[,3]*HaploEstimate[,3]/HaploEstimate$freq20 ## (1,0) is contributing to (2,0)
n[,26]=HaplotypeFreq_panel_testindiv$freq21*HaploEstimate[,3]*HaploEstimate[,4]/HaploEstimate$freq21 ## (1,0) is contributing to (2,1)
n[,27]=0
n[,28]=0
n[,29]=0
n[,30]=0
n[,31]=0
n[,32]=HaplotypeFreq_panel_testindiv$freq11*HaploEstimate[,4]*HaploEstimate[,1]/HaploEstimate$freq11 ## (1,1) is contributing to (1,1)
n[,33]=HaplotypeFreq_panel_testindiv$freq12*HaploEstimate[,4]*HaploEstimate[,2]/HaploEstimate$freq12 ## (1,1) is contributing to (1,2)
n[,34]=0
n[,35]=HaplotypeFreq_panel_testindiv$freq21*HaploEstimate[,4]*HaploEstimate[,3]/HaploEstimate$freq21 ## (1,1) is contributing to (2,1)
n[,36]=HaplotypeFreq_panel_testindiv$freq22*HaploEstimate[,4]*HaploEstimate[,4]/HaploEstimate$freq22 ## (1,1) is contributing to (2,2)

### Maximize the expected value of haplotype.  For each haplotype configuration, sum up the probability of that configuration over all diploid configurations.
d=data.frame(1:nrow(HaploEstimate))
d[,1]=(n[,1]+n[,1])*HaplotypeFreq_panel_testindiv$freq00+(n[,2]+n[,2])*HaplotypeFreq_panel_testindiv$freq01+(n[,4]+n[,4])*HaplotypeFreq_panel_testindiv$freq10+(n[,5]+n[,5])*HaplotypeFreq_panel_testindiv$freq11
d[,2]=(n[,11]+n[,11])*HaplotypeFreq_panel_testindiv$freq01+(n[,12]+n[,12])*HaplotypeFreq_panel_testindiv$freq02+(n[,14]+n[,14])*HaplotypeFreq_panel_testindiv$freq11+(n[,15]+n[,15])*HaplotypeFreq_panel_testindiv$freq12
d[,3]=(n[,22]+n[,22])*HaplotypeFreq_panel_testindiv$freq10+(n[,23]+n[,23])*HaplotypeFreq_panel_testindiv$freq11+(n[,25]+n[,25])*HaplotypeFreq_panel_testindiv$freq20+(n[,26]+n[,26])*HaplotypeFreq_panel_testindiv$freq21
d[,4]=(n[,32]+n[,32])*HaplotypeFreq_panel_testindiv$freq11+(n[,33]+n[,33])*HaplotypeFreq_panel_testindiv$freq12+(n[,35]+n[,35])*HaplotypeFreq_panel_testindiv$freq21+(n[,36]+n[,36])*HaplotypeFreq_panel_testindiv$freq22

### Get the haplotype distributions.
h=data.frame(1:nrow(HaploEstimate))
h[,1]=d[,1]/(d[,1]+d[,2]+d[,3]+d[,4])
h[,2]=d[,2]/(d[,1]+d[,2]+d[,3]+d[,4])
h[,3]=d[,3]/(d[,1]+d[,2]+d[,3]+d[,4])
h[,4]=d[,4]/(d[,1]+d[,2]+d[,3]+d[,4])

total1=h[,1]-HaploEstimate[,1]
total2=h[,2]-HaploEstimate[,2]
total3=h[,3]-HaploEstimate[,3]
total4=h[,4]-HaploEstimate[,4]
total5=(sum(total1))^2+(sum(total2))^2+(sum(total3))^2+(sum(total4))^2
}


####### Solve for alpha:

#### First get p, the backround probability of unrelated "haplotypes".
p=data.frame(1:nrow(HaploEstimate))
## freq of 0 at first SNP
p[,1]=(2*HaplotypeFreq_panel_testindiv$freq00+2*HaplotypeFreq_panel_testindiv$freq01+2*HaplotypeFreq_panel_testindiv$freq02+HaplotypeFreq_panel_testindiv$freq10+HaplotypeFreq_panel_testindiv$freq11+HaplotypeFreq_panel_testindiv$freq12)/4
## freq of 1 at first SNP
p[,2]=(2*HaplotypeFreq_panel_testindiv$freq20+2*HaplotypeFreq_panel_testindiv$freq21+2*HaplotypeFreq_panel_testindiv$freq22+HaplotypeFreq_panel_testindiv$freq10+HaplotypeFreq_panel_testindiv$freq11+HaplotypeFreq_panel_testindiv$freq12)/4
## freq of 0 at 2nd SNP
p[,3]=(2*HaplotypeFreq_panel_testindiv$freq00+2*HaplotypeFreq_panel_testindiv$freq10+2*HaplotypeFreq_panel_testindiv$freq20+HaplotypeFreq_panel_testindiv$freq01+HaplotypeFreq_panel_testindiv$freq11+HaplotypeFreq_panel_testindiv$freq21)/4
## freq of 1 at 2nd SNP
p[,4]=(2*HaplotypeFreq_panel_testindiv$freq02+2*HaplotypeFreq_panel_testindiv$freq12+2*HaplotypeFreq_panel_testindiv$freq22+HaplotypeFreq_panel_testindiv$freq01+HaplotypeFreq_panel_testindiv$freq11+HaplotypeFreq_panel_testindiv$freq21)/4

HaplotypeFreq_panel_testindiv$exp00=p[,1]*p[,3] ##random probability of 00
HaplotypeFreq_panel_testindiv$exp01=p[,1]*p[,4] ##random probability of 01
HaplotypeFreq_panel_testindiv$exp10=p[,2]*p[,3] ##random probability of 10
HaplotypeFreq_panel_testindiv$exp11=p[,2]*p[,4] ##random probability of 11

### Normalize HaplotypeFreq_panel_testindiv$exp00, etc. values
p[,5]=HaplotypeFreq_panel_testindiv$exp00/(HaplotypeFreq_panel_testindiv$exp00+HaplotypeFreq_panel_testindiv$exp01+HaplotypeFreq_panel_testindiv$exp10+HaplotypeFreq_panel_testindiv$exp11)
p[,6]=HaplotypeFreq_panel_testindiv$exp01/(HaplotypeFreq_panel_testindiv$exp00+HaplotypeFreq_panel_testindiv$exp01+HaplotypeFreq_panel_testindiv$exp10+HaplotypeFreq_panel_testindiv$exp11)
p[,7]=HaplotypeFreq_panel_testindiv$exp10/(HaplotypeFreq_panel_testindiv$exp00+HaplotypeFreq_panel_testindiv$exp01+HaplotypeFreq_panel_testindiv$exp10+HaplotypeFreq_panel_testindiv$exp11)
p[,8]=HaplotypeFreq_panel_testindiv$exp11/(HaplotypeFreq_panel_testindiv$exp00+HaplotypeFreq_panel_testindiv$exp01+HaplotypeFreq_panel_testindiv$exp10+HaplotypeFreq_panel_testindiv$exp11)
HaplotypeFreq_panel_testindiv$exp00=p[,5]
HaplotypeFreq_panel_testindiv$exp01=p[,6]
HaplotypeFreq_panel_testindiv$exp10=p[,7]
HaplotypeFreq_panel_testindiv$exp11=p[,8]

### Calculate q, the expected haplotype distribution of the test individual.
HaplotypeFreq_panel_testindiv[,1] = h[,1]/2 + HaplotypeFreq_panel_testindiv$exp00/2
HaplotypeFreq_panel_testindiv[,2] = h[,2]/2 + HaplotypeFreq_panel_testindiv$exp01/2
HaplotypeFreq_panel_testindiv[,3] = h[,3]/2 + HaplotypeFreq_panel_testindiv$exp10/2
HaplotypeFreq_panel_testindiv[,4] = h[,4]/2 + HaplotypeFreq_panel_testindiv$exp11/2

### Save original files to regenerate them for jackknife.
Haplotypes_testindiv_full=Haplotypes_testindiv
HaplotypeFreq_panel_testindiv_full=HaplotypeFreq_panel_testindiv

### Make a function to sum up all LOD scores for a particular alpha (a) value.
getvalues <- function(a,Haplotypes_testindiv,HaplotypeFreq_panel_testindiv){
num1=sum(Haplotypes_testindiv[,1]*log10(((1-2*a+2*(a^2))*HaplotypeFreq_panel_testindiv[,1]+2*a*(1-a)*HaplotypeFreq_panel_testindiv$exp00)/HaplotypeFreq_panel_testindiv[,1]))
num2=sum(Haplotypes_testindiv[,2]*log10(((1-2*a+2*(a^2))*HaplotypeFreq_panel_testindiv[,2]+2*a*(1-a)*HaplotypeFreq_panel_testindiv$exp01)/HaplotypeFreq_panel_testindiv[,2]))
num3=sum(Haplotypes_testindiv[,3]*log10(((1-2*a+2*(a^2))*HaplotypeFreq_panel_testindiv[,3]+2*a*(1-a)*HaplotypeFreq_panel_testindiv$exp10)/HaplotypeFreq_panel_testindiv[,3]))
num4=sum(Haplotypes_testindiv[,4]*log10(((1-2*a+2*(a^2))*HaplotypeFreq_panel_testindiv[,4]+2*a*(1-a)*HaplotypeFreq_panel_testindiv$exp11)/HaplotypeFreq_panel_testindiv[,4]))
return(num1+num2+num3+num4)
}

## Table2 for damage ratio
Table3=data.frame(1:22)
Table3[,2]=0

#### Plot LOD scores to solve for a.
#### Make a dataframe with all alpha from 0 to 0.5 
Table=data.frame(0:6000)
Table[,1]=seq(-0.1,0.5,by=0.0001)

####### Get jackknife scores (loop through each chromosome).:
for(f in 1:22){
### Regenerate original file
Haplotypes_testindiv=Haplotypes_testindiv_full
HaplotypeFreq_panel_testindiv=HaplotypeFreq_panel_testindiv_full

### Calculate Haplotype frequencies for each chromosome
Haplotypes_testindiv=data.frame(Haplotypes_testindiv[LDfile[,8]==f,])
HaplotypeFreq_panel_testindiv=HaplotypeFreq_panel_testindiv[LDfile[,8]==f,]

### Get the sum of all readpairs for each chromosome.
Table3[f,2]=sum(Haplotypes_testindiv[,1])+sum(Haplotypes_testindiv[,2])+sum(Haplotypes_testindiv[,3])+sum(Haplotypes_testindiv[,4])

### Sum it across all values of a and add it to a new column.
Table[,(f+1)]=sapply(Table[,1],getvalues,Haplotypes_testindiv=Haplotypes_testindiv,HaplotypeFreq_panel_testindiv=HaplotypeFreq_panel_testindiv)
}

# Total readpairs minus readpairs for the specific chromosome (for jackknife)
Table3[,2]=sum(Table3[,2])-Table3[,2]
write.table(Table3,file=paste(args[3],"_fullreadpairs.jin",sep=""),sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)

Table4=Table
Table4[,1]=NULL
# Table5 stores final values of LOD scores.
Table5=data.frame(0:6000)
Table5[,1]=seq(-0.1,0.5,by=0.0001)
Table5[,2]=rowSums(Table4)
# First round is with no chromosomes remove
top=Table5[which.max(Table5[,2]),]
write.table(Table5,file=paste(args[3],"_",args[4],"panel_1ind_noChr","0","_alphatable_fullreads.txt",sep=""),sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(top,file=paste(args[3],"_",args[4],"panel_1ind_noChr","0","_topscore_fullreads.txt",sep=""),sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)

### Remove chromosomes one at a time for jackknife.
for(f in 1:22){
temp=Table4
temp[,f]=NULL
Table5[,2]=rowSums(temp)
top=Table5[which.max(Table5[,2]),]
write.table(Table5,file=paste(args[3],"_",args[4],"panel_1ind_noChr",f,"_alphatable_fullreads.txt",sep=""),sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(top,file=paste(args[3],"_",args[4],"panel_1ind_noChr",f,"_topscore_fullreads.txt",sep=""),sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
}
