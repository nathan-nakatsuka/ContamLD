"""
#Script allows user to downsample BAM files.

#Usage: bash downsample_merge.sh target contam coverage contam_coverage alpha outfile
"""

#!/bin/bash

target=$1
contam=$2
coverage=$3
contam_coverage=$4
alpha=$5
outfile=$6

###Test parameters
#contam="/n/data1/hms/genetics/reich/1000Genomes/amh_samples/amorim2018_6thCentury_barbarian_samples/B-wgs/B-fix/samples/SZ43/hg19/SZ43.mapped.rmdupse_adna_v2.md_no1kgvar.1240k.bam"
#target="/n/data1/hms/genetics/reich/1000Genomes/amh_samples/krzewinska2018_viking_samples/A-round1/B-fix/samples/vik_84001/hg19/84001_GAATCTC_merged.160227.hs37d5.fa.cons.90perc_rg.bam"
#coverage=3.682839
#alpha=0.01
#outfile=vik_84001.SG_SZ43.SG_0.01


#determine coverage to downsample to
rate_adjust=$(echo "$alpha/(1-$alpha)*$coverage/$contam_coverage" | bc -l)

if (( $(echo "$rate_adjust >  0" |bc -l) )); then

	##downsample contam file to specified value
	samtools view -s ${rate_adjust} -b ${contam} > temp_${outfile}


	###Consider Change the Readgroups in files before merging###

	##merge bam files
	samtools merge temp_${outfile}_unsorted ${target} temp_${outfile}


else

	#copy bam to temp_${outfile}
	cp ${target} temp_${outfile}_unsorted

fi

##sort bam files
samtools sort -o ${outfile}.bam temp_${outfile}_unsorted
samtools index ${outfile}.bam


rm -rf temp_${outfile}
rm -rf temp_${outfile}_unsorted





