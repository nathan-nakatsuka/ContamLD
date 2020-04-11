"""
#Script allows user to simulate contamination in BAM files.

#Note: This script has many features hardcoded and cannot be run by users unless they replace the hardcoded elements with their own paths to the relevant files and software.
"""
#!/bin/bash

###To run ###########################
#Note that this file follows the same workflow as bam_contam.sh except for the way in which the damage restriction pulldown works
#bash bam_contam_ancientcontam.sh ${file}
#####################################

### Define Contaminant Individual 

#contaminant info
contam_ID="I10895"
contam_bam_loc="filtbam/outbam/${contam_ID}.out.bam"
contam_coverage="1.087059"
contam_RG="Aenar_1240k_plus_half_S10895.E1.L1_20180402_HLYM2BGX5_1:Aenar_1240k_plus_half_S10895.E1.L1_20180402_HLYM2BGX5_2:Aenar_1240k_plus_half_S10895.E1.L1_20180402_HLYM2BGX5_3:Aenar_1240k_plus_half_S10895.E1.L1_20180402_HLYM2BGX5_4:Aenar_1240k_plus_half_S10895.E1.L1_20180513_HMLVNBGX5_1:Aenar_1240k_plus_half_S10895.E1.L1_20180513_HMLVNBGX5_2:Aenar_1240k_plus_half_S10895.E1.L1_20180513_HMLVNBGX5_3:Aenar_1240k_plus_half_S10895.E1.L1_20180513_HMLVNBGX5_4"

### 1) Read in individuals to produce simulated contamination data for

#Take Reich lab anno file as input. Input file should contain at minimum, the following info. 
#ind name at column 1
#pointers to the original bam in column 4 & 5
#bam read groups defined in column 6
#coverage on the 1240k SNPs in column 12

file=$1

###########Set up loop part 1##########
#iterate through entire input file and run contamination on every individual

len=$(wc $file | awk '{print $1}') 

for line in $(seq 1 ${len})
do

coverage=$(awk -F "\t" -v var=${line} 'NR==var {print $12}' ${file})
target_ID=$(awk -F "\t" -v var=${line} 'NR==var {print $1}' ${file})


### 2) create filtered bams (i.e. Bams that use the same filtering parameters as applied to the 1240k pulldown, so that the coverage values are consistent)  

ind=$(awk -F "\t" -v var=${line} 'NR==var {print $1}' ${file})
bam_loc1=$(awk -F "\t" -v var=${line} 'NR==var {print $4}' ${file})
bam_loc2=$(awk -F "\t" -v var=${line} 'NR==var {print $5}' ${file})

#determine which bam format is used in the input file and define bam_loc pointer accordingly
if [ $bam_loc2 != "." ] && [ $bam_loc2 != ".." ]; then
	bam_loc="${bam_loc1}../${bam_loc2}bam"
else
	bam_loc=$bam_loc1
fi

#make the filtered bams
bash mkfiltbam.sh ${ind} ${bam_loc}
/home/np29/o2bin/pulldown -p filtbam/parfiles/${ind}.par 


#### START LOOP FOR CONTAMINATION RATES ####

for alpha in 0.000 0.005 0.010 0.015 0.020 0.025 0.030 0.040 0.050 0.060 0.070 0.080 0.090 0.100 0.150 0.200 0.300 0.400 0.500
do

#determine what level the contaminant individual should be downsampled to to make desired contamination rate (alpha)
rate_adjust=$(echo "$alpha/(1-$alpha)*$coverage/$contam_coverage" | bc -l)
target="outfiles/${target_ID}_${contam_ID}_${alpha}"
outfile="outfiles/${target_ID}_${contam_ID}_${alpha}"
outbam_filt="filtbam/outbam/${target_ID}.out.bam"

### 3) generate contaminated bams (if possible to achieve desired contamination rate for specified individuals)

if (( $(echo "$rate_adjust < 1" |bc -l) )); then
	bash downsample_merge.sh $outbam_filt $contam_bam_loc $coverage $contam_coverage $alpha $outfile
else
	echo "Not possible: Target- ${target_ID} Contaminant- ${contam_ID} Alpha- ${alpha}" >> ${file}_not_possible
fi


### 4) Pulldown contaminated individuals on 1240k datasets


if (( $(echo "$rate_adjust < 1" |bc -l) )); then
	target_RG=$(awk -F "\t" -v var=${line} 'NR==var {print $6}' ${file})

	ind="${target_ID}_${contam_ID}_${alpha}"
	bam_loc="/n/groups/reich/eadaoin/contamLD/simulations/v30.2/outfiles/${ind}.bam"

	bash mkdb.sh ${ind} ${bam_loc} ${target_RG} ${contam_RG}
	bash mkpar.sh ${ind}
	/n/groups/reich/matt/pipeline/static/pulldown -p pulldown/parfiles/${ind}.par > pulldown/outfiles/${ind}.out

fi


### 5-alternative) Set up damage restriction pulldown, in this script damaged reads also come from the contaminant individual 

#make downsampled bam with damage rate 0.05 (note that the samtools view -s parameter was giving weird results when 0 was used as the seed (i.e. the number before the decimal), so 1.05 was used to change the seed)
rate=1.05

if (( $(echo "$rate_adjust < 1" |bc -l) )); then

	target_RG=$(awk -F "\t" -v var=${line} 'NR==var {print $6}' ${file})

	ind="${target_ID}_${contam_ID}_${alpha}"
	bam="/n/groups/reich/eadaoin/contamLD/simulations/v30.2/outfiles/${ind}"

	samtools view -s ${rate} -b ${bam}.bam > restrict_${outfile}_ancientcontam.restrict.bam
	samtools index restrict_${outfile}_ancientcontam.restrict.bam

	#run dr pulldown

	restrict_bam_loc=/n/groups/reich/eadaoin/contamLD/simulations/v30.2/restrict_${outfile}_ancientcontam.restrict.bam

	bash mkdb_damage_ancientcontam.sh ${ind} ${restrict_bam_loc} ${target_RG} ${contam_RG}
	bash mkpar_damage.sh ${ind}_ancientcontam
	/home/np29/o2bin/pulldown -p pulldown/parfiles/${ind}_ancientcontam.dam.par > pulldown/outfiles/${ind}_ancientcontam.dam.out

fi


### 5b) Alternatively, if you want to set up damage restriction pulldown without data from contaminant (i.e. regular version)

#make downsampled bam

#rate=$(awk -F "\t" -v var=${line} 'NR==var {print $15}' ${file})
#rate=1.05

#if (( $(echo "$rate_adjust < 1" |bc -l) )); then

#        target_RG=$(awk -F "\t" -v var=${line} 'NR==var {print $6}' ${file})

#        ind="${target_ID}_${contam_ID}_${alpha}"
#        bam="/n/groups/reich/eadaoin/contamLD/simulations/v30.2/outfiles/${ind}"

#        samtools view -s ${rate} -b ${bam}.bam > restrict_${outfile}_ancientcontam.restrict.bam
#        samtools index restrict_${outfile}_ancientcontam.restrict.bam

         #run dr pulldown

#        restrict_bam_loc=/n/groups/reich/eadaoin/contamLD/simulations/v30.2/restrict_${outfile}_ancientcontam.restrict.bam

#        bash mkdb_damage.sh ${ind} ${restrict_bam_loc} ${target_RG}
#        bash mkpar_damage.sh ${ind}
#        /home/np29/o2bin/pulldown -p pulldown/parfiles/${ind}.dam.par > pulldown/outfiles/${ind}.dam.out

#fi


### 6 Run ANGSD


#if (( $(echo "$rate_adjust < $contam_coverage" |bc -l) )); then
#
#	#clip bams
#	java -jar /n/groups/reich/matt/pipeline/static/adnascreen-1.5.0-SNAPSHOT.jar softclip -b -n 2 -i outfiles/${ind}.bam -o clipped_bams/${ind}.bam
#	samtools index clipped_bams/${ind}.bam
#
#	
#	#run ANGSD
#	/home/mym11/external/angsd/angsd -i clipped_bams/${ind}.bam  -r X:5000000-154900000 -doCounts 1  -iCounts 1 -minMapQ 10 -minQ 20 -out ANGSD/${ind}.angsd
#	Rscript /home/mym11/external/angsd/R/contamination.R mapFile="/home/mym11/external/angsd/RES/map100.chrX.gz" hapFile="/home/mym11/external/angsd/RES/hapMapCeuXlift.map.gz" countFile="ANGSD/${ind}.angsd.icnts.gz" mc.cores=24 > ANGSD/${ind}.angsd.out
#
#fi



### 7) Delete Stuff

if (( $(echo "$rate_adjust < 1" |bc -l) )); then

	if [ "${target_ID}" != "${contam_ID}" ]; then
		
		ind="${target_ID}_${contam_ID}_${alpha}"

		rm -rf outfiles/${ind}.bam
		rm -rf outfiles/${ind}.bam.bai

		rm -rf clipped_bams/${ind}.bam
		rm -rf clipped_bams/${ind}.bam.bai
		
		rm -rf restrict_outfiles/${ind}.restrict.bam
		rm -rf restrict_outfiles/${ind}.restrict.bam.bai

	fi
		
fi

done

if [ "${target_ID}" != "${contam_ID}" ]; then

	rm -rf filtbam/outbam/${target_ID}.out.bam
	rm -rf filtbam/outbam/${target_ID}.out.bam.bai
fi

done

