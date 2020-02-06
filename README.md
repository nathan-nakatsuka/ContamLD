# ContamLD

ContamLD is a software is designed to estimate autosomal contamination in ancient DNA samples.

**Citation:**  Nakatsuka, N.\*; Harney, E.\*; Mallick, S..; Mah, M.; Patterson, N.; Reich, D. "Estimation of Ancient Nuclear DNA Contamination Using Breakdown of Linkage Disequilibrium." BioRxiv.<br/>
<br/>
**Contact:** Nathan Nakatsuka: nathan_nakatsuka@hms.harvard.edu
<br/>
<br/>
## <p>Steps for use:</p>
### <p>Section 1)  Pre-processing steps</p>
#### Part 1)  Download panels or prepare your own panels.<br/>
**Step 1)** Download panels from https://reichdata.hms.harvard.edu/pub/datasets/release/contamLD<br/>
Note: In most cases you should download the 1240K panels. If you have low coverage (<0.5X) whole-genome shotgun sequences, then you can try the SG_panels for improved power at the expense of significantly increased running time and memory requirements.<br/>
<br/>
**Step 2)** Put the panels in the "panels" folder where the "helperdir" folder also is located (referred to as "directory_orig" below).

Note: If you have a SNP set that is very different than the 1240k SNP set or whole-genome shotgun set then follow steps in PreProcessing folder to make your own panel.


#### Part 2) Pull down reads onto SNP set.<br/>
ContamLD requires files with all reads and with only damaged reads in the following format (tab delimited): <br/>	
SNP_ID	Chrom\*	Position\*	REF	ALT	::	IndName*	REF_allelecount\*	ALT_allelecount\*<br/>
<br/>
Chrom is chromosome number. Position is the position of the SNP on the chromosome in Hg19 coordinates. REF_allelecount and ALT_allelecount refer to the number of reads mapping to the reference or alternative allele, respectively (based on Hg19). \* indicates the necessary columns (the other columns can be filled with place holders, but the Chrom, Position, IndName, REF_allelecount, and ALT_allelecount must be in the 2nd, 3rd, 7th, 8th, and 9th columns, respectively). The files should have reads corresponding to the 1240K.snp or SG.snp files in the "PreProcessing" folder, or the snp files corresponding to the panels prepared by the user.<br/>	
<br/>	
**Step 1)** To obtain damaged reads, use PMDtools (https://github.com/pontussk/PMDtools) with a PMDscore threshold of 3.<br/>
<br/>
**Step 2)** The read count information can be obtained in any of the following ways:<br/>
Note:  Each readdepth file must contain only a single individual for Section 2.<br/>
<br/>
*Option 1)*  Use samtools mpileup with the appropriate base and mapping quality cutoffs (e.g. 30). <br/>
<br/>
*Option 2)* If you have eigenstrat files and are unable to pull down read information from bams, use the eig2readdepth.py script in the PreProcessing folder to transform eigenstrat files to readdepth files.<br/>
(Note: this has less power than the read based method because it ignores reads that map to the same site)<br/>
-Put files in the format: Prefix.snp, Prefix.ind, Prefix.geno and damaged reads: Prefix_dam.snp, Prefix_dam.ind, Prefix_dam.geno<br/>
Use this file for eigenstrat format files in pseudo-haploid format (one read chosen to represent the genotype, either 0=ALT or 2=REF; no heterozygotes). If you have diploid data with heterozygotes, use -d flag.
```python
python eig2readdepth.py [-d] Prefix
```
<br/>
*Option 3)* (future) Before the end of 2020 we plan to release a pulldown program that will allow users to automatically generate the readdepth files required for ContamLD (as well as genotype files for other purposes).<br/>
<br/>
**Step 3)**  Name the files IndName_All.readdepth and IndName_dam.readdepth (IndName is the name of that particular individual. IndName_All.readdepth is the file corresponding to all reads, and IndName_dam.readdepth is the file corresponding to only damaged reads for that individual.).<br/>


#### Part 3) Determine what panel the target individual is genetically closest to.<br/>
Follow HowtoDetermineClosestPanel_Outgroupf3.txt instructions in the "PreProcessing" folder to run outgroupf3 statistics to determine which panel is genetically closest to the target individual.<br/>
Note: Guessing on this step is okay as long as the sample is within continental ancestry variation of the 1000 Genomes population.
<br/>

#### Part 4) Create a file with the names of all individuals and their corresponding panels.<br/>
Create file called "GroupName_inds.txt" (where GroupName is the name of your collection of individuals) in the following format, where the 1000Genomes population is determined from Section 1 Part 3, and put it in same directory as the readdepth files:<br/>
IndName_1 1000Genomes_Pop_closesttoIndName_1<br/>
IndName_2 1000Genomes_Pop_closesttoIndName_2<br/>
IndName_3 1000Genomes_Pop_closesttoIndName_3<br/>
<br/>
<br/>

### <p>Section 2)  Run Contamination Estimate</p>
After you have readdepth files, the panels (placed in the panels folder), and the GroupName_inds.txt file, run ContamLD.<br/>
Note: Run this with 3 cores if possible.<br/>
-In the following notation: "directory_orig" is the directory with helperdir and panels folders are; "directory_files" is the directory where your .readdepth and GroupName_inds.txt are.<br/>

Run the following:<br/>
```
cd directory_files
mkdir -p directories
cd directory_orig

#!/bin/bash
while read IndName panel; do
bash directory_orig/helperdir/ContamLD_RunningScript.txt directory_orig directory_files ${IndName} ${panel}
done < directory_files/GroupName_inds.txt
```
<br/>

### <p>Section 3) Post-processing</p>
After ContamLD is run, the final values and standard errors are obtained with this script.
Note: The script will automatically do both the damage correction and the external correction version. Set "External_Correction_Value" to be the external correction score of an uncontaminated individual of the same group as your target individual. If you do not have this, set the score to 0. Panel_Type is the type of panel: 1240K, SG, or your own.<br/>
Note: The first time this script is run, sometimes it has an error because the files are not yet finished before they are needed for another script. If this happens, re-run the script. If it does not work the second time, then something went wrong upstream of this.<br/>
```
cd directory_orig
bash directory_orig/helperdir/Post_Processing_New.txt directory_orig directory_files GroupName_inds.txt External_Correction_Value Panel_Type
```


**Examining the output of ContamLD:**<br/>
The files of interest will be in the directory "directory_files". The damage restriction corrected estimates will be named "FinalContamScoresK_damageratio_Panel_Type_GroupName_inds.txt". The external correction estimates will be named "FinalContamScores_ExtCorr_Panel_Type_GroupName_inds.txt"<br/>
If the warning "Model_Misspecified" shows up, this usually means the coverage is very low and the estimate might not be reliable. The warning "Very_High_Contamination" sometimes comes up even if there is low contamination if there is contamination from another ancient source or the panel is very diverged from the ancestry of the target individual (e.g. this occurs with most Native American individuals and some African groups).


**Example to test:**<br/>
There is an example to test in the folder "exampledir". These are 2 samples (I7210 and I7278) with 0.04 contamination from an ancient West Eurasian (I10895).<br/>
If you run these samples with the Example_inds.txt file, you should get approximately the following estimates:<br/>
I7210:  0.039<br/>
I7278:  0.040

The files that should be generated are in the "directories" folder in the "exampledir" folder. To test if ContamLD is running properly for you, move these files to another location and see if the files you generate are approximately the same as these files.


