# ContamLD

ContamLD is a software is designed to estimate autosomal contamination in ancient DNA samples.

**Citation:**  Nakatsuka, N.\*; Harney, E.\*; Mallick, S..; Mah, M.; Patterson, N.; Reich, D. "Estimation of Ancient Nuclear DNA Contamination Using Breakdown of Linkage Disequilibrium." BioRxiv.

### <p>Steps for use:</p>
#### <p>Section 1)  Pre-processing steps:</p>
##### Part 1)  Download panels<br/>
	Step 1) Download panels from https://reichdata.hms.harvard.edu/pub/datasets/release/contamLD<br/>
		Note: In most cases you should download the 1240K panels. If you have low coverage (<0.5X) whole-genome shotgun sequences, then you can try the SG_panels for improved power at the expense of significantly increased running time and memory requirements.<br/>
	Step 2) Put the panels in the same folder (referred to as "directory_orig" below) that you have the helperdir folder in.

Note: If you have a SNP set that is very different than the 1240k SNP set or whole-genome shotgun set then follow steps in PreProcessing folder to make your own panel.


##### Part 2) Pull down reads onto SNP set.<br/>
	Step 1)  Obtain individual readdepth files for damaged or undamaged reads for each sample.  (Could also use genotype call information)<br/>
	Note:  This should be done using a pulldown on bam files. (See separate note on pulldown; you will need to grep out each sample from the full readdepth file).<br/>
	Name them Prefix_All.readdepth and Prefix_dam.readdepth
	(Prefix is the name of your file, same as Sample_ID below)

Optional) If you have eigenstrat files and are unable to pull down read information from bams, use the eig2readdepth.py script in the PreProcessing folder to transform eigenstrat files to readdepth files:
(Note: this has less power than the read based method because it ignores reads that map to the same site)<br/>
Put files in the format: Prefix.snp, Prefix.ind, Prefix.geno and damaged reads: Prefix_dam.snp, Prefix_dam.ind, Prefix_dam.geno<br/>
Use this file for eigenstrat format files in pseudo-haploid format (one read chosen to represent the genotype either 0=ALT or 2=REF; no heterozygotes). If you have diploid data with heterozygotes, use -d flag.
```python
python eig2readdepth.py [-d] Prefix
```


##### Part 3) Determine what panel the target individual is genetically closest to:<br/>
Step 1) Use outgroupf3.R script in PreProcessing folder to run outgroupf3 statistics to determine which panel is genetically closest to the target individual.<br/>
Note: Guessing on this step is okay as long as the sample is within continental ancestry variation of the 1000 Genomes population.

#### <p>Section 2:  Run Contamination Estimate:</p>
Note: Run this with 3 cores if possible.<br/>
-In the following notation: directory_orig is the directory with helperdir and panels; directory_files is the directory where your .readdepth and Prefix_inds.txt are.<br/>

Step 1) Run the following, where directory_files is the directory your files are in, directory_orig is the directory where the helperdir and panels are:<br/>
```
cd directory_files
mkdir -p directories
cd directory_orig

#!/bin/bash
while read sampleID panel; do
bash ./helperdir/Contamination_Analysis.txt directory_orig directory_files ${sampleID} ${panel}
done < directory_files/Prefix_inds.txt
```

#### <p>Section 3) Post-processing</p>
Note: The script will automatically do both the damage correction and the external correction version. Set External Correction value to the external correction score of on an uncontaminated individual of the same group as your target individual.<br/>
Note: The first time this script is run, sometimes it has an error because the files are not yet finished before they are needed for another script. If this happens, re-run the script.<br/>
```
cd directory_orig
bash ./helperdir/Post_Processing_New.txt directory_orig directory_files Prefix_inds.txt External_Correction_Value
```


**Examining the output of ContamLD:**<br/>
The file of interest will be in the directory "directory_files" and will be named "FinalContamScoresK_damageratio_Prefix_inds.txt"<br/>
If the warning "Model_Misspecified" shows up, this usually means the coverage is very low and the estimate might not be reliable. The warning "Very_High_Contamination" sometimes comes up even if there is low contamination if there is contamination from another ancient source or the panel is very diverged from the ancestry of the target individual (e.g. this occurs with most Native American individuals and some African groups).


**Example to test:**<br/>
There is an example to test in the folder exampledir. These are 2 samples (I7210 and I7278) with 0.04 contamination from an ancient West Eurasian (I10895).<br/>
If you run these samples with the Example_inds.txt file, you should get approximately the following estimates:<br/>
I7210:  0.039<br/>
I7278:  0.040

The files that should be generated are in the directories folder in the exampledir. To test if it is running properly for you, move these files to another location and see if the files you generate are approximately the same as these files.


