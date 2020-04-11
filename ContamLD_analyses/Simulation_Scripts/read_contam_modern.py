"""
Simulate contamination in a modern sample using a modern contaminant and use an ancient file to decide how many reads are created per snp

run using the following commands:

python read_contam_modern.py target_ind contam_ind ancient_ind output_directory

use the following file formats
contamfile = "%s.geno" %contamind
targetfile = "%s.geno" %targetind
ancientfile = "%s_pull.out.readdepth" %ancient_ind
ancientdamagefile = "%s_pull_dam.out.readdepth"	%ancient_ind
"""

import numpy as np
import random
import sys

#define input parameters
targetind = sys.argv[1]
contamind = sys.argv[2]
ancientind = sys.argv[3]
outdir = sys.argv[4]

#define contamination rates to simulate
contam_rate_list=[0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15]

ancientfile = "%s_pull.out.readdepth" %ancientind
ancientdamagefile = "%s_pull_dam.out.readdepth" %ancientind
contamfile = "%s.geno" %contamind
targetfile = "%s.geno" %targetind

ancient = np.loadtxt(fname = ancientfile, dtype = 'string')
contam = np.loadtxt(fname = contamfile, dtype = 'string')
damage = np.loadtxt(fname = ancientdamagefile, dtype = 'string')
target = np.loadtxt(fname = targetfile, dtype = 'string')

#make a standard damage file to use for all contamination rates
damoutfile = "%s/%s_%s_damaged.readdepth" %(outdir, targetind, contamind)
out_dam = []

#iterate through snps in ancient data
pos=-1
for row in ancient:
	pos += 1
	ref_dam = 0
	alt_dam = 0
        # determine number of damaged reads
        dam_ref = int(damage[pos][7])
        dam_alt = int(damage[pos][8])
        dam_tot = dam_ref+dam_alt
        #add damage reads to totals
        for num in range(dam_tot):
        	#if ind is homozygous at site, then assign the corresponding allele
                 if target[pos] == '0':
                        alt_dam += 1
                 elif target[pos] == '2':
                        ref_dam += 1
                 #if ind is a het at site, assign ref or alt allele with 0.5 frequency
                 elif target[pos] == '1':
                        rand_ref = random.random()
                        if rand_ref <= .5:
                                ref_dam += 1
                        else:
                                alt_dam += 1
	out_dam.append("%s  %s  %s	%s	%s	%s	%s	%s	%s" %(row[0], row[1], row[2], row[3], row[4], row[5], row[6], ref_dam, alt_dam))

#write damage data to file
file = open(damoutfile, 'w')
for row in out_dam:
        file.write('%s\n' %row)
file.close()

#iterate through all the desired contamination rates and make the non damage restricted versions of the data

#load newly created damaged file for reference
damagefile = np.loadtxt(fname = damoutfile, dtype = 'string')

for contam_rate in contam_rate_list:
	#name outfile
	outfile = "%s/%s_%s_%s_all.readdepth" %(outdir, targetind, contamind, contam_rate)
	contam_rate=contam_rate/float(1.-contam_rate)
	#reset output data (out) and position to read in geno file
	pos = -1
	out = []
	out_dam = []
	#iterate through snps in ancient data
	for row in ancient:
		pos += 1
		#set counts of ref and alt alleles to zero, for both damaged and total
		ref = 0
		alt = 0
		ref_dam = 0
		alt_dam = 0
		# determine number of damaged reads
		dam_ref = int(damagefile[pos][7])
		dam_alt = int(damagefile[pos][8])
		dam_tot = dam_ref+dam_alt
		#determine the total number of undamaged reads
		undam_ref = int(row[7]) - dam_ref
		undam_alt = int(row[8])	- dam_alt
		undam_tot = undam_ref + undam_alt
		#simulate calls for 'sum' number of reads
		for num in range(undam_tot):
			cont_rand = random.random()
			#read is contaminated
			if cont_rand > 1-contam_rate:
				#if ind is homozygous at site, then assign the corresponding allele
    				if contam[pos] == '0':
    					alt += 1
    				elif contam[pos] == '2':
    					ref += 1
				#if ind is a het at site, assigne ref or alt allele with 0.5 frequency
    				elif contam[pos] == '1':
    					rand_ref = random.random()
    					if rand_ref <= .5:
    						ref += 1
    					else:
    						alt += 1
			#read is not contaminated
			else:
			        #if ind is homozygous at site, then assign the corresponding allele
                        	if target[pos] == '0':
                                        alt += 1
                                elif target[pos] == '2':
                                        ref += 1
                                #if ind is a het at site, assign ref or alt allele with 0.5 frequency
                                elif target[pos] == '1':
                                     	rand_ref = random.random()
                                        if rand_ref <= .5:
                                                ref += 1
                                        else:
                                             	alt += 1
		#add damage reads to totals
		ref += dam_ref
		alt += dam_alt			

		out.append("%s	%s	%s	%s	%s	%s	%s	%s	%s" %(row[0], row[1], row[2], row[3], row[4], row[5], row[6], ref, alt))

	
	#write data to file
	file = open(outfile, 'w')
	for row in out:
		file.write('%s\n' %row)
	file.close()



