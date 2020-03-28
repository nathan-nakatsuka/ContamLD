"""
Script takes data from an ancient individual and contaminates the readdepth file with data from another modern individual

Takes 4 inputs
1) name of ancient readdepth file - full data 
2) name of ancient readdepth file - damage restricted
3) name of contamination geno file 
4) output file prefix 
"""

import numpy as np
import random
import sys

contam_rate_list = [0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15]

ancientfile = sys.argv[1]
ancientdamagefile = sys.argv[2]
contamfile = sys.argv[3]
outfilename = sys.argv[4]

file = open('%s.log' %outfilename, 'w')
file.write('ancientfile = %s\n' %ancientfile)
file.write('ancientdamagefile = %s\n' %ancientdamagefile)
file.write('contamfile = %s\n' %contamfile)
file.write('outfilename = %s\n' %outfilename)
file.close()

ancient = np.loadtxt(fname = ancientfile, dtype = 'string')
contam = np.loadtxt(fname = contamfile, dtype = 'string')
damage = np.loadtxt(fname = ancientdamagefile, dtype = 'string')

for contam_rate in contam_rate_list:
	#name outfile
	contam_rate2=str(contam_rate)
	outfile = "%s_%s_contam_auto.out" %(outfilename, contam_rate2)
	contam_rate=contam_rate/float(1.-contam_rate)
	#reset output data (out) and position to read in geno file
	pos = -1
	out = []

	#iterate through snps in ancient data
	for row in ancient:
		pos += 1
		#set counts of ref and alt alleles to zero
		ref = 0
		alt = 0
		# determine number of damaged reads
		dam_ref = int(damage[pos][7])
		dam_alt = int(damage[pos][8])
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
    				elif contam[pos] =='2':
    					ref += 1
				#if ind is a het at site, assign ref or alt allele with 0.5 frequency (there are no 9s in this file)
    				else:
    					rand_ref = random.random()
    					if rand_ref <= .5:
    						ref += 1
    					else:
    						alt += 1
		#add damage reads to totals
		ref += (undam_ref + dam_ref)
		alt += (undam_alt + dam_alt)
		out.append("%s	%s	%s	%s	%s	%s	%s	%s	%s" %(row[0], row[1], row[2], row[3], row[4], row[5], row[6], ref, alt))

	#write data to file
	file = open(outfile, 'w')
	for row in out:
		file.write('%s\n' %row)
	file.close()
