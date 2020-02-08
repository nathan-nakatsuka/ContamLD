#!/bin/python
"""
This tool converts data in eigenstrat format in to readdepth format. Note that when possible, 
we recommend instead generating readdepth files directly from the bams.

Use the -d flag when using diploid data. Otherwise pseudohaploid data is assumed.

usage: eig2readdepth.py [-d] filename 
filename is the prefix to .snp, .geno, .ind files.
"""
import numpy as np
import argparse
import sys


parser = argparse.ArgumentParser(description="Generates an approximation of a readdepth file from eigenstrat data")

#add options to argparser
parser.add_argument(dest="filename", type=str)
parser.add_argument('-d', action="store_true", dest="diploid")


try:
	options=parser.parse_args()
except:
	parser.print_help()
	sys.exit(0)

filename=options.filename
diploid=options.diploid


#load data
snp = np.genfromtxt('%s.snp' %filename, dtype = str)
geno = np.genfromtxt('%s.geno' %filename, dtype = str)
ind = np.genfromtxt('%s.ind' %filename, dtype = str)

#check to make sure that snp and geno files are the same length
if len(snp) != len(geno):
	print('ERROR: snp and geno files are different lengths')
	sys.exit()

#open output file
file = open('%s.readdepth' %filename, 'w')

#check to see if data is pseudohaploid or diploid
if diploid == True:
	PH_adjust = 2
else:
	PH_adjust = 1

#determine number of individuals in dataset
if len(ind.shape) == 2:
	ind_len, junk = ind.shape
elif len(ind.shape) ==1:
	ind_len = 1
else:
	print('ERROR: check ind file')
	sys.exit(0)

 
#generate readdepth information
for num in range(len(snp)):
	for ind_pos in range(0,ind_len):
		alt_count =0 
		ref_count = 0
		if geno[num][ind_pos] == '0':
			alt_count += 1*PH_adjust
		elif geno[num][ind_pos] == '2':
			ref_count += 1*PH_adjust
		elif geno[num][ind_pos] == '1':
			if diploid == True:
				alt_count += 1
				ref_count += 1
			else:
				print('Error: geno file contains a 1 at line %s for ind %s, is not pseudohapoid' %(num, ind[ind_pos][0]))
				sys.exit()
		elif geno[num][ind_pos] != '9':
			print('Error: geno file contains unexpected value at line %s' %num)
			sys.exit()
		file.write('%s %s       %s %s %s :: %s  %s  %s\n' %(snp[num][0], snp[num][1], snp[num][3], snp[num][4], snp[num][5], ind[ind_pos][0], ref_count, alt_count)) 

file.close()
