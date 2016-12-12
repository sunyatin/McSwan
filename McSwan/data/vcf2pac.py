#!/usr/bin/python
__author__ = "Remi Tournebize 14042016"
__email__ = "remi.tournebize at gmail dot com"
import argparse
import os
import sys

## NOTA BENE (14042016):
# any SNP with at least one missing genotype is SKIPPED
# works only for biallelic SNPs (the script skips automatically the non-biallelic SNPs)
# PAC = Population (Derived) Allele Count
# now handles SPECTRUM FOLDING

### ATTENTION! le sorted iterkeys fait de facon alphabetique !! pas bon !!

# remove monomorphic SNPs

####################################################################
####################################################################
####################################################################
# ARGUMENT PARSING [UNIX compatible]

parser = argparse.ArgumentParser()

parser.add_argument("-vcf", '--fvcf', type=str, required=True, action='store',
         help='path to input VCF file')
parser.add_argument("-o", '--fpac', type=str, required=True, action='store',
         help='path to output PAC file')
parser.add_argument("-p", '--pop_sample_file', type=str, required=True, action='store',
         help='path to file specifying assignment of samples to populations\n Format of the file: samples on each line with tab-separated fields (the two first ones are mandatory, the third optional), note that POPULATION_ID should be an integer identifying the population and respecting the order of islands specified in the MS command:\n SAMPLE_NAME POPULATION_ID COMMENTS')
parser.add_argument("-chr", '--chrom', type=str, required=True, action='store',
         help='chromosome to analyse (must match the name in VCF)')
parser.add_argument("-fold", action='store_true',
         help='add this switch if you want to fold the multidimensional joint site-frequency-spectrum')
parser.add_argument("-minQ", '--minQual', type=float, default=10.0, action='store',
         help='SNPs with QUAL < minQual will be skipped')

args = parser.parse_args()

fvcf = args.fvcf
fpac = args.fpac
pop_sample_file = args.pop_sample_file
chrom = args.chrom
minQual = args.minQual
doFold = args.fold

# prompt to user
print "> Input file: "+fvcf
print "> Minimum QUAL: "+str(minQual)
print "> Chromosome: "+chrom
if doFold==False:
	print "> SFS type: Unfolded"
else:
	print "> SFS type: Folded"
#print "Populations: "+str(pop_sample_names)

####################################################################
####################################################################
####################################################################
# SCRIPT

def get_ind_AC(GT):
	# GT always first subfield in sample_data (VCF format 4.0)
	
	#GT = ind.split(':')[0]
	#AC = 9 # by default, 9=NA
	#if '.' not in GT:
	#	GT = GT.replace('/','').replace('|','')
	#	AC = sum(int(x) for x in GT)
	# return AC
	
	# because filter on biallelism is already done, extract simply first and second characters:
	
	if GT[0] == '.':
		return -1
	else:
		return int(GT[0])+int(GT[2])

def fold(AC, pop_n_gametes, threshold):
	sAC = sum(AC)
	if sAC > threshold:
		return [[pop_n_gametes[i]-AC[i] for i in range(len(AC))]]
	elif sAC == threshold:
		return [AC, [pop_n_gametes[i]-AC[i] for i in range(len(AC))]]
	else:
		return [AC]

# linarr[3]: REF
# linarr[4]: ALT
# linarr[5]: QUAL

if os.path.isfile(pop_sample_file) is False: sys.exit("Population file does not exist.")

# get population structure
# store into a dictionary:  pops = { popID => [sampleIDs] }
pops = {}
with open(pop_sample_file, 'r') as fp:
	for line in fp:
		line = line.rstrip()
		linarr = line.split('\t')
		if linarr[1].isdigit() == False: sys.exit("POPULATION_ID at line "+str(line)+" is not an integer")
		if linarr[1] not in pops.keys(): pops[linarr[1]] = set()
		pops[linarr[1]].add(linarr[0])
# convert to list of lists
# pop_sample_names = [ [sam01, sam02, ...], [sam11, sam12, ..] ]
pop_sample_names = []
print "> Populations:\n"
popIDs = sorted([int(x) for x in pops])
pop_n_gametes = []
for key in popIDs:
	key = str(key)
	pop_sample_names.append( list(pops[key]) )
	pop_n_gametes.append( 2 * len(pops[key]) )
	print "|___ Population "+str(key)+" with "+str(len(pops[key]))+" samples:"
	print list(pops[key])
	print "\n"
threshold = 1. * sum(pop_n_gametes) / 2

# check existence of input file
if os.path.isfile(fvcf) is False: sys.exit("Input file does not exist.")

# open output file
fout = open(fpac, 'w')
fout.write('# Input VCF: '+fvcf+'\n')
fout.write('# Chromosome: '+chrom+'\n')
fout.write('# Min. QUAL: '+str(minQual)+'\n')
if doFold == True:
	fout.write('# Minor Allele Counts => folded SFS\n')
else:
	fout.write('# Derived Allele Counts => unfolded SFS\n')

# loop over VCF file
pop_sample_idx = []
with open(fvcf, 'r') as f:
	for line in f:
		# if header
		if line.startswith('#'):
			# if #CHROM line (mandatory), store sample index
			if line.startswith('#CHROM'):
				line = line.rstrip()
				linarr = line.split('\t')
				for p in pop_sample_names:
					pop_sample_idx.append([linarr.index(x) for x in p])
		elif line.startswith(chrom):
				line = line.rstrip()
				linarr = line.split('\t')
				# if strictly biallelic and with sufficient quality support
				if len(linarr[3]) == 1 and len(linarr[4]) == 1 and float(linarr[5]) >= minQual:
					pac = []
					for p in pop_sample_idx:
						ac = [get_ind_AC(x) for x in [linarr[x] for x in p]]
						if -1 in ac:
							continue
						sum_ac = sum(ac)
						pac.append(sum_ac)
					if doFold:
						pac = fold(pac, pop_n_gametes, threshold)
						total_ac = sum(pac[0])
						if len(pac)==1:
							# output SNP if NO monomorphic across populations
							if total_ac != 0:
								fout.write(linarr[1]+' 1 '+'ac.'+'.'.join([str(x) for x in pac[0]])+'\n')
						else:
							if total_ac != 0:
								fout.write(linarr[1]+' 0.5 '+'ac.'+'.'.join([str(x) for x in pac[0]])+'\n')
								fout.write(linarr[1]+' 0.5 '+'ac.'+'.'.join([str(x) for x in pac[1]])+'\n')
					else:
						total_ac = sum(pac)
						# output SNP if NO monomorphic across populations
						if total_ac != 0:
							fout.write(linarr[1]+' 1 '+'ac.'+'.'.join([str(x) for x in pac])+'\n')

fout.close()

print "VCF reading: done."



