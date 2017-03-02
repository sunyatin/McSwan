"""
Created on Fri Feb 19 2016
@author: Remi Tournebize
"""
### MODIFIED TO HANDLE >> temp/temp.txt with multiple ".ms" lines 
### MODIFIED TO HANDLE downsampling 02 March 2017

## added functionality to get POSITIONS
## WORKS ONLY IF SINGLE SIMULATION @input file

import argparse
import os.path
import time
import sys
start_time = time.time()

# convert MS or MSMS output to multidimensional SFS

# v3: incorporates island merging capacity

# input format
# output format

# only strict biallelic SNP encoded as 0 or 1
# no NA allowed in segsites
# haplotypes in rows

###########################################
###########################################
###########################################
# functions

def get_islands(linarr):
	if '-N' in linarr:
		msms = linarr.index('-N')
		linarr.pop(msms)
		linarr.pop(msms)
	digits = [x for x in linarr if x.isdigit()]
	ngam = int(digits[0])
	if '-I' in linarr:
		islstart = linarr.index('-I')
		nisl = int(linarr[islstart+1])
		islands = []
		for i in range(nisl):
			islands = islands + [int(linarr[islstart+1+i+1])]
	else:
		islands = [ngam]
	return(islands)

def rep_lengthout(x, l):
	if type(x) is not list: x = [x]
	n = l/len(x)
	y = []
	for i in range(n): y = y + x
	return(y)

def rep_each(x, e):
	y = []
	for i in x:
		for h in range(e): y = y + [i]
	return(y)

def prod(x):
	p= 1
	for n in x:
		p *= n
	return p

def multisfs_template(isl):
	l = prod([x+1 for x in isl]) # number of items in multiSFS
	ni = len(isl)
	for d in list(reversed(range(1,ni+1))):
		if d == ni:
			r = rep_lengthout(range(isl[d-1]+1), l)
			s = [str(x) for x in r]
		else:
			ii = prod([x+1 for x in isl[d:(ni+1)]])
			r = rep_lengthout(rep_each(range(isl[d-1]+1), ii), l)
			for i in range(len(r)):
				s[i] = str(r[i])+'.'+s[i]
	return(s)

def get_multisfs(segsites, islands, merge, multisfs_template, pos, downsizes):
	thresh = 1. * sum(islands) / 2

	if len(segsites) == 0:
		SFS = [0] * len(multisfs_template)
	else:

		nseg = len(segsites[0])

		prev_idx = 0
		ds = []
		for i in range(len(islands)):
			if islands[i] == 0:
				s = [0] * nseg
			else:
				idx = range(prev_idx, prev_idx + downsizes[i])
				#prev_idx = idx[-1] + 1
				prev_idx = sum(islands[0:i]) - 1
				s = segsites[idx[0]]
				for i in range(1, len(idx)):
					s = [x + y for x, y in zip(s, segsites[idx[i]])]
			ds.append(s)

		# merge, if requested
		if len(merge) > 0:
			isl2rem = set()
			for m in merge:
				isl_from = m[0] - 1
				isl_to = m[1] - 1
				ds[isl_to] = [x + y for x, y in zip(ds[isl_to], ds[isl_from])]
				isl2rem.add(isl_from)
			for x in sorted(isl2rem, reverse=True):
				del ds[x]

		# get keys
		#keys = [str(x) for x in ds[0]]
		#for i in range(1, len(ds)):
		#	keys = [x+'.'+str(y) for x, y in zip(keys, ds[i])]

		# fold, if requested
		if doFold == True:
			keys = ds[0]
			for i in range(1, len(ds)):
				keys = [[x,y] for x, y in zip(keys, ds[i])]
			for k in range(len(keys)):
				sumAC = sum(keys[k])
				if sumAC == sum(islands) or sumAC == 0:
					return None
				if sumAC > thresh:
					keys[k] = [islands[w]-keys[k][w] for w in range(len(keys[k]))]
					keys[k] = '.'.join([str(x) for x in keys[k]])
					keys[k] = pos[k]+' 1 '+'ac.'+keys[k]
				elif sumAC == thresh:
					fkey = [islands[w]-keys[k][w] for w in range(len(keys[k]))]
					fkey = '.'.join([str(x) for x in fkey])
					fkey = pos[k]+' 0.5 '+'ac.'+fkey
					keys[k] = '.'.join([str(x) for x in keys[k]])
					# combine as string
					keys[k] = pos[k]+' 0.5 '+'ac.'+keys[k]
					keys[k] = keys[k]+'\n'+fkey					
				else:
					keys[k] = '.'.join([str(x) for x in keys[k]])
					keys[k] = pos[k]+' 1 '+'ac.'+keys[k]
		# else, unfold
		else:
			# get keys
			keys = [str(x) for x in ds[0]]
			for i in range(1, len(ds)):
				keys = [x+'.'+str(y) for x, y in zip(keys, ds[i])]
			# add pos
			for i in range(len(keys)):
				keys[i] = pos[i]+' 1 '+'ac.'+keys[i]

	return(keys)

###########################################
###########################################
###########################################
# arguments
parser = argparse.ArgumentParser(description='options to multisfs computation')
parser.add_argument('-i', '--inputFile', type=str, required=True,
                    action="store", help="input file path")
parser.add_argument('-o', '--outputFile', type=str, required=True,
                    action="store", help="output file path")
parser.add_argument('-d', '--downsample', required=True, type=int, nargs="+", action="store", help="array of downsampled size for each population")
parser.add_argument('-m', '--merge', type=int, nargs="+", default=-1, required=False,
                    action="store", help="if multisfs needs to be done on merged island, specify island merging")
parser.add_argument('-L', '--chrSize', type=float, required=True,
                    action="store", help="size of the chromosome fragment")
parser.add_argument('--fold', required=False,
                    action="store_true", help="add this option to fold the multisfs")
parser.add_argument('-v', '--verbose', required=False,
                    action="store_true", help="add this option to increase software verbosity")

args = parser.parse_args()
fpath = args.inputFile
fout = args.outputFile
merge_sp = args.merge
doFold = args.fold
verbose = args.verbose
chrSize = args.chrSize
downsizes = args.downsample

if doFold == True:
	print "Folding:\ttrue"
else:
	print "Folding:\tfalse"
	
print "Downsampling to sizes: "+str(downsizes)+"\n"

if merge_sp != -1:
	sys.exit("The current implementation cannot handle merging (-m option) while downsampling.\n")

###########################################
###########################################
###########################################
# script

ff = open(fout, 'w')
ff.write('# Input VCF: '+fpath+'\n')
ff.write('# Chromosome size: '+str(chrSize)+'\n')
ff.write('# \n')
if doFold == True:
	ff.write('# Minor Allele Counts => folded SFS\n')
else:
	ff.write('# Derived Allele Counts => unfolded SFS\n')

with open(fpath) as f:

	# get info from ms command
	line = f.readline()
	line = line.rstrip()
	linarr = line.split(' ')
	islands = get_islands(linarr)

	# get info about island merging
	if merge_sp != -1:
		nmerge = len(merge_sp) / 2
		merge = []
		isl2rem = set()
		merged_islands = list(islands) # deep copy
		for i in range(nmerge):
			isl_from = int(merge_sp[i])
			isl_to = int(merge_sp[i+1])
			merge.append([isl_from, isl_to])
			isl2rem.add(isl_from-1)
			merged_islands[isl_to-1] = merged_islands[isl_to-1] + merged_islands[isl_from-1]
			if verbose: print "merge: "+str(isl_from)+" => "+str(isl_to)
		for x in sorted(isl2rem, reverse=True):
			del merged_islands[x]
	else:
		merged_islands = list(islands)
		merge = []
	
	# since downsampling
	if len(merged_islands) != len(downsizes): 
		sys.exit("The number of islands in the MS command must exactly match the number the downsampled sizes (-d option)\n.")
	for i in range(len(merged_islands)):
		if merged_islands[i] < downsizes[i]:
			sys.exit('The downsampled size must be inferior to the original island sizes specified in the MS command.\n')

	# print header for out
	multisfs = multisfs_template(downsizes)
	#ff.write('\t'.join(['ac.'+str(x) for x in multisfs])+'\n')

	segsites = []
	pos = []
	read = False
	first = True
	for line in f:

		if read == True and '-t' in line:
			read = False

		if line.startswith('//'):
			read = True
			if first == False:
				sfs = get_multisfs(segsites, islands, merge, multisfs, pos, downsizes)
				if sfs != None:
					ff.write('\n'.join(str(x) for x in sfs)+'\n')
				# reinitialize
				first = False
				segsites = []

		if read == True and line[0].isdigit():
			line = line.rstrip()
			segsites.append([int(x) for x in list(line)])
			first = False
			
		if read == True and line[0] == 'p':
			line = line.rstrip()
			pos = line.split(' ')
			pos.pop(0)
			pos = [str(int(float(p)*chrSize)) for p in pos]

# output of final rep
sfs = get_multisfs(segsites, islands, merge, multisfs, pos, downsizes)
if sfs != None:
	ff.write('\n'.join(str(x) for x in sfs)+'\n')

ff.close()

if verbose: print("--- %s seconds ---" % (time.time() - start_time))



