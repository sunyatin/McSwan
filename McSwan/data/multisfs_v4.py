"""
Created on Fri Feb 19 2016
@author: Remi Tournebize
"""
### MODIFIED TO HANDLE >> temp/temp.txt with multiple ".ms" lines 

import argparse
import os.path
import time
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

def get_multisfs(segsites, islands, merge, multisfs_template):
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
				idx = range(prev_idx, prev_idx + islands[i])
				prev_idx = idx[-1] + 1
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

		# fold, if requested
		if doFold == True:
			keys = ds[0]
			for i in range(1, len(ds)):
				keys = [[x,y] for x, y in zip(keys, ds[i])]
			sfs = {key: 0 for key in multisfs_template}
			for key in keys:
				sumAC = sum(key)
				if sumAC > thresh:
					key = [islands[w]-key[w] for w in range(len(key))]
					key = '.'.join([str(x) for x in key])
					sfs[key] = sfs[key] + 1
				elif sumAC == thresh:
					fkey = [islands[w]-key[w] for w in range(len(key))]
					key = '.'.join([str(x) for x in key])
					fkey = '.'.join([str(x) for x in fkey])
					sfs[key] = sfs[key] + .5
					sfs[fkey] = sfs[fkey] + .5
				else:
					key = '.'.join([str(x) for x in key])
					sfs[key] = sfs[key] + 1
			SFS = [sfs[x] for x in multisfs_template]
			# always set the monomorphic cell at ZERO
			SFS[0] = 0
		# else, unfold
		else:
			# get keys
			keys = [str(x) for x in ds[0]]
			for i in range(1, len(ds)):
				keys = [x+'.'+str(y) for x, y in zip(keys, ds[i])]
			# return sfs
			sfs = {key: 0 for key in multisfs_template}
			for key in keys:
				sfs[key] = sfs[key] + 1
			SFS = [sfs[x] for x in multisfs_template]
			# always set the monomorphic cell at ZERO
			SFS[0] = 0

	return(SFS)

###########################################
###########################################
###########################################
# arguments
parser = argparse.ArgumentParser(description='options to multisfs computation')
parser.add_argument('-i', '--inputFile', type=str, required=True,
                    action="store", help="input file path")
parser.add_argument('-o', '--outputFile', type=str, required=True,
                    action="store", help="output file path")
parser.add_argument('-m', '--merge', type=int, nargs="+", default=-1, required=False,
                    action="store", help="if multisfs needs to be done on merged island, specify island merging as a space-separated list of pairs, eg IslFrom1 IslDest1 IslFrom2 IslDest2 .. IslFromN IslDestN")
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

if verbose == True:
	if doFold == True:
		print "Folding:\ttrue"
	else:
		print "Folding:\tfalse"

###########################################
###########################################
###########################################
# script

ff = open(fout, 'w')

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

	# print header for out
	multisfs = multisfs_template(merged_islands)
	ff.write('\t'.join(['ac.'+str(x) for x in multisfs])+'\n')

	segsites = []
	read = False
	first = True
	for line in f:

		if read == True and '-t' in line:
			read = False

		if line.startswith('//'):
			read = True
			if first == False:
				sfs = get_multisfs(segsites, islands, merge, multisfs)
				ff.write('\t'.join(str(x) for x in sfs)+'\n')
				# reinitialize
				first = False
				segsites = []

		if read == True and line[0].isdigit():
			line = line.rstrip()
			segsites.append([int(x) for x in list(line)])
			first = False

# output of final rep
sfs = get_multisfs(segsites, islands, merge, multisfs)
ff.write('\t'.join(str(x) for x in sfs)+'\n')

ff.close()

if verbose: print("--- %s seconds ---" % (time.time() - start_time))



