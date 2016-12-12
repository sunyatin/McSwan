"""
Created on Aug 19 2016
@author: Remi Tournebize
@description: In R, system command output is split with \n's when a line has more than 8095 characters, this script reconstruct the original line-ity
"""

import argparse
import os.path

###########################################
###########################################
###########################################
# arguments
parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', '--inputFile', type=str, required=True,
                    action="store", help="input file path")
parser.add_argument('-o', '--outputFile', type=str, required=True,
                    action="store", help="output file path")

args = parser.parse_args()
fpath = args.inputFile
fout = args.outputFile

###########################################
###########################################
###########################################
# script

ff = open(fout, 'w')

with open(fpath) as f:

	# get info from ms command
	line = f.readline()
	line = line.rstrip('\r\n')
	ff.write(line)
	
	while 1:
		line = f.readline()
		line = line.rstrip('\r\n')
		if not line:
			ff.write('\n\n')
			break
		else:
			ff.write(line)

	readPos = False
	nSeg = 0
	for line in f:
		line = line.rstrip('\r\n')
		
		if len(line) == 0:
			ff.write('\n')
			continue

		if line[0] == 's':
			linarr = line.split(' ')
			nSegsites = int(linarr[1])
			ff.write(line+'\n')
			continue
		elif line[0] == 'p':
			readPos = True

		if readPos == True:
			if '.' in line:
				ff.write(line)
			else:
				readPos = False
				ff.write('\n')
		
		if readPos == False:
			if line[0].isdigit():
				if nSegsites < 8095:
					ff.write(line+'\n')
				else:
					nchr = len(line)
					if nSeg + nchr < nSegsites:
						ff.write(line)
						nSeg = nSeg + nchr
					else:
						ff.write(line+'\n')
						nSeg = 0

			else:
				ff.write(line+'\n')

ff.write('\n')
ff.close()




