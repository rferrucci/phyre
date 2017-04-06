#!/usr/bin/env python

"""
CODENAME:	  PhyRe
DESCRIPTION:  

Copyright (c) 2010-2017 Ronald R. Ferrucci, Federico Plazzi, and Marco Passamonti..

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following
conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

"""
import argparse
import sys
from PhyReMethods import *

"""p = permutations for confidence intervals, d1 and d2 are range for number of 
species for funnel plot. parameter: m = AvTD, v = VarTD, e = euler, b = AvTD and VarTd. 
ci = confidence intervals b = batch file. l = user-defined path lengths
defaults are chosen when option is not indicated:
d1, d2 = 30, 70
permutations = 1000
confidence intervals = yes
user-defined path lengths = no
batch file = no
missing data = no
outfile = phyre.out
"""

def main():
	"""reads command line and obtains parameters for program
	"""
	parser = argparse.ArgumentParser(prog="PhyRe", prefix_chars='-+')

	parser.add_argument('samplefile', type=argparse.FileType('r', encoding="ISO-8859-1"),
		default=sys.stdin)
	parser.add_argument('popfile', type=argparse.FileType('r', encoding="ISO-8859-1"),
		default=sys.stdin)

	parser.add_argument("-d1", "--dimenesion1", dest="d1", nargs="?", default=30, type=int)
	parser.add_argument("-d2", "--dimenesion2", dest="d2", nargs="?", default=70, type=int)

	###-----------------options-------------------------###

	ci = 'y'; b = 'n'; l = 'n'
	batch = b; pathlengths = l; missing = 'n'

	parser.add_argument("-p", "--permutations", dest="p", type=int, nargs="?", default=1000)
	parser.add_argument("-ci", "--confidence", dest="ci", nargs="?", default='y')
	parser.add_argument("-b", "--batch", dest="b", nargs="?", default='n')
	parser.add_argument("-l", "--pathlengths", dest="l", nargs="?", default='n')
	parser.add_argument("-m", "--missingdata", dest="m", nargs="?", default='n')
	parser.add_argument("-o", "--outfile", dest="out")

	args = parser.parse_args()
	return args
	
if __name__ == "__main__":
	args = main()
	sample, population, taxon = getDataStructures(args)
	coef, pathLengths, taxonpopN = PathLength(args, population, taxon)
	atd,taxonN, Taxon = ATDmean(sample, population, taxon, coef)
	vtd = ATDvariance(taxonN,sample,atd, taxon, coef)
	Eresults = euler(sample,atd, taxonN, taxon, Taxon, coef)

	f = args.samplefile.name
	
	results = {}; popStats = {}
	popStats['taxon'] = taxon
	popStats['taxonPopN'] = taxonpopN
	popStats['pathlengths'] = pathLengths
	
	results[f] = {}
	results[f]['atd'] = atd
	results[f]['vtd'] = vtd
	results[f]['euler'] = Eresults
	results[f]['N'] = taxonN
	results[f]['n'] = len(sample)
	results[f]['TaxonN'] = {t: len(Taxon[t]) for t in taxon}
	printResults(args, results, popStats, population, taxon, coef)
	#o.write ("---------------------------------------------------")
