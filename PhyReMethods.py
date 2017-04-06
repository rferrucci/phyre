#!/usr/bin/env python

"""
CODENAME:	  PhyReMethods
DESCRIPTION: Methods as required by other PhyRe scripts

Copyright (c) 2010-2017 Ronald R. Ferrucci, Federico Plazzi, and Marco Passamonti.

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

import sys
import argparse
from re import *
import random

def getDataStructures(args):
	""" returns list of samples, master list from the population at large, and list of taxon included"""
	popfile=args.popfile; samplefile=args.samplefile

	#if batch processing, take a list of files
	if args.b == 'y': 
		Files = []
		for i in open(samplefile):	
				j = i.strip()
				Files.append(j)
	f = popfile.readlines()
	
	header = f[0] #get header information
	
	sample = [s.split()[0] for s in samplefile.readlines()]
	sample = set(sample)

	f = [f[i].split() for i in range(1, len(f)) if f[i].split() != []]

	# taxon is list of taxon levels for which we have information
	taxon = [header.split()[i] for i in range(len(header.split()))[1:]]
	population = getPopulationData(args, f, taxon)
	species = population.keys()
	
	# test whether species are found in master list, if not program outputs missing species then dies
	try:
		missing = [s for s in sample if s not in species]
		if len(missing) != 0: raise Exception(', '.join(missing) + " not in master list " + popfile.name)
	except Exception as e:
		print(e)
		quit()
	else: return sample, population, taxon

def getPopulationData(args, f, taxon):
	""" get full taxonomic data from dataset """
	population = {f[i][0]: {} for i in range(len(f))}

	for line in f:
		species = line[0]
		if args.m == 'y':
			mtax = ''
			population[species][taxon[0]] = line[1]
			for i in range(1, len(taxon)):
				if x[Index[t]] == '/': population[species][taxon[i]] = population[species][taxon[i-1]]
				else: population[species][taxon[i]] = line[i+1]

		else:
			population[species] = {taxon[i]: line[i+1] for i in range(len(taxon))}
	return population
		
def StandPathLength(taxon):
	Ln = 100/2*(len(taxon) - 1)
	return Ln

def PathLength(args, population, taxon):
	""" get coefficients for path lengths """
	coef = {}; pathLengths= {}; taxonN ={}
	X = {}
	
	Taxon = {t: {} for t in taxon}
	sample = population.keys()
	for t in taxon:
		# at each taxon level, making list of taxon present
		Taxon[t] = [population[i][t] for i in sample if population[i][t] != '/']
		
	#count of taxon present at each level, add one for upper level
	n = [len(set([i for i in Taxon[t]])) for t in taxon]
	n.insert(0,1.0)

	raw = [1 - n[i]/n[i + 1] for i in range(len(n)-1)]

	s = sum(raw)

	adjco = [i*100/s for i in raw]
	coef = {}; pathLengths = {}
	for i in range(len(taxon)):
		t = taxon[i]
		coef[t] = sum(adjco[i:])
		pathLengths[t] = adjco[i]
		
	return coef,  pathLengths, {taxon[i]: n[i+1] for i in range(len(taxon))}

def ATDmean(sample, data, taxon, coef):
	""" calculate Average Taxonomic Distinctiveness """
	N = len(sample)
	
	taxonN = {}; AvTD = 0; n = 0
	#Taxon are counts of taxa at each level, taxonN are numbers of pairwise differences
	#at each level, with n being the accumlation of pairwise differences at that level. the difference
	#between n and TaxonN is the number of species that are in different taxa in that level
	#but not in upper levels

	Taxon = {t: {} for t in taxon}

	for t in taxon:
		x = [data[i][t] for i in sample]
		Taxon[t] = {i: x.count(i) for i in set(x)}
	
	for t in taxon:
		taxonN[t] = sum([Taxon[t][i] * Taxon[t][j] for i in Taxon[t] for j in Taxon[t] if i != j])
		n = taxonN[t] - n
		AvTD += (n * coef[t]) 
		n = taxonN[t]

	#print sample
	AvTD /= (N * (N - 1))
	return AvTD,taxonN, Taxon

def ATDvariance(taxonN, sample, atd, taxon, coef):
	""" calculate Variation in Taxonomic Distinctiveness """
	n = 0; vtd=0;
	for t in taxon:
		n = taxonN[t] - n
		vtd += n * coef[t]**2 
		n = taxonN[t]

	N = len(sample)
	n = N * (N - 1)
	
	vtd = (vtd - ((atd*n)**2)/n)/n

	return vtd
	
def euler(sample, atd, TaxonN, taxon, Taxon, coef):
	"""Von Euler's index of imbalance"""
	sample = list(sample)
	n = len(sample)
	TDmin = 0; N = 0
	
	for t in taxon:
		k = len(Taxon[t])
		TDmin += coef[t] * (((k-1)*(n-k +1)* 2+ (k-1)*(k-2))-N)
		N += ((k-1)*(n-k +1)* 2 + (k-1)*(k-2))-N
	
	TDmin /=(n * (n-1))
	taxon.reverse()
	TaxMax = {t: [] for t in taxon}
	
	taxonN = {}

	for t in taxon:
		j = 0
		n = len(Taxon[t])
		if taxon.index(t) == 0:	
			TaxMax[t] = [[i] for i in sample]
			TaxMax[taxon[taxon.index(t) + 1]] = [[] for i in range(len(Taxon[taxon[taxon.index(t) + 1]]))]
			for i in range(n):
				# cycle through list of taxa to be placed into higher taxa
				TaxMax[taxon[taxon.index(t) + 1]][j].append(sample[i])
				# j is incremented by 1 as we cycle through higher taxa
				j += 1
				if j == len(Taxon[taxon[taxon.index(t) + 1]]):
					# reset j to zero when last of higher taxa has been reached to start the cycle over again
					j =0
		elif taxon.index(t) == len(taxon) -1:
			continue
		else:
			TaxMax[taxon[taxon.index(t) + 1]] = [[] for i in range(len(Taxon[taxon[taxon.index(t) + 1]]))]
			
			for i in range(n):
				# cycle through list of taxa to be placed into higher taxa
				TaxMax[taxon[taxon.index(t) + 1]][j].extend(TaxMax[t][i])

				# j is incremented by 1 as we cycle through higher taxa
				j += 1
				if j == len(Taxon[taxon[taxon.index(t) + 1]]):
					# reset j to zero when last of higher taxa has been reached to start the cycle over again
					j =0

	taxon.reverse(); TDmax = 0; n = 0; N = len(sample)
	
	for t in taxon:
		taxonN[t] = sum([len(TaxMax[t][i]) * len(TaxMax[t][j]) for i in range(len(TaxMax[t])) for j in range(len(TaxMax[t])) if i != j])
		n = taxonN[t] - n
		TDmax += n * coef[t]
		n = taxonN[t]

	TDmax /= (N * (N-1))
	EI = (TDmax-atd)/(TDmax-TDmin)

	Eresults = {'EI':EI, 'TDmin':TDmin,'TDmax':TDmax}
	return Eresults

def printResults(args, results, popStats, population, taxon, coef):
	if args.out is not None and args.out == 'y': args.out.split('.') + '.out'
	else: outfile = args.samplefile.name.split('.')[0] + '.out'
	print("Writing to " + outfile + "\n")
	o = open(outfile,'w')	
	o.write("Output from Average Taxonomic Distinctness\n\n")
	o.write ("Number of taxa and path lengths for each taxonomic level:")

	for t in popStats['taxon']:
		l = '%-10s\t%d\t%.4f' %(t,popStats['taxonPopN'][t],popStats['pathlengths'][t])
		o.write (l)
	
	o.write ("\n")
	for f in results:
		print ("---------------------------------------------------")
		print ("Results for sample: ", f,'\n')
		print ("Dimension for this sample is", results[f]['n'], '\n\n')
		print ("Number of taxa and pairwise comparisons	at each taxon level:")
			
		n = 0
		for t in popStats['taxon']:
			N = results[f]['N'][t] - n
			print('{0:20s}\t{1:d}\t{2:d}'.format(t,results[f]['TaxonN'][t],N))
			
		print ("\nNumber of pairwise comparisons is for pairs that differ \
at each level excluding comparisons that differ at upper levels")
		print ("\n")
		  
		print ("Average taxonomic distinctness	    = %.4f" % results[f]['atd'])
		print ("Variation in taxonomic distinctness = %.4f" % results[f]['vtd'])
		print ("Minimum taxonomic distinctness	    = %.4f" % results[f]['euler']['TDmin'])
		print ("Maximum taxonomic distinctness	    = %.4f" % results[f]['euler']['TDmax'])
		print ("von Euler's index of imbalance	    = %.4f" % results[f]['euler']['EI'])
		print ('\n')

	o.close()
	if args.ci == 'y':
		ConfienceIntervals(args, outfile, population, taxon, coef)
		
def ConfienceIntervals(args, outfile, population, taxon, coef):		
	output = outfile.split('_')[0] + '_funnel.out'
	p = args.p; d1 = args.d1; d2 = args.d2
	o = open(output,'w')

	print ("""Confidence limits for average taxonomic distinctness and variation \
in taxonomic distinctness limits are lower 95% limit for AvTD and upper 95% limit \
for VarTD\n""")

	print ("Number of permutations for confidence limits =", p, '\n')

	pop = list(population.keys())
	up =[]; lo=[]; means=[]

	dims = [d for d in range(d1, d2+1)]

	for d in range(d1, d2 + 1):	
		AvTDci = []; VarTDci = []
		for j in range(p):
			rsamp = random.sample(pop,d)
			atd,taxonN, Taxon = ATDmean(rsamp, population, taxon, coef); AvTDci.append(atd)
			vtd = ATDvariance(taxonN,rsamp,atd, taxon, coef); VarTDci.append(vtd)

		AvTDci.sort()
		VarTDci.sort()

		AvTD = AvTDci[int(.05 * p)], sum(AvTDci)/p, AvTDci[int(.95 * p)], max(AvTDci)
		VarTD = min(VarTDci), VarTDci[int(.05 * p)],sum(VarTDci)/p,VarTDci[int(.95 * p)] 
		
		l = '{0:10s} {1:s}  {2:s}  {3:s}  {4:s}  {5:s}  {6:s}  {7:s}  {8:s}' \
		 .format("dimension", "AvTD05%", "AvTDmean","AvTD95%", "AvTDup", 
		"VarTDlow","VarTD05%", "VarTDmean", "VarTD95%")
		o.write(l)
		l = '{0:<10d} {1:6.4f}  {2:6.4f}  {3:6.4f}  {4:6.4f} {5:6.4f} {6:6.4f} {7:6.4f} {8:6.4f}' \
		.format(d, AvTD[0], AvTD[1], AvTD[2], AvTD[3], VarTD[0], VarTD[1], VarTD[2], VarTD[3])
		#o.write ('%i		 %6.4f	 %6.4f	 %6.4f	 %6.4f	 %6.4f	 %6.4f	 %6.4f	 %6.4f' \
		o.write(l)
		
	