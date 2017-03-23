#!/usr/bin/env python

import sys
import os
from random import *
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(prog="PhyloSample", prefix_chars='-+',
	description='choose mode of sampling', )
	
parser.add_argument('popfile', type=argparse.FileType('r'),
	default=sys.stdin)
parser.add_argument('outfile')
	
parser.add_argument("-m", "--mode", dest="mode")

###-------these are used if the mode is automated-------###
parser.add_argument("-n", "--number", dest="number", type=int)	#number of new files
parser.add_argument("-t", "--taxon", dest="taxon")	#taxon level
parser.add_argument("-s", "--split", dest="split",type=int)	#split
parser.add_argument("-e", "--merge", dest="merge",type=int)	#merge
parser.add_argument("-v", "--move", dest="move", type=int)	#move

args = parser.parse_args()

popfile=args.popfile
outfile=args.outfile
mode = args.mode

f = popfile.readlines()
header = f[0]

# ta is list of taxon levels for which we have information
ta = [taxon for taxon in header.split()]
ta.remove('Taxon:')


def getDataStructures():
	"""Obtains the necessary data structures for use in this program:
	1) TA = dictionary of taxonomic ranks containing taxonomic dictionaries of 
		taxon found at that level (key) with lists of next lower rank taxon
	2) species = similar to TA, but with lists of species for each taxon group
	3) master = dictionary of all species with dictionaries of taxonomic ranks
	4) changes = dictionary of changes made in the program
	"""
	# TA is dictionary of taxon levels composed of dictionaries of taxon groups at that level with 
	# sub-taxon falling under each taxon group
	TA = {taxon: defaultdict(list) for taxon in ta[:-1]}

	# species contains dictionary of each taxon level as TA, but with each taxon group showing a list of species
	# rather than sub-taxon groups
	species = {taxon: defaultdict(list) for taxon in ta[:-1]}
	
	# dictionary af species with taxon listings for each taxon level
	master = {f[i].split()[0]: {ta[j]: f[i].split()[j+1] for j in range(len(ta))} for i in range(1,len(f))};

	#keeps track of the changes, i.e. what sort of change and what taxon group was involved
	#simulatedData lists the new simulated data
	changes = {taxon: {} for taxon in ta}

	#form dictionary of taxon with list of sub taxon at each taxon level

	for i in range(1,len(f)):
		x = f[i].split();
		x = x[1:]
		for j in range(len(ta)-1):
			t = ta[j]
			if x[j] == '/': continue
			TA[t][x[j]].append(x[j + 1])
			species[t][x[j]].append(x[-1])

	#original lists will contain duplicates, so eliminate them
	for taxon in TA:
		TA[taxon] = {t: set(TA[taxon][t]) for t in TA[taxon]}

	return TA, species, master, changes

def MergeTaxa(simPop, ch, t, k=2):
	"""merges two randomly selected taxon at rank t, makes necessary changes
	to master dictionary, and outputs dictionary of changes made plus list
	of species to be output to simulated dataset. All species found within
	the two taxon will be the contained within the new taxon
	"""
	while True:
		if len(TA[t]) < k: 
			#if there are less than k number of taxon at a particular level, merge is not possible
			message = "there are not enough originaltaxon at the level to merge"
			break

		else:
			S = sample(TA[t].keys(),k)
			taxon = '{0:s}-{1:s}'.format(S[0],S[1])
			n = len(species[t][S[0]]) + len(species[t][S[1]])
			desc = '{0:d} species moved to composite taxon {1:s}'.format(n, taxon)
			
			ch[t][S[0]] = ['merge', desc]; ch[t][S[1]] = ['merge', desc]
			simPop.extend(species[t][S[0]])
			simPop.extend(species[t][S[1]])
			
			
			for s in species[t][S[0]]:
				master[s][t] = taxon

			for s in species[t][S[1]]:
				master[s][t] = taxon

			del TA[t][S[0]]; del TA[t][S[1]]

			message = "Two taxon at %s level merged, original taxon removed from dataset" %t
			break
	
	return simPop, ch, message

def SplitTaxa(simPop, ch,t):
	"""splits one randomly selected taxon at rank into two new taxon, makes
	necessary changes to master dictionary, and outputs dictionary of changes 
	made plus list of populations to be output to simulated dataset. all species
	from original species will now be divided into two new species
	"""
	while True:
		if len(TA[t].keys()) == 0:
			message = "there are no original taxon left at that level to split"
			break
		else:
			S = sample(TA[t].keys(),1)
			T = S[0]
			S = list(TA[t][T])
			n=0

			for j in range(0,len(S),2):
				taxon='{0:s}_1'.format(S[0])
				if ta.index(t) == len(ta) - 2:
					simPop.extend(species[t][S[j]])	
					n+=len(species[t][S[j]])
					for s in species[t][S[j]]:
						master[s][t] = taxon
				else:
					simPop.extend(species[ta[ta.index(t)+1]][S[j]])	
					for s in species[ta[ta.index(t)+1]][S[j]]:
						n+=1
						master[s][t] = taxon
								
			n=0
			for j in range(1,len(S),2):
				taxon='{0:s}_2'.format(S[0])
				n+=len(species[ta[ta.index(t)+1]][S[j]])
				if ta.index(t) == len(ta) - 2:
					for s in species[t][S[j]]:
						master[s][t] = taxon
				else:
					n+=len(species[ta[ta.index(t)+1]][S[j]])
					for s in species[ta[ta.index(t)+1]][S[j]]:
						master[s][t] = taxon

			desc = '{0:d} species split between {1:s}_1 and {1:s}_2'.format(n, taxon)	

			ch[t][T] = ['split', desc];
			del TA[t][T];			
			message = "Single taxon at %s level split, original taxon removed from dataset" % t
			break
	return simPop, ch,message

def MoveTaxa(simPop, ch, t):
	"""moves one randomly selected taxon at rank t+1 from a randomly selected
	taxon at rank t into another randomly selected taxon at rank t, makes 
	necessary changes to master dictionary, and outputs dictionary of changes 
	made plus list of populations to be output to simulated dataset. all species
	from the subtaxon at rank t+1 will be moved form taxon 1 to taxon 2
	"""
	#if ta.index(t) == 0: print ("cannot move taxa within uppermost level")
	#this function moves a subtaxon group within upper taxon level
	while True:
		if len(TA[t].keys()) < 2:
			message = "less than two taxon at this level, no separate source and destination populations"
			break			
		else:
			T = sample(TA[t].keys(),2)
			if TA[t][T[0]] == '/': continue
			if len(TA[t][T[0]]) == 0:
				message = "there are no original taxon left at that level to move"
				break
			S = sample(TA[t][T[0]],1)[0]
			simPop.extend(species[t][T[0]])
			if ta.index(t) == len(ta) - 2:				
				n = len(species[ta[ta.index(t)]][S])
				for s in species[ta[ta.index(t)]][S]:
					master[s][t]=T[1]
					species[t][T[0]].remove(s)

			else:
				n = len(species[ta[ta.index(t) + 1]][S])
				for s in species[ta[ta.index(t) + 1]][S]:
					master[s][t]=T[1]
					species[t][T[0]].remove(s)
			
			TA[t][T[0]].remove(S)
			ch[t][S] = ['move', '{0:d} species from {1:s} to {2:s}'.format(n, T[0], T[1])]
			message = "Single taxon at %s level moved to another upper taxon" % t
			break
		
	return simPop,ch, message

def taxonSimulation(changes):
	"""controls simulation. prompts user to taxon change to perform, taxon tank.
	and number of repetitions. Continues until user chooses to end simulations.
	"""
	t = "".join([str(ta.index(t)) + " for " + t + " " for t in ta[0:len(ta) -1]])
	simPop=[]
	while True:
		prompt = 'what function would you like to perform? input m for move, e for merge, or s for split: '
		function = str(input(prompt))

		prompt = 'please input taxon level to manipulate: %s ' %t
		taxon = ta[int(input(prompt))]
		
		prompt = 'please input number of repetitions: ' 
		n = int(input(prompt))
		
		if function == 'm':
			for r in range(n): 
				simPop, changes,message = MoveTaxa(simPop, changes,taxon)
				print(message + '\n')
		elif function == 'e':
			#prompt = 'how many taxon would you like to merge '
			#K = int(input(prompt))
			for r in range(n): 
				simPop, changes, message = MergeTaxa(simPop, changes,taxon, k=2)
				print(message + '\n')	
		elif function == 's':
			for r in range(n):
				simPop, changes, message = SplitTaxa(simPop, changes,taxon)
				print(message + '\n')
		else: print ('that is not a legal option, please input again m for move, e for merge, or s for split: ')

		prompt = 'would you like to continue? please input y for yes or n for no: '
		i = (input(prompt))
		
		if i == 'n':
			break
	return simPop,changes

def printChanges(master, simPop, changes, out, dir='.'):
	"""prints out dataset with simulated species data and file of changes.
	"""
	out = dir + '/' + out.split(".")[0]
	o = open(out + '.sim','w')
	
	o.write("\t".join(header.split()))

	for m in simPop:
		dat = [m]
		dat.extend([master[m][t] for t in ta])

		l = '\t'.join(dat)
		o.write(l)
		o.write('\n')
	
	o.close()
	
	o = open(out + '.changes','w')
	l = 'taxon level\ttaxon\tchange\tdescription\n'
	o.write(l)
	
	for i in range(len(ta)-1):
		t=ta[i]
		for c in changes[t]:
			taxon = c
			change = changes[t][c][0]
			desc = changes[t][c][1]
			l = '{0:s}\t{1:s}\t{2:s}\t{3:s}\n'.format(t,taxon,change,desc)
			o.write(l)
	o.close()
	
if mode == 'i':
	"""if mode is interactive, that is user will be asked for changes to 
	be performed
	"""
	doctext = "PhyloSample will ask a number of prompts, beginning with the function you wish to use---merge, \
	split, or move, followed by the taxon level to be manipulted. Choices of taxon for manipulation include all \
	the species level. When you choose to merge or split taxon, taxon at the level chosen will be changed; however \
	if you choose to move a taxon it is actually a taxon at the taxon level below the chosen one that will be moved. \
	Your chosen taxon level is the one within which a taxon will moved. For example, choosing 'genus' will cause a \
	species to be moved from one taxon to another. Finally, you will be asked to choose number of reps to be \
	performed. Following simulations, you will be asked whether you want to continue with simulations. Choosing \
	'yes' allows you to choose other functions to be performed and taxon levels to be maniuplated. Choosing 'no' \
	outputs new simulated file along with a file of changes performed."

	print (doctext)
	TA, species, master, changes = getDataStructures()
	simulatedPopulation, changes = taxonSimulation(changes)

	print ("printing new master list and changes file\n")
	printChanges(master, simulatedPopulation, changes, outfile)

if mode == 'a':
	"""if mode is automated. Used to perform more than one simulation at once.
	limitation is that you can only perform at one taxon rank for each 
	simulation.
	"""
	fdir = '%sSimulations' % outfile.split()[0]
	#fc = '%sMasterChanges' % outfile
	try:
		os.mkdir(f)
		os.mkdir(fc)
	except:
		print ('directory exists')

	Pdict = {}
	Pdict['n'] = args.number
	Pdict['t'] = args.taxon
	Pdict['s'] = args.split
	Pdict['e'] = args.merge
	Pdict['v'] = args.move

	for n in range(Pdict['n']):
		TA, species, master, changes = getDataStructures()
		taxon =Pdict['t']
		simPop=[]

		for i in range(max(Pdict['v'],Pdict['e'],Pdict['s'])):
			if Pdict['v'] <= i:	simPop, changes, message = MoveTaxa(simPop, changes,taxon)
			if Pdict['e'] <= i:	simPop, changes, message = MergeTaxa(simPop, changes,taxon, k=2)
			if Pdict['s'] <= i:	simPop, changes, message = SplitTaxa(simPop, changes,taxon)
		outfile = outfile.split('.')[0]
		out= '{0:s}-{1:d}'.format(outfile,n+1)

		printChanges(master, simPop, changes, out, dir=fdir)

