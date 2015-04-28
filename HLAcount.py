#!/usr/bin/env python
import sys
import random

def allelicCount(infile, digits, perm=False):
	'''
	count all alleles
	return the count of each allele in case and control
	'''
	case = {}      # counts for each allele
	control = {}
	np = {}         # total non-NA alleles for each gene
	nc = {}

	f = open(infile) # read all phenotype
	pht = []
	for line in f:
		line = line.rstrip()
		pht.append(line.split()[1])
	f.close()
	if perm:
		random.shuffle(pht)
	l = 0     # line index
	nn = len(pht)

	f = open(infile)
	for line in f:
		line = line.rstrip()
		alleles = line.split()
		for i in range(2,len(alleles),2):
			j = i + 1
			if alleles[i] != 'NA' and alleles[j] != 'NA':
				names1 = alleles[i].split(":")
				names2 = alleles[j].split(":")
				if digits == 6:
					if len(names1) < 3 or len(names2) < 3:
						sys.exit("--digits 6 requires at least 6 digits resolution genotype!")
					a4d1 = names1[0] + ":" + names1[1] + ":" + names1[2]
					a4d2 = names2[0] + ":" + names2[1] + ":" + names2[2]
				if digits == 4:
					if len(names1) < 2 or len(names2) < 2:
						sys.exit("--digits 4 requires at least 4 digits resolution genotype!")
					a4d1 = names1[0] + ":" + names1[1]
					a4d2 = names2[0] + ":" + names2[1]
				elif digits == 2:
					if len(names1) < 1 or len(names2) < 1:
						sys.exit("--digits 2 requires at least 2 digits resolution genotype!")
					a4d1 = names1[0]
					a4d2 = names2[0]

				if pht[l] == "2":
					if a4d1.split('*')[0] in np:
						np[a4d1.split('*')[0]] += 1
					else:
						np[a4d1.split('*')[0]] = 1

					if a4d1 in case:
						case[a4d1] += 1
					else:
						case[a4d1] = 1

					if a4d1 != a4d2:
						if a4d2 in case:
							case[a4d2] += 1
						else:
							case[a4d2] = 1
						

				elif pht[l] == "1":
					if a4d1.split('*')[0] in nc:
						nc[a4d1.split('*')[0]] += 1
					else:
						nc[a4d1.split('*')[0]] = 1

					if a4d1 in control:
						control[a4d1] += 1
					else:
						control[a4d1] = 1

					if a4d1 != a4d2:
						if a4d2 in control:
							control[a4d2] += 1
						else:
							control[a4d2] = 1
		l += 1
	f.close()
	return case, control, np, nc, nn

def domCount(infile, digits, perm=False):
	'''
	count all alleles
	return the count of each allele in case and control
	'''
	case = {}     # counts for each allele
	control = {}
	nc = {}
	np = {}

	f = open(infile) # read all phenotype
	pht = []
	for line in f:
		line = line.rstrip()
		pht.append(line.split()[1])
	f.close()
	if perm:
		random.shuffle(pht)
	l = 0     # line index
	nn = len(pht)

	f = open(infile)
	for line in f:
		line = line.rstrip()
		alleles = line.split()
		for i in range(2,len(alleles),2):
			j = i + 1
			if alleles[i] != 'NA' and alleles[j] != 'NA':
				names1 = alleles[i].split(":")
				names2 = alleles[j].split(":")
				if digits == 6:
					if len(names1) < 3 or len(names2) < 3:
						sys.exit("--digits 6 requires at least 6 digits resolution genotype!")
					a4d1 = names1[0] + ":" + names1[1] + ":" + names1[2]
					a4d2 = names2[0] + ":" + names2[1] + ":" + names2[2]
				if digits == 4:
					if len(names1) < 2 or len(names2) < 2:
						sys.exit("--digits 4 requires at least 4 digits resolution genotype!")
					a4d1 = names1[0] + ":" + names1[1]
					a4d2 = names2[0] + ":" + names2[1]
				elif digits == 2:
					if len(names1) < 1 or len(names2) < 1:
						sys.exit("--digits 2 requires at least 2 digits resolution genotype!")
					a4d1 = names1[0]
					a4d2 = names2[0]


				if pht[l] == "2":
					if a4d1.split('*')[0] in np:
						np[a4d1.split('*')[0]] += 1
					else:
						np[a4d1.split('*')[0]] = 1

					if a4d1 in case:
						case[a4d1] += 1
					else:
						case[a4d1] = 1

					if a4d1 != a4d2:
						if a4d2 in case:
							case[a4d2] += 1
						else:
							case[a4d2] = 1
						

				elif pht[l] == "1":
					if a4d1.split('*')[0] in nc:
						nc[a4d1.split('*')[0]] += 1
					else:
						nc[a4d1.split('*')[0]] = 1

					if a4d1 in control:
						control[a4d1] += 1
					else:
						control[a4d1] = 1

					if a4d1 != a4d2:
						if a4d2 in control:
							control[a4d2] += 1
						else:
							control[a4d2] = 1
		l += 1
	f.close()
	return case, control, np, nc, nn

def recCount(infile, digits, perm=False):
	'''
	count all alleles
	return the count of each allele in case and control
	'''
	case = {}     # counts for each allele
	control = {}
	nc = {}
	np = {}

	f = open(infile) # read all phenotype
	pht = []
	for line in f:
		line = line.rstrip()
		pht.append(line.split()[1])
	f.close()
	if perm:
		random.shuffle(pht)
	l = 0     # line index
	nn = len(pht)

	f = open(infile)
	for line in f:
		line = line.rstrip()
		alleles = line.split()
		for i in range(2,len(alleles),2):
			j = i + 1
			if alleles[i] != 'NA' and alleles[j] != 'NA':
				names1 = alleles[i].split(":")
				names2 = alleles[j].split(":")
				if digits == 6:
					if len(names1) < 3 or len(names2) < 3:
						sys.exit("--digits 6 requires at least 6 digits resolution genotype!")
					a4d1 = names1[0] + ":" + names1[1] + ":" + names1[2]
					a4d2 = names2[0] + ":" + names2[1] + ":" + names2[2]
				if digits == 4:
					if len(names1) < 2 or len(names2) < 2:
						sys.exit("--digits 4 requires at least 4 digits resolution genotype!")
					a4d1 = names1[0] + ":" + names1[1]
					a4d2 = names2[0] + ":" + names2[1]
				elif digits == 2:
					if len(names1) < 1 or len(names2) < 1:
						sys.exit("--digits 2 requires at least 2 digits resolution genotype!")
					a4d1 = names1[0]
					a4d2 = names2[0]


				if pht[l] == "2":
					if a4d1.split('*')[0] in np:
						np[a4d1.split('*')[0]] += 1
					else:
						np[a4d1.split('*')[0]] = 1

					if a4d1 == a4d2:
						if a4d2 in case:
							case[a4d2] += 1
						else:
							case[a4d2] = 1

				elif pht[l] == "1":
					if a4d1.split('*')[0] in nc:
						nc[a4d1.split('*')[0]] += 1
					else:
						nc[a4d1.split('*')[0]] = 1

					if a4d1 == a4d2:
						if a4d2 in control:
							control[a4d2] += 1
						else:
							control[a4d2] = 1
		l += 1
	f.close()
	return case, control, np, nc, nn
def hlaFreq(caseAlleles, ctrlAlleles, np, nc, nn):
	'''
	Frequency of all alleles
	return the frequency of each allele in case, control, and both
	'''
	temp = []
	temp.extend(caseAlleles.keys())
	temp.extend(ctrlAlleles.keys())
	alleles = sorted(set(temp))
	freq = {}
	for allele in alleles:
		ff = []
		if allele in caseAlleles:
			ff.append(1.0 * caseAlleles[allele] / np[allele.split('*')[0]])
		else:
			ff.append(0)
		if allele in ctrlAlleles:	
			ff.append(1.0 * ctrlAlleles[allele] / nc[allele.split('*')[0]])
		else:
			ff.append(0)
		if allele in caseAlleles and allele in ctrlAlleles:
			ff.append(1.0 * (ctrlAlleles[allele] + caseAlleles[allele]) / (np[allele.split('*')[0]] + nc[allele.split('*')[0]]))
		elif allele in caseAlleles:
			ff.append(1.0 * caseAlleles[allele] / np[allele.split('*')[0]])
		else:
			ff.append(1.0 * ctrlAlleles[allele] / nc[allele.split('*')[0]])
		freq[allele] = ff
	return freq, alleles
def alleleCount(infile, digits):
	'''
	count all alleles for quantitative traits
	return the count of each allele in case and control
	'''
	allelesN = {}
	genesN = {}
	N = 0

	f = open(infile)
	for line in f:
		line = line.rstrip()
		alleles = line.split()
		N += 1
		for i in range(2,len(alleles),2):
			j = i + 1
			if alleles[i] != 'NA' and alleles[j] != 'NA':
				names1 = alleles[i].split(":")
				names2 = alleles[j].split(":")
				if digits == 6:
					if len(names1) < 3 or len(names2) < 3:
						sys.exit("--digits 6 requires at least 6 digits resolution genotype!")
					a4d1 = names1[0] + ":" + names1[1] + ":" + names1[2]
					a4d2 = names2[0] + ":" + names2[1] + ":" + names2[2]
				if digits == 4:
					if len(names1) < 2 or len(names2) < 2:
						sys.exit("--digits 4 requires at least 4 digits resolution genotype!")
					a4d1 = names1[0] + ":" + names1[1]
					a4d2 = names2[0] + ":" + names2[1]
				elif digits == 2:
					if len(names1) < 1 or len(names2) < 1:
						sys.exit("--digits 2 requires at least 2 digits resolution genotype!")
					a4d1 = names1[0]
					a4d2 = names2[0]
				if a4d1.split('*')[0] in genesN:
					genesN[a4d1.split('*')[0]] += 1
				else:
					genesN[a4d1.split('*')[0]] = 1
				if a4d2.split('*')[0] in genesN:
					genesN[a4d2.split('*')[0]] += 1
				else:
					genesN[a4d2.split('*')[0]] = 1
				if a4d1 in allelesN:
					allelesN[a4d1] += 1
				else:
					allelesN[a4d1] = 1
				if a4d2 in allelesN:
					allelesN[a4d2] += 1
				else:
					allelesN[a4d2] = 1
					
	f.close()
	return allelesN, genesN, N
def quantTrait(infile):
	f = open(infile)
	for line in f:
		line = line.rstrip()
		alleles = line.split()
		if alleles[1] != '1' and alleles[1] != '2' and alleles[1] != 'NA':
			return True
	f.close()
	return False
