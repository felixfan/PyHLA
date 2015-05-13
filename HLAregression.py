#!/usr/bin/env python

import pandas as pd
import statsmodels.formula.api as smf
import math
import os
import sys
import time
import string
import random
import HLAassoc

def getAlleles(infile, digits):
	'''
	get all alleles for each gene
	return a dictionary, keys are the gene names, values are the alleles
	return second dictionary, keys are the gene names, values are the start column
	'''
	geneAlleles = {}
	geneCol = {}
	f = open(infile)
	for line in f:
		line = line.rstrip()
		alleles = line.split()
		for i in range(2,len(alleles)):
			if alleles[i] != 'NA':
				names = alleles[i].split(":")
				gene = names[0].split("*")
				if i not in geneCol:
					geneCol[i] = gene[0]
				if gene[0] not in geneAlleles:
					geneAlleles[gene[0]] = []
				if digits == 4:
					temp = names[0] + ":" + names[1]
					if temp not in geneAlleles[gene[0]]:
						geneAlleles[gene[0]].append(temp)
				elif digits == 2:
					temp = names[0]
					if temp not in geneAlleles[gene[0]]:
						geneAlleles[gene[0]].append(temp)
	f.close()
	return geneAlleles, geneCol
def allelicRecode(infile, digits, test):
	'''
	allele dosage coding
	Assume A*01:01 is the test allele, then A*01:01 A*01:01 is code as 2
	A*01:01 A*01:02 is code as 1, and A*01:02 A*01:03 is code as 0
	return a list contains the alleles
	return a dictionary contains the coding for alleles
	'''
	geneAlleles, geneCol = getAlleles(infile, digits)
	header = ['IID','PHT']
	ans = {}
	f = open(infile)
	for line in f:
		line = line.rstrip()
		alleles = line.split()
		if test == 'logistic':
			ans[alleles[0]] = [alleles[0],int(alleles[1])-1]
		elif test == 'linear':
			ans[alleles[0]] = [alleles[0],alleles[1]]
		for i in range(2,len(alleles),2):
			j = i + 1
			if alleles[i] != 'NA' and alleles[j] != 'NA':
				gene1 = alleles[i].split('*')[0]
				gene2 = alleles[j].split('*')[0]
				if gene1 == gene2:
					if digits == 4:
						allele1 = alleles[i].split(':')[0] + ':' + alleles[i].split(':')[1]
						allele2 = alleles[j].split(':')[0] + ':' + alleles[j].split(':')[1]
					else:
						allele1 = alleles[i].split(':')[0]
						allele2 = alleles[j].split(':')[0]
					gAlleles = sorted(geneAlleles[gene1])
					for ga in gAlleles:
						if ga not in header:
							header.append(ga)
						if allele1 == ga and allele2 == ga:
							ans[alleles[0]].append(2)
						elif allele1 == ga or allele2 == ga:
							ans[alleles[0]].append(1)
						else:
							ans[alleles[0]].append(0)
				else:
					sys.exit("input format is wrong!")
			else:
				for gg in geneAlleles[geneCol[i]]:
					ans[alleles[0]].append('NA')
	return ans,header
def writeRecode(infile, digits, test):
	'''
	write coding to a temp file
	return the temp file name
	'''
	tmp = time.strftime("%H%M%S%d%b%Y")
	tmp = tmp + '.txt'
	f = open(tmp,'w')
	ans, header = allelicRecode(infile, digits, test)
	for i in header:
		i = string.replace(i, '*', '_') # change A*01:01 to A_01_01
		i = string.replace(i, ':', '_')
		f.write("%12s" % i,)
	f.write('\n')
	for i in ans:
		for j in ans[i]:
			f.write("%12s" % j,)
		f.write('\n')
	f.close()
	return tmp
##########################################
def regressionLogistic(infile, digits, freq, adjust = 'FDR', exclude=None, covfile=None, covname=None, perm=None, seed=None, test='logistic'):
	'''
	logistitic regression
	output: dictionary, key: allele, value: statistic includes count, freq, p and OR
	'''
	### geno
	tfile = writeRecode(infile, digits, test)
	geno = pd.read_csv(tfile,delim_whitespace= True, header = 0)
	os.remove(tfile)
	alleles = list(geno.columns.values)[2:]
	### exclude alleles
	excludeAlleles =[]
	if exclude is not None:
		ef = open(exclude)
		for line in ef:
			line = line.strip()
			excludeAlleles.append(line)
		ef.close()
	### covar file
	if covfile is not None:
		cov = pd.read_csv(covfile,delim_whitespace= True, header = 0)
		covindex = list(cov.columns.values)[1:]
		### covar name
		if covname is None:                           # default: use all covariants
			covname = covindex
		else:
			covname = covname.split(',')
	### assoc
	usedAllele = []
	assoc = {}
	for allele in alleles:
		if allele not in excludeAlleles:
			n1 = int(geno[allele].where(geno['PHT'] == 1).sum(axis=0)) # case
			n2 = int(geno[allele].where(geno['PHT'] == 1).count()) * 2 # case
			n3 = int(geno[allele].where(geno['PHT'] == 0).sum(axis=0)) # control
			n4 = int(geno[allele].where(geno['PHT'] == 0).count()) * 2 # control
			f1 = 1.0 * n1 / n2
			f2 = 1.0 * n3 / n4
			f12 = 1.0 * (n1 + n3) / (n2 + n4)
			if f12 > freq:
				usedAllele.append(allele)
				n2 -= n1
				n4 -= n3
				myformula = 'PHT ~ ' + allele
				if covfile is None:
					mydata = geno.ix[:, ['IID', 'PHT', allele]]
				else:
					for name in covname:
						if name in covindex:
							myformula = myformula + ' + ' + name
						else:
							print 'can not find covariant name ' + '"' + name + '" in covariant file'
							sys.exit()
					geno9 = geno.ix[:, ['IID', 'PHT', allele]]
					mydata = pd.merge(geno9, cov, on='IID', how='inner')
				try:
					lr = smf.logit(formula = myformula, data = mydata).fit(maxiter=100, disp=False)
					p = lr.pvalues[1]
					try:
						OR = math.exp(lr.params[1])
					except:
						OR = 'NA'
					try:
						L95 = math.exp(lr.conf_int()[0][1])
					except:
						L95 = 'NA'
					try:
						U95 = math.exp(lr.conf_int()[1][1])
					except:
						U95 = 'NA'
				except:
					p = 'NA'
					OR = 'NA'
					L95 = 'NA'
					U95 = 'NA'
				aname = allele.split('_')
				nname = aname[0] + '*' + aname[1]
				if digits == 4:
					nname = nname + ':' + aname[2]
				elif digits == 6:
					nname = nname + ':' + aname[2] + ':' + aname[3]
				ss = []
				ss.append(n1)
				ss.append(n2)
				ss.append(n3)
				ss.append(n4)
				ss.append(f1)
				ss.append(f2)
				ss.append(f12)
				if p != 'NA':
					ss.append(p)
				else:
					ss.append('NA')	
				if OR != 'NA':
					ss.append(OR)
				else:
					ss.append('NA')
				if L95 != 'NA':
					ss.append(L95)
				else:
					ss.append('NA')
				if U95 != 'NA':
					ss.append(U95)
				else:
					ss.append('NA')
				assoc[nname] = ss
	### adjust
	genes = []
	for a in assoc:
		genes.append(a.split('*')[0])
	genes = set(genes)
	for g in sorted(genes):      ### GENE BY GENE
		ps = [] # p value
		ns = [] # allele name
		for a in assoc:
			if a.startswith(g) and assoc[a][7] != 'NA':  # p value at 7 col start from 0
				ps.append(assoc[a][7])
				ns.append(a)
		cp = HLAassoc.adjustP(ps,adjust)
		for i in range(len(ns)):
			if assoc[ns[i]][7] != 'NA':
				assoc[ns[i]].append(cp[i])
			else:
				assoc[ns[i]].append('NA')
	### perm
	if perm is None:
		return assoc
	else:
		random.seed(seed)
		permP = {}
		permN = {}
		permNL = {}
		permNA = {}
		for a in assoc:
			permNA[a] = 0
			permNL[a] = 0
			permN[a] = 0
		if perm > 10:
			pf = perm / 10
		else:
			pf = 2
		for i in range(perm):
			if i % pf == 1:
				print 'permutation {}/{} ...'.format(i, perm)
			xxx = geno['PHT'].values.flatten()
			random.shuffle(xxx)
			geno['PHT'] = xxx
			for allele in alleles:
				aname = allele.split('_')
				nname = aname[0] + '*' + aname[1]
				if digits == 4:
					nname = nname + ':' + aname[2]
				elif digits == 6:
					nname = nname + ':' + aname[2] + ':' + aname[3]
				if nname in assoc:
					myformula = 'PHT ~ ' + allele
					if covfile is None:
						mydata = geno.ix[:, ['IID', 'PHT', allele]]
					else:
						for name in covname:
							myformula = myformula + ' + ' + name
						geno9 = geno.ix[:, ['IID', 'PHT', allele]]
						mydata = pd.merge(geno9, cov, on='IID', how='inner')
					try:
						lr = smf.logit(formula = myformula, data = geno).fit(maxiter=100, disp=False)
						p = lr.pvalues[1]	
					except:
						p = 'NA'
					if p == 'NA':
						permNA[nname] += 1
					else:
						if assoc[nname][7] == 'NA':
							permNA[nname] += 1
						else:
							if p < assoc[nname][7]:
								permN[nname] += 1
							else:
								permNL[nname] += 1
		for a in assoc:
			if assoc[a][7] == 'NA':
				permP[a] = 'NA'
			else:
				if permNA[a] == perm:
					permP[a] = 'NA'
				else:
					permP[a] = 1.0 * (permN[a] + 1) / (perm + 1 - permNA[a])
	return assoc, permP, permN, permNA
##################################################################
def regressionLinear(infile, digits, freq, adjust = 'FDR', exclude=None, covfile=None, covname=None, perm=None, seed=None, test='linear'):
	'''
	linear regression with covariants
	output: dictionary, key: allele, value: statistic includes freq, p and beta
	'''
	### geno
	tfile = writeRecode(infile, digits, test)
	geno = pd.read_csv(tfile,delim_whitespace= True, header = 0)
	os.remove(tfile)
	alleles = list(geno.columns.values)[2:]
	### exclude alleles
	excludeAlleles =[]
	if exclude is not None:
		ef = open(exclude)
		for line in ef:
			line = line.strip()
			excludeAlleles.append(line)
		ef.close()
	### covar file
	if covfile is not None:
		cov = pd.read_csv(covfile,delim_whitespace= True, header = 0)
		covindex = list(cov.columns.values)[1:]
		### covar name
		if covname is None:                           # default: use all covariants
			covname = covindex
		else:
			covname = covname.split(',')
	### assoc
	usedAllele = []
	assoc = {}
	for allele in alleles:
		if allele not in excludeAlleles:
			n1 = int(geno[allele].sum(axis=0)) # number of allele
			n2 = int(geno[allele].count()) * 2 # total allele
			f12 = 1.0 * n1 / n2
			if f12 > freq:
				usedAllele.append(allele)
				myformula = 'PHT ~ ' + allele
				if covfile is None:
					mydata = geno.ix[:, ['IID', 'PHT', allele]]
				else:
					geno9 = geno.ix[:, ['IID', 'PHT', allele]]
					mydata = pd.merge(geno9, cov, on='IID', how='inner')
					for name in covname:
						if name in covindex:
							myformula = myformula + ' + ' + name
						else:
							print 'can not find covariant name ' + '"' + name + '" in covariant file'
							sys.exit()
				try:
					lr = smf.ols(formula = myformula, data = mydata).fit(maxiter=100, disp=False)
					p = lr.pvalues[1]
					beta = lr.params[1]
					L95 = lr.conf_int()[0][1]
					U95 = lr.conf_int()[1][1]
				except:
					p = 'NA'
					beta = 'NA'
					L95 = 'NA'
					U95 = 'NA'
				aname = allele.split('_')
				nname = aname[0] + '*' + aname[1]
				if digits == 4:
					nname = nname + ':' + aname[2]
				elif digits == 6:
					nname = nname + ':' + aname[2] + ':' + aname[3]
				ss = []
				ss.append(f12)
				if p != 'NA':
					ss.append(p)
				else:
					ss.append('NA')
				if beta != 'NA':
					ss.append(beta)
				else:
					ss.append('NA')
				if L95 != 'NA':
					ss.append(L95)
				else:
					ss.append('NA')
				if U95 != 'NA':
					ss.append(U95)
				else:
					ss.append('NA')
				assoc[nname] = ss
	### adjust
	genes = []
	for a in assoc:
		genes.append(a.split('*')[0])
	genes = set(genes)
	for g in sorted(genes):      ### GENE BY GENE
		ps = [] # p value
		ns = [] # allele name
		for a in assoc:
			if a.startswith(g) and assoc[a][1] != 'NA':  # p value at 1 col start from 0
				ps.append(assoc[a][1])
				ns.append(a)
		cp = HLAassoc.adjustP(ps,adjust)
		for i in range(len(ns)):
			if assoc[ns[i]][1] != 'NA':
				assoc[ns[i]].append(cp[i])
			else:
				assoc[ns[i]].append('NA')

	### perm
	if perm is None:
		return assoc
	else:
		random.seed(seed)
		permP = {}
		permN = {}
		permNL = {}
		permNA = {}
		for a in assoc:
			permNA[a] = 0
			permNL[a] = 0
			permN[a] = 0
		if perm > 10:
			pf = perm / 10
		else:
			pf = 2
		for i in range(perm):
			if i % pf == 1:
				print 'permutation {}/{} ...'.format(i, perm)
			xxx = geno['PHT'].values.flatten()
			random.shuffle(xxx)
			geno['PHT'] = xxx
			for allele in alleles:
				aname = allele.split('_')
				nname = aname[0] + '*' + aname[1]
				if digits == 4:
					nname = nname + ':' + aname[2]
				elif digits == 6:
					nname = nname + ':' + aname[2] + ':' + aname[3]
				if nname in assoc:
					myformula = 'PHT ~ ' + allele
					if covfile is None:
						mydata = geno.ix[:, ['IID', 'PHT', allele]]
					else:
						for name in covname:
							myformula = myformula + ' + ' + name
						geno9 = geno.ix[:, ['IID', 'PHT', allele]]
						mydata = pd.merge(geno9, cov, on='IID', how='inner')
					try:
						lr = smf.ols(formula = myformula, data = mydata).fit(maxiter=100, disp=False)
						p = lr.pvalues[1]
					except:
						p = 'NA'
					if p == 'NA':
						permNA[nname] += 1
					else:
						if assoc[nname] == 'NA':
							permNA[nname] += 1
						else:
							if p < assoc[nname][1]:
								permN[nname] += 1
							else:
								permNL[nname] += 1
		for a in assoc:
			if assoc[a] == 'NA':
				permP[a] = 'NA'
			else:
				if permNA[a] == perm:
					permP[a] = 'NA'
				else:
					permP[a] = 1.0 * (permN[a] + 1) / (perm + 1 - permNA[a])
	return assoc, permP, permN, permNA
