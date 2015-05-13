#!/usr/bin/env python

from __future__ import division
import scipy.stats
import math
import random
import HLAcount

def adjustP(pvalues, method = "Benjamini-Hochberg"):                
	"""                                                                                                   
	correct p-values for multiple testing
	methods: Bonferroni, Bonferroni-Holm or Holm, Benjamini-Hochberg or FDR
	"""
	n = len(pvalues)
	cp = [1]*n
	if method == "Bonferroni":
		cp = map(lambda x:min(x*n,1.0), pvalues)
	elif method == "Bonferroni-Holm" or method == "Holm":
		values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
		values = sorted(values)
		for rank, vals in enumerate(values):
			pvalue, i = vals
			cp[i] = (n-rank) * pvalue
		for rank, vals in enumerate(values):
			pvalue, i = vals                                                      
			if rank > 0:
				cp[i] = min(1.0, max(cp[i], cp[j]))
			else:
				cp[i] = min(1.0, cp[i])
			j = i
	elif method == "Benjamini-Hochberg" or method == "FDR":
		values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
		values = sorted(values,reverse=True)
		for rank, vals in enumerate(values):
			pvalue, i = vals
			cp[i] = n * pvalue / (n-rank)
		for rank, vals in enumerate(values):
			pvalue, i = vals
			if rank > 0:
				cp[i] = min(1.0, min(cp[i], cp[j]))
			else:
				cp[i] = min(1.0, cp[i])
			j = i
	elif method == "Benjamini-Yekutieli" or method == "FDR_BY":
		q = 0
		for i in range(1,n+1):
			q += 1.0 / i
		values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
		values = sorted(values,reverse=True)
		for rank, vals in enumerate(values):
			pvalue, i = vals
			cp[i] = q * pvalue * n/(n-rank)
		for rank, vals in enumerate(values):
			pvalue, i = vals
			if rank > 0:
				cp[i] = min(1.0, min(cp[i], cp[j]))
			else:
				cp[i] = min(1.0, cp[i])
			j = i
	return cp
def assocADRChiFisher(infile, digit, freq, test='chisq', model = 'allelic', adjust='FDR', exclude=None, perm=None, seed=None):
	'''
	Association Analysis for Allelic, Dominant or Recessive Model
	Pearson's Chi-squared test or Fisher exact test
	return allele counts, frequency, [chi-square, df,] p, and OR, adjustedP, permutationP
	'''
	if model == 'allelic':
		caseAlleles, ctrlAlleles, np, nc, nn = HLAcount.allelicCount(infile, digit)
	elif model == 'dom':
		caseAlleles, ctrlAlleles, np, nc, nn = HLAcount.domCount(infile, digit)
	elif model == 'rec':
		caseAlleles, ctrlAlleles, np, nc, nn = HLAcount.recCount(infile, digit)
	alleleFreq, alleles = HLAcount.hlaFreq(caseAlleles, ctrlAlleles, np, nc, nn)

	excludeAlleles =[]
	if exclude is not None:
		ef = open(exclude)
		for line in ef:
			line = line.strip()
			excludeAlleles.append(line)
		ef.close()

	usedAllele = []
	assoc = {}
	for allele in caseAlleles.keys():
		if allele in ctrlAlleles:
			if allele not in excludeAlleles:
				if (alleleFreq[allele][2]) > freq:
					usedAllele.append(allele)
					n1 = caseAlleles[allele]
					n2 = np[allele.split('*')[0]] - n1
					n3 = ctrlAlleles[allele]
					n4 = nc[allele.split('*')[0]] - n3
					data = [[n1, n2], [n3, n4]]
					if test == "chisq":
						chi2, p, dof, expected = scipy.stats.chi2_contingency(data)
					OR, pvalue = scipy.stats.fisher_exact(data)
					se = math.sqrt(1.0/n1  + 1.0/n2 +  1.0/n3 + 1.0/n4)
					l95 = math.exp(math.log(OR) - 1.96 * se)
					u95 = math.exp(math.log(OR) + 1.96 * se)
					ss = []
					ss.append(n1)
					ss.append(n2)
					ss.append(n3)
					ss.append(n4)
					ss.append(alleleFreq[allele][0])
					ss.append(alleleFreq[allele][1])
					ss.append(alleleFreq[allele][2])	
					if test == "chisq":
						ss.append(p)
						ss.append(chi2)
						ss.append(dof)
					elif test == "fisher":
						ss.append(pvalue)
					ss.append(OR)
					ss.append(l95)
					ss.append(u95)
					assoc[allele] = ss
	### adjust
	genes = []
	for g in np.keys():
		if g in nc.keys():
			genes.append(g)
	for g in sorted(genes):      ### GENE BY GENE
		ps = []
		ns = []
		for a in assoc:
			if a.startswith(g) and assoc[a][7] != 'NA':  # p value at 7 col start from 0
				ps.append(assoc[a][7])
				ns.append(a)
		cp = adjustP(ps,adjust)
		for i in range(len(ns)):
			if assoc[ns[i]][7] != 'NA':
				assoc[ns[i]].append(cp[i])
			else:
				assoc[ns[i]].append('NA')

	if perm is None:
		return assoc
	else:
		random.seed(seed)
		permP = {}       # perm p value
		permN = {}      # perm p < orig p
		permNL = {}    # perm p > orig p
		permNA = {}    # perm NA
		for a in assoc:
			permNA[a] = 0
			permNL[a] = 0
			permN[a] = 0
		if perm > 10:
			pf = perm / 10
		else:
			pf = 2
		pn = 0   # effect perm number
		while True:
			if model == 'allelic':
				case9, ctrl9, np9, nc9, nn9 = HLAcount.allelicCount(infile,digit, True)
			elif model == 'dom':
				case9, ctrl9, np9, nc9, nn9 = HLAcount.domCount(infile,digit, True)
			elif model == 'rec':
				case9, ctrl9, np9, nc9, nn9 = HLAcount.recCount(infile,digit, True)
			ca = []   # current alleles
			for a in case9:
				if a in ctrl9:
					ca.append(a)
			if set(usedAllele) <= set(ca):  # current alleles contain used allele
	        		for a in usedAllele:
					n1 = case9[a]
					n2 = np9[a.split('*')[0]] - n1
					n3 = ctrl9[a]
					n4 = nc9[a.split('*')[0]] - n3
					data = [[n1, n2], [n3, n4]]
					if test == "chisq":
						chi2, p, dof, expected = scipy.stats.chi2_contingency(data)
					elif test == 'fisher':
						OR, p = scipy.stats.fisher_exact(data)
					if not isinstance(p, float):
						permNA[a] += 1
					else:
						if assoc[a][7] == 'NA':
							permNA[a] += 1
						else:
							if p < assoc[a][7]:
								permN[a] += 1
							else:
								permNL[a] += 1
				pn += 1
				if pn % pf == 1:
					print 'permutation {}/{} ...'.format(pn, perm)
			if pn == perm:
				break
		for a in assoc:
			if assoc[a][7] == 'NA':
				permP[a] = 'NA'
			else:
				if permNA[a] == perm:
					permP[a] = 'NA'
				else:
					permP[a] = 1.0 * (permN[a] + 1) / (perm + 1 - permNA[a])
		return assoc, permP, permN, permNA, permNL
def assocRaw(infile, digit, freq, exclude=None, perm=None, seed=None):
	'''
	Association Analysis (2 x m)
	Pearson's Chi-squared test
	return chi-square, df, p
	'''
	caseAlleles, ctrlAlleles, np, nc, nn = HLAcount.allelicCount(infile, digit)
	alleleFreq, alleles = HLAcount.hlaFreq(caseAlleles, ctrlAlleles, np, nc, nn)
	assoc = {}
	usedAllele = []

	excludeAlleles =[]
	if exclude is not None:
		ef = open(exclude)
		for line in ef:
			line = line.strip()
			excludeAlleles.append(line)
		ef.close()

	### genes
	gene = []  # get all genes name
	for a in np:
		if a in nc:
			gene.append(a)

	for g in gene:
		n1 = []
		n2 = []
		for a in caseAlleles:
			if a in ctrlAlleles:
				if a not in excludeAlleles:
					if a.startswith(g):
						if alleleFreq[a][2] > freq:
							n1.append(caseAlleles[a])
							n2.append(ctrlAlleles[a])
							usedAllele.append(a)
		data = [n1, n2]
		chi2, p, dof, expected = scipy.stats.chi2_contingency(data)
		ss = []
		if not isinstance(chi2, float):
			ss.append('NA')
		else:
			ss.append(chi2)
		if not isinstance(dof, int):
			ss.append('NA')
		else:
			ss.append(dof)
		if not isinstance(p, float):
			ss.append('NA')
		else:
			ss.append(p)
		assoc[g] = ss
	if perm is None:
		return assoc, usedAllele
	else:
		random.seed(seed)
		permP = {}       # perm p value
		permN = {}      # perm p < orig p
		permNL = {}    # perm p > orig p
		permNA = {}    # perm NA
		for a in assoc:
			permNA[a] = 0
			permNL[a] = 0
			permN[a] = 0
		if perm > 10:
			pf = perm / 10
		else:
			pf = 2
		pn = 0   # effect perm number
		while True:
			case9, ctrl9, np9, nc9, nn9 = HLAcount.allelicCount(infile,digit, True)
			ca = []   # current alleles
			for a in case9:
				if a in ctrl9:
					ca.append(a)
			if set(usedAllele) <= set(ca):  # current alleles contain used allele
	        		for g in gene:
					n1 = []
					n2 = []
					for a in case9:
						if a in ctrl9:
							if a.startswith(g):
								if a in usedAllele:
									n1.append(case9[a])
									n2.append(ctrl9[a])

					data = [n1, n2]
					chi2, p, dof, expected = scipy.stats.chi2_contingency(data)
					if not isinstance(p, float):
						permNA[g] += 1
					else:
						if assoc[g][2] == 'NA':
							permNA[g] += 1
						else:
							if p < assoc[g][2]:
								permN[g] += 1
							else:
								permNL[g] += 1
				pn += 1
				if pn % pf == 1:
					print 'permutation {}/{} ...'.format(pn, perm)
			if pn == perm:
				break
		for a in assoc:
			if assoc[a][2] == 'NA':
				permP[a] = 'NA'
			else:
				if permNA[a] == perm:
					permP[a] = 'NA'
				else:
					permP[a] = 1.0 * (permN[a] + 1) / (perm + 1 - permNA[a])
		return assoc, usedAllele, permP, permN, permNA, permNL
def assocScoreU(infile, digit, freq, exclude=None, perm=None, seed=None):
	'''
	Association Analysis (2 x m)
	Score test
	return score test U
	'''
	caseAlleles, ctrlAlleles, np, nc, nn = HLAcount.allelicCount(infile, digit)
	alleleFreq, alleles = HLAcount.hlaFreq(caseAlleles, ctrlAlleles, np, nc, nn)
	assoc = {}
	usedAllele = []

	excludeAlleles =[]
	if exclude is not None:
		ef = open(exclude)
		for line in ef:
			line = line.strip()
			excludeAlleles.append(line)
		ef.close()

	### genes
	gene = []  # get all genes name
	for a in np:
		if a in nc:
			gene.append(a)

	for g in gene:
		u = 0
		n1 = np[g]
		for a in caseAlleles:
			if a in ctrlAlleles:
				if a not in excludeAlleles:
					if a.startswith(g):
						if alleleFreq[a][2] > freq:
							usedAllele.append(a)
							u = u + (caseAlleles[a] - n1 * alleleFreq[a][2]) ** 2 / alleleFreq[a][2] - (caseAlleles[a] - n1 * alleleFreq[a][2]) / alleleFreq[a][2]
		if not isinstance(u, float):
			u = 'NA'
		assoc[g] = u

	if perm is None:
		return assoc, usedAllele
	else:
		random.seed(seed)
		permP = {}       # perm p value
		permN = {}      # perm p < orig p
		permNL = {}    # perm p > orig p
		permNA = {}    # perm NA
		for a in assoc:
			permNA[a] = 0
			permNL[a] = 0
			permN[a] = 0
		if perm > 10:
			pf = perm / 10
		else:
			pf = 2
		pn = 0   # effect perm number
		while True:
			case9, ctrl9, np9, nc9, nn9 = HLAcount.allelicCount(infile,digit, True)
			alleleFreq9, alleles9 = HLAcount.hlaFreq(caseAlleles, ctrlAlleles, np, nc, nn)
			ca = []   # current alleles
			for a in case9:
				if a in ctrl9:
					ca.append(a)
			if set(usedAllele) <= set(ca):  # current alleles contain used allele
	        		for g in gene:
					u = 0
					n1 = np[g]
					for a in case9:
						if a in ctrl9:
							if a.startswith(g):
								if a in usedAllele:
									u = u + (case9[a] - n1 * alleleFreq9[a][2]) ** 2 / alleleFreq9[a][2] - (case9[a] - n1 * alleleFreq9[a][2]) / alleleFreq9[a][2]
					if not isinstance(u, float):
						permNA[g] += 1
					else:
						if assoc[g] == 'NA':
							permNA[g] += 1
						else:
							if u > assoc[g]:
								permN[g] += 1
							else:
								permNL[g] += 1
				pn += 1
				if pn % pf == 1:
					print 'permutation {}/{} ...'.format(pn, perm)
			if pn == perm:
				break
		for a in assoc:
			if assoc[a] == 'NA':
				permP[a] = 'NA'
			else:
				if permNA[a] == perm:
					permP[a] = 'NA'
				else:
					permP[a] = 1.0 * (permN[a] + 1) / (perm + 1 - permNA[a])
		return assoc, usedAllele, permP, permN, permNA, permNL
#################################################################
def assocDelta(infile, digit, freq=0.05, adjust='FDR', exclude=None, perm=None, seed=None):
	popCase, popCtrl, popP, popC, popN = HLAcount.domCount(infile, digit)
	caseAlleles99, ctrlAlleles99, np99, nc99, nn99 = HLAcount.allelicCount(infile, digit)
	alleleFreq, alleles = HLAcount.hlaFreq(caseAlleles99, ctrlAlleles99, np99, nc99, nn99)

	assoc ={}
	usedAllele = []

	excludeAlleles =[]
	if exclude is not None:
		ef = open(exclude)
		for line in ef:
			line = line.strip()
			excludeAlleles.append(line)
		ef.close()

	for allele in alleles:
		if allele not in excludeAlleles:
			if alleleFreq[allele] > freq:
				usedAllele.append(allele)
				np = popP[allele.split('*')[0]]
				nc = popC[allele.split('*')[0]]
				if allele in popCase:
					n1 = popCase[allele]	
				else:
					n1 = 0
				tfp = 1.0 * n1 / np
				n2 = np - n1
				if allele in popCtrl:
					n3 =  popCtrl[allele]
				else:
					n3 = 0
				tfc = 1.0 * n3 / nc
				n4 = nc - n3
				delta = tfp - tfc
				data = [[n1, n2], [n3, n4]]
				OR, p = scipy.stats.fisher_exact(data)
				OR = (n1 + 0.5) * (n4 + 0.5) / (n2 + 0.5) / (n3 + 0.5)
				tmp = []
				tmp.append(delta)
				tmp.append(p)
				tmp.append(OR)
				assoc[allele] = tmp

	### adjust
	genes = []
	for g in popP.keys():
		if g in popC.keys():
			genes.append(g)
	for g in sorted(genes):      ### GENE BY GENE
		ps = []
		ns = []
		for a in assoc:
			if a.startswith(g) and assoc[a][1] != 'NA':  # p value at 7 col start from 0
				ps.append(assoc[a][1])
				ns.append(a)
		cp = adjustP(ps,adjust)
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
		permN = {}      # perm delta > orig delta
		for a in assoc:
			permN[a] = 0
		if perm > 10:
			pf = perm / 10
		else:
			pf = 2
		pn = 0   # effect perm number
		while True:
			case9, ctrl9, np9, nc9, nn9 = HLAcount.domCount(infile,digit, True)
        		for allele in usedAllele:
        			np = np9[allele.split('*')[0]]
				nc = nc9[allele.split('*')[0]]
				if allele in case9:
					n1 = case9[allele]	
				else:
					n1 = 0
				tfp = 1.0 * n1 / np
				n2 = np - n1
				if allele in ctrl9:
					n3 =  ctrl9[allele]
				else:
					n3 = 0
				tfc = 1.0 * n3 / nc
				n4 = nc - n3
				delta = tfp - tfc
				if delta > assoc[allele][0]:
					permN[allele] += 1
			pn += 1
			if pn % pf == 1:
				print 'permutation {}/{} ...'.format(pn, perm)
			if pn == perm:
				break
		for a in assoc:
			permP = 1.0 * (1 + permN[a]) / (1 + perm)
			assoc[a].append(permP)
		return assoc