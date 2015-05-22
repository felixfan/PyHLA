#!/usr/bin/env python

import numpy as np
import pandas as pd
from collections import Counter
import scipy.stats
import sys

def keyDicts(dict1, dict2):
	k1 = dict1.keys()
	k1.extend(dict2.keys())
	return set(k1)
def getGenes(alleles):
	genes = []
	for allele in alleles:
		genes.append(allele.split('*')[0])
	return sorted(set(genes))
def readAAseq(aafile):
	'''
	read allele sequence
	'''
	seq = {}
	f = open(aafile)
	for i in f:
		if i.strip():
			i = i.strip()
			arr = i.split()
			if len(arr) == 2:
				seq[arr[0]] = arr[1]
	f.close()
	return seq
def readGeno(infile):
	'''
	read geno data
	return allele for each individual, number of individual with each allele
	'''
	case = {}         # number of cases take one allele
	ctrl = {}            # number of ctrls take one allele
	caseGeno = []  # alleles for each case
	ctrlGeno = []     # alleles for each ctrl
	ncase = 0       # number of case
	nctrl = 0         # number of ctrl
	f = open(infile)
	for i in f:
		if i.strip():
			i = i.strip()
			arr = i.split()
			if arr[1] == '1':
				nctrl += 1
				ctrlGeno.append(arr[2:])
				for j in range(2, len(arr), 2):
					k = j + 1
					if arr[j] in ctrl:
						ctrl[arr[j]] += 1
					else:
						ctrl[arr[j]] = 1
					if arr[k] != arr[j]:
						if arr[k] in ctrl:
							ctrl[arr[k]] += 1
						else:
							ctrl[arr[k]] = 1
			elif arr[1] == '2':
				ncase += 1
				caseGeno.append(arr[2:])
				for j in range(2, len(arr), 2):
					k = j + 1
					if arr[j] in case:
						case[arr[j]] += 1
					else:
						case[arr[j]] = 1
					if arr[k] != arr[j]:
						if arr[k] in case:
							case[arr[k]] += 1
						else:
							case[arr[k]] = 1
	f.close()
	return case, ctrl, caseGeno, ctrlGeno, ncase, nctrl
def consensusSeq(seqs):
	'''
	input: list of amino acid sequences
	output: one consensus sequences (len = shortest one)
	'''
	conSeq = ''
	if len(seqs) == 1:
		conSeq = seqs[0]
	else:
		for i in range(0, len(seqs[0])):
			flag = 1
			lenflag = 0               # end of shortest seq
			for j in range(1, len(seqs)):
				if i == len(seqs[j]):
					lenflag = 1
					break
				elif seqs[0][i] == '*':
					conSeq += '*'
					flag = 0
					break
				elif seqs[0][i] != seqs[j][i]:
					conSeq += '*'
					flag = 0
					break
			if lenflag == 1:
				break
			elif flag == 1:
				conSeq += seqs[0][i]
	return conSeq
def getSeq(alleles, seq, consensus=True):
	'''
	find the consensus seq for each allele
	'''
	aseq = {}
	for allele in alleles:
		if allele in seq:
			aseq[allele] = seq[allele]
		else:
			if consensus:
				tmp = []
				for k in sorted(seq.keys()):
					if allele in k:
						tmp.append(seq[k])
				if len(tmp):
					aseq[allele] = consensusSeq(tmp)
				else:
					print 'no sequence is avaiable for: %s' % allele
			else: # use the first allele seq
				for k in sorted(seq.keys()):
					if allele in k:
						aseq[allele] = seq[k]
						break
				if allele not in aseq:
					print 'no sequence is avaiable for: %s' % allele
	return aseq
def aaAlign(case, ctrl, seq):
	'''
	amino acid alignment
	'''
	aln = {}
	alleles = keyDicts(case, ctrl)
	genes = getGenes(alleles)
	for gene in genes:
		tmp = []
		for allele in alleles:
			if allele.startswith(gene):
				if allele in seq:
					tmp.append([allele, seq[allele]])
		aln[gene] = tmp
	return aln
def convertID(oldID, allID = 'Allelelist_history.txt'):
	'''
	convert allele id to Release 3.20.0, 2015-04-17
	'''
	ids = {}
	f = open(allID)
	for line in f:
		if line.strip():
			line = line.strip()
			arr = line.split()
			for i in range(1, len(arr)):
				if arr[i] != 'NA':
					ids[arr[i]] = arr[0]
			if arr[1] != 'NA':
				ids[arr[0]] = arr[1]
	f.close()
	newid = {}
	for i in oldID:
		newid[i] = 'NA'
		if i in ids:
			if ids[i] in ids:
				newid[i] = ids[ids[i]]
	return newid
def deltaCal(case, ctrl, ncase, nctrl):
	'''
	calculate dalta
	'''
	alleles = keyDicts(case, ctrl)
	delta = {}
	for allele in alleles:
		if allele not in case:
			case[allele] = 0
		if allele not in ctrl:
			ctrl[allele] = 0
		delta[allele] = 1.0 * case[allele] / ncase - 1.0 * ctrl[allele] / nctrl
	return delta
def aaCount(Geno, seq, gene, pos):
	'''
	count aa at a locus
	'''
	nn = {}
	gpa = {} #gene_pos_aa
	for item in Geno:
		for j in range(0,len(item),2):
			k = j + 1
			if item[j].startswith(gene) and item[k].startswith(gene):
				lastaa = ''
				if item[j] in seq: # exact same id
					if len(seq[item[j]]) > pos:
						if seq[item[j]][pos] in nn:
							nn[seq[item[j]][pos]] += 1
							gpa[seq[item[j]][pos]].append(item[j])
						else:
							nn[seq[item[j]][pos]] = 1
							gpa[seq[item[j]][pos]] = [item[j]]
						lastaa = seq[item[j]][pos]
				if item[j] != item[k]:
					if item[k] in seq:
						if len(seq[item[k]]) > pos:
							if lastaa != seq[item[k]][pos]:
								if seq[item[k]][pos] in nn:
									nn[seq[item[k]][pos]] += 1
									gpa[seq[item[k]][pos]].append(item[k])
								else:
									nn[seq[item[k]][pos]] = 1
									gpa[seq[item[k]][pos]] = [item[k]]
	return nn, gpa
def aaAssoc(case, ctrl, caseGeno, ctrlGeno, ncase, nctrl, seq, test='fisher'):
	'''
	amino acid association
	'''
	assoc = {}
	alleles = keyDicts(case,ctrl)
	genes = getGenes(alleles)
	aln = aaAlign(case, ctrl, seq)
	for gene in genes:
		df = pd.DataFrame(aln[gene],columns=['allele','seq'])
		df12 = df['seq'].apply(lambda x: pd.Series([i for i in list(x)]))
		col12 = df12.shape[1]
		for i in range(0, col12):
			t12 = Counter(df12[i])
			t12.pop(np.nan, None)
			t12.pop('*', None)
			keys = t12.keys()
			if len(t12) > 1:
				nn, gpaP = aaCount(caseGeno, seq, gene, i)
				nnc, gpaC = aaCount(ctrlGeno, seq, gene, i)
				for key in keys:
					if key in nn:
						n1 = nn[key]
					else:
						n1 = 0
					n2 = ncase - n1
					if key in nnc:
						n3 = nnc[key]
					else:
						n3 = 0
					n4 = nctrl-n3
					data = [[n1,n2],[n3,n4]]
					if test == 'fisher':
						try:
							OR, pvalue = scipy.stats.fisher_exact(data)
						except:
							pvalue = 'NA'
					elif test == 'chisq':
						try:
							chi2, pvalue, dof, expected = scipy.stats.chi2_contingency(data)
						except:
							pvalue = 'NA'
					else:
						sys.exit("only 'fisher' or 'chisq' can be used for amino acid association!")
					OR = (n1+0.5) * (n4+0.5) / (n2+0.5) / (n3+0.5)
					# ALLELES HAS THIS AA
					awr = []
					if key in gpaP:
						awr.extend(gpaP[key])
					if key in gpaC:
						awr.extend(gpaC[key])
					awr = sorted(set(awr))
					assoc[(gene, i+1, key)]=[n1, n2, n3, n4, pvalue, OR, awr]
	return assoc
############################################################
def printAA(aln):
	'''
	print amino acid, 10 aa in each block, 5 blocks per line
	'''
	for g in sorted(aln.keys()):
		aadict = aln[g]
		length = len(aadict[0][1])
		for i in range(0, int(round((length+24)/50.0))):
			for k in aadict:
				for r in range(5):
					s = i*50+r*10
					e = s + 10
					if r == 0:
						print "%-20s%-11s" % (k[0], k[1][s:e]),
					elif r == 4:
						print "%-11s" % k[1][s:e],
						if e < length:
							print "%-5d" % e
						else:
							print "%-5d" % length
					else:
						print "%-11s" % k[1][s:e],
			print
def writeAA(aln, outfile):
	fw = open(outfile, 'w')
	for g in sorted(aln.keys()):
		aadict = aln[g]
		length = len(aadict[0][1])
		for i in range(0, int(round((length+24)/50.0))):
			for k in aadict:
				for r in range(5):
					s = i*50+r*10
					e = s + 10
					if r == 0:
						fw.write("%-20s%-11s" % (k[0], k[1][s:e]))
					elif r == 4:
						fw.write("%-11s" % k[1][s:e])
						if e < length:
							fw.write("%-5d\n" % e)
						else:
							fw.write("%-5d\n" % length)
					else:
						fw.write("%-11s" % k[1][s:e])
			fw.write('\n')
############################################################
def printAAA(assoc):
	print '%-20s' % 'ID',
	for h in ('A_case', 'B_case', 'A_ctrl', 'B_ctrl'):
		print '%8s' % h,
	print '%10s' % 'P',
	print '%8s' % 'OR',
	print '\t%s' % 'ACR'
	for k in sorted(assoc.keys()):
		print "%-20s" % (k[0] + '_' + str(k[1]) + '_' + k[2]),
		for i in range(4):
			print "%8d" % assoc[k][i],
		if assoc[k][4] == 'NA':
			print "%10s" % 'NA',
		elif assoc[k][4] > 0.001:
			print "%10.5f" % assoc[k][4],
		else:
			print "%10.2e" % assoc[k][4],
		print "%8.2f" % assoc[k][5],
		##
		awrs = ''
		si = 0
		for it in assoc[k][6]:
			if si > 0:
				awrs += (',' + it)
			else:
				awrs = it
			si += 1
		print '\t%s' % awrs
def writeAAA(assoc, outfile):
	fw = open(outfile, 'w')
	fw.write('%-20s' % 'ID')
	for h in ('A_case', 'B_case', 'A_ctrl', 'B_ctrl'):
		fw.write('%8s' % h)
	fw.write('%10s' % 'P')
	fw.write('%8s' % 'OR')
	fw.write('\t%s\n' % 'ACR')
	for k in sorted(assoc.keys()):
		fw.write("%-20s" % (k[0] + '_' + str(k[1]) + '_' + k[2]))
		for i in range(4):
			fw.write("%8d" % assoc[k][i])
		if assoc[k][4] == 'NA':
			fw.write("%10s" % assoc[k][4])
		elif assoc[k][4] > 0.001:
			fw.write("%10.5f" % assoc[k][4])
		else:
			fw.write("%10.2e" % assoc[k][4])
		fw.write("%8.2f" % assoc[k][5])
		##
		awrs = ''
		si = 0
		for it in assoc[k][6]:
			if si > 0:
				awrs += (',' + it)
			else:
				awrs = it
			si += 1
		fw.write("\t%s\n" % awrs)
	fw.close()
############################################################
