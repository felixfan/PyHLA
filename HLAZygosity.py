#!/usr/bin/env python

import sys
import HLAInteraction

def countAAzyg(Geno, seq, gene, pos, aa):
	'''
	count of homozygous, heterozygous, and absent of amino acid
	'''
	pp, pn, nn = (0,0,0)
	for item in Geno:
		flag1 = 0
		flag2 = 0
		for j in range(0,len(item),2):
			k = j + 1
			if item[j].startswith(gene) and item[k].startswith(gene):
				if item[j] in seq: 
					if len(seq[item[j]]) > pos:
						if seq[item[j]][pos] == aa:
							flag1 = 1
				if item[j] != item[k]:
					if item[k] in seq:
						if len(seq[item[k]]) > pos:
							if seq[item[k]][pos] ==aa:
								flag2 =1
				else:
					flag2 = flag1
				if flag1 + flag2 == 2:
					pp += 1
				elif flag1 + flag2 == 1:
					pn += 1
				else:
					nn += 1
				break
	return [pp, pn, nn]
def threeTests(x, y, test = 'fisher'):
	'''
	zygosity test
	'''
	n1 = n2 = n3 = n4 = 0
	# test 1
	n1 = x[0]
	n2 = x[1]
	n3 = y[0]
	n4 = y[1]
	data = [[n1, n2], [n3, n4]]
	p1 = HLAInteraction.simpleTest(data, test)
	# test2
	n1 = x[2]
	n2 = x[1]
	n3 = y[2]
	n4 = y[1]
	data = [[n1, n2], [n3, n4]]
	p2 = HLAInteraction.simpleTest(data, test)
	# test 3
	n1 = x[0]
	n2 = x[2]
	n3 = y[0]
	n4 = y[2]
	data = [[n1, n2], [n3, n4]]
	p3 = HLAInteraction.simpleTest(data, test)
	return [p1,p2,p3]
def zygosityAA(keys, caseGeno, ctrlGeno, myseq, test):
	ans = {}
	for k in keys:
		gene = k[0]
		pos = k[1] - 1
		aa = k[2]
		x = countAAzyg(caseGeno, myseq, gene, pos, aa)
		y = countAAzyg(ctrlGeno, myseq, gene, pos, aa)
		ps = threeTests(x, y, test)
		ans[k] = ps
	return ans
####################################################
def countAlleleZyg(Geno, allele):
	'''
	count of homozygous, heterozygous, and absent of alleles
	'''
	pp, pn, nn = (0,0,0)
	for item in Geno:
		for j in range(0,len(item),2):
			k = j + 1
			if item[j].startswith(allele.split('*')[0]) and item[k].startswith(allele.split('*')[0]):
				flag1 = 0
				flag2 = 0
				if item[j] == allele:
					flag1 = 1
				if item[k] == allele:
					flag2 = 1
				if flag1 + flag2 == 2:
					pp += 1
				elif flag1 + flag2 == 1:
					pn += 1
				else:
					nn += 1
				break
	return [pp, pn, nn]
def zygosityAllele(keys, caseGeno, ctrlGeno, test):
	ans = {}
	for k in keys:
		x = countAlleleZyg(caseGeno, k)
		y = countAlleleZyg(ctrlGeno, k)
		ps = threeTests(x, y, test)
		ans[k] = ps
	return ans
####################################################
