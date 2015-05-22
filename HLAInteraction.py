#!/usr/bin/env python

import scipy.stats

def simpleTest(data, test = 'fisher'):
	'''
	help function
	'''
	pvalue = 'NA'
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
	return pvalue
def tenTests(x, y, test = 'fisher'):
	'''
	strongest association
	'''
	n1 = n2 = n3 = n4 = 0
	# test 1
	n1 = x[0] + x[1]
	n2 = x[2] + x[3]
	n3 = y[0] + y[1]
	n4 = y[2] + y[3]
	data = [[n1, n2], [n3, n4]]
	p1 = simpleTest(data, test)
	# test2
	n1 = x[0] + x[2]
	n2 = x[1] + x[3]
	n3 = y[0] + y[2]
	n4 = y[1] + y[3]
	data = [[n1, n2], [n3, n4]]
	p2 = simpleTest(data, test)
	# test 3
	n1 = x[0]
	n2 = x[2]
	n3 = y[0]
	n4 = y[2]
	data = [[n1, n2], [n3, n4]]
	p3 = simpleTest(data, test)
	# test 4
	n1 = x[1]
	n2 = x[3]
	n3 = y[1]
	n4 = y[3]
	data = [[n1, n2], [n3, n4]]
	p4 = simpleTest(data, test)
	# test5
	n1 = x[0]
	n2 = x[1]
	n3 = y[0]
	n4 = y[1]
	data = [[n1, n2], [n3, n4]]
	p5 = simpleTest(data, test)
	# test6
	n1 = x[2]
	n2 = x[3]
	n3 = y[2]
	n4 = y[3]
	data = [[n1, n2], [n3, n4]]
	p6 = simpleTest(data, test)
	# test7
	n1 = x[1]
	n2 = x[2]
	n3 = y[1]
	n4 = y[2]
	data = [[n1, n2], [n3, n4]]
	p7 = simpleTest(data, test)
	# test8
	n1 = x[0]
	n2 = x[3]
	n3 = y[0]
	n4 = y[3]
	data = [[n1, n2], [n3, n4]]
	p8 = simpleTest(data, test)
	# test9
	n1 = x[0]
	n2 = x[1]
	n3 = x[2]
	n4 = x[3]
	data = [[n1, n2], [n3, n4]]
	p9 = simpleTest(data, test)
	# test10
	n1 = y[0]
	n2 = y[1]
	n3 = y[2]
	n4 = y[3]
	data = [[n1, n2], [n3, n4]]
	p10 = simpleTest(data, test)
	return [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10]
def readAlleleZygInteract(infile, digits):
	caseGeno = []
	ctrlGeno = []
	f = open(infile)
	for line in f:
		tpg = []
		tcg = []
		line = line.strip()
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
			else:
				a4d1 = 'NA'
				a4d2 = 'NA'
			if alleles[1] == '1':
				tcg.append(a4d1)
				tcg.append(a4d2)
			elif alleles[1] == '2':
				tpg.append(a4d1)
				tpg.append(a4d2)
		if tcg:
			ctrlGeno.append(tcg)
		if tpg:
			caseGeno.append(tpg)
	return caseGeno, ctrlGeno
################################################################
def countAA(Geno, seq, gene1, gene2, pos1, pos2, aa1, aa2):
	'''
	count of factor one and factor two: ++, +-, -+, --
	'''
	pp, pn, np, nn = (0,0,0,0)
	for item in Geno:
		flag1 = 0
		flag2 = 0
		for j in range(0,len(item),2):
			k = j + 1
			if item[j].startswith(gene1) and item[k].startswith(gene1):
				if item[j] in seq: 
					if len(seq[item[j]]) > pos1:
						if seq[item[j]][pos1] == aa1:
							flag1 = 1
				if item[j] != item[k]:
					if item[k] in seq:
						if len(seq[item[k]]) > pos1:
							if seq[item[k]][pos1] ==aa1:
								flag1 =1
			if item[j].startswith(gene2) and item[k].startswith(gene2):
				if item[j] in seq: 
					if len(seq[item[j]]) > pos2:
						if seq[item[j]][pos2] == aa2:
							flag2 = 1
				if item[j] != item[k]:
					if item[k] in seq:
						if len(seq[item[k]]) > pos2:
							if seq[item[k]][pos2] ==aa2:
								flag2 =1
		if flag1 == 1 and flag2 == 1:
			pp += 1
		elif flag1 ==1 and flag2 == 0:
			pn += 1
		elif flag1 == 0 and flag2 == 1:
			np += 1
		else:
			nn += 1
	return [pp, pn, np, nn]
def interactAA(keys, caseGeno, ctrlGeno, myseq, test):
	ans = {}
	for k1 in range(0, len(keys)):
		for k2 in range(k1+1,len(keys)):
			if keys[k1][0] != keys[k2][0] or keys[k1][1] != keys[k2][1]: # differ gene or differ pos
				gene1 = keys[k1][0]
				gene2 = keys[k2][0]
				pos1 = keys[k1][1] - 1
				pos2 = keys[k2][1] - 1
				aa1 = keys[k1][2]
				aa2 = keys[k2][2]
				x = countAA(caseGeno, myseq, gene1, gene2, pos1, pos2, aa1, aa2)
				y = countAA(ctrlGeno, myseq, gene1, gene2, pos1, pos2, aa1, aa2)
				ps = tenTests(x, y, test)
				newkey = keys[k1] + keys[k2]
				ans[newkey] = ps
	return ans
####################################################
def countAlleleInteract(Geno, allele1, allele2):
	'''
	count of factor one and factor two: ++, +-, -+, --
	'''
	pp, pn, np, nn = (0,0,0,0)
	for item in Geno:
		flag1 = 0
		flag2 = 0
		for j in range(0,len(item)):
			if item[j] == allele1:
				flag1 = 1
			if item[j] == allele2:
				flag2 =1
		if flag1 == 1 and flag2 == 1:
			pp += 1
		elif flag1 ==1 and flag2 == 0:
			pn += 1
		elif flag1 == 0 and flag2 == 1:
			np += 1
		else:
			nn += 1
	return [pp, pn, np, nn]
def interactAllele(keys, caseGeno, ctrlGeno, test):
	ans = {}
	for k1 in range(0, len(keys)):
		for k2 in range(k1+1,len(keys)):
			if keys[k1].split('*')[0] != keys[k2].split('*')[0]: # differ gene
				allele1 = keys[k1]
				allele2 = keys[k2]
				x = countAlleleInteract(caseGeno, allele1, allele2)
				y = countAlleleInteract(ctrlGeno, allele1, allele2)
				ps = tenTests(x, y, test)
				newkey = (keys[k1], keys[k2])
				ans[newkey] = ps
	return ans
####################################################
def printInteract(assoc, level):
	header = ('ID1', 'ID2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9','P10')
	for h in header:
		if h == 'ID1' or h == 'ID2':
			print "%-12s" % h,
		else:
			print "%12s" % h,
	print

	for k in sorted(assoc.keys()):
		if level == 'residue':
			print "%-12s" % (k[0]+'_'+str(k[1])+'_'+k[2]), 
			print "%-12s" % (k[3]+'_'+str(k[4])+'_'+k[5]),
		else:
			print "%-12s"  % k[0], 
			print "%-12s" % k[1],
		for i in range(2, len(assoc[k])):
			if assoc[k][i] == 'NA':
				print "%12s" % "NA",
			elif assoc[k][i] > 0.001:
				print "%12.4f" % assoc[k][i],
			else:
				print "%12.2e" % assoc[k][i],
		print
def writeInteract(assoc, level, outfile):
	fw = open(outfile, 'w')
	header = ('ID1', 'ID2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9','P10')
	for h in header:
		if h == 'ID1' or h == 'ID2':
			fw.write("%-12s" % h)
		else:
			fw.write("%12s" % h)
	fw.write('\n')

	for k in sorted(assoc.keys()):
		if level == 'residue':
			fw.write("%-12s" % (k[0]+'_'+str(k[1])+'_'+k[2]))
			fw.write("%-12s" % (k[3]+'_'+str(k[4])+'_'+k[5]))
		else:
			fw.write("%-12s"  % k[0])
			fw.write("%-12s" % k[1])
		for i in range(2, len(assoc[k])):
			if assoc[k][i] == 'NA':
				fw.write("%12s" % "NA")
			elif assoc[k][i] > 0.001:
				fw.write("%12.4f" % assoc[k][i])
			else:
				fw.write("%12.2e" % assoc[k][i])
		fw.write('\n')
	fw.close()