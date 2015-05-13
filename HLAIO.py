#!/usr/bin/env python

def printSummary(alleles, freq, caseAlleles, ctrlAlleles, np, nc, nn, popCase, popCtrl, popP, popC):
	print 'Sample size: %d'  % nn[0]
	print 'Number of cases: %d'  % nn[1]
	print 'Number of controls: %d'  % nn[2]
	print
	print 'Gene level summary'
	print '------------------------------------------------------------------'
	print '%12s%12s%12s%12s' % ('Gene', 'CaseCount', 'CtrlCount', 'TotalCount')
	for g in sorted(np.keys()):
		print '%12s%12d%12d%12d' % (g, np[g], nc[g], np[g] + nc[g])
	print
	print 'Allele level summary'
	print '------------------------------------------------------------------'
	print "%20s" % 'Allele',
	for a in ('CaseCount',  'CtrlCount', 'TotalCount', 'CaseFreq', 'CtrlFreq',  'TotalFreq'):
		print "%12s" % a,
	print
	for allele in alleles:
		print  "%20s" % allele,
		ttt = 0
		if allele in caseAlleles:
			print "%12d" % caseAlleles[allele], 
			ttt = caseAlleles[allele]
		else:
			print "%12d" % 0,
		if allele in ctrlAlleles:	
			print  "%12d" % ctrlAlleles[allele],
			print "%12d" % (ctrlAlleles[allele] + ttt),
		else:
			print "%12d" % 0,
			print "%12d" % ttt,
		for f in freq[allele]:
			print "%12.4f" % f,
		print
	print '\nPopulation level summary'
	print '------------------------------------------------------------------'
	print "%20s" % 'Allele',
	for a in ('popCaseCount',  'popCaseFreq', 'popCtrlCount', 'popCtrlFreq'):
		print "%14s" % a,
	print
	for allele in alleles:
		print  "%20s" % allele,
		if allele in popCase:
			print "%14d" % popCase[allele],
			print "%14.4f" % (1.0 * popCase[allele] / popP[allele.split('*')[0]]),
		else:
			print "%14d" % 0,
			print "%14d" % 0,
		if allele in popCtrl:
			print "%14d" % popCtrl[allele],
			print "%14.4f" % (1.0 * popCtrl[allele] / popC[allele.split('*')[0]])
		else:
			print "%14d" % 0,
			print "%14d" % 0
def writeSummary(alleles, freq, caseAlleles, ctrlAlleles, np, nc, nn, OUTFILE, popCase, popCtrl, popP, popC):
	fp = open(OUTFILE, 'w')
	fp.write('Sample size: %d\n'  % nn[0])
	fp.write('Number of cases: %d\n'  % nn[1])
	fp.write('Number of controls: %d\n\n'  % nn[2])

	fp.write('Gene level summary\n')
	fp.write('------------------------------------------------------------------\n')
	fp.write('%20s%12s%12s%12s\n' % ('Gene', 'CaseCount', 'CtrlCount', 'TotalCount'))
	for g in sorted(np.keys()):
		fp.write('%20s%12d%12d%12d\n' % (g, np[g], nc[g], np[g] + nc[g]))
	fp.write('\n')
	fp.write('Allele level summary\n')
	fp.write('------------------------------------------------------------------\n')
	fp.write('%20s%12s%12s%12s%12s%12s%12s\n' % ('Allele', 'CaseCount',  'CtrlCount', 'TotalCount', 'CaseFreq', 'CtrlFreq',  'TotalFreq'))
	for allele in alleles:
		fp.write("%20s" % allele)
		ttt = 0
		if allele in caseAlleles:
			fp.write("%12d" % caseAlleles[allele])
			ttt = caseAlleles[allele]
		else:
			fp.write("%12d" % 0)
		if allele in ctrlAlleles:
			fp.write("%12d" % ctrlAlleles[allele])
			fp.write("%12d" % (ctrlAlleles[allele] + ttt))
		else:
			fp.write("%12d" % 0)
			fp.write("%12d" % ttt)
		for f in freq[allele]:
			fp.write("%12.4f" % f)
		fp.write('\n')
	fp.write('\nPopulation level summary\n')
	fp.write('------------------------------------------------------------------\n')
	fp.write("%20s" % 'Allele')
	for a in ('popCaseCount',  'popCaseFreq', 'popCtrlCount', 'popCtrlFreq'):
		fp.write("%14s" % a)
	fp.write('\n')
	for allele in alleles:
		fp.write("%20s" % allele)
		if allele in popCase:
			fp.write("%14d" % popCase[allele])
			fp.write("%14.4f" % (1.0 * popCase[allele] / popP[allele.split('*')[0]]))
		else:
			fp.write("%14d" % 0)
			fp.write("%14d" % 0)
		if allele in popCtrl:
			fp.write("%14d" % popCtrl[allele])
			fp.write("%14.4f\n" % (1.0 * popCtrl[allele] / popC[allele.split('*')[0]]))
		else:
			fp.write("%14d" % 0)
			fp.write("%14d\n" % 0)
	fp.close()
###########################################################################
def printSummaryQuant(alleles, genes, n):
	print 'Sample size: %d'  % n
	print
	print 'Gene level summary'
	print '------------------------------------------------------------------'
	print '%12s%12s' % ('Gene',  'TotalCount')
	for g in sorted(genes.keys()):
		print '%12s%12d' % (g, genes[g])
	print
	print 'Allele level summary'
	print '------------------------------------------------------------------'
	print "%20s" % 'Allele',
	for a in ('TotalCount', 'TotalFreq'):
		print "%12s" % a,
	print
	for allele in sorted(alleles.keys()):
		print "%20s" % allele,
		print "%12d" % alleles[allele],
		print "%12.4f" % (1.0 * alleles[allele] / genes[allele.split('*')[0]])
def writeSummaryQuant(alleles, genes, n, outfile):
	fp = open(outfile, 'w')
	fp.write('Sample size: %d\n\n'  % n)
	fp.write('Gene level summary\n')
	fp.write('------------------------------------------------------------------\n')
	fp.write('%12s%12s\n' % ('Gene', 'TotalCount'))
	for g in sorted(genes.keys()):
		fp.write('%12s%12d\n' % (g, genes[g]))
	fp.write('\n')
	fp.write('Allele level summary\n')
	fp.write('------------------------------------------------------------------\n')
	fp.write('%20s%12s%12s\n' % ('Allele', 'TotalCount', 'TotalFreq'))
	for allele in sorted(alleles.keys()):
		fp.write("%20s" % allele)
		fp.write("%12d" % alleles[allele])
		fp.write("%12.4f\n" % (1.0 * alleles[allele] / genes[allele.split('*')[0]]))
	fp.close()
###########################################################################
def printAssocChiFisher(assoc, test, permP=None, permN=None, permNA=None):
	print "%20s" % 'Allele',
	for a in ("A_case","B_case","A_ctrl","B_ctrl","F_case","F_ctrl","Freq"):
		print "%8s" % a,
	if test == 'chisq':
		print "%10s" % 'P_Chisq',
		print "%8s" % 'Chisq',
		print "%4s" % 'DF',
	else:
		print "%10s" % 'P_FET',
	for a in ("OR","L95","U95"):
		print "%8s" % a,
	if permP is None:
		print "%10s" % 'P_adj'
	else:
		print "%10s" % 'P_adj',
		print "%10s" % 'P_perm',
		print "%8s" % 'PermN',
		print "%8s" % 'permNA'
	for a in sorted(assoc.keys()):
		print "%20s" % a,
		for i in range(4):
			print "%8d" % assoc[a][i],
		for i in range(4,7):
			print "%8.4f" % assoc[a][i],
		if assoc[a][7] == 'NA':
			print "%10s" % 'NA',
		elif assoc[a][7] > 0.001:
			print "%10.4f" % assoc[a][7],
		else:
			print "%10.2e" % assoc[a][7],
		if test =='chisq':
			print "%8.4f" % assoc[a][8],
			print "%4d" % assoc[a][9],
			for i in range(10, 13):
				print "%8.4f" % assoc[a][i],
			if assoc[a][13] == 'NA':
				print "%10s" % 'NA',
			elif assoc[a][13] > 0.001:
				print "%10.4f" % assoc[a][13],
			else:
				print "%10.2e" % assoc[a][13],
		else:
			for i in range(8, 11):
				print "%8.4f" % assoc[a][i],
			if assoc[a][11] == 'NA':
				print "%10s" % 'NA',
			elif assoc[a][11] > 0.001:
				print "%10.4f" % assoc[a][11],
			else:
				print "%10.2e" % assoc[a][11],
		if permP is None:
			print
		else:
			if permP[a] > 0.001:
				print "%10.4f" % permP[a],
			else:
				print "%10.2e" % permP[a],
			print "%8d" % permN[a],
			print "%8d" % permNA[a]
def writeAssocChiFisher(assoc, test, outfile, permP=None, permN=None, permNA=None):
	fp = open(outfile, 'w')
	fp.write("%20s" % 'Allele')
	for a in ("A_case","B_case","A_ctrl","B_ctrl","F_case","F_ctrl","Freq"):
		fp.write("%8s" % a)
	if test == 'chisq':
		fp.write("%10s" % 'P_Chisq')
		fp.write("%8s" % 'Chisq')
		fp.write("%4s" % 'DF')
	else:
		fp.write("%10s" % 'P_FET')
	for a in ("OR","L95","U95"):
		fp.write("%8s" % a)
	if permP is None:
		fp.write("%10s\n" % 'P_adj')
	else:
		fp.write("%10s" % 'P_adj')
		fp.write("%10s" % 'P_perm')
		fp.write("%8s" % 'PermN')
		fp.write("%8s\n" % 'permNA')
	for a in sorted(assoc.keys()):
		fp.write("%20s" % a)
		for i in range(4):
			fp.write("%8d" % assoc[a][i])
		for i in range(4,7):
			fp.write("%8.4f" % assoc[a][i])
		if assoc[a][7] == 'NA':
			fp.write("%10s" % 'NA')
		elif assoc[a][7] > 0.001:
			fp.write("%10.4f" % assoc[a][7])
		else:
			fp.write("%10.2e" % assoc[a][7])
		if test =='chisq':
			fp.write("%8.4f" % assoc[a][8])
			fp.write("%4d" % assoc[a][9])
			for i in range(10, 13):
				fp.write("%8.4f" % assoc[a][i])
			if assoc[a][13] == 'NA':
				fp.write("%10s" % 'NA')
			elif assoc[a][13] > 0.001:
				fp.write("%10.4f" % assoc[a][13])
			else:
				fp.write("%10.2e" % assoc[a][13])
		else:
			for i in range(8, 11):
				fp.write("%8.4f" % assoc[a][i])
			if assoc[a][11] == 'NA':
				fp.write("%10s" % 'NA')
			elif assoc[a][11] > 0.001:
				fp.write("%10.4f" % assoc[a][11])
			else:
				fp.write("%10.2e" % assoc[a][11])
		if permP is None:
			fp.write('\n')
		else:
			if permP[a] > 0.001:
				fp.write("%10.4f" % permP[a])
			else:
				fp.write("%10.2e" % permP[a])
			fp.write("%8d" % permN[a])
			fp.write("%8d\n" % permNA[a])
	fp.close()
###########################################################################
def printAssocRaw(assoc, alleles, permP=None, permN=None, permNA=None):
	for h in ('Gene', 'Chisq',  'DF', 'P_raw'):
		print "%10s" % h,
	if permP is not None:
		for h in ('P_perm', 'PermN', 'permNA'):
			print "%10s" % h,
	print
	for a in assoc:
		print "%10s" % a,
		print "%10.4f" % assoc[a][0],
		print "%10d" % assoc[a][1],
		if assoc[a][2] == 'NA':
				print "%10s" % 'NA',
		elif assoc[a][2] > 0.001:
			print "%10.4f" % assoc[a][2],
		else:
			print "%10.2e" % assoc[a][2],
		if permP is None:
			print
		else:
			if permP[a] > 0.001:
				print "%10.4f" % permP[a],
			else:
				print "%10.2e" % permP[a],
			print "%10d" % permN[a],
			print "%10d" % permNA[a]
	### alleles
	print '-----------------------------------------------------------------------------'
	print 'Alleles were used:'
	print '-----------------------------------------------------------------------------'
	g = ''
	i = 0
	for a in sorted(alleles):
		if g == '':
			g = a.split('*')[0]
			print "%20s" % a,
			i +=1
		elif g != a.split('*')[0]:
			g = a.split('*')[0]
			print
			print 
			print "%20s" % a,
			i = 1
		else:
			print "%20s" % a,
			i += 1
			if i % 5 == 0:
				print
				i = 0
	print
def writeAssocRaw(assoc, alleles, outfile, permP=None, permN=None, permNA=None):
	fp = open(outfile, 'w')
	for h in ('Gene', 'Chisq',  'DF', 'P_raw'):
		fp.write("%10s" % h)
	if permP is not None:
		for h in ('P_perm', 'PermN', 'permNA'):
			fp.write("%10s" % h)
	fp.write('\n')
	for a in assoc:
		fp.write("%10s" % a)
		fp.write("%10.4f" % assoc[a][0])
		fp.write("%10d" % assoc[a][1])
		if assoc[a][2] == 'NA':
			fp.write("%10s" % 'NA')
		elif assoc[a][2] > 0.001:
			fp.write("%10.4f" % assoc[a][2])
		else:
			fp.write("%10.2e" % assoc[a][2])
		if permP is None:
			fp.write('\n')
		else:
			if permP[a] > 0.001:
				fp.write("%10.4f" % permP[a])
			else:
				fp.write("%10.2e" % permP[a])
			fp.write("%10d" % permN[a])
			fp.write("%10d\n" % permNA[a])
	### alleles
	fp.write('-----------------------------------------------------------------------------\n')
	fp.write('Alleles were used:\n')
	fp.write('-----------------------------------------------------------------------------\n')
	g = ''
	i = 0
	for a in sorted(alleles):
		if g == '':
			g = a.split('*')[0]
			fp.write("%20s" % a)
			i +=1
		elif g != a.split('*')[0]:
			g = a.split('*')[0]
			fp.write('\n\n')
			fp.write("%20s" % a)
			i = 1
		else:
			fp.write("%20s" % a)
			i += 1
			if i % 5 == 0:
				fp.write('\n')
				i = 0
	fp.write('\n')
	fp.close()
###########################################################################
def printAssocScore(assoc, alleles, permP=None, permN=None, permNA=None):
	for h in ('Gene', 'U'):
		print "%10s" % h,
	if permP is not None:
		for h in ('P_perm', 'PermN', 'permNA'):
			print "%10s" % h,
	print
	for a in assoc:
		print "%10s" % a,
		print "%10.4f" % assoc[a],
		if permP is None:
			print
		else:
			if permP[a] > 0.001:
				print "%10.4f" % permP[a],
			else:
				print "%10.2e" % permP[a],
			print "%10s" % permN[a],
			print "%10s" % permNA[a]
	### alleles
	print '-----------------------------------------------------------------------------'
	print 'Alleles were used:'
	print '-----------------------------------------------------------------------------'
	g = ''
	i = 0
	for a in sorted(alleles):
		if g == '':
			g = a.split('*')[0]
			print "%20s" % a,
			i +=1
		elif g != a.split('*')[0]:
			g = a.split('*')[0]
			print
			print 
			print "%20s" % a,
			i = 1
		else:
			print "%20s" % a,
			i += 1
			if i % 5 == 0:
				print
				i = 0
	print
def writeAssocScore(assoc, alleles, outfile, permP=None, permN=None, permNA=None):
	fp = open(outfile, 'w')
	for h in ('Gene', 'U'):
		fp.write("%12s" % h)
	if permP is not None:
		for h in ('P_perm', 'PermN', 'permNA'):
			fp.write("%12s" % h)
	fp.write('\n')
	for a in assoc:
		fp.write("%12s" % a)
		fp.write("%12.4f" % assoc[a])
		if permP is None:
			fp.write('\n')
		else:
			if permP[a] > 0.001:
				fp.write("%12.4f" % permP[a])
			else:
				fp.write("%12.2e" % permP[a])
			fp.write("%12d" % permN[a])
			fp.write("%12d\n" % permNA[a])
	### alleles
	fp.write('-----------------------------------------------------------------------------\n')
	fp.write('Alleles were used:\n')
	fp.write('-----------------------------------------------------------------------------\n')
	g = ''
	i = 0
	for a in sorted(alleles):
		if g == '':
			g = a.split('*')[0]
			fp.write("%20s" % a)
			i +=1
		elif g != a.split('*')[0]:
			g = a.split('*')[0]
			fp.write('\n\n')
			fp.write("%20s" % a)
			i = 1
		else:
			fp.write("%20s" % a)
			i += 1
			if i % 5 == 0:
				fp.write('\n')
				i = 0
	fp.write('\n')
	fp.close()
###########################################################################
def printLogistic(assoc, permP=None, permN=None, permNA=None):
	print "%20s" % 'Allele',
	for a in ("A_case","B_case","A_ctrl","B_ctrl","F_case","F_ctrl","Freq"):
		print "%8s" % a,
	print "%10s" % 'P_Logit',
	for a in ("OR","L95","U95"):
		print "%8s" % a,
	if permP is None:
		print "%10s" % 'P_adj'
	else:
		print "%10s" % 'P_adj',
		print "%10s" % 'P_perm',
		print "%8s" % 'PermN',
		print "%8s" % 'permNA'
	for a in sorted(assoc.keys()):
		print "%20s" % a,
		for i in range(4):
			print "%8d" % assoc[a][i],
		for i in range(4,7):
			print "%8.4f" % assoc[a][i],
		if assoc[a][7] == 'NA':
			print "%10s" % 'NA',
		elif assoc[a][7] > 0.001:
			print "%10.4f" % assoc[a][7],
		else:
			print "%10.2e" % assoc[a][7],
		for i in range(8, 11):
				print "%8.4f" % assoc[a][i],
		if assoc[a][11] == 'NA':
			print "%10s" % 'NA',
		elif assoc[a][11] > 0.001:
			print "%10.4f" % assoc[a][11],
		else:
			print "%10.2e" % assoc[a][11],
		if permP is None:
			print
		else:
			if permP[a] == 'NA':
				print "%10s" % 'NA',
			elif permP[a] > 0.001:
				print "%10.4f" % permP[a],
			else:
				print "%10.2e" % permP[a],
			print "%8d" % permN[a],
			print "%8d" % permNA[a]
def writeLogistic(assoc, outfile, permP=None, permN=None, permNA=None):
	fp = open(outfile, 'w')
	fp.write("%20s" % 'Allele')
	for a in ("A_case","B_case","A_ctrl","B_ctrl","F_case","F_ctrl","Freq"):
		fp.write("%8s" % a)
	fp.write("%10s" % 'P_Logit')
	for a in ("OR","L95","U95"):
		fp.write("%8s" % a)
	if permP is None:
		fp.write("%10s\n" % 'P_adj')
	else:
		fp.write("%10s" % 'P_adj')
		fp.write("%10s" % 'P_perm')
		fp.write("%8s" % 'PermN')
		fp.write("%8s\n" % 'permNA')
	for a in sorted(assoc.keys()):
		fp.write("%20s" % a)
		for i in range(4):
			fp.write("%8d" % assoc[a][i])
		for i in range(4,7):
			fp.write("%8.4f" % assoc[a][i])
		if assoc[a][7] == 'NA':
			fp.write("%10s" % 'NA')
		elif assoc[a][7] > 0.001:
			fp.write("%10.4f" % assoc[a][7])
		else:
			fp.write("%10.2e" % assoc[a][7])
		for i in range(8, 11):
				fp.write("%8.4f" % assoc[a][i])
		if assoc[a][11] == 'NA':
			fp.write("%10s" % 'NA')
		elif assoc[a][11] > 0.001:
			fp.write("%10.4f" % assoc[a][11])
		else:
			fp.write("%10.2e" % assoc[a][11])
		if permP is None:
			fp.write('\n')
		else:
			if permP[a] == 'NA':
				fp.write("%10s" % 'NA')
			elif permP[a] > 0.001:
				fp.write("%10.4f" % permP[a])
			else:
				fp.write("%10.2e" % permP[a])
			fp.write("%8d" % permN[a])
			fp.write("%8d\n" % permNA[a])
###########################################################################
def printLinear(assoc, permP=None, permN=None, permNA=None):
	print "%20s" % 'Allele',
	print "%8s" % "Freq",
	print "%10s" % 'P_Linear',
	for a in ("beta","L95","U95"):
		print "%8s" % a,
	if permP is None:
		print "%10s" % 'P_adj'
	else:
		print "%10s" % 'P_adj',
		print "%10s" % 'P_perm',
		print "%8s" % 'PermN',
		print "%8s" % 'permNA'
	for a in sorted(assoc.keys()):
		print "%20s" % a,
		print "%8.4f" % assoc[a][0],
		if assoc[a][1] == 'NA':
			print "%10s" % 'NA',
		elif assoc[a][1] > 0.001:
			print "%10.4f" % assoc[a][1],
		else:
			print "%10.2e" % assoc[a][1],
		for i in range(2, 5):
				print "%8.4f" % assoc[a][i],
		if assoc[a][5] == 'NA':
			print "%10s" % 'NA',
		elif assoc[a][5] > 0.001:
			print "%10.4f" % assoc[a][5],
		else:
			print "%10.2e" % assoc[a][5],
		if permP is None:
			print
		else:
			if permP[a] > 0.001:
				print "%10.4f" % permP[a],
			else:
				print "%10.2e" % permP[a],
			print "%8d" % permN[a],
			print "%8d" % permNA[a]
def writeLinear(assoc, outfile, permP=None, permN=None, permNA=None):
	fp = open(outfile, 'w')
	fp.write("%20s" % 'Allele')
	fp.write("%8s" % "Freq")
	fp.write("%10s" % 'P_Linear')
	for a in ("beta","L95","U95"):
		fp.write("%8s" % a)
	if permP is None:
		fp.write("%10s\n" % 'P_adj')
	else:
		fp.write("%10s" % 'P_adj')
		fp.write("%10s" % 'P_perm')
		fp.write("%8s" % 'PermN')
		fp.write("%8s\n" % 'permNA')
	for a in sorted(assoc.keys()):
		fp.write("%20s" % a)
		fp.write("%8.4f" % assoc[a][0])
		if assoc[a][1] == 'NA':
			fp.write("%10s" % 'NA')
		elif assoc[a][1] > 0.001:
			fp.write("%10.4f" % assoc[a][1])
		else:
			fp.write("%10.2e" % assoc[a][1])
		for i in range(2, 5):
				fp.write("%8.4f" % assoc[a][i])
		if assoc[a][5] == 'NA':
			fp.write("%10s" % 'NA')
		elif assoc[a][5] > 0.001:
			fp.write("%10.4f" % assoc[a][5])
		else:
			fp.write("%10.2e" % assoc[a][5])
		if permP is None:
			fp.write('\n')
		else:
			if permP[a] > 0.001:
				fp.write("%10.4f" % permP[a])
			else:
				fp.write("%10.2e" % permP[a])
			fp.write("%8d" % permN[a])
			fp.write("%8d\n" % permNA[a])
###########################################################################
def printAssocDelta(assoc, perm=None):
	print "%20s" % 'Allele',
	for a in ("Delta","P_FET","OR","P_adj"):
		print "%10s" % a,
	if perm is None:
		print
	else:
		print "%10s" % "P_perm"
	for a in sorted(assoc.keys()):
		print "%20s" % a,
		print "%10.4f" % assoc[a][0],
		if assoc[a][1] == 'NA':
			print "%10s" % 'NA',
		elif assoc[a][1] > 0.001:
			print "%10.4f" % assoc[a][1],
		else:
			print "%10.2e" % assoc[a][1],
		print "%10.4f" % assoc[a][2],
		if assoc[a][3] > 0.001:
			print "%10.4f" % assoc[a][3],
		else:
			print "%10.2e" % assoc[a][3],
		if perm is None:
			print
		else:
			if assoc[a][4] > 0.001:
				print "%10.4f" % assoc[a][4]
			else:
				print "%10.2e" % assoc[a][4]
def writeAssocDelta(assoc, outfile, perm=None):
	fp = open(outfile, 'w')
	fp.write("%20s" % 'Allele')
	for a in ("Delta","P_FET","OR","P_adj"):
		fp.write("%10s" % a)
	if perm is None:
		fp.write('\n')
	else:
		fp.write("%10s\n" % "P_perm")
	for a in sorted(assoc.keys()):
		fp.write("%20s" % a)
		fp.write("%10.4f" % assoc[a][0])
		if assoc[a][1] == 'NA':
			fp.write("%10s" % 'NA')
		elif assoc[a][1] > 0.001:
			fp.write("%10.4f" % assoc[a][1])
		else:
			fp.write("%10.2e" % assoc[a][1])
		fp.write("%10.4f" % assoc[a][2])
		if assoc[a][3] == 'NA':
			fp.write("%10s" % 'NA')
		elif assoc[a][3] > 0.001:
			fp.write("%10.4f" % assoc[a][3])
		else:
			fp.write("%10.2e" % assoc[a][3])
		if perm is None:
			fp.write('\n')
		else:
			if assoc[a][4] > 0.001:
				fp.write("%10.4f\n" % assoc[a][4])
			else:
				fp.write("%10.2e\n" % assoc[a][4])
###########################################################################