#!/usr/bin/env python

def printSummary(alleles, freq, caseAlleles, ctrlAlleles, np, nc, nn):
	print 'Sample size: %d'  % nn
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
def writeSummary(alleles, freq, caseAlleles, ctrlAlleles, np, nc, nn, OUTFILE):
	fp = open(OUTFILE, 'w')
	fp.write('Sample size: %d\n\n'  % nn)
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
		print "%6s" % a,
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
		if assoc[a][7] > 0.001:
			print "%10.4f" % assoc[a][7],
		else:
			print "%10.2e" % assoc[a][7],
		if test =='chisq':
			print "%8.4f" % assoc[a][8],
			print "%4d" % assoc[a][9],
			for i in range(10, 13):
				print "%6.4f" % assoc[a][i],
			if assoc[a][13] > 0.001:
				print "%10.4f" % assoc[a][13],
			else:
				print "%10.2e" % assoc[a][13],
		else:
			for i in range(8, 11):
				print "%6.4f" % assoc[a][i],
			if assoc[a][11] > 0.001:
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
		fp.write("%6s" % a)
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
		if assoc[a][7] > 0.001:
			fp.write("%10.4f" % assoc[a][7])
		else:
			fp.write("%10.2e" % assoc[a][7])
		if test =='chisq':
			fp.write("%8.4f" % assoc[a][8])
			fp.write("%4d" % assoc[a][9])
			for i in range(10, 13):
				fp.write("%6.4f" % assoc[a][i])
			if assoc[a][13] > 0.001:
				fp.write("%10.4f" % assoc[a][13])
			else:
				fp.write("%10.2e" % assoc[a][13])
		else:
			for i in range(8, 11):
				fp.write("%6.4f" % assoc[a][i])
			if assoc[a][11] > 0.001:
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
		if assoc[a][2] > 0.001:
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
		if assoc[a][2] > 0.001:
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

###########################################################################
