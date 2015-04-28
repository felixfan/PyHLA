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
###########################################################################
