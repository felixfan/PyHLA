#!/usr/bin/env python

import argparse
import time
import sys
import HLAcount
import HLAassoc
import HLAIO
import HLAregression
import HLAAA
import HLAZygosity
import HLAInteraction

strattime = time.time()
###################################################
parser = argparse.ArgumentParser(description='Python for HLA analysis', prog="PyHLA.py")
parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0.0')
parser.add_argument('-V', '--print', help='print output to screen', action='store_true')
parser.add_argument('-i', '--file', help='input file', required=True, type=str)
parser.add_argument('-o', '--out', help='output file', default='output.txt')
parser.add_argument('-d', '--digit', help='digit to test, default 4', default=4, type=int, choices=[2,4,6])
### summary
parser.add_argument('-s', '--summary', help='data summary', action='store_true')
### association analysis
parser.add_argument('-a', '--assoc', help='association analysis', action='store_true')
parser.add_argument('-m', '--model', help='genetic model, default allelic', default='allelic', type=str, choices=['allelic','dom','rec'])
parser.add_argument('-t', '--test', help='statistical test method, default fisher', default='fisher', type=str, choices=['fisher','chisq','logistic','linear','raw','score', 'delta'])
parser.add_argument('-f', '--freq', help='minimal frequency, default 0', default=0, type=float)
parser.add_argument('-j', '--adjust', help='p value correction, default FDR', default='FDR', type=str,choices=['FDR','FDR_BY','Bonferroni','Holm'])
parser.add_argument('-e', '--exclude', help='exclude alleles file', type=str)
parser.add_argument('-c', '--covar', help='covariants file', type=str)
parser.add_argument('-n', '--covarname', help='select a particular subset of covariates', type=str)
parser.add_argument('-p', '--perm', help='number of permutation', type=int)
parser.add_argument('-r', '--seed', help='random seed', type=int)
### amino acid association analysis
parser.add_argument('-A', '--assocAA', help='amino acid association analysis', action='store_true')
parser.add_argument('-u', '--consensus', help='use the sonsensus amino acid senquence', action='store_true')
### amino acid alignment
parser.add_argument('-w', '--align', help='amino acid senquence alignment', action='store_true')
### zygosity
parser.add_argument('-z', '--zygosity', help='zygosity test', action='store_true')
parser.add_argument('-L', '--level', help='level to test', default='residue', type=str, choices=['residue', 'allele'])
### Interaction
parser.add_argument('-I', '--interaction', help='Interaction test', action='store_true')
###################################################
aafile = 'aa.aln.txt'

args = vars(parser.parse_args())

INFILE = args['file']
OUTFILE = args['out']
DIGIT = args['digit']
PRINT =  args['print']

SUMMARY = args['summary']

ASSOC= args['assoc']
TEST = args['test']
MODEL = args['model']
FREQ = args['freq']
ADJUST = args['adjust']
SEED = args['seed'] if 'seed' in args else None
PERM = args['perm'] if 'perm' in args else None
COVFILE = args['covar'] if 'covar' in args else None
COVNAME = args['covarname'] if 'covarname' in args else None
EXCLUDE = args['exclude'] if 'exclude' in args else None

AAA = args['assocAA']
CONSENSUS = args['consensus']

ALN = args['align']

ZYG = args['zygosity']
LEVEL = args['level']

INT = args['interaction']
###################################################
print "@-------------------------------------------------------------@"
print "|       PyHLA       |     v1.0.0      |      28 May 2015      |"
print "|-------------------------------------------------------------|"
print "|  (C) 2015 Felix Yanhui Fan, GNU General Public License, v2  |"
print "|-------------------------------------------------------------|"
print "|    For documentation, citation & bug-report instructions:   |"
print "|          http://felixfan.github.io/PyHLA                    |"
print "@-------------------------------------------------------------@"
print "\n\tOptions in effect:"
print "\t--file", INFILE
if PRINT:
	print "\t--print"
if SUMMARY:
	print "\t--digit", DIGIT
	print "\t--summary"
elif ASSOC:
	print "\t--digit", DIGIT
	print "\t--assoc"
	print "\t--test", TEST
	if TEST == 'chisq' or TEST == 'fisher':
		print "\t--model", MODEL
	elif TEST == 'logistic' or 'linear':
		if COVFILE:
			print "\t--covar", COVFILE
			if COVNAME:
				print "\t--covarname", COVNAME
	print "\t--freq", FREQ
	if EXCLUDE:
		print "\t--exclude", EXCLUDE
	if TEST != 'raw' and TEST != 'score':
		print "\t--adjust", ADJUST
	if PERM:
		print "\t--perm", PERM
		if SEED:
			print "\t--seed", SEED
elif AAA:
	print "\t--assocAA"
	print "\t--test", TEST
	if CONSENSUS:
		print "\t--consensus"
elif ALN:
	print "\t--align"
	if CONSENSUS:
		print "\t--consensus"
elif ZYG or INT:
	if ZYG:
		print "\t--zygosity"
	elif INT:
		print "\t--interaction"
	print "\t--test", TEST
	print "\t--level", LEVEL
	if LEVEL == 'allele':
		print "\t--digit", DIGIT
		print "\t--freq", FREQ
	else:
		if CONSENSUS:
			print "\t--consensus"
print "\t--out", OUTFILE
print
###################################################
if SUMMARY:
	if HLAcount.quantTrait(INFILE):
		alleles, genes, n = HLAcount.alleleCount(INFILE, DIGIT)
		if PRINT:
			HLAIO.printSummaryQuant(alleles, genes, n)
		HLAIO.writeSummaryQuant(alleles, genes, n, OUTFILE)
	else:
		caseAlleles, ctrlAlleles, np, nc, nn = HLAcount.allelicCount(INFILE, DIGIT)
		freq, alleles = HLAcount.hlaFreq(caseAlleles, ctrlAlleles, np, nc, nn[0])
		popCase, popCtrl, popP, popC, popN = HLAcount.domCount(INFILE, DIGIT)
		if PRINT:
			HLAIO.printSummary(alleles, freq, caseAlleles, ctrlAlleles, np, nc, nn, popCase, popCtrl, popP, popC)
		HLAIO.writeSummary(alleles, freq, caseAlleles, ctrlAlleles, np, nc, nn, OUTFILE, popCase, popCtrl, popP, popC)
elif ASSOC:
	if HLAcount.quantTrait(INFILE):
		if TEST != 'linear':
			sys.exit("quantitative trait was detected, only linear regression can be used for association analysis!")
	else:
		if TEST =='linear':
			sys.exit("case-control trait was detected, linear regression can not be used for association analysis!")
	if TEST == 'chisq' or TEST == 'fisher':
		if PERM is None:
			assoc = HLAassoc.assocADRChiFisher(INFILE, DIGIT, FREQ, TEST,MODEL,ADJUST,EXCLUDE, PERM, SEED)
			if PRINT:
				HLAIO.printAssocChiFisher(assoc, TEST)
			HLAIO.writeAssocChiFisher(assoc, TEST, OUTFILE)
		else:
			assoc, permP, permN, permNA, permNL= HLAassoc.assocADRChiFisher(INFILE, DIGIT, FREQ, TEST,MODEL,ADJUST,EXCLUDE, PERM, SEED)
			if PRINT:
				HLAIO.printAssocChiFisher(assoc, TEST, permP, permN, permNA)
			HLAIO.writeAssocChiFisher(assoc, TEST, OUTFILE, permP, permN, permNA)
	elif TEST == 'raw':
		if PERM is None:
			assoc, alleles = HLAassoc.assocRaw(INFILE, DIGIT, FREQ,EXCLUDE, PERM, SEED)
			if PRINT:
				HLAIO.printAssocRaw(assoc,alleles)
			HLAIO.writeAssocRaw(assoc, alleles, OUTFILE)
		else:
			assoc, alleles, permP, permN, permNA, permNL = HLAassoc.assocRaw(INFILE, DIGIT, FREQ,EXCLUDE, PERM, SEED)
			if PRINT:
				HLAIO.printAssocRaw(assoc, alleles, permP, permN, permNA)
			HLAIO.writeAssocRaw(assoc, alleles, OUTFILE, permP, permN, permNA)
	elif TEST == 'score':
		if PERM is None:
			assoc, alleles = HLAassoc.assocScoreU(INFILE, DIGIT, FREQ,EXCLUDE, PERM, SEED)
			if PRINT:
				HLAIO.printAssocScore(assoc, alleles)
			HLAIO.writeAssocScore(assoc, alleles, OUTFILE)
		else:
			assoc, alleles, permP, permN, permNA, permNL = HLAassoc.assocScoreU(INFILE, DIGIT, FREQ,EXCLUDE, PERM, SEED)
			if PRINT:
				HLAIO.printAssocScore(assoc, alleles, permP, permN, permNA)
			HLAIO.writeAssocScore(assoc, alleles, OUTFILE, permP, permN, permNA)
	elif TEST == 'logistic':
		if PERM is None:
			assoc = HLAregression.regressionLogistic(INFILE, DIGIT, FREQ, ADJUST, EXCLUDE, COVFILE, COVNAME, PERM, SEED, TEST)
			if PRINT:
				HLAIO.printLogistic(assoc)
			HLAIO.writeLogistic(assoc, OUTFILE)
		else:
			assoc, permP, permN, permNA = HLAregression.regressionLogistic(INFILE, DIGIT, FREQ, ADJUST, EXCLUDE, COVFILE, COVNAME, PERM, SEED, TEST)
			if PRINT:
				HLAIO.printLogistic(assoc, permP, permN, permNA)
			HLAIO.writeLogistic(assoc, OUTFILE, permP, permN, permNA)
	elif TEST == 'linear':
		if PERM is None:
			assoc = HLAregression.regressionLinear(INFILE, DIGIT, FREQ, ADJUST, EXCLUDE, COVFILE, COVNAME, PERM, SEED, TEST)
			if PRINT:
				HLAIO.printLinear(assoc)
			HLAIO.writeLinear(assoc, OUTFILE)
		else:
			assoc, permP, permN, permNA = HLAregression.regressionLinear(INFILE, DIGIT, FREQ, ADJUST, EXCLUDE, COVFILE, COVNAME, PERM, SEED, TEST)
			if PRINT:
				HLAIO.printLinear(assoc, permP, permN, permNA)
			HLAIO.writeLinear(assoc, OUTFILE, permP, permN, permNA)
	elif TEST == 'delta':
		assoc = HLAassoc.assocDelta(INFILE, DIGIT, FREQ, ADJUST, EXCLUDE, PERM, SEED)
		if PRINT:
			HLAIO.printAssocDelta(assoc, PERM)
		HLAIO.writeAssocDelta(assoc, OUTFILE, PERM)
elif AAA or ALN:
	if HLAcount.quantTrait(INFILE):
		sys.exit("quantitative trait was detected, only case-control data can be used for amino acid analysis!")
	case, ctrl, caseGeno, ctrlGeno, ncase, nctrl = HLAAA.readGeno(INFILE)
	seq = HLAAA.readAAseq(aafile)
	alleles = HLAAA.keyDicts(case, ctrl)
	myseq = HLAAA.getSeq(alleles, seq, CONSENSUS)
	if AAA:
		assoc = HLAAA.aaAssoc(case, ctrl, caseGeno, ctrlGeno, ncase, nctrl, myseq, TEST)
		if PRINT:
			HLAAA.printAAA(assoc)
		HLAAA.writeAAA(assoc, OUTFILE)
	elif ALN:
		myaln = HLAAA.aaAlign(case, ctrl, myseq)
		if PRINT:
			HLAAA.printAA(myaln)
		HLAAA.writeAA(myaln, OUTFILE)
elif ZYG or INT:
	if TEST != 'fisher' and TEST != 'chisq':
		sys.exit("only 'fisher' and 'chisq' test can be used!")
	if LEVEL == 'residue':
		case, ctrl, caseGeno, ctrlGeno, ncase, nctrl = HLAAA.readGeno(INFILE)
		seq = HLAAA.readAAseq(aafile)
		alleles = HLAAA.keyDicts(case, ctrl)
		myseq = HLAAA.getSeq(alleles, seq, CONSENSUS)
		assoc = HLAAA.aaAssoc(case, ctrl, caseGeno, ctrlGeno, ncase, nctrl, myseq, TEST)
		sig = {}
		for k in assoc:
			if assoc[k][4] != 'NA' and assoc[k][4] < 0.05:
				sig[k] = 1
		keys = sorted(sig.keys())
		if ZYG: # zygosity test
			ans = HLAZygosity.zygosityAA(keys, caseGeno, ctrlGeno, myseq, TEST)
			if PRINT:
				HLAZygosity.printZygosity(ans, LEVEL)
			HLAZygosity.writeZygosity(ans, LEVEL,OUTFILE)
		else:    # interaction test
			ans = HLAInteraction.interactAA(keys, caseGeno, ctrlGeno, myseq, TEST)
			if PRINT:
				HLAInteraction.printInteract(ans, LEVEL)
			HLAInteraction.writeInteract(ans, LEVEL, OUTFILE)
	else:
		assoc = HLAassoc.assocADRChiFisher(INFILE, DIGIT, FREQ, TEST)
		caseGeno, ctrlGeno = HLAInteraction.readAlleleZygInteract(INFILE, DIGIT)
		sig = {}
		for k in assoc:
			if assoc[k][7] != 'NA' and assoc[k][7] < 0.05:
				sig[k] = 1
		keys = sorted(sig.keys())
		if  ZYG:
			ans = HLAZygosity.zygosityAllele(keys, caseGeno, ctrlGeno, TEST)
			if PRINT:
				HLAZygosity.printZygosity(ans, LEVEL)
			HLAZygosity.writeZygosity(ans, LEVEL,OUTFILE)
		else:
			ans = HLAInteraction.interactAllele(keys, caseGeno, ctrlGeno, TEST)
			if PRINT:
				HLAInteraction.printInteract(ans, LEVEL)
			HLAInteraction.writeInteract(ans, LEVEL, OUTFILE)
else:
	c = ('--summary', '--assoc', '--assocAA', '--align', '--zygosity', '--interaction')
	sys.exit("one of the following option must be used: \n\n%s\n%s\n%s\n%s\n%s\n%s\n" % c)
###################################################
usedtime = time.time() - strattime
print "Time used:",
if usedtime >=60:
	ts = int(usedtime) % 60
	usedtime = int(usedtime) / 60
	tm = int(usedtime) % 60
	usedtime = int(usedtime) / 60
	th = int(usedtime) % 60
	if th > 0:
		print "%d hours"  % th,
		print "%d minutes"  % tm,
	elif tm > 0:
		print "%d minutes"  % tm,
else:
	ts = usedtime
print '%.2f seconds' % ts
print "Finished at ",
print time.strftime("%H:%M:%S %d %b %Y")
