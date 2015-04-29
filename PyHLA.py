#!/usr/bin/env python

import argparse
import HLAcount
import HLAassoc
import HLAIO
import HLAregression

###################################################
parser = argparse.ArgumentParser(description='Python for HLA analysis', prog="PyHLA.py")
parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.4')
parser.add_argument('-V', '--print', help='print output to screen', action='store_true')
parser.add_argument('-i', '--infile', help='input file', required=True, type=str)
parser.add_argument('-o', '--out', help='output file', default='output.txt')
parser.add_argument('-d', '--digits', help='digits to test, default 4', default=4, type=int, choices=[2,4,6])
### summary
parser.add_argument('-s', '--summary', help='data summary', action='store_true')
### quality control
parser.add_argument('-q', '--qc', help='quality control', action='store_true')
### association analysis
parser.add_argument('-a', '--assoc', help='association analysis', action='store_true')
parser.add_argument('-m', '--model', help='genetic model, default allelic', default='allelic', type=str, choices=['allelic','dom','rec'])
parser.add_argument('-t', '--test', help='statistical test method, default chisq', default='chisq', type=str, choices=['chisq','fisher','logistic','linear','raw','score'])
parser.add_argument('-f', '--freq', help='minimal frequency, default 0.05', default=0.05, type=float)
parser.add_argument('-j', '--adjust', help='p value correction, default FDR', default='FDR', type=str,choices=['FDR','FDR_BY','Bonferroni','Holm'])
parser.add_argument('-e', '--exclude', help='exclude alleles file', type=str)
parser.add_argument('-c', '--covar', help='covariants file', type=str)
parser.add_argument('-n', '--covarname', help='select a particular subset of covariates', type=str)
parser.add_argument('-p', '--perm', help='number of permutation', type=int)
parser.add_argument('-r', '--seed', help='random seed', type=int)
### annotation
parser.add_argument('-A', '--annotation', help='annotation', action='store_true')
###################################################
args = vars(parser.parse_args())

INFILE = args['infile']
OUTFILE = args['out']
DIGIT = args['digits']
PRINT =  args['print']

SUMMARY = args['summary']

QC = args['qc']

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

ANNOT = args['annotation']
###################################################
if SUMMARY:
	if HLAcount.quantTrait(INFILE):
		alleles, genes, n = HLAcount.alleleCount(INFILE, DIGIT)
		if PRINT:
			HLAIO.printSummaryQuant(alleles, genes, n)
		HLAIO.writeSummaryQuant(alleles, genes, n, OUTFILE)
	else:
		caseAlleles, ctrlAlleles, np, nc, nn = HLAcount.allelicCount(INFILE, DIGIT)
		freq, alleles = HLAcount.hlaFreq(caseAlleles, ctrlAlleles, np, nc, nn)
		if PRINT:
			HLAIO.printSummary(alleles, freq, caseAlleles, ctrlAlleles, np, nc, nn)
		HLAIO.writeSummary(alleles, freq, caseAlleles, ctrlAlleles, np, nc, nn, OUTFILE)
elif QC:
	pass
elif ASSOC:
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
		pass
elif ANNOT:
	pass
else:
	pass
###################################################

