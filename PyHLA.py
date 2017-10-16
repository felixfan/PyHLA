#coding=utf-8
#!/usr/bin/env python

from __future__ import division
import statsmodels.formula.api as smf
from collections import Counter
import numpy as np
import pandas as pd
import scipy.stats
import argparse
import time
import sys
import random
import math
import os
import string
###################### HLAcount ############################################################
def allelicCount(infile, digits, perm=False):
    '''
    count all alleles
    return the count of each allele in case and control
    '''
    case = {}      # counts for each allele
    control = {}
    np = {}         # total non-NA alleles for each gene
    nc = {}
    f = open(infile) # read all phenotype
    pht = []
    for line in f:
        line = line.rstrip()
        pht.append(line.split()[1])
    f.close()
    if perm:
        random.shuffle(pht)
    l = 0     # line index
    nn = []
    nn.append(len(pht))
    nnp = 0
    nnc = 0
    for p in pht:
        if p == '1':
            nnc += 1
        elif p == '2':
            nnp += 1
    nn.append(nnp)
    nn.append(nnc)
    f = open(infile)
    for line in f:
        line = line.rstrip()
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
                if pht[l] == "2":
                    if a4d1.split('*')[0] in np:
                        np[a4d1.split('*')[0]] += 1
                    else:
                        np[a4d1.split('*')[0]] = 1
                    if a4d1 in case:
                        case[a4d1] += 1
                    else:
                        case[a4d1] = 1
                    ##
                    if a4d2.split("*")[0] in np:
                        np[a4d2.split("*")[0]] += 1
                    else:
                        np[a4d2.split("*")[0]] = 1
                    if a4d2 in case:
                        case[a4d2] += 1
                    else:
                        case[a4d2] = 1
                elif pht[l] == "1":
                    if a4d1.split('*')[0] in nc:
                        nc[a4d1.split('*')[0]] += 1
                    else:
                        nc[a4d1.split('*')[0]] = 1

                    if a4d1 in control:
                        control[a4d1] += 1
                    else:
                        control[a4d1] = 1
                    ##
                    if a4d2.split("*")[0] in nc:
                        nc[a4d2.split("*")[0]] += 1
                    else:
                        nc[a4d2.split("*")[0]] = 1
                    if a4d2 in control:
                        control[a4d2] += 1
                    else:
                        control[a4d2] = 1
        l += 1
    f.close()
    return case, control, np, nc, nn
def domCount(infile, digits, perm=False):
    '''
    count all alleles
    return the count of each allele in case and control
    '''
    case = {}     # counts for each allele
    control = {}
    nc = {}
    np = {}
    f = open(infile) # read all phenotype
    pht = []
    for line in f:
        line = line.rstrip()
        pht.append(line.split()[1])
    f.close()
    if perm:
        random.shuffle(pht)
    l = 0     # line index
    nn = len(pht)
    f = open(infile)
    for line in f:
        line = line.rstrip()
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
                if pht[l] == "2":
                    if a4d1.split('*')[0] in np:
                        np[a4d1.split('*')[0]] += 1
                    else:
                        np[a4d1.split('*')[0]] = 1
                    if a4d1 in case:
                        case[a4d1] += 1
                    else:
                        case[a4d1] = 1
                    if a4d1 != a4d2:
                        if a4d2 in case:
                            case[a4d2] += 1
                        else:
                            case[a4d2] = 1
                elif pht[l] == "1":
                    if a4d1.split('*')[0] in nc:
                        nc[a4d1.split('*')[0]] += 1
                    else:
                        nc[a4d1.split('*')[0]] = 1
                    if a4d1 in control:
                        control[a4d1] += 1
                    else:
                        control[a4d1] = 1
                    if a4d1 != a4d2:
                        if a4d2 in control:
                            control[a4d2] += 1
                        else:
                            control[a4d2] = 1
        l += 1
    f.close()
    return case, control, np, nc, nn
def recCount(infile, digits, perm=False):
    '''
    count all alleles
    return the count of each allele in case and control
    '''
    case = {}     # counts for each allele
    control = {}
    nc = {}
    np = {}
    f = open(infile) # read all phenotype
    pht = []
    for line in f:
        line = line.rstrip()
        pht.append(line.split()[1])
    f.close()
    if perm:
        random.shuffle(pht)
    l = 0     # line index
    nn = len(pht)
    f = open(infile)
    for line in f:
        line = line.rstrip()
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
                if pht[l] == "2":
                    if a4d1.split('*')[0] in np:
                        np[a4d1.split('*')[0]] += 1
                    else:
                        np[a4d1.split('*')[0]] = 1
                    if a4d1 == a4d2:
                        if a4d2 in case:
                            case[a4d2] += 1
                        else:
                            case[a4d2] = 1
                elif pht[l] == "1":
                    if a4d1.split('*')[0] in nc:
                        nc[a4d1.split('*')[0]] += 1
                    else:
                        nc[a4d1.split('*')[0]] = 1
                    if a4d1 == a4d2:
                        if a4d2 in control:
                            control[a4d2] += 1
                        else:
                            control[a4d2] = 1
        l += 1
    f.close()
    return case, control, np, nc, nn
def hlaFreq(caseAlleles, ctrlAlleles, np, nc, nn):
    '''
    Frequency of all alleles
    return the frequency of each allele in case, control, and both
    '''
    temp = []
    temp.extend(caseAlleles.keys())
    temp.extend(ctrlAlleles.keys())
    alleles = sorted(set(temp))
    if 'NA' in alleles:
        alleles.remove('NA')
    freq = {}
    for allele in alleles:
        ff = []
        if allele in caseAlleles:
            ff.append(1.0 * caseAlleles[allele] / np[allele.split('*')[0]])
        else:
            ff.append(0)
        if allele in ctrlAlleles:
            ff.append(1.0 * ctrlAlleles[allele] / nc[allele.split('*')[0]])
        else:
            ff.append(0)
        if allele in caseAlleles and allele in ctrlAlleles:
            ff.append(1.0 * (ctrlAlleles[allele] + caseAlleles[allele]) / (np[allele.split('*')[0]] + nc[allele.split('*')[0]]))
        elif allele in caseAlleles:
            ff.append(1.0 * caseAlleles[allele] / (np[allele.split('*')[0]] + nc[allele.split('*')[0]]))
        else:
            ff.append(1.0 * ctrlAlleles[allele] / (np[allele.split('*')[0]] + nc[allele.split('*')[0]]))
        freq[allele] = ff
    return freq, alleles
def alleleCount(infile, digits):
    '''
    count all alleles for quantitative traits
    return the count of each allele in case and control
    '''
    allelesN = {}
    genesN = {}
    N = 0
    f = open(infile)
    for line in f:
        line = line.rstrip()
        alleles = line.split()
        N += 1
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
                if a4d1.split('*')[0] in genesN:
                    genesN[a4d1.split('*')[0]] += 1
                else:
                    genesN[a4d1.split('*')[0]] = 1
                if a4d2.split('*')[0] in genesN:
                    genesN[a4d2.split('*')[0]] += 1
                else:
                    genesN[a4d2.split('*')[0]] = 1
                if a4d1 in allelesN:
                    allelesN[a4d1] += 1
                else:
                    allelesN[a4d1] = 1
                if a4d2 in allelesN:
                    allelesN[a4d2] += 1
                else:
                    allelesN[a4d2] = 1
    f.close()
    return allelesN, genesN, N
def quantTrait(infile):
    f = open(infile)
    for line in f:
        line = line.rstrip()
        alleles = line.split()
        if alleles[1] != '1' and alleles[1] != '2' and alleles[1] != 'NA':
            return True
    f.close()
    return False
######################### HLAIO ############################################################
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
                if assoc[a][i] == 'NA':
                    print "%8s" % 'NA',
                else:
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
                if assoc[a][i] == 'NA':
                    fp.write("%8s" % 'NA')
                else:
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
################  HLAassoc   ###############################################################
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
        caseAlleles, ctrlAlleles, np, nc, nn = allelicCount(infile, digit)
    elif model == 'dom':
        caseAlleles, ctrlAlleles, np, nc, nn = domCount(infile, digit)
    elif model == 'rec':
        caseAlleles, ctrlAlleles, np, nc, nn = recCount(infile, digit)
    alleleFreq, alleles = hlaFreq(caseAlleles, ctrlAlleles, np, nc, nn)

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
            if a.startswith(g):
                if assoc[a][7] != 'NA':  # p value at 7 col start from 0
                    ps.append(assoc[a][7])
                    ns.append(a)
                else:
                    assoc[a].append('NA')
        cp = adjustP(ps,adjust)
        for i in range(len(ns)):
            assoc[ns[i]].append(cp[i])
    if perm is None:
        return assoc
    ### perm
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
                case9, ctrl9, np9, nc9, nn9 = allelicCount(infile,digit, True)
            elif model == 'dom':
                case9, ctrl9, np9, nc9, nn9 = domCount(infile,digit, True)
            elif model == 'rec':
                case9, ctrl9, np9, nc9, nn9 = recCount(infile,digit, True)
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
    caseAlleles, ctrlAlleles, np, nc, nn = allelicCount(infile, digit)
    alleleFreq, alleles = hlaFreq(caseAlleles, ctrlAlleles, np, nc, nn)
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
            case9, ctrl9, np9, nc9, nn9 = allelicCount(infile,digit, True)
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
    caseAlleles, ctrlAlleles, np, nc, nn = allelicCount(infile, digit)
    alleleFreq, alleles = hlaFreq(caseAlleles, ctrlAlleles, np, nc, nn)
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
            case9, ctrl9, np9, nc9, nn9 = allelicCount(infile,digit, True)
            alleleFreq9, alleles9 = hlaFreq(caseAlleles, ctrlAlleles, np, nc, nn)
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
def assocDelta(infile, digit, freq=0.05, adjust='FDR', exclude=None, perm=None, seed=None):
    popCase, popCtrl, popP, popC, popN = domCount(infile, digit)
    caseAlleles99, ctrlAlleles99, np99, nc99, nn99 = allelicCount(infile, digit)
    alleleFreq, alleles = hlaFreq(caseAlleles99, ctrlAlleles99, np99, nc99, nn99)

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
            if a.startswith(g):
                if  assoc[a][1] != 'NA':  # p value at 7 col start from 0
                    ps.append(assoc[a][1])
                    ns.append(a)
                else:
                    assoc[a].append('NA')
        cp = adjustP(ps,adjust)
        for i in range(len(ns)):
            assoc[ns[i]].append(cp[i])
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
            case9, ctrl9, np9, nc9, nn9 = domCount(infile,digit, True)
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
################  HLAregression ############################################################
def getAlleles(infile, digits):
    '''
    get all alleles for each gene
    return a dictionary, keys are the gene names, values are the alleles
    return second dictionary, keys are the start column, values are the gene names
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
                    if len(names)>=2:
                        temp = names[0] + ":" + names[1]
                        if temp not in geneAlleles[gene[0]]:
                            geneAlleles[gene[0]].append(temp)
                    else:
                        sys.exit('please use a lower digits!')
                elif digits == 2:
                    if len(names)>=1:
                        temp = names[0]
                        if temp not in geneAlleles[gene[0]]:
                            geneAlleles[gene[0]].append(temp)
                    else:
                        sys.exit('please use a lower digits!')
                elif digits == 6:
                    if len(names)>=3:
                        temp = names[0] + ":" + names[1] + ":" + names[2]
                        if temp not in geneAlleles[gene[0]]:
                            geneAlleles[gene[0]].append(temp)
                    else:
                        sys.exit('please use a lower digits!')
    f.close()
    return geneAlleles, geneCol
def allelicRecode(infile, digits, test, model):
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
                    elif digits == 2:
                        allele1 = alleles[i].split(':')[0]
                        allele2 = alleles[j].split(':')[0]
                    elif digits == 6:
                        allele1 = alleles[i].split(':')[0] + ':' + alleles[i].split(':')[1] + ':' + alleles[i].split(':')[2]
                        allele2 = alleles[j].split(':')[0] + ':' + alleles[j].split(':')[1] + ':' + alleles[j].split(':')[2]
                    gAlleles = sorted(geneAlleles[gene1])
                    for ga in gAlleles:
                        if ga not in header:
                            header.append(ga)
                        if model == 'additive':
                            if allele1 == ga and allele2 == ga:
                                ans[alleles[0]].append(2)
                            elif allele1 == ga or allele2 == ga:
                                ans[alleles[0]].append(1)
                            else:
                                ans[alleles[0]].append(0)
                        elif model == 'dom':
                            if allele1 == ga and allele2 == ga:
                                ans[alleles[0]].append(1)
                            elif allele1 == ga or allele2 == ga:
                                ans[alleles[0]].append(1)
                            else:
                                ans[alleles[0]].append(0)
                        elif model =='rec':
                            if allele1 == ga and allele2 == ga:
                                ans[alleles[0]].append(1)
                            elif allele1 == ga or allele2 == ga:
                                ans[alleles[0]].append(0)
                            else:
                                ans[alleles[0]].append(0)
                else:
                    sys.exit("input format is wrong!")
            else:
                for gg in geneAlleles[geneCol[i]]:
                    ans[alleles[0]].append('NA')
    return ans,header
def writeRecode(infile, digits, test, model):
    '''
    write coding to a temp file
    return the temp file name
    '''
    tmp = time.strftime("%H%M%S%d%b%Y")
    tmp = tmp + '.txt'
    f = open(tmp,'w')
    ans, header = allelicRecode(infile, digits, test, model)
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
def regressionLogistic(infile, digits, freq, model = 'additive', adjust = 'FDR', exclude=None, covfile=None, covname=None, perm=None, seed=None, test='logistic'):
    '''
    logistitic regression
    output: dictionary, key: allele, value: statistic includes count, freq, p and OR
    '''
    ### geno
    tfile = writeRecode(infile, digits, test, model)
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
            if a.startswith(g):
                if assoc[a][7] != 'NA':  # p value at 7 col start from 0
                    ps.append(assoc[a][7])
                    ns.append(a)
                else:
                    assoc[a].append('NA')
        cp = adjustP(ps,adjust)
        for i in range(len(ns)):
            assoc[ns[i]].append(cp[i])
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
                        lr = smf.logit(formula = myformula, data = mydata).fit(maxiter=100, disp=False)
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
def regressionLinear(infile, digits, freq, model = 'additive', adjust = 'FDR', exclude=None, covfile=None, covname=None, perm=None, seed=None, test='linear'):
    '''
    linear regression with covariants
    output: dictionary, key: allele, value: statistic includes freq, p and beta
    '''
    ### geno
    tfile = writeRecode(infile, digits, test, model)
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
            if a.startswith(g):
                if assoc[a][1] != 'NA':  # p value at 1 col start from 0
                    ps.append(assoc[a][1])
                    ns.append(a)
                else:
                    assoc[a].append('NA')
        cp = adjustP(ps,adjust)
        for i in range(len(ns)):
            assoc[ns[i]].append(cp[i])
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
###################  HLAAA  ################################################################
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
    if 'NA' in alleles:
        alleles.remove('NA')
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
    if 'NA' in alleles:
        alleles.remove('NA')
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
###################   HLAInteraction   #####################################################
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
    OR1 = (n1+0.5) * (n4+0.5) / (n2+0.5) / (n3+0.5)
    p1 = simpleTest(data, test)
    # test2
    n1 = x[0] + x[2]
    n2 = x[1] + x[3]
    n3 = y[0] + y[2]
    n4 = y[1] + y[3]
    data = [[n1, n2], [n3, n4]]
    OR2 = (n1+0.5) * (n4+0.5) / (n2+0.5) / (n3+0.5)
    p2 = simpleTest(data, test)
    # test 3
    n1 = x[0]
    n2 = x[2]
    n3 = y[0]
    n4 = y[2]
    data = [[n1, n2], [n3, n4]]
    OR3 = (n1+0.5) * (n4+0.5) / (n2+0.5) / (n3+0.5)
    p3 = simpleTest(data, test)
    # test 4
    n1 = x[1]
    n2 = x[3]
    n3 = y[1]
    n4 = y[3]
    data = [[n1, n2], [n3, n4]]
    OR4 = (n1+0.5) * (n4+0.5) / (n2+0.5) / (n3+0.5)
    p4 = simpleTest(data, test)
    # test5
    n1 = x[0]
    n2 = x[1]
    n3 = y[0]
    n4 = y[1]
    data = [[n1, n2], [n3, n4]]
    OR5 = (n1+0.5) * (n4+0.5) / (n2+0.5) / (n3+0.5)
    p5 = simpleTest(data, test)
    # test6
    n1 = x[2]
    n2 = x[3]
    n3 = y[2]
    n4 = y[3]
    data = [[n1, n2], [n3, n4]]
    OR6 = (n1+0.5) * (n4+0.5) / (n2+0.5) / (n3+0.5)
    p6 = simpleTest(data, test)
    # test7
    n1 = x[1]
    n2 = x[2]
    n3 = y[1]
    n4 = y[2]
    data = [[n1, n2], [n3, n4]]
    OR7 = (n1+0.5) * (n4+0.5) / (n2+0.5) / (n3+0.5)
    p7 = simpleTest(data, test)
    # test8
    n1 = x[0]
    n2 = x[3]
    n3 = y[0]
    n4 = y[3]
    data = [[n1, n2], [n3, n4]]
    OR8 = (n1+0.5) * (n4+0.5) / (n2+0.5) / (n3+0.5)
    p8 = simpleTest(data, test)
    # test9
    n1 = x[0]
    n2 = x[1]
    n3 = x[2]
    n4 = x[3]
    data = [[n1, n2], [n3, n4]]
    OR9 = (n1+0.5) * (n4+0.5) / (n2+0.5) / (n3+0.5)
    p9 = simpleTest(data, test)
    # test10
    n1 = y[0]
    n2 = y[1]
    n3 = y[2]
    n4 = y[3]
    data = [[n1, n2], [n3, n4]]
    OR10 = (n1+0.5) * (n4+0.5) / (n2+0.5) / (n3+0.5)
    p10 = simpleTest(data, test)
    return [p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,OR1,OR2,OR3,OR4,OR5,OR6,OR7,OR8,OR9,OR10]
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
def printInteract(assoc, level):
    header = ('ID1', 'ID2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9','P10')
    for h in header:
        if h == 'ID1' or h == 'ID2':
            print "%-12s" % h,
        else:
            print "%12s" % h,
    header2 = ('OR3', 'OR4', 'OR5', 'OR6', 'OR7', 'OR8', 'OR9','OR10')
    for h in header2:
        print "%10s" % h,
    print

    for k in sorted(assoc.keys()):
        if level == 'residue':
            print "%-12s" % (k[0]+'_'+str(k[1])+'_'+k[2]),
            print "%-12s" % (k[3]+'_'+str(k[4])+'_'+k[5]),
        else:
            print "%-12s"  % k[0],
            print "%-12s" % k[1],
        for i in range(2, 10):
            if assoc[k][i] == 'NA':
                print "%12s" % "NA",
            elif assoc[k][i] > 0.001:
                print "%12.4f" % assoc[k][i],
            else:
                print "%12.2e" % assoc[k][i],
        for i in range(12, len(assoc[k])):
            if assoc[k][i] > 100:
                print "%10.0f" %  assoc[k][i],
            else:
                print "%10.2f" %  assoc[k][i],
        print
def writeInteract(assoc, level, outfile):
    fw = open(outfile, 'w')
    header = ('ID1', 'ID2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9','P10')
    for h in header:
        if h == 'ID1' or h == 'ID2':
            fw.write("%-12s" % h)
        else:
            fw.write("%12s" % h)
    header2 = ('OR3', 'OR4', 'OR5', 'OR6', 'OR7', 'OR8', 'OR9','OR10')
    for h in header2:
        fw.write("%10s" % h)
    fw.write('\n')

    for k in sorted(assoc.keys()):
        if level == 'residue':
            fw.write("%-12s" % (k[0]+'_'+str(k[1])+'_'+k[2]))
            fw.write("%-12s" % (k[3]+'_'+str(k[4])+'_'+k[5]))
        else:
            fw.write("%-12s"  % k[0])
            fw.write("%-12s" % k[1])
        for i in range(2, 10):
            if assoc[k][i] == 'NA':
                fw.write("%12s" % "NA")
            elif assoc[k][i] > 0.001:
                fw.write("%12.4f" % assoc[k][i])
            else:
                fw.write("%12.2e" % assoc[k][i])
        for i in range(12, len(assoc[k])):
            if assoc[k][i] > 100:
                fw.write("%10.0f" %  assoc[k][i])
            else:
                fw.write("%10.2f" %  assoc[k][i])
        fw.write('\n')
    fw.close()
##############   HLAZygosity   #############################################################
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
    OR1 = (n1+0.5) * (n4+0.5) / (n2+0.5) / (n3+0.5)
    p1 = simpleTest(data, test)
    # test2
    n1 = x[2]
    n2 = x[1]
    n3 = y[2]
    n4 = y[1]
    data = [[n1, n2], [n3, n4]]
    OR2 = (n1+0.5) * (n4+0.5) / (n2+0.5) / (n3+0.5)
    p2 = simpleTest(data, test)
    # test 3
    n1 = x[0]
    n2 = x[2]
    n3 = y[0]
    n4 = y[2]
    data = [[n1, n2], [n3, n4]]
    OR3 = (n1+0.5) * (n4+0.5) / (n2+0.5) / (n3+0.5)
    p3 = simpleTest(data, test)
    return [p1, p2, p3, OR1, OR2, OR3]
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
def printZygosity(ans, level):
    print '%-20s' % 'ID',
    for h in ('Hom_P', 'Het_P', 'Zyg_P'):
        print '%12s' % h,
    for h in ('Hom_OR', 'Het_OR', 'Zyg_OR'):
        print '%8s' % h,
    print
    for k in sorted(ans.keys()):
        if level == 'residue':
            print '%-20s' % (k[0] + '_' + str(k[1]) + '_' + k[2]),
        else:
            print '%-20s' % k,
        for i in range(3):
            if ans[k][i] == 'NA':
                print "%12s" % 'NA',
            elif ans[k][i] > 0.001:
                print "%12.4f" % ans[k][i],
            else:
                print "%12.2e" % ans[k][i],
        for i in range(3,6):
            print "%8.4f" % ans[k][i],
        print
def writeZygosity(ans, level, outfile):
    fw = open(outfile, 'w')
    fw.write('%-20s' % 'ID')
    for h in ('Hom_P', 'Het_P', 'Zyg_P'):
        fw.write('%12s' % h)
    for h in ('Hom_OR', 'Het_OR', 'Zyg_OR'):
        fw.write('%8s' % h)
    fw.write('\n')
    for k in sorted(ans.keys()):
        if level == 'residue':
            fw.write('%-20s' % (k[0] + '_' + str(k[1]) + '_' + k[2]))
        else:
            fw.write('%-20s' % k)
        for i in range(3):
            if ans[k][i] == 'NA':
                fw.write("%12s" % 'NA')
            elif ans[k][i] > 0.001:
                fw.write("%12.4f" % ans[k][i])
            else:
                fw.write("%12.2e" % ans[k][i])
        for i in range(3,6):
            fw.write("%8.4f" % ans[k][i])
        fw.write('\n')
    fw.close()
###################### arguments ######################################################################
strattime = time.time()
parser = argparse.ArgumentParser(description='Python for HLA analysis', prog="PyHLA.py")
parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0.0')
parser.add_argument('-V', '--print', help='print output to screen', action='store_true')
parser.add_argument('-i', '--input', help='input file', required=True, type=str)
parser.add_argument('-o', '--out', help='output file', default='output.txt')
parser.add_argument('-d', '--digit', help='digit to test, default 4', default=4, type=int, choices=[2,4,6])
### summary
parser.add_argument('-s', '--summary', help='data summary', action='store_true')
### association analysis
parser.add_argument('-a', '--assoc', help='association analysis', action='store_true')
parser.add_argument('-m', '--model', help='genetic model, default allelic', default='allelic', type=str, choices=['allelic','dom','rec','additive'])
parser.add_argument('-t', '--test', help='statistical test method, default fisher', default='fisher', type=str, choices=['fisher','chisq','logistic','linear'])
parser.add_argument('-f', '--freq', help='minimal frequency, default 0', default=0.05, type=float)
parser.add_argument('-j', '--adjust', help='p value correction, default FDR', default='FDR', type=str,choices=['FDR','FDR_BY','Bonferroni','Holm'])
parser.add_argument('-e', '--exclude', help='exclude alleles file', type=str)
parser.add_argument('-c', '--covar', help='covariants file', type=str)
parser.add_argument('-n', '--covar-name', help='select a particular subset of covariates', type=str)
parser.add_argument('-p', '--perm', help='number of permutation', type=int)
parser.add_argument('-r', '--seed', help='random seed', type=int)
### amino acid association analysis
parser.add_argument('-A', '--assoc-AA', help='amino acid association analysis', action='store_true')
parser.add_argument('-u', '--consensus', help='use the sonsensus amino acid senquence', action='store_true')
### amino acid alignment
parser.add_argument('-w', '--align', help='amino acid senquence alignment', action='store_true')
### zygosity
parser.add_argument('-z', '--zygosity', help='zygosity test', action='store_true')
parser.add_argument('-L', '--level', help='level to test', default='residue', type=str, choices=['residue', 'allele'])
### Interaction
parser.add_argument('-I', '--interaction', help='Interaction test', action='store_true')
##################### parser arguments  ###########################################################
aafile = 'aa.aln.txt'

args = vars(parser.parse_args())

INFILE = args['input'] if 'input' in args else None
OUTFILE = args['out'] if 'out' in args else None
DIGIT = args['digit'] if 'digit' in args else None
PRINT =  args['print'] if 'print' in args else None

SUMMARY = args['summary'] if 'summary' in args else None

ASSOC= args['assoc'] if 'assoc' in args else None
TEST = args['test'] if 'test' in args else None
MODEL = args['model'] if 'model' in args else None
FREQ = args['freq'] if 'freq' in args else None
ADJUST = args['adjust'] if 'adjust' in args else None
SEED = args['seed'] if 'seed' in args else None
PERM = args['perm'] if 'perm' in args else None
COVFILE = args['covar'] if 'covar' in args else None
COVNAME = args['covar_name'] if 'covar_name' in args else None
EXCLUDE = args['exclude'] if 'exclude' in args else None

AAA = args['assoc_AA'] if 'assoc_AA' in args else None
CONSENSUS = args['consensus'] if 'consensus' in args else None

ALN = args['align'] if 'align' in args else None

ZYG = args['zygosity'] if 'zygosity' in args else None
LEVEL = args['level'] if 'level' in args else None

INT = args['interaction'] if 'interaction' in args else None
##########################  log infor  ##################################################################
print "@-------------------------------------------------------------@"
print "|       PyHLA       |     v1.1.1      |      16 Oct 2017      |"
print "|-------------------------------------------------------------|"
print "|  (C) 2017 Felix Yanhui Fan, GNU General Public License, v2  |"
print "|-------------------------------------------------------------|"
print "|    For documentation, citation & bug-report instructions:   |"
print "|          http://felixfan.github.io/PyHLA                    |"
print "@-------------------------------------------------------------@"
print "\n\tOptions in effect:"
print "\t--input", INFILE
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
        mm = ['allelic', 'dom', 'rec']
        if MODEL in mm:
            print "\t--model", MODEL
        else:
            sys.exit('model should be one of {}'.format(mm))
    elif TEST == 'logistic' or 'linear':
        mm = ['additive', 'dom', 'rec']
        if MODEL in mm:
            print "\t--model", MODEL
        else:
            sys.exit('model should be one of {}'.format(mm))
        if COVFILE:
            print "\t--covar", COVFILE
            if COVNAME:
                print "\t--covar-name", COVNAME
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
    print "\t--assoc-AA"
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
################################ run program ############################################################
if SUMMARY:
    if quantTrait(INFILE):
        alleles, genes, n = alleleCount(INFILE, DIGIT)
        if PRINT:
            printSummaryQuant(alleles, genes, n)
        writeSummaryQuant(alleles, genes, n, OUTFILE)
    else:
        caseAlleles, ctrlAlleles, np, nc, nn = allelicCount(INFILE, DIGIT)
        freq, alleles = hlaFreq(caseAlleles, ctrlAlleles, np, nc, nn[0])
        popCase, popCtrl, popP, popC, popN = domCount(INFILE, DIGIT)
        if PRINT:
            printSummary(alleles, freq, caseAlleles, ctrlAlleles, np, nc, nn, popCase, popCtrl, popP, popC)
        writeSummary(alleles, freq, caseAlleles, ctrlAlleles, np, nc, nn, OUTFILE, popCase, popCtrl, popP, popC)
elif ASSOC:
    if quantTrait(INFILE):
        if TEST != 'linear':
            sys.exit("quantitative trait was detected, only linear regression can be used for association analysis!")
    else:
        if TEST =='linear':
            sys.exit("case-control trait was detected, linear regression can not be used for association analysis!")
    if TEST == 'chisq' or TEST == 'fisher':
        if PERM is None:
            assoc = assocADRChiFisher(INFILE, DIGIT, FREQ, TEST,MODEL,ADJUST,EXCLUDE, PERM, SEED)
            if PRINT:
                printAssocChiFisher(assoc, TEST)
            writeAssocChiFisher(assoc, TEST, OUTFILE)
        else:
            assoc, permP, permN, permNA, permNL= assocADRChiFisher(INFILE, DIGIT, FREQ, TEST,MODEL,ADJUST,EXCLUDE, PERM, SEED)
            if PRINT:
                printAssocChiFisher(assoc, TEST, permP, permN, permNA)
            writeAssocChiFisher(assoc, TEST, OUTFILE, permP, permN, permNA)
    elif TEST == 'logistic':
        if PERM is None:
            assoc = regressionLogistic(INFILE, DIGIT, FREQ, MODEL, ADJUST, EXCLUDE, COVFILE, COVNAME, PERM, SEED, TEST)
            if PRINT:
                printLogistic(assoc)
            writeLogistic(assoc, OUTFILE)
        else:
            assoc, permP, permN, permNA = regressionLogistic(INFILE, DIGIT, FREQ, MODEL, ADJUST, EXCLUDE, COVFILE, COVNAME, PERM, SEED, TEST)
            if PRINT:
                printLogistic(assoc, permP, permN, permNA)
            writeLogistic(assoc, OUTFILE, permP, permN, permNA)
    elif TEST == 'linear':
        if PERM is None:
            assoc = regressionLinear(INFILE, DIGIT, FREQ, MODEL, ADJUST, EXCLUDE, COVFILE, COVNAME, PERM, SEED, TEST)
            if PRINT:
                printLinear(assoc)
            writeLinear(assoc, OUTFILE)
        else:
            assoc, permP, permN, permNA = regressionLinear(INFILE, DIGIT, FREQ, MODEL, ADJUST, EXCLUDE, COVFILE, COVNAME, PERM, SEED, TEST)
            if PRINT:
                printLinear(assoc, permP, permN, permNA)
            writeLinear(assoc, OUTFILE, permP, permN, permNA)
elif AAA or ALN:
    if quantTrait(INFILE):
        sys.exit("quantitative trait was detected, only case-control data can be used for amino acid analysis!")
    case, ctrl, caseGeno, ctrlGeno, ncase, nctrl = readGeno(INFILE)
    seq = readAAseq(aafile)
    alleles = keyDicts(case, ctrl)
    if 'NA' in alleles:
        alleles.remove('NA')
    myseq = getSeq(alleles, seq, CONSENSUS)
    if AAA:
        assoc = aaAssoc(case, ctrl, caseGeno, ctrlGeno, ncase, nctrl, myseq, TEST)
        if PRINT:
            printAAA(assoc)
        writeAAA(assoc, OUTFILE)
    elif ALN:
        myaln = aaAlign(case, ctrl, myseq)
        if PRINT:
            printAA(myaln)
        writeAA(myaln, OUTFILE)
elif ZYG or INT:
    if TEST != 'fisher' and TEST != 'chisq':
        sys.exit("only 'fisher' and 'chisq' test can be used!")
    if LEVEL == 'residue':
        case, ctrl, caseGeno, ctrlGeno, ncase, nctrl = readGeno(INFILE)
        seq = readAAseq(aafile)
        alleles = keyDicts(case, ctrl)
        if 'NA' in alleles:
            alleles.remove('NA')
        myseq = getSeq(alleles, seq, CONSENSUS)
        assoc = aaAssoc(case, ctrl, caseGeno, ctrlGeno, ncase, nctrl, myseq, TEST)
        sig = {}
        for k in assoc:
            if assoc[k][4] != 'NA' and assoc[k][4] < 0.05:
                sig[k] = 1
        keys = sorted(sig.keys())
        if ZYG: # zygosity test
            ans = zygosityAA(keys, caseGeno, ctrlGeno, myseq, TEST)
            if PRINT:
                printZygosity(ans, LEVEL)
            writeZygosity(ans, LEVEL,OUTFILE)
        else:    # interaction test
            ans = interactAA(keys, caseGeno, ctrlGeno, myseq, TEST)
            if PRINT:
                printInteract(ans, LEVEL)
            writeInteract(ans, LEVEL, OUTFILE)
    else:
        assoc = assocADRChiFisher(INFILE, DIGIT, FREQ, TEST)
        caseGeno, ctrlGeno = readAlleleZygInteract(INFILE, DIGIT)
        sig = {}
        for k in assoc:
            if assoc[k][7] != 'NA' and assoc[k][7] < 0.05:
                sig[k] = 1
        keys = sorted(sig.keys())
        if  ZYG:
            ans = zygosityAllele(keys, caseGeno, ctrlGeno, TEST)
            if PRINT:
                printZygosity(ans, LEVEL)
            writeZygosity(ans, LEVEL,OUTFILE)
        else:
            ans = interactAllele(keys, caseGeno, ctrlGeno, TEST)
            if PRINT:
                printInteract(ans, LEVEL)
            writeInteract(ans, LEVEL, OUTFILE)
else:
    c = ('--summary', '--assoc', '--assoc-AA', '--align', '--zygosity', '--interaction')
    sys.exit("one of the following option must be used: \n\n%s\n%s\n%s\n%s\n%s\n%s\n" % c)
###############################  time used #############################################################
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
