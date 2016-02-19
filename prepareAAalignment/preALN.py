#!/usr/bin/env python

import re

def oneline(infile, outfile):
	pattern = re.compile(r'^(\w+)(\*){1}(\d+)(\:?)')
	f = open(infile)
	geno = {}
	ref = ''
	flag = False
	for i in f:
		i = i.strip()
		if i and pattern.search(i):
			fs = i.split()
			if fs[0] not in geno:
				geno[fs[0]]=''
			for j in range(1, len(fs)):
				geno[fs[0]] += fs[j]
			if not flag: # the first allele is reference
				ref = fs[0]
				flag = True
	f.close()

	f = open(outfile, 'w')
	keys = sorted(geno.keys())
	for a in keys:
		if a == ref:
			f.write(a)
			f.write('\t')
			f.write(geno[a])
			f.write('\n')
		else:
			f.write(a)
			f.write('\t')
			tmp = geno[a]
			for k in range(len(tmp)):
				if tmp[k] == '*':
					f.write('*')
				elif tmp[k] == '-':
					f.write(geno[ref][k])
				else:
					f.write(tmp[k])
			if len(geno[ref]) > len(tmp):
				for k in range(len(tmp),len(geno[ref])):
					f.write('*')
			f.write('\n')

def catfiles(fileList, outfile):
	fw = open(outfile, 'w')
	for f in fileList:
		fr = open(f)
		for r in fr:
			fw.write(r)
		fr.close()
	fw.close()

if __name__ == '__main__':
	infiles = ['A_prot.txt', 'B_prot.txt', 'C_prot.txt', 'DMA_prot.txt', 'DMB_prot.txt', 'DOA_prot.txt', 'DOB_prot.txt','DPA_prot.txt', 'DPB_prot.txt', 'DQA_prot.txt', 'DQB_prot.txt', 'DRA_prot.txt','DRB_prot.txt']
	outfiles = ['A.aln', 'B.aln', 'C.aln', 'DMA.aln','DMB.aln','DOA.aln','DOB.aln','DPA.aln', 'DPB.aln', 'DQA.aln', 'DQB.aln', 'DRA.aln','DRB.aln']
	for i in range(0, len(infiles)):
		oneline(infiles[i], outfiles[i])
	catfiles(outfiles, "aa.aln.txt")
	