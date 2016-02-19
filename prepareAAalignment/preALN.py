#!/usr/bin/env python

import re

def oneline(infile, outfile):
	pattern = re.compile(r'^(\w+)(\*){1}(\d+)(\:?)')
	f = open(infile)
	geno = {}
	for i in f:
		i = i.strip()
		if i and pattern.search(i):
			fs = i.split()
			if fs[0] not in geno:
				geno[fs[0]]=''
			for j in range(1, len(fs)):
				geno[fs[0]] += fs[j]
	f.close()

	f = open(outfile, 'w')
	keys = sorted(geno.keys())
	index = 0
	ref = ''
	for a in keys:
		index += 1
		if index == 1:
			f.write(a)
			f.write('\t')
			f.write(geno[a])
			f.write('\n')
			ref = geno[a]
		else:
			f.write(a)
			f.write('\t')
			tmp = geno[a]
			for k in range(len(tmp)):
				if tmp[k] == '*':
					f.write('*')
				elif tmp[k] == '-':
					f.write(ref[k])
				else:
					f.write(tmp[k])
			if len(ref) > len(tmp):
				for k in range(len(tmp),len(ref)):
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
	