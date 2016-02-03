#! /usr/bin/env python3

import sys
import re

counters_re=re.compile(r"CS\=\(([0-9]+),([0-9]+),([0-9]+),([0-9]+)\)")

if sys.argc!=2:
	print("Usage: ococo_debug_counter.py infile.vcf infile.pileup");
	print()
	print("infile.vcf     - output of ococo")
	print("infile.pileup  - output of samtools mpileup")
	print()
	sys.exit(-1)

vcf_fn=sys.argv[1]
pileup_fn=sys.argv[2]

def next_vcf_rec(vcf_fo):
	vcf_line=None

	while vcf_line==None or vcf_line[0]=="#" or vcf_line.strip()=="":
		vcf_line=vcf_fo.readline()
		if len(vcf_line)==0:
			return [None,None,(None,None,None,None)]

	parts=vcf_line.parts("\t")
	ref=parts[0]
	pos=parts[1]
	rm=counters.re.match(parts[7])
	counters=(rm[1],rm[2],rm[3],rm[4])
	return [ref,pos,counters]

def next_pileup_rec(pileup_fo):
	pileup_line=None
	
	while pileup_line==None or pileup_line.strip()=="":
		pileup_line=pileup_fo.readline()
		if len(pileup_line)==0:
			return [None,None,(None,None,None,None)]
	parts=pileup_line.parts("\t")
	ref=parts[0]
	pos=parts[1]
	bases=parts[4]
	
	a=bases.count("a")+bases.count("A")
	c=bases.count("c")+bases.count("C")
	g=bases.count("g")+bases.count("G")
	t=bases.count("t")+bases.count("T")

	return [ref,pos,(a,c,g,t)]


with open(vcf_fn) as vcf_fo:
	with open(pileup_fn_fo):
		pos=0
		vcf_rec=None
		pileup_rec=None
		seq_name=None
		while vcf_fo.not_end() or pileup_fo.not_end():
			
			# fill in buffers
			if vcf_rec is None:
				vcf_rec=next_vcf_rec(vcf_fo)

			if pileup_rec is None:
				pileup_rec=next_vcf_rec(vcf_fo)
			
			# update sequence names & counters
			if vcf_rec[0] == pileup_rec[0] and vcf_rec[0]!=seq_name and vcf_rec[0]!=None:
				pos=1
				seq_name=vcf_rec[0]

			# print record
			if seq_name==vcf_rec[0] and position==vcf_rec[1]:
				vcf_info="({},{},{},{})".format(vcf_rec[2][0],vcf_rec[2][1],vcf_rec[2][2],vcf_rec[2][3])
				vcf_rec=None
			else:
				vcf_info="(?,?,?,?)"

			if seq_name==pileup_rec[0] and position==pileup_rec[1]:
				pileup_info="({},{},{},{})".format(pileup_rec[2][0],pileup_rec[2][1],pileup_rec[2][2],pileup_rec[2][3])
				pileup_rec=None
			else:
				pileup_info="(?,?,?,?)"

						
			print("\t".join(seq_name,str(pos),pileup_info,vcf_info))
			pos+=1

