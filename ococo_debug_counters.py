#! /usr/bin/env python3

import sys
import re

counters_re=re.compile(r"CS=([0-9]+),([0-9]+),([0-9]+),([0-9]+)")

# ../../ococo/ococo_debug_counters.py <(bcftools view ./exp.1.02__Borrelia__0.07-ococo/2_alignments.dyn/5_vcf/00000.vcf.gz ) <(samtools mpileup ./exp.1.01__Borrelia__0.07/2_alignments.dyn/3.2_sorted_bam/00000.bam)

if len(sys.argv)!=3:
	print()
	print("Usage: ococo_debug_counter.py infile.vcf infile.pileup");
	print()
	print("\tinfile.vcf     - output of ococo")
	print("\tinfile.pileup  - output of samtools mpileup")
	print()
	sys.exit(-1)

vcf_fn=sys.argv[1]
pileup_fn=sys.argv[2]

def next_vcf_rec(vcf_fo):
	vcf_line=None

	while vcf_line==None or vcf_line[0]=="#" or vcf_line.strip()=="":
		vcf_line=vcf_fo.readline()
		#print(vcf_line)
		if len(vcf_line)==0:
			return [None,None,(None,None,None,None)]

	parts=vcf_line.split("\t")
	ref=parts[0]
	pos=int(parts[1])
	rm=counters_re.match(parts[7])
	counters=(rm.group(1),rm.group(2),rm.group(3),rm.group(4))
	return [ref,pos,counters]

def next_pileup_rec(pileup_fo):
	pileup_line=None
	
	while pileup_line==None or pileup_line.strip()=="":
		pileup_line=pileup_fo.readline()
		if len(pileup_line)==0:
			return None
	parts=pileup_line.split("\t")
	ref=parts[0]
	pos=int(parts[1])
	bases=parts[4]
	
	a=bases.count("a")+bases.count("A")
	c=bases.count("c")+bases.count("C")
	g=bases.count("g")+bases.count("G")
	t=bases.count("t")+bases.count("T")

	return [ref,pos,(a,c,g,t)]


with open(vcf_fn) as vcf_fo:
	with open(pileup_fn) as pileup_fo:
		pos=0
		vcf_rec=[]
		pileup_rec=[]
		seq_name=None

		while 1:
			
			# fill in buffers
			if vcf_rec == [] or vcf_rec is None:
				vcf_rec=next_vcf_rec(vcf_fo)

			if pileup_rec == [] or pileup_rec is None:
				pileup_rec=next_pileup_rec(pileup_fo)

			if vcf_rec is None and pileup_rec is None:
				break

			if seq_name is None:
				seq_name=vcf_rec[0]

			
			# print record
			eq=None

			if pos==vcf_rec[1]:
				vcf_info="({},{},{},{})".format(vcf_rec[2][0],vcf_rec[2][1],vcf_rec[2][2],vcf_rec[2][3])
				vcf_rec=None
			else:
				vcf_info="( , , , )"
				eq="?"

			if pos==pileup_rec[1]:
				pileup_info="({},{},{},{})".format(pileup_rec[2][0],pileup_rec[2][1],pileup_rec[2][2],pileup_rec[2][3])
				pileup_rec=None
			else:
				pileup_info="( , , , )"
				eq="?"

			if eq is None:
				eq = "=" if pileup_info==vcf_info else "X"

			print("\t".join([seq_name,str(pos),eq,pileup_info,vcf_info]))
			pos+=1

