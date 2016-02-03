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
	
	#a=bases.count("a")+bases.count("A")
	#c=bases.count("c")+bases.count("C")
	#g=bases.count("g")+bases.count("G")
	#t=bases.count("t")+bases.count("T")

	dict_freq= {
		"a": 0,
		"c": 0,
		"g": 0,
		"t": 0,
		"n": 0,
	}

	number=0
	number_string=None

	circ=False
	for x in bases:
		if circ==True:
			circ=False
			continue
		if x=="^":
			circ=True
		elif x in "+-":
			number_string=""
			assert number==0
		elif x in "0123456789":
			assert number_string is not None, "'{}'".format(pileup_line.strip())
			number_string+=x
		elif x in "ACGTNacgtn":
			if number_string is not None:
				number=int(number_string)

			if number==0:
				dict_freq[x.lower()]+=1
			else:
				number-=1
				number_string=None

	assert number==0, "'{}'".format(pileup_line.strip())

	return [ref,pos,(dict_freq["a"],dict_freq["c"],dict_freq["g"],dict_freq["t"])]


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

			if vcf_rec is not None and pos==vcf_rec[1]:
				vcf_info="({},{},{},{})".format(vcf_rec[2][0],vcf_rec[2][1],vcf_rec[2][2],vcf_rec[2][3])
				vcf_rec=None
			else:
				vcf_info="( , , , )"
				eq="?"

			if pileup_rec is not None and pos==pileup_rec[1]:
				pileup_info="({},{},{},{})".format(pileup_rec[2][0],pileup_rec[2][1],pileup_rec[2][2],pileup_rec[2][3])
				pileup_rec=None
			else:
				pileup_info="( , , , )"
				eq="?"

			if eq is None:
				eq = "=" if pileup_info==vcf_info else "X"

			print("\t".join([seq_name,str(pos),eq,pileup_info,vcf_info]))
			pos+=1
