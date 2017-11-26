.PHONY: all help clean unaligned aligned1 aligned2 import_export

SHELL=/usr/bin/env bash -eo pipefail

.SECONDARY:

all: unaligned aligned1 aligned2 noparams import_export

OCOCO=../ococo

noparams:
	$(OCOCO) -h
	$(OCOCO) -v


unaligned:
	$(OCOCO) \
		-m batch \
		-i data/alignment_AA_unm.sam \
		-f data/fasta_NN.fa \
		-F output/fasta_NN.fa \
		-V -
	diff output/fasta_NN.fa data/fasta_NN.fa


aligned1:
	$(OCOCO) \
		-m batch \
		-i data/alignment_A_2.sam \
		-f data/fasta_NN.fa \
		-F output/fasta_NA.fa \
		-V -
	diff output/fasta_NA.fa data/fasta_NA.fa


aligned2:
	$(OCOCO) \
		-m batch \
		-i data/alignment_AA_1.sam \
		-f data/fasta_NN.fa \
		-F output/fasta_AA.fa \
		-V -
	diff output/fasta_AA.fa data/fasta_AA.fa


import_export:
	rm -f output/stats.ococo
	$(OCOCO) \
		-m batch \
		-i data/alignment_C_2.sam \
		-f data/fasta_NN.fa \
		-F output/fasta_NA.fa \
		-S output/stats.ococo \
		-V -
	diff output/fasta_NA.fa data/fasta_NC.fa
	#
	#
	$(OCOCO) \
		-m batch \
		-i data/alignment_A_2.sam \
		-F output/fasta_NA.fa \
		-s output/stats.ococo \
		-S output/stats.ococo \
		-V -
	$(OCOCO) \
		-m batch \
		-i data/alignment_A_2.sam \
		-F output/fasta_NA.fa \
		-s output/stats.ococo \
		-S output/stats.ococo \
		-V -
	diff output/fasta_NA.fa data/fasta_NA.fa


help: ## Print help message
	@echo "$$(grep -hE '^\S+:.*##' $(MAKEFILE_LIST) | sed -e 's/:.*##\s*/:/' -e 's/^\(.\+\):\(.*\)/\\x1b[36m\1\\x1b[m:\2/' | column -c2 -t -s : | sort)"


clean: ## Clean
	rm -f output/*

