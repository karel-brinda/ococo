.PHONY: all help clean unaligned aligned1 aligned2 import_export

SHELL=/usr/bin/env bash -eo pipefail

.SECONDARY:

all: unaligned aligned_AA aligned_A_A noparams import_export

OCOCO=../ococo

noparams:
	$(OCOCO) -h
	$(OCOCO) -v


# NN
# ++
#      AA
#      AA
# ==
# NN
unaligned:
	@printf "\n\n\n===== $@\n\n\n"
	$(OCOCO) \
		-m batch \
		-f data/fasta_NN.fa \
		-i data/alignment_AA_unm.sam \
		-F output/fasta_NN.fa \
		-V -
	diff output/fasta_NN.fa data/fasta_NN.fa


# NN
# ++
#  A
#  A
# ==
# NA
aligned_A_A:
	@printf "\n\n\n===== $@\n\n\n"
	$(OCOCO) \
		-m batch \
		-f data/fasta_NN.fa \
		-i data/alignment_A_2.sam \
		-F output/fasta_NA.fa \
		-V -
	diff output/fasta_NA.fa data/fasta_NA.fa


# NN
# ++
# AA
# ==
# AA
aligned_AA:
	@printf "\n\n\n===== $@\n\n\n"
	$(OCOCO) \
		-m batch \
		-f data/fasta_NN.fa \
		-i data/alignment_AA_1.sam \
		-F output/fasta_AA.fa \
		-V -
	diff output/fasta_AA.fa data/fasta_AA.fa


# NN
# ++
#  C
#    save&load
#  A
#  A
#    save&load
#  A
#  A
# ==
# NA
import_export:
	@printf "\n\n\n===== $@\n\n\n"
	rm -f output/stats.ococo
	$(OCOCO) \
		-m batch \
		-f data/fasta_NN.fa \
		-i data/alignment_C_2.sam \
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

