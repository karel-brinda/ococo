.PHONY: all help clean

SHELL=/usr/bin/env bash -eo pipefail

.SECONDARY:

all:
	ococo \
		-m batch \
		-i ../data/alignment_AA_unm.sam \
		-f ../data/fasta_NN.fa \
		-F _fasta_NN.fa \
		-V - \

	diff _fasta_NN.fa ../data/fasta_NN.fa


help: ## Print help message
	@echo "$$(grep -hE '^\S+:.*##' $(MAKEFILE_LIST) | sed -e 's/:.*##\s*/:/' -e 's/^\(.\+\):\(.*\)/\\x1b[36m\1\\x1b[m:\2/' | column -c2 -t -s : | sort)"

clean: ## Clean


