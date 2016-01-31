#! /usr/bin/env bash

set -eu
set -o xtrace
set -o pipefail

for strat in "majority" "stochastic";
do
	../../ococo \
		-m batch \
		-i data/alignment_AA_unm.sam \
		-f data/fasta_NN.fa \
		-c output/fasta_NN.fa \
		-S $strat \
		-v - \

	diff output/fasta_NN.fa data/fasta_NN.fa

	echo
	echo "==============================="
	echo
done
