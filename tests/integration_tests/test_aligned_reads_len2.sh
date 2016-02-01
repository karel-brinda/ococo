#! /usr/bin/env bash

set -eu
set -o xtrace
set -o pipefail

for strat in "majority" "stochastic";
do
	../../ococo \
		-m batch \
		-i data/alignment_AA_1.sam \
		-f data/fasta_NN.fa \
		-F output/fasta_AA.fa \
		-t $strat \
		-v - \

	diff output/fasta_AA.fa data/fasta_AA.fa

	echo
	echo "==============================="
	echo
done
