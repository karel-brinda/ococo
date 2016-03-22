##! /usr/bin/env bash

set -eu
set -o xtrace
set -o pipefail

for strat in "majority" "stochastic";
do
	../../ococo \
		-m batch \
		-i data/alignment_A_2.sam \
		-f data/fasta_NN.fa \
		-F output/fasta_NA.fa \
		-t $strat \
		-V - \

	diff output/fasta_NA.fa data/fasta_NA.fa

	echo
	echo "==============================="
	echo
done