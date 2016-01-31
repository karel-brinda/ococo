##! /usr/bin/env bash

set -eux
set -o pipefail

../../ococo \
	-m batch \
	-i data/alignments_A_2.sam \
	-f data/fasta_NN.fa \
	-c output/fasta_NA.fa \

diff output/fasta_NA.fa data/fasta_NA.fa
