#! /usr/bin/env bash

set -eux
set -o pipefail

ococo \
	-m batch \
	-i data/alignments_AA_1.sam \
	-f data/fasta_NN.fa \
	-c output/fasta_AA.fa \

diff output/fasta_AA.fa data/fasta_AA.fa
