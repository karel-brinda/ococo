#! /usr/bin/env bash

set -eux
set -o pipefail

../../ococo \
	-m batch \
	-i data/alignments_AA_unm.sam \
	-f data/fasta_NN.fa \
	-c output/fasta_NN.fa \

