##! /usr/bin/env bash

set -eu
set -o pipefail
set -o xtrace

rm -f output/stats.ococo

../../ococo \
	-m batch \
	-i data/alignment_C_2.sam \
	-f data/fasta_NN.fa \
	-F output/fasta_NA.fa \
	-S output/stats.ococo \
	-t majority \
	-v - \

diff output/fasta_NA.fa data/fasta_NC.fa

echo
echo "==============================="
echo

../../ococo \
	-m batch \
	-i data/alignment_A_2.sam \
	-F output/fasta_NA.fa \
	-s output/stats.ococo \
	-S output/stats.ococo \
	-t majority \
	-v - \

echo
echo "==============================="
echo

../../ococo \
	-m batch \
	-i data/alignment_A_2.sam \
	-F output/fasta_NA.fa \
	-s output/stats.ococo \
	-S output/stats.ococo \
	-t majority \
	-v - \

echo
echo "==============================="
echo

diff output/fasta_NA.fa data/fasta_NA.fa
