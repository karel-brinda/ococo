#!/usr/bin/env bats

export PATH=$PATH:../..

@test "Test of ococo without parameters" {
	run ococo

	[ "$status" -ne 0 ]

}

@test "Test of unaligned reads" {
	rm -f output/*.fa
	run ococo \
		-m batch \
		-i data/alignments_AA_unm.sam \
		-f data/fasta_NN.fa \
		-c output/fasta_NN.fa

	[ "$status" -eq 0 ]

	diff output/fasta_NN.fa data/fasta_NN.fa
}

@test "Test aligned reads of length 1" {
	rm -f output/*.fa
	run ococo \
		-m batch \
		-i data/alignments_A_2.sam \
		-f data/fasta_NN.fa \
		-c output/fasta_NA.fa

	[ "$status" -eq 0 ]

	diff output/fasta_NA.fa data/fasta_NA.fa
}

@test "Test aligned reads of length 2" {
	rm -f output/*.fa
	run ococo \
		-m batch \
		-i data/alignments_AA_1.sam \
		-f data/fasta_NN.fa \
		-c output/fasta_AA.fa

	[ "$status" -eq 0 ]

	diff output/fasta_AA.fa data/fasta_AA.fa
}
