#!/usr/bin/env bats

export PATH=$PATH:..

@test "Test of ococo without parameters" {

	run ./test_no_parameters

	[ "$status" -ne 0 ]

}

@test "Test of unaligned reads" {
	rm -f output/*.fa

	run ./test_unaligned_reads.sh

	[ "$status" -eq 0 ]
}

@test "Test aligned reads of length 1" {
	rm -f output/*.fa

	run ./test_aligned_reads_len1.sh

	[ "$status" -eq 0 ]
}

@test "Test aligned reads of length 2" {
	rm -f output/*.fa

	run ./test_aligned_reads_len2.sh

	[ "$status" -eq 0 ]
}

@test "Test of import/export" {
	rm -f output/*.fa

	run ./test_stats_import_export.sh

	[ "$status" -eq 0 ]
}
