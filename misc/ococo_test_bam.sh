#! /usr/bin/env bash

set -e
set -o xtrace
set -o pipefail


if [ $# -ne 3 ]; then
	echo
	echo "illegal number of parameters, 3 are expected ($# are provided):"
	echo
	echo "	- BAM file (possibly unsorted)"
	echo "	- minimum mapping quality"
	echo "	- minimum base quality"
	echo
	exit 1
fi

script_dir=`dirname "$0"`
tmp_dir=$script_dir/"tmp"
unsorted_bam=$1
sorted_bam=$tmp_dir/bam.bam
min_mq=$2
min_bq=3
vcf=$tmp_dir/vcf.vcf
pileup=$tmp_dir/pileup.pileup

# create tmp dir
rm -fr $tmp_dir
mkdir -p $tmp_dir

# sort reads (BAM)
samtools sort -T ${unsorted_bam}.tmp -o $sorted_bam $unsorted_bam

# index bam
samtools index $sorted_bam

# pileup
samtools mpileup --min-BQ $min_bq --min-MQ $min_mq $sorted_bam > $pileup

# vcf
ococo -i $sorted_bam -v $vcf --min-BQ $min_bq --min-MQ $min_mq

# compare 'ococo output' vs. 'samtools mpileup output'
$script_dir/_ococo_debug_compare_vcf_pileup.py $vcf $pileup
