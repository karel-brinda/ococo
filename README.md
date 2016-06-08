# OCOCO - Online Consensus Caller

[![Build Status](https://travis-ci.org/karel-brinda/ococo.svg?branch=master)](https://travis-ci.org/karel-brinda/ococo)
[![Arxiv](https://img.shields.io/badge/arXiv-1605.09070-green.svg?style=flat)](https://arxiv.org/abs/1605.09070)

Welcome to OCOCO, an online consensus caller.

## Prerequisities

* G++ 4.9+ or equivalent
* CMake (http://cmake.org/)
* Boost (http://www.boost.org/)
* Git (https://git-scm.com/)

## Getting started

```bash
git clone https://github.com/karel-brinda/ococo
cd ococo && cmake . && make
./ococo -i test.bam -f test.fa -m real-time --vcf-cons -
```

## Command line parameters

```
Generic options:
  -v [ --version ]                      Print version and exit.
  -h [ --help ]                         Print this message and exit.

Input options:
  -i [ --input ] arg                    Input SAM/BAM file (- for standard
                                        input).
  -f [ --fasta-ref ] arg                Initial FASTA reference (if not
                                        provided, sequence of N's is considered
                                        as the reference).
  -s [ --stats-in ] arg                 Input statistics.

Output options:
  -F [ --fasta-cons ] arg               FASTA file with consensus.
  -S [ --stats-out ] arg                Outputs statistics.
  -V [ --vcf-cons ] arg                 VCF file with updates of consensus (-
                                        for standard output).
  -P [ --pileup ] arg                   Truncated pileup (- for standard
                                        output).
  --verbose                             Verbose mode.

Parameters of consensus calling:
  -x [ --counters ] arg (=ococo16)      Counters configuration:
                                         - ococo16 (3 bits per counter)
                                         - ococo32 (7 bits per counter)
                                         - ococo64 (15 bits per counter)
  -m [ --mode ] arg (=batch)            Mode: real-time / batch.
  -t [ --strategy ] arg (=majority)     Strategy for updates: no-updates /
                                        majority / stochastic.
  -a [ --allow-amb ]                    Allow updates to ambiguous nucleotides.
  -q [ --min-MQ ] arg (=1)              Skip alignments with mapping quality
                                        smaller than INT.
  -Q [ --min-BQ ] arg (=13)             Skip bases with base quality smaller
                                        than INT.
  -w [ --ref-weight ] arg (=0)          Initial counter value for nucleotides
                                        from the reference.
  -c [ --min-coverage ] arg (=2)        Minimum coverage required for update.
  -M [ --majority-threshold ] arg (=0.59999999999999998)
                                        Majority threshold.
```

## Citing OCOCO

* K. Brinda, V. Boeva, G. Kucherov. Dynamic read mapping and online consensus calling for better variant detection. arXiv:1605.09070v1 [q-bio.GN], 2016. http://arxiv.org/abs/1605.09070
