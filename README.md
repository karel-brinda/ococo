# OCOCO - Online Consensus Caller

[![Build Status](https://travis-ci.org/karel-brinda/ococo.svg?branch=master)](https://travis-ci.org/karel-brinda/ococo)
[![Arxiv](https://img.shields.io/badge/arXiv-1605.09070-green.svg?style=flat)](https://arxiv.org/abs/1605.09070)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](https://anaconda.org/bioconda/ococo)

Welcome to OCOCO, an online consensus caller.

## Prerequisities

* GCC 4.8+ or equivalent
* ZLib

## Getting started

```bash
git clone --recursive https://github.com/karel-brinda/ococo
cd ococo && make -j
./ococo -i test.bam -f test.fa -m real-time --vcf-cons -
```

**Installation:** ``make install``

### Alternative ways of installation

Using Conda
```
conda install ococo
```

Using Brew
```
brew install ococo
```

## Command line parameters

```
Usage:   ococo -i <SAM/BAM file> [other options]

Input options:
  -i, --input FILE      input SAM/BAM file (- for standard input)
  -f, --fasta-ref FILE  initial FASTA reference (otherwise seq of N's is used)
  -s, --stats-in FILE   input statistics

Output options:
  -F, --fasta-cons FILE FASTA file with consensus
  -S, --stats-out FILE  output statistics
  -V, --vcf-cons FILE   VCF file with updates of consensus (- for standard output)
  -P, --pileup FILE     truncated pileup (- for standard output)
  --verbose             verbose mode (report every update of a counter)

Parameters for consensus calling:
  -x, --counters STR    counter configuration: [ococo16]
                           - ococo16 (3b/counter, 16b/position)
                           - ococo32 (7b/counter, 32b/position)
                           - ococo64 (15b/counter, 64b/position)
  -m, --mode STR        mode: [batch]
                           - real-time (updates reported immediately)
                           - batch (updates reported after end of algn stream)
  -t, --strategy STR    strategy for updates: [majority]
                           - majority (update to majority base)
                           - stochastic (update to stochastically drawn base)
                           - no-updates (no updates, only counters updated)
  -q, --min-MQ INT      skip alignments with mapping quality smaller than INT [1]
  -Q, --min-BQ INT      skip bases with base quality smaller than INT [13]
  -w, --ref-weight INT  initial counter value for nucleotides from ref [0]
  -c, --min-cov INT     minimum coverage required for update [2]
  -M, --maj-thres FLOAT majority threshold [0.51]

Examples:
   ococo -i test.bam -f test.fa -m real-time -V -
   ococo -x ococo64 -i test.bam -f test.fa -P - -V variants.vcf
```

## Issues

Please use [Github issues](https://github.com/karel-brinda/ococo/issues).


## Changelog

See [Releases](https://github.com/karel-brinda/ococo/releases).


## Licence

[MIT](https://github.com/karel-brinda/ococo/blob/master/LICENSE)


## Author

[Karel Brinda](http://brinda.cz) <kbrinda@hsph.harvard.edu>


## Citing OCOCO

* K. Brinda, V. Boeva, G. Kucherov. **Dynamic read mapping and online consensus calling for better variant detection.** arXiv:1605.09070v1 [q-bio.GN], 2016. http://arxiv.org/abs/1605.09070
