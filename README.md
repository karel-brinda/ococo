# OCOCO - the first online variant and consensus caller

[![Build Status](https://travis-ci.org/karel-brinda/ococo.svg?branch=master)](https://travis-ci.org/karel-brinda/ococo)
[![Arxiv](https://img.shields.io/badge/arXiv-1712.01146-green.svg?style=flat)](https://arxiv.org/abs/1712.01146)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](https://anaconda.org/bioconda/ococo)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1066531.svg)](https://doi.org/10.5281/zenodo.1066531)


## Abstract

**Motivation:** Identifying genomic variants is an essential step for connecting genotype and
phenotype. The usual approach consists of statistical inference of variants
from alignments of sequencing reads. State-of-the-art variant callers can
resolve a wide range of different variant types with high accuracy. However,
they require that all read alignments be available from the beginning of
variant calling and be sorted by coordinates. Sorting is computationally
expensive, both memory- and speed-wise, and the resulting pipelines suffer from
storing and retrieving large alignments files from external memory. Therefore,
there is interest in developing methods for resource-efficient variant calling.

**Results:** We present Ococo, the first program capable of inferring variants in a
real-time, as read alignments are fed in. Ococo inputs unsorted alignments from
a stream and infers single-nucleotide variants, together with a genomic
consensus, using statistics stored in compact several-bit counters. Ococo
provides a fast and memory-efficient alternative to the usual variant calling.
It is particularly advantageous when reads are sequenced or mapped
progressively, or when available computational resources are at a premium.


## Getting started

```bash
git clone --recursive https://github.com/karel-brinda/ococo
cd ococo && make -j
./ococo -i test.bam -f test.fa --vcf-cons -
```

## Citation

> Brinda K, Boeva V, Kucherov G. **Ococo: an online variant and consensus
> caller.** arXiv:1712.01146 [q-bio.GN], 2017. https://arxiv.org/abs/1712.01146


## Installation

### From Bioconda

```
conda install -c bioconda ococo
```

### Building from source

**Prerequisities**

* GCC 4.8+ or equivalent
* ZLib

**Compilation:** ``make``

**Installation:** ``make install``


## Command line parameters

<!---
USAGE-BEGIN
-->
```Usage:   ococo -i <SAM/BAM file> [other options]

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

Note:
   For more details, see the manual page 'man ./ococo.1'.

```
<!---
USAGE-END
-->


## Issues

Please use [Github issues](https://github.com/karel-brinda/ococo/issues).


## Changelog

See [Releases](https://github.com/karel-brinda/ococo/releases).


## Licence

[MIT](https://github.com/karel-brinda/ococo/blob/master/LICENSE)


## Author

[Karel Brinda](http://brinda.cz) \<kbrinda@hsph.harvard.edu\>


