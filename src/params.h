/* The MIT License

   Copyright (c) 2016-2019 Karel Brinda (kbrinda@hsph.harvard.edu)

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#pragma once

#include "misc.h"
#include "types.h"

#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>

#include <htslib/faidx.h>
#include <htslib/khash.h>
#include <htslib/kseq.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>

/****************************
 *** Consensus parameters ***
 ****************************/

namespace ococo {

const int default_c   = 2;
const float default_M = 0.51;
const int default_w   = 0;
const int default_q   = 1;
const int default_Q   = 13;
const int default_C   = -1;

enum mode_t { BATCH, REALTIME };

enum strategy_t {
    NO_UPDATES,
    STOCHASTIC,
    STOCHASTIC_AMB,
    MAJORITY,
    MAJORITY_AMB,
    count
};

enum counter_configuration_t {
    OCOCO16,
    OCOCO32,
    OCOCO64,
};

struct params_t {
    bool correctly_initialized;
    int return_code;

    std::string command;

    /*
     * Counter parameters
     */
    counter_configuration_t counter_configuration;
    std::string counters_str;
    std::string counters_str_descr;
    int32_t stats_bits_per_position;
    int32_t stats_bits_per_nucleotide;

    /*
     * Input parameters
     */
    std::string in_sam_fn;
    std::string in_fasta_fn;
    std::string in_stats_fn;

    /*
     * Output parameters
     */
    bool verbose;

    std::string out_sam_fn;
    std::string out_vcf_fn;
    std::string out_fasta_fn;
    std::string out_stats_fn;
    std::string out_pileup_fn;
    std::string out_log_fn;

    /* Filter alignments when coverage is greater than */
    int coverage_filter;

    /*
     * Files
     */
    FILE *out_vcf_file;
    FILE *out_pileup_file;
    FILE *out_fasta_file;
    FILE *out_log_file;
    samFile *out_sam_file;

    samFile *in_sam_file;

    /*
     * Consensus calling parameters
     */

    mode_t mode;
    strategy_t strategy;

    /* minimum mapping quality for update */
    int32_t min_mapq;

    /* minimum base quality for update */
    int32_t min_baseq;

    /* initial values for counters corresponding to ref */
    int32_t init_ref_weight;

    /* minimum coverage for update (does not include init_ref_weight */
    int32_t min_coverage;

    /* threshold for having majority */
    double majority_threshold;

    /* auxiliary */
    std::string strategy_str;
    std::string mode_str;
    int64_t n_upd;

    params_t();

    params_t(int argc, const char **argv);

    void parse_commandline(int argc, const char **argv);

    void print_help();

    void init_default_values();
};
}  // namespace ococo
