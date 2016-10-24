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

    /* filter alignments when coverage is greater than */
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

    /*
     * Array of consensus calling functions
     */
    char (*cons_alg[strategy_t::count])(const pos_stats_uncompr_t &psu,
                                        const params_t &params);

    params_t();

    params_t(int argc, const char **argv);

    ~params_t();

    void parse_commandline(int argc, const char **argv);

    void print_help();

    void init_default_values();
};
}
