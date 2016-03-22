#pragma once

#include "ococo_types.h"
#include "ococo_misc.h"

#include <cstdio>
#include <iostream>
#include <string>

#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <htslib/faidx.h>
#include <htslib/khash.h>
#include <htslib/kseq.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>
#include <htslib/kseq.h>

/****************************
 *** Consensus parameters ***
 ****************************/

namespace ococo {

enum mode_t { BATCH, REALTIME };

enum strategy_t {
    NO_UPDATES,
    STOCHASTIC,
    STOCHASTIC_AMB,
    MAJORITY,
    MAJORITY_AMB,
    count
};


struct params_t {
    bool correctly_initialized;

    mode_t mode;
    strategy_t strategy;

    std::string command;

    /*
        Input parameters
    */
    std::string sam_fn;
    std::string fasta_in_fn;
    std::string stats_in_fn;


    /*
        Output parameters
    */
    std::string vcf_fn;
    std::string fasta_out_fn;
    std::string stats_out_fn;
    std::string pileup_fn;

    /*
        Files
    */

    FILE *vcf_file;
    FILE *pileup_file;
    FILE *fasta_out_file;
    samFile *sam_file;


    /*
        Consensus calling parameters
    */

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

    /*
        Array of consensus calling functions
    */
    char (*cons_alg[strategy_t::count])(const pos_stats_uncompr_t &psu,
                                        const params_t &params);

    params_t();

    params_t(int argc, const char *argv[]);

    ~params_t();
    
    void parse_commandline(int argc, const char *argv[]);

    void init_default_values();
};

}
