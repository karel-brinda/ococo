#pragma once

#include "ococo_types.h"

#include <string>

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


struct consensus_params_t {
    mode_t mode;
    strategy_t strategy;

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

    FILE *vcf_file = nullptr;
    FILE *pileup_file = nullptr;
    FILE *fasta_out_file = nullptr;


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
                                        const consensus_params_t &params);

    consensus_params_t();
};

}
