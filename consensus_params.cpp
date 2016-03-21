#include "consensus_params.h"
#include "consensus_functions.h"

/****************************
 *** Consensus parameters ***
 ****************************/

ococo::consensus_params_t::consensus_params_t()
    : mode(BATCH), strategy(MAJORITY), min_mapq(1), min_baseq(13),
      init_ref_weight(0), min_coverage(2), majority_threshold(0.60) {
    cons_alg[strategy_t::NO_UPDATES] = &cons_call_no_updates;
    cons_alg[strategy_t::STOCHASTIC] = &cons_call_stoch;
    cons_alg[strategy_t::STOCHASTIC_AMB] = &cons_call_stoch_amb;
    cons_alg[strategy_t::MAJORITY] = &cons_call_maj;
    cons_alg[strategy_t::MAJORITY_AMB] = &cons_call_maj_amb;
}
