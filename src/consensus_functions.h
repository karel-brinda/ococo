#pragma once

#include "ococo_types.h"
#include "consensus_params.h"

#include "cassert"
#include "cmath"

namespace ococo {

inline char cons_call_no_updates(const pos_stats_uncompr_t &psu,
                                 const consensus_params_t &params) {
    return nt16_nt256[psu.nt16];
}

inline char cons_call_stoch(const pos_stats_uncompr_t &psu,
                            const consensus_params_t &params) {
    if (psu.sum == 0) {
        return nt16_nt256[psu.nt16];
    }

    if (psu.nt16 != nt256_nt16[static_cast<uint32_t>('N')]) {
        if (psu.sum < params.min_coverage + params.init_ref_weight) {
            return nt16_nt256[psu.nt16];
        }
    }

    const int32_t prefsum[] = {
        psu.counters[0], psu.counters[0] + psu.counters[1],
        psu.counters[0] + psu.counters[1] + psu.counters[2],
        psu.counters[0] + psu.counters[1] + psu.counters[2] + psu.counters[3]};

    assert(prefsum[3] == psu.sum);

    const int32_t rn = rand() % psu.sum;
    for (int32_t i = 0; i < 4; i++) {
        if (rn < prefsum[i]) {
            return nt4_nt256[i];
        }
    }

    return 'n';
}

inline char cons_call_stoch_amb(const pos_stats_uncompr_t &psu,
                                const consensus_params_t &params) {
    if (psu.sum == 0) {
        return nt16_nt256[psu.nt16];
    }

    if (psu.nt16 != nt256_nt16[static_cast<uint32_t>('N')]) {
        if (psu.sum < params.min_coverage + params.init_ref_weight) {
            return nt16_nt256[psu.nt16];
        }
    }

    nt16_t nucl_nt16 = nt256_nt16[static_cast<uint32_t>('N')];

    while (nucl_nt16 == nt256_nt16[static_cast<uint32_t>('N')]) {
        nucl_nt16 = 0;
        for (int32_t i = 0; i < 4; i++) {
            const int32_t rn = rand() % psu.sum;

            if (rn < psu.counters[i]) {
                nucl_nt16 |= nt4_nt16[i];
            }
        }
    }

    return nt16_nt256[nucl_nt16];
}

inline char cons_call_maj(const pos_stats_uncompr_t &psu,
                          const consensus_params_t &params) {
    if (psu.sum == 0) {
        return nt16_nt256[psu.nt16];
    }

    if (psu.nt16 != nt256_nt16[static_cast<uint32_t>('N')]) {
        if (psu.sum < params.min_coverage + params.init_ref_weight) {
            return nt16_nt256[psu.nt16];
        }
    }

    char nucl_nt256 = nt16_nt256[psu.nt16];

    int32_t required_min =
        static_cast<int32_t>(ceil(params.majority_threshold * psu.sum));
    int32_t max = 0;
    for (int32_t i = 0; i < 4; i++) {
        if (psu.counters[i] >= required_min) {
            if (psu.counters[i] > max) {
                max = psu.counters[i];
                nucl_nt256 = nt4_nt256[i];
            }
        }
    }

    return nucl_nt256;
}

inline char cons_call_maj_amb(const pos_stats_uncompr_t &psu,
                              const consensus_params_t &params) {
    if (psu.sum == 0) {
        return nt16_nt256[psu.nt16];
    }

    if (psu.nt16 != nt256_nt16[static_cast<uint32_t>('N')]) {
        if (psu.sum < params.min_coverage + params.init_ref_weight) {
            return nt16_nt256[psu.nt16];
        }
    }

    char nucl_nt16 = psu.nt16;

    int32_t required_min =
        static_cast<int32_t>(round(params.majority_threshold * psu.sum));
    int32_t max = 0;
    for (int32_t i = 0; i < 4; i++) {
        if (psu.counters[i] >= required_min) {
            if (psu.counters[i] > max) {
                max = psu.counters[i];
                nucl_nt16 = nt4_nt16[i];
            } else if (psu.counters[i] >= max) {
                nucl_nt16 |= nt4_nt16[i];
            }
        }
    }

    return nt16_nt256[static_cast<int32_t>(nucl_nt16)];
}

}
