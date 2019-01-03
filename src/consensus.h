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

#include "counters.h"
#include "params.h"
#include "types.h"

#include <cassert>
#include <cmath>

namespace ococo {

/*
 * Call consensus using the majority rule.
 */
inline char cons_call_maj(const pos_stats_uncompr_t &psu,
                          const params_t &params) {
    char cons = nt16_nt256[psu.nt16];  // initial consensus

    /* Has sufficiently many alignments been collected? */
    if (psu.sum >= params.min_coverage_upd) {
        /* Calculate the minimal required counter value for an
         * update. */
        int32_t required_min =
            static_cast<int32_t>(ceil(params.majority_threshold * psu.sum));

        /* Find the maximal counter with such a value. */
        int32_t max = 0;
        for (int32_t i = 0; i < 4; i++) {
            if (psu.counters[i] >= required_min) {
                if (psu.counters[i] > max) {
                    max  = psu.counters[i];
                    cons = nt4_nt256[i];
                }
            }
        }
    }

    return cons;
}

}  // namespace ococo
