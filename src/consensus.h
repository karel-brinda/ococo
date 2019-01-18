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
#include "types.h"

#include <cmath>

namespace ococo {

/*
 * Call consensus using the majority rule.
 */
inline char cons_call_maj(const PosStats &ps, int32_t min_coverage_upd,
                          float majority_threshold) {
    char cons = nt16_nt256[ps.nt16_];  // initial consensus

    /* Has sufficiently many alignments been collected? */
    if (ps.sum_ >= min_coverage_upd) {
        /* Calculate the minimal required counter value for an
         * update. */
        int32_t required_min = ceil(majority_threshold * ps.sum_);

        /* Find the maximal counter with such a value. */
        int32_t max = 0;
        for (int32_t i = 0; i < 4; i++) {
            if (ps.counters_[i] >= required_min) {
                if (ps.counters_[i] > max) {
                    max  = ps.counters_[i];
                    cons = nt4_nt256[i];
                }
            }
        }
    }

    return cons;
}

}  // namespace ococo
