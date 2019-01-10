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

#include <cassert>
#include <iomanip>
#include <iostream>

#include "misc.h"
#include "types.h"

namespace ococo {

template <typename T>
constexpr int counter_size() {
    return (8 * sizeof(T) - 4) / 4;
}

/*! @struct
    @abstract          Structure for pileup statistics for 1 position.
    @field nt16        Consensus base.
                       0x00 and 0x0f = N
                       0x01, 0x02, 0x04, 0x08 = bases
                       other values: error
    @field counters    Nucleotide counters.
    @field sum         Sum of the nucleotide counters.
    @field bitshifted  Already bit-shifted? (i.e., in-exact).
*/
struct PosStats {
    nt16_t nt16;
    int32_t counters[4];
    int32_t sum;
    bool bitshifted;

    PosStats() : nt16(0), counters{0, 0, 0, 0}, sum(0), bitshifted(false) {}

    inline void increment(nt4_t nt4) {
        counters[nt4]++;
        sum++;
    }

    inline void normalize(int nbits) {
        // std::cerr << "     " << nbits << "\n";
        int32_t mask =
            (counters[0] | counters[1] | counters[2] | counters[3]) >> nbits;
        int32_t shifts = 0;
        while (mask > 0) {
            shifts++;
            mask >>= 1;
        }

        for (int i = 0; i < 4; i++) {
            counters[i] >>= shifts;
        }

        sum        = counters[0] + counters[1] + counters[2] + counters[3];
        bitshifted = (shifts > 0) ? true : bitshifted;
    }

    template <typename T>
    inline void pull(T psc) {
        // std::cerr << "     " << __PRETTY_FUNCTION__ << " " << psc << "\n";

        const int C = counter_size<T>();

        // 1. reference base(s) (before correction)
        nt16 = psc & right_full_mask<T>(4);
        psc >>= 4;

        // 2. are the values exact?
        int nones = bitsset_table256[nt16];
        assert(nones != 2);
        if (nones == 1) {
            bitshifted = false;
        } else {
            if (nones == 3) {
                bitshifted = true;
                // if not exact, invert base bits
                nt16 ^= right_full_mask<T>(4);
            }
        }

        // 3. count of individual nucleotides and the sum
        sum = 0;
        for (int32_t i = 3; i >= 0; i--) {
            counters[i] = psc & right_full_mask<T>(C);
            sum += counters[i];
            psc >>= C;
        }
    }

    template <typename T>
    inline void push(T &psc) {
        // std::cerr << "     " << __PRETTY_FUNCTION__ << "\n";

        const int C = counter_size<T>();
        psc         = 0;

        // 1. bitshift counters if necessary
        normalize(C);

        // 2. incorporate counters
        for (int32_t i = 0; i < 4; i++) {
            psc <<= C;
            psc |= counters[i] & right_full_mask<T>(C);
        }

        // 3. incorporate ref base
        psc <<= 4;
        psc |= nt16 & right_full_mask<T>(4);

        // 4. if not exact, invert the base bits
        if (bitshifted) {
            psc ^= right_full_mask<T>(4);
        }
    }
};

}  // namespace ococo
