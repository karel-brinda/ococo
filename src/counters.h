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

#include <iomanip>
#include <iostream>
#include <sstream>

#include "misc.h"
#include "types.h"

using namespace std;

namespace ococo {

template <typename T>
constexpr int counter_size() {
    return (8 * sizeof(T) - 4) / 8;
}

/*! @struct
    @abstract          Structure for uncompressed pileup statistics for
                       1 position.
    @field nt16        Consensus base.
                       0x00 and 0x0f = N
                       0x01, 0x02, 0x04, 0x08 = bases
                       other values: error
    @field counters    Nucleotide counters.
    @field sum         Sum of the nucleotide counters.
    @field bitshifted  Already bit-shifted? (i.e., in-exact).
*/
struct pos_stats_uncompr_t {
    nt16_t nt16;
    int32_t counters[4];
    int32_t sum;
    bool bitshifted;

    pos_stats_uncompr_t()
        : nt16(0), counters{0, 0, 0, 0}, sum(0), bitshifted(false) {}

    void increment(nt4_t nt4) {
        counters[nt4]++;
        sum++;
    }

    template <typename T>
    void decompress(T psc) {
        // std::cerr << "     " << __PRETTY_FUNCTION__ << " " << psc <<
        // std::endl;

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
    void compress(T psc) {
        // std::cerr << "     " << __PRETTY_FUNCTION__ << " " << psc <<
        // std::endl;

        // todo: bitshift before compression

        const int C = counter_size<T>();
        psc         = 0;

        // remove if you want to support ambiguous nucleotides
        assert(bitsset_table256[nt16] != 2);

        // 1. incorporate counters
        for (int32_t i = 0; i < 4; i++) {
            psc <<= C;
            psc |= counters[i] & right_full_mask<T>(C);
        }

        // 2. incorporate ref base
        psc <<= 4;
        psc |= nt16 & right_full_mask<T>(4);

        // 3. if not exact, invert the base bits
        if (bitshifted) {
            psc ^= right_full_mask<T>(4);
        }
    }

    void bitshift(int n) {
        if (n > 0) {
            counters[0] >>= n;
            counters[1] >>= n;
            counters[2] >>= n;
            counters[3] >>= n;

            sum        = counters[0] + counters[1] + counters[2] + counters[3];
            bitshifted = true;
        }
    }
};

/*
string _pos_stats_uncompr(pos_stats_uncompr_t psu) {
    stringstream ss;
    ss << showbase << internal << setfill('0');
    ss << "[" << nt16_nt256[psu.nt16] << "]"
       << "(" << psu.counters[0] << "," << psu.counters[1] << ","
       << psu.counters[2] << "," << psu.counters[3] << ")";

    return ss.str();
}

template <typename T>
string _pos_stats_compr(T psc) {
    stringstream ss;
    ss << showbase << internal << setfill('0');
    ss << psc;

    return ss.str();
}
*/
template <typename T>
void _print_pos_stats(T psc) {
    // pos_stats_uncompr_t psu;
    // psu.decompress(psc);
    // cerr << _pos_stats_compr(psc) << " " << _pos_stats_uncompr(psu) << endl;
}

}  // namespace ococo
