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

#include "counters.h"
#include "misc.h"
#include "types.h"

namespace ococo {

/*
   Auxiliary debugging functions for printing position statistics before and
   after compression.
*/

std::string __pos_stats_uncompr(pos_stats_uncompr_t psu) {
    std::stringstream ss;
    ss << std::showbase << std::internal << std::setfill('0');
    ss << "[" << nt16_nt256[psu.nt16] << "]"
       << "(" << psu.counters[0] << "," << psu.counters[1] << ","
       << psu.counters[2] << "," << psu.counters[3] << ")";

    return ss.str();
}

template <typename T>
std::string __pos_stats_compr(T psc) {
    std::stringstream ss;
    ss << std::showbase << std::internal << std::setfill('0');
    ss << psc;

    return ss.str();
}

template <typename T>
void _print_pos_stats(T psc) {
    pos_stats_uncompr_t psu;
    psu.decompress(psc);
    std::cerr << __pos_stats_compr(psc) << " " << __pos_stats_uncompr(psu)
              << "\n";
}

template <typename T>
void _print_pos_stats(const pos_stats_uncompr_t &psu) {
    T psc;
    psu.compress<T>(psc);
    std::cerr << __pos_stats_compr(psc) << " " << __pos_stats_uncompr(psu)
              << "\n";
}

}  // namespace ococo

// template <typename T>
// void stats_t<T>::debug_print_counters() const {
//     for (int seqid = 0; seqid < n_seqs; seqid++) {
//         fprintf(stderr, "%s\n", seq_name[seqid]);
//         for (int64_t pos = 0; pos < seq_len[seqid]; pos++) {
//             fprintf(stderr, "%8" PRId64 " %04x \n", pos,
//             seq_stats[seqid][pos]);
//         }
//     }
// }

// template <typename T>
// std::string stats_t<T>::debug_str_counters(int32_t seqid, int64_t pos) const
// {
//     pos_stats_uncompr_t psu;
//     psu.decompress(seq_stats[seqid][pos]);
//     std::stringstream ss;
//     ss << "[" << nt16_nt256[psu.nt16] << "]"
//        << "(" << psu.counters[0] << "," << psu.counters[1] << ","
//        << psu.counters[2] << "," << psu.counters[3] << ")";
//     return ss.str();
// }
