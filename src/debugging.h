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

template <typename T>
void _print_pos_stats(T psc) {
    pos_stats_uncompr_t psu;
    psu.decompress(psc);
    cerr << _pos_stats_compr(psc) << " " << _pos_stats_uncompr(psu) << endl;
}

}  // namespace ococo
