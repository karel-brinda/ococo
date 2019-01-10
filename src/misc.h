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

#include <sys/resource.h>
#include <sys/time.h>

namespace ococo {

double realtime() {
    struct timeval tp;
    // struct timezone tzp;
    // gettimeofday(&tp, &tzp);
    gettimeofday(&tp, nullptr);
    return tp.tv_sec + tp.tv_usec * 1e-6;
}

double cputime() {
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);

    // todo: check also memory
    // std::cerr << r.ru_maxrss << std::endl;
    return r.ru_utime.tv_sec + r.ru_stime.tv_sec +
           1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

/*
 * Get a right full mask (right n bits set to 1)
 *
 * T - type
 * size - number of 1's
 */
template <typename T>
constexpr T right_full_mask(int size) {
    return (size == 0) ? 0
                       : (((static_cast<T>(0x1) << (size - 1)) - 1) << 1) | 1;
}
}  // namespace ococo
