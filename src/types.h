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

#include <unistd.h>

#ifdef _LIBCPP_VERSION
#include <cinttypes>
#else
#include <tr1/cinttypes>
#endif

namespace ococo {

const int fasta_line_l  = 50;
const int stats_delim_l = 10;

typedef uint8_t nt4_t;
typedef uint8_t nt16_t;
typedef uint8_t nt256_t;

/******************
 *                *
 *   Structures   *
 *                *
 ******************/

/*************
 *** Ococo ***
 *************/

/*! @struct
    @abstract           Structure for metadata for 1 sequence.
    @field seq_active   Active.
    @field seq_len      Sequence length.
    @field seq_name     Name of the sequence.
    @field seq_comment  Comment of the sequence.
*/
struct single_seq_serial_t {
    bool seq_active;
    int64_t seq_len;
    char seq_name[1000];
    char seq_comment[1000];
};

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
};

/**************************
 *** Translation tables ***
 **************************/

// clang-format off

static const uint8_t nt256_nt4[] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

static const uint8_t nt16_nt4[] = {4, 0, 1, 4, 2, 4, 4, 4,
                                   3, 4, 4, 4, 4, 4, 4, 4};

static const uint8_t nt256_nt16[] = {
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,      15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,      15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,      15, 15,
    1,  2,  4,  8,  15, 15, 15, 15, 15, 15, 15, 15, 15, 0 /*=*/, 15, 15,
    15, 1,  14, 2,  13, 15, 15, 4,  11, 15, 15, 12, 15, 3,       15, 15,
    15, 15, 5,  6,  8,  15, 7,  9,  15, 10, 15, 15, 15, 15,      15, 15,
    15, 1,  14, 2,  13, 15, 15, 4,  11, 15, 15, 12, 15, 3,       15, 15,
    15, 15, 5,  6,  8,  15, 7,  9,  15, 10, 15, 15, 15, 15,      15, 15,

    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,      15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,      15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,      15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,      15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,      15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,      15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,      15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,      15, 15};

// clang-format on

static const uint8_t nt16_nt256[] =
    "NACMGRSVTWYHKDBN";  // modified, the first character is usually '='

static const uint8_t nt4_nt256[] = "ACGTN";

static const uint8_t nt4_nt16[] = {1, 2, 4, 8, 15};
}  // namespace ococo

/*
 * Number of bits in an integer.
 */
static const uint8_t bitsset_table256[256] = {
#define B2(n) n, n + 1, n + 1, n + 2
#define B4(n) B2(n), B2(n + 1), B2(n + 1), B2(n + 2)
#define B6(n) B4(n), B4(n + 1), B4(n + 1), B4(n + 2)
    B6(0), B6(1), B6(1), B6(2)};
