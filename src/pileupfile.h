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

#include <array>
#include <cstdio>
#include <string>

#include "counters.h"
#include "misc.h"

namespace ococo {

const int PILEUP_MAX_DEPTH = 20000;

struct PileupFile {
    std::string fn;
    FILE *file;

    std::array<char, PILEUP_MAX_DEPTH + 1> bases;
    std::array<char, PILEUP_MAX_DEPTH + 1> qualities;

    PileupFile(std::string fn) : fn(fn) {
        if (fn.size() > 0) {
            ococo::info("Opening the Pileup stream ('%s').\n", fn.c_str());

            if (fn == std::string("-")) {
                file = stdout;
            } else {
                file = fopen(fn.c_str(), "w+");
                if (file == nullptr) {
                    ococo::fatal_error(
                        "Problem with opening the Pileup file '%s'.\n",
                        fn.c_str());
                    // todo:
                    // correctly_initialized = false;
                    return;
                }
            }
        }
    }

    ~PileupFile() {
        if (file != nullptr && fn != "") {
            int error_code = fclose(file);
            if (error_code != 0) {
                ococo::error("Output Pileup file could not be closed.\n");
                // todo:
                // return_code = -1;
            }
        }
    }

    int print_position(const std::string &seq_name, int64_t pos,
                       const pos_stats_uncompr_t &psu) {
        bool overflow = false;

        if (psu.sum >= PILEUP_MAX_DEPTH) {
            ococo::warning(
                "Too high coverage at position %" PRId64
                " in '%s'. "  //
                "Pileup does not support coverage higher than %" PRId32
                "."            //
                " A=%" PRId32  //
                " A=%" PRId32  //
                " G=%" PRId32  //
                " T=%" PRId32  //
                "\n",
                pos, seq_name.c_str(), PILEUP_MAX_DEPTH, psu.counters[0],
                psu.counters[1], psu.counters[2], psu.counters[3]);
        }

        char ref_nt256 = nt16_nt256[psu.nt16];

        if (psu.sum == 0) {
            return 0;
        }

        int32_t j = 0;

        // todo: check that j does not exceeds buffer size
        for (int32_t nt4 = 0; nt4 < 4; nt4++) {
            const char filling_char =
                nt4_nt16[nt4] == psu.nt16 ? '.' : nt4_nt256[nt4];
            for (int32_t i = 0; i < psu.counters[nt4] && j < PILEUP_MAX_DEPTH;
                 i++, j++) {
                bases[j]     = filling_char;
                qualities[j] = '~';
            }
        }

        bases[j]     = '\0';
        qualities[j] = '\0';

        fprintf(file, "%s\t%" PRId64 "\t%c\t%" PRId32 "\t%s\t%s\n",
                seq_name.c_str(), pos + 1, ref_nt256,
                overflow ? PILEUP_MAX_DEPTH : psu.sum, bases.data(),
                qualities.data());

        return 0;
    }
};

}  // namespace ococo