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
    std::string fn_;
    FILE *file_;

    std::array<char, PILEUP_MAX_DEPTH + 1> bases_;
    std::array<char, PILEUP_MAX_DEPTH + 1> qualities_;

    PileupFile(std::string fn) : fn_(fn), file_(nullptr) {
        if (fn.size()) {
            info("Opening the Pileup stream ('%s').\n", fn.c_str());

            if (fn == std::string("-")) {
                file_ = stdout;
            } else {
                file_ = fopen(fn.c_str(), "w+");
                if (file_ == nullptr) {
                    fatal_error("Problem with opening the Pileup file '%s'.\n",
                                fn.c_str());
                }
            }
        }
    }

    ~PileupFile() {
        if (file_ != nullptr && fn_ != "-") {
            int error_code = fclose(file_);
            if (error_code != 0) {
                fatal_error("Output Pileup file could not be closed.\n");
            }
        }
    }

    void print_position(const std::string &seq_name, int64_t pos,
                        const PosStats &ps) {
        if (file_ == nullptr) {
            return;
        }

        if (ps.sum_ == 0) {
            return;
        }

        bool overflow = false;

        if (ps.sum_ >= PILEUP_MAX_DEPTH) {
            warning("Too high coverage at position %" PRId64
                    " in '%s'. "  //
                    "Pileup does not support coverage higher than %" PRId32
                    "."            //
                    " A=%" PRId32  //
                    " A=%" PRId32  //
                    " G=%" PRId32  //
                    " T=%" PRId32  //
                    "\n",
                    pos, seq_name.c_str(), PILEUP_MAX_DEPTH, ps.counters_[0],
                    ps.counters_[1], ps.counters_[2], ps.counters_[3]);
            overflow = true;
        }

        char ref_nt256 = nt16_nt256[ps.nt16_];

        int32_t j = 0;

        // todo: check that j does not exceeds buffer size
        for (int32_t nt4 = 0; nt4 < 4; nt4++) {
            const char filling_char =
                nt4_nt16[nt4] == ps.nt16_ ? '.' : nt4_nt256[nt4];
            for (int32_t i = 0; i < ps.counters_[nt4] && j < PILEUP_MAX_DEPTH;
                 i++, j++) {
                bases_[j]     = filling_char;
                qualities_[j] = '~';
            }
        }

        bases_[j]     = '\0';
        qualities_[j] = '\0';

        fprintf(file_, "%s\t%" PRId64 "\t%c\t%" PRId32 "\t%s\t%s\n",
                seq_name.c_str(), pos + 1, ref_nt256,
                overflow ? PILEUP_MAX_DEPTH : ps.sum_, bases_.data(),
                qualities_.data());
    }
};

}  // namespace ococo
