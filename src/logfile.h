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

#include <cstdio>
#include <string>

#include "types.h"

namespace ococo {

struct LogFile {
    std::string fn;
    FILE *file;

    LogFile(std::string fn) : fn(fn) {
        if (!fn.empty()) {
            info("Opening the log file ('%s').\n", fn.c_str());

            file = fopen(fn.c_str(), "w+");
            if (file == nullptr) {
                fatal_error("Problem with opening the log file '%s'.\n",
                            fn.c_str());
            }
        }
    }

    ~LogFile() {
        if (file != nullptr) {
            int error_code = fclose(file);
            if (error_code != 0) {
                fatal_error("Output log file '%s' could not be closed.\n",
                            fn.c_str());
            }
        }
    }

    void print(int64_t i_reads, const char *rname, int64_t nupd) {
        if (file != nullptr) {
            fprintf(file, "%" PRId64 "\t%s\t%" PRId64 "\n", i_reads, rname,
                    nupd);
        }
    }
};

}  // namespace ococo
