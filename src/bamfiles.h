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

#include <htslib/faidx.h>
#include <htslib/khash.h>
#include <htslib/kseq.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>

#include "counters.h"
#include "io.h"
#include "misc.h"

namespace ococo {

struct BamFiles {
    std::string fn_in, fn_out;
    samFile *file_in, *file_out;

    bam1_t *b;
    bam_hdr_t *header;

    /* SAM/BAM variables */
    char *rname;
    uint8_t *seq;
    uint8_t *qual;
    uint32_t *cigar;
    int32_t n_cigar;
    int32_t seqid;
    int64_t mapping_pos;
    int32_t mapq;
    int32_t flags;

    BamFiles(std::string fn_in, std::string fn_out)
        : fn_in(fn_in), fn_out(fn_out), b(bam_init1()), header(nullptr) {
        file_in = sam_open(fn_in.c_str(), "r");

        if (file_in == nullptr) {
            fatal_error("Problem with opening the input SAM/BAM file ('%s').\n",
                        fn_in.c_str());
            // todo:
            // correctly_initialized = false;
            return;
        }

        if ((header = sam_hdr_read(file_in)) == nullptr) {
            fatal_error("SAM/BAM headers are missing or corrupted.\n");
            // todo:
            // correctly_initialized = false;
            return;
        }

        if (fn_out.size() > 0) {
            file_out = sam_open(fn_out.c_str(), "w");
            if (file_out == nullptr) {
                fatal_error(
                    "Problem with opening the output SAM/BAM file ('%s').\n",
                    fn_out.c_str());
                // todo
                // correctly_initialized = false;
                return;
            }

            int error_code = sam_hdr_write(file_out, header);
            if (error_code != 0) {
                error("Construction of the SAM header failed (error %d)",
                      error_code);
                // todo
                // return_code = EXIT_FAILURE;
                return;
            }
        }
    }

    ~BamFiles() {
        bam_destroy1(b);
        bam_hdr_destroy(header);

        if (file_in != nullptr) {
            int error_code = sam_close(file_in);
            if (error_code != 0) {
                error("Input SAM file could not be closed.\n");
                // todo:
                // return_code = -1;
            }
        }

        if (file_out != nullptr) {
            int error_code = sam_close(file_out);
            if (error_code != 0) {
                error("Output SAM file could not be closed.\n");
                // todo:
                // return_code = -1;
            }
        }
    }

    int read_alignment() {
        int return_value = sam_read1(file_in, header, b);

        rname       = bam_get_qname(b);
        seq         = bam_get_seq(b);
        qual        = bam_get_qual(b);
        cigar       = bam_get_cigar(b);
        n_cigar     = b->core.n_cigar;
        seqid       = b->core.tid;
        mapping_pos = b->core.pos;
        mapq        = b->core.qual;
        flags       = b->core.flag;
        return return_value;
    }

    int print_alignment() {
        if (file_out != nullptr) {
            return sam_write1(file_out, header, b);
        }
        return 0;
    }
};

}  // namespace ococo