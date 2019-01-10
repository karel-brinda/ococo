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
#include <vector>

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
    std::string fn_in_, fn_out_;
    samFile *file_in_, *file_out_;

    bam1_t *b_;
    bam_hdr_t *header_;

    /* SAM/BAM variables */
    char *rname_;
    uint8_t *seq_;
    uint8_t *qual_;
    uint32_t *cigar_;
    int32_t n_cigar_;
    int32_t seqid_;
    int64_t mapping_pos_;
    int32_t mapq_;
    int32_t flags_;

    BamFiles(std::string fn_in, std::string fn_out)
        : fn_in_(fn_in),
          fn_out_(fn_out),
          file_in_(nullptr),
          file_out_(nullptr),
          b_(bam_init1()),
          header_(nullptr) {
        file_in_ = sam_open(fn_in.c_str(), "r");
        if (file_in_ == nullptr) {
            fatal_error("Problem with opening the input SAM/BAM file ('%s').\n",
                        fn_in_.c_str());
        }

        header_ = sam_hdr_read(file_in_);
        if (header_ == nullptr) {
            fatal_error("SAM/BAM headers are missing or corrupted.\n");
        }

        if (!fn_out_.empty()) {
            file_out_ = sam_open(fn_out_.c_str(), "w");
            if (file_out_ == nullptr) {
                fatal_error(
                    "Problem with opening the output SAM/BAM file ('%s').\n",
                    fn_out.c_str());
            }

            int error_code = sam_hdr_write(file_out_, header_);
            if (error_code != 0) {
                fatal_error("Construction of the SAM header failed (error %d)",
                            error_code);
            }
        }
    }

    ~BamFiles() {
        bam_destroy1(b_);
        bam_hdr_destroy(header_);

        if (file_in_ != nullptr) {
            int error_code = sam_close(file_in_);
            if (error_code != 0) {
                fatal_error("Input SAM file could not be closed.\n");
            }
        }

        if (file_out_ != nullptr) {
            int error_code = sam_close(file_out_);
            if (error_code != 0) {
                fatal_error("Output SAM file could not be closed.\n");
            }
        }
    }

    std::vector<int64_t> get_refseq_lens() {
        int n_seqs = header_->n_targets;
        std::vector<int64_t> lens(n_seqs);
        for (int i = 0; i < n_seqs; i++) {
            lens[i] = header_->target_len[i];
        }
        return lens;
    }

    std::vector<std::string> get_refseq_names() {
        int n_seqs = header_->n_targets;
        std::vector<std::string> names(n_seqs);
        for (int i = 0; i < n_seqs; i++) {
            names[i] = std::string(header_->target_name[i]);
        }
        return names;
    }

    int read_alignment() {
        int return_value = sam_read1(file_in_, header_, b_);

        rname_       = bam_get_qname(b_);
        seq_         = bam_get_seq(b_);
        qual_        = bam_get_qual(b_);
        cigar_       = bam_get_cigar(b_);
        n_cigar_     = b_->core.n_cigar;
        seqid_       = b_->core.tid;
        mapping_pos_ = b_->core.pos;
        mapq_        = b_->core.qual;
        flags_       = b_->core.flag;

        return return_value;
    }

    void print_alignment() {
        if (file_out_ != nullptr) {
            int error_code = sam_write1(file_out_, header_, b_);
            if (error_code != 0) {
                fatal_error("Writing SAM failed (error %d)", error_code);
            }
        }
    }
};

}  // namespace ococo
