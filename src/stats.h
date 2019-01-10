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

#include <cassert>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#include <zlib.h>

#include <htslib/faidx.h>
#include <htslib/khash.h>
#include <htslib/kseq.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>

#include "consensus.h"
#include "counters.h"
#include "io.h"
#include "pileupfile.h"
#include "types.h"
#include "vcffile.h"

/***********************
 *** Main statistics ***
 ***********************/

namespace ococo {

KSEQ_INIT(gzFile, gzread);

template <typename T>
struct Stats {
    int32_t n_seqs;
    std::vector<bool> seq_active;
    std::vector<int64_t> seq_len;
    std::vector<std::string> seq_name;
    std::vector<std::string> seq_comment;
    std::vector<std::vector<T>> seq_stats;

    Params params;

    Stats(Params params, const std::vector<std::string> &seq_name,
          const std::vector<int64_t> &seq_len)
        : n_seqs(seq_len.size()),
          seq_active(std::vector<bool>(n_seqs, true)),
          seq_len(seq_len),
          seq_name(seq_name),
          seq_comment(std::vector<std::string>(n_seqs)),
          seq_stats(std::vector<std::vector<T>>(n_seqs)),
          params(params) {
        for (int seqid = 0; seqid < n_seqs; seqid++) {
            seq_stats[seqid].resize(seq_len[seqid]);
        }
    }

    /*******
     * I/O *
     *******/

    void import_stats(const std::string &stats_fn) {
        int error_code = 0;

        FILE *fo = fopen(stats_fn.c_str(), "r");
        if (fo == nullptr) {
            fatal_error("File with statistics could not be opened ('%s').\n",
                        stats_fn.c_str());
        }

        /* number of seqs */
        int32_t n_seqs_loaded;
        fread(&n_seqs_loaded, sizeof(int32_t), 1, fo);

        if (n_seqs_loaded != n_seqs) {
            fatal_error(
                "Numbers of sequences in stats and SAM/BAM do not correspond "
                "%" PRId32 "!=%" PRId32 ").\n",
                n_seqs_loaded, n_seqs);
        }

        for (int seqid = 0; seqid < n_seqs; seqid++) {
            /* sequence */

            single_seq_serial_t seq_ser;
            fread(&seq_ser, sizeof(single_seq_serial_t), 1, fo);

            if (seq_ser.seq_active != seq_active[seqid]) {
                fatal_error(
                    "Active sequences in stats and SAM/BAM do not correspond "
                    "(seqid %" PRId32 ").\n",
                    seqid);
            }

            if (seq_ser.seq_len != seq_len[seqid]) {
                fatal_error(
                    "Sequence lengths in stats and SAM/BAM do not correspond "
                    "(seqid %" PRId32 ", %" PRId64 "!=%" PRId64 ").\n",
                    seqid, seq_ser.seq_len, seq_len[seqid]);
            }

            if (seq_name[seqid].compare(seq_ser.seq_name) != 0) {
                fatal_error(
                    "Sequence names in stats and SAM/BAM do not correspond "
                    "(seqid %" PRId32 ", '%s'!='%s').\n",
                    seqid, seq_ser.seq_name, seq_name[seqid].c_str());
            }

            fread(&(seq_stats[seqid][0]), sizeof(T), seq_len[seqid], fo);
        }

        error_code = fclose(fo);
        if (error_code != 0) {
            fatal_error("File with statistics could not be closed ('%s').\n",
                        stats_fn.c_str());
        }
    }

    void export_stats(const std::string &stats_fn) const {
        FILE *fo = fopen(stats_fn.c_str(), "w+");
        if (fo == nullptr) {
            fatal_error("File with statistics could not be opened ('%s').\n",
                        stats_fn.c_str());
        }

        /* number of seqs */
        fwrite(&n_seqs, sizeof(int32_t), 1, fo);

        for (int seqid = 0; seqid < n_seqs; seqid++) {
            /* sequence */
            single_seq_serial_t seq_ser = {0};
            seq_ser.seq_active          = seq_active[seqid];
            seq_ser.seq_len             = seq_len[seqid];
            strncpy(seq_ser.seq_name, seq_name[seqid].c_str(), 999);
            strncpy(seq_ser.seq_name, seq_name[seqid].c_str(), 999);
            uint64_t written = 0;
            written += fwrite(&seq_ser, sizeof(single_seq_serial_t), 1, fo);
            written +=
                fwrite(&(seq_stats[seqid][0]), sizeof(T), seq_len[seqid], fo);
            if (written != 1 + static_cast<uint64_t>(seq_len[seqid])) {
                fatal_error(
                    "Problem with writting to the file with statistics "
                    "('%s').\n",
                    stats_fn.c_str());
            }
        }

        int error_code = fclose(fo);
        if (error_code != 0) {
            fatal_error("File with statistics could not be closed ('%s').\n",
                        stats_fn.c_str());
        }
    }

    // Call consensus probabilistically.
    void call_consensus(const VcfFile &vcf_file, PileupFile &pileup_file) {
        PosStats psu;

        for (int32_t seqid = 0; seqid < n_seqs; seqid++) {
            for (int64_t pos = 0; pos < seq_len[seqid]; pos++) {
                psu.decompress(seq_stats[seqid][pos]);
                call_consensus_position(vcf_file, pileup_file, seqid, pos, psu);
                psu.compress(seq_stats[seqid][pos]);
            }
        }
    }

    void call_consensus_position(const VcfFile &vcf_file,
                                 PileupFile &pileup_file, int32_t seqid,
                                 int64_t pos, PosStats &psu) {
        const char old_base_nt256 = nt16_nt256[psu.nt16];
        const char new_base_nt256 = cons_call_maj(psu, params.min_coverage_upd,
                                                  params.majority_threshold);

        if (old_base_nt256 != new_base_nt256) {
            params.n_upd += 1;
            psu.nt16 = nt256_nt16[static_cast<int16_t>(new_base_nt256)];
        }

        if (old_base_nt256 != new_base_nt256 || params.verbose) {
            vcf_file.print_substitution(seq_name[seqid], pos, old_base_nt256,
                                        new_base_nt256, psu);
        }

        pileup_file.print_position(seq_name[seqid], pos, psu);
    }

    // Load header and data from a FASTA file and initialize statistics.
    void load_fasta(const std::string &fasta_fn) {
        gzFile fp;
        kseq_t *seq;
        int l;
        fp  = gzopen(fasta_fn.c_str(), "r");
        seq = kseq_init(fp);

        if (fp == nullptr) {
            fatal_error("The FASTA file '%s' could not be opened.\n",
                        fasta_fn.c_str());
        }

        for (int seqid = 0; (l = kseq_read(seq)) >= 0; seqid++) {
            if (seq_name[seqid].compare(seq->name.s) != 0) {
                fatal_error(
                    "Sequence names in BAM/SAM and in FASTA do not correspond "
                    "('%s'!='%s').\n",
                    seq_name[seqid].c_str(), seq->name.s);
            }

            if (seq_len[seqid] != static_cast<int64_t>(seq->seq.l)) {
                fatal_error(
                    "Sequence lengths in BAM/SAM and in FASTA do not "
                    "correspond "
                    "(%" PRId64 "!=%" PRId64 ").\n",
                    static_cast<int64_t>(seq->seq.l),
                    static_cast<int64_t>(seq_len[seqid]));
            }

            if (seq->comment.l && seq_comment[seqid].empty()) {
                seq_comment[seqid] = std::string(seq->comment.s);
            }

            for (int64_t pos = 0; pos < static_cast<int64_t>(seq->seq.l);
                 pos++) {
                assert(seq_stats[seqid][pos] == 0);
                PosStats psu;
                psu.nt16 = nt256_nt16[static_cast<int32_t>(seq->seq.s[pos])];
                psu.compress(seq_stats[seqid][pos]);
            }
        }
        kseq_destroy(seq);
        gzclose(fp);
    }

    void save_fasta(const std::string &fasta_fn) const {
        FILE *fasta_file = nullptr;
        fasta_file       = fopen(fasta_fn.c_str(), "w+");
        if (fasta_file == nullptr) {
            fatal_error("Problem with opening the FASTA file: '%s'.\n",
                        fasta_fn.c_str());
        }

        char fasta_buffer[fasta_line_l];
        PosStats psu;

        for (int s = 0; s < n_seqs; s++) {
            // printf("%s\n",seq_name[s]);
            if (!seq_comment[s].empty()) {
                fprintf(fasta_file, ">%s %s\n", seq_name[s].c_str(),
                        seq_comment[s].c_str());
            } else {
                fprintf(fasta_file, ">%s\n", seq_name[s].c_str());
            }

            for (int64_t i = 0, j = 0; i < seq_len[s]; i++, j++) {
                psu.decompress(seq_stats[s][i]);
                fasta_buffer[j] = nt16_nt256[psu.nt16];

                if (j == fasta_line_l - 1 || i == seq_len[s] - 1) {
                    fwrite(fasta_buffer, 1, j + 1, fasta_file);
                    fwrite("\n", 1, 1, fasta_file);
                    j = -1;
                }
            }
        }

        int error_code = fclose(fasta_file);
        if (error_code != 0) {
            fatal_error("File with consensus could not be closed ('%s').\n",
                        fasta_fn.c_str());
        }
    }

    /***********************
     * Debuging & checking *
     ***********************/

    // Check if a BAM header corresponds to the stats.
    bool check_headers_bam_hdr(const bam_hdr_t *h) const {
        // todo: is it used anywhere?
        assert(h != nullptr);
        for (int32_t seqid = 0; seqid < n_seqs; seqid++) {
            if (seq_len[seqid] != static_cast<int64_t>(h->target_len[seqid])) {
                return false;
            }
            if (seq_name[seqid].compare(h->target_name[seqid]) != 0) {
                return false;
            }
        }

        return true;
    }
};

}  // namespace ococo
