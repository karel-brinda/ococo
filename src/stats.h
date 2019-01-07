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
#include "vcffile.h"

/***********************
 *** Main statistics ***
 ***********************/

namespace ococo {

KSEQ_INIT(gzFile, gzread);

template <typename T>
struct stats_t {
    int32_t n_seqs;
    bool *seq_active;
    int64_t *seq_len;
    std::string *seq_name;
    std::string *seq_comment;
    T **seq_stats;

    params_t *params;

    // stats_t();
    stats_t(params_t *params, bam_hdr_t &h);
    ~stats_t();

    /*******
     * I/O *
     *******/

    int import_stats(const std::string &stats_fn);
    int export_stats(const std::string &stats_fn) const;

    // Call consensus probabilistically.
    int call_consensus(const VcfFile &vcf_file, PileupFile &pileup_file);
    int call_consensus_position(const VcfFile &vcf_file,
                                PileupFile &pileup_file, int32_t seqid,
                                int64_t pos, pos_stats_uncompr_t &psu);

    // Load header and data from a FASTA file and initialize statistics.
    int load_fasta(const std::string &fasta_fn);
    int save_fasta(const std::string &fasta_fn) const;

    /***********************
     * Debuging & checking *
     ***********************/

    // Check if a BAM header corresponds to the stats.
    bool check_headers_bam_hdr(const bam_hdr_t &h) const;

    void debug_print_counters() const;

    std::string debug_str_counters(int32_t seqid, int64_t pos) const;
};

template <typename T>
stats_t<T>::stats_t(params_t *params, bam_hdr_t &h)
    : n_seqs(h.n_targets),
      seq_active(new (std::nothrow) bool[n_seqs]()),
      seq_len(new (std::nothrow) int64_t[n_seqs]()),
      seq_name(new (std::nothrow) std::string[n_seqs]()),
      seq_comment(new (std::nothrow) std::string[n_seqs]()),
      seq_stats(new (std::nothrow) T *[n_seqs]()),
      params(params) {
    for (int seqid = 0; seqid < n_seqs; seqid++) {
        seq_len[seqid]    = h.target_len[seqid];
        seq_active[seqid] = true;
        seq_name[seqid]   = std::string(h.target_name[seqid]);

        seq_stats[seqid] = new (std::nothrow) T[seq_len[seqid]]();
    }
}

template <typename T>
stats_t<T>::~stats_t() {
    if (seq_stats != nullptr) {
        for (int32_t seqid = 0; seqid < n_seqs; seqid++) {
            delete[] seq_stats[seqid];
        }
    }
    delete[] seq_active;
    delete[] seq_len;
    delete[] seq_name;
    delete[] seq_comment;
    delete[] seq_stats;
}

template <typename T>
int stats_t<T>::load_fasta(const std::string &fasta_fn) {
    gzFile fp;
    kseq_t *seq;
    int l;
    fp  = gzopen(fasta_fn.c_str(), "r");
    seq = kseq_init(fp);

    if (fp == nullptr) {
        error("File '%s' could not be opened.\n", fasta_fn.c_str());
        return -1;
    }

    for (int seqid = 0; (l = kseq_read(seq)) >= 0; seqid++) {
        if (seq_name[seqid].compare(seq->name.s) != 0) {
            error(
                "Sequence names in BAM/SAM and in FASTA do not correspond "
                "('%s'!='%s').\n",
                seq_name[seqid].c_str(), seq->name.s);
            return -1;
        }

        if (seq_len[seqid] != static_cast<int64_t>(seq->seq.l)) {
            error(
                "Sequence lengths in BAM/SAM and in FASTA do not correspond "
                "(%" PRId64 "!=%" PRId64 ").\n",
                static_cast<int64_t>(seq->seq.l),
                static_cast<int64_t>(seq_len[seqid]));
            return -1;
        }

        if (seq->comment.l && seq_comment[seqid].empty()) {
            seq_comment[seqid] = std::string(seq->comment.s);
        }

        for (int64_t pos = 0; pos < static_cast<int64_t>(seq->seq.l); pos++) {
            assert(seq_stats[seqid][pos] == 0);
            pos_stats_uncompr_t psu;
            psu.nt16 = nt256_nt16[static_cast<int32_t>(seq->seq.s[pos])];
            psu.compress(seq_stats[seqid][pos]);
        }
    }
    kseq_destroy(seq);
    gzclose(fp);
    return 0;
}

template <typename T>
int stats_t<T>::save_fasta(const std::string &fasta_fn) const {
    FILE *fasta_file = nullptr;
    fasta_file       = fopen(fasta_fn.c_str(), "w+");
    if (fasta_file == nullptr) {
        error("Problem with opening the FASTA file: '%s'.\n", fasta_fn.c_str());
        return -1;
    }

    char fasta_buffer[fasta_line_l];
    pos_stats_uncompr_t psu;

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
        error("File with consensus could not be closed ('%s').\n",
              fasta_fn.c_str());
        return -1;
    }

    return 0;
}

template <typename T>
bool stats_t<T>::check_headers_bam_hdr(const bam_hdr_t &h) const {
    for (int32_t seqid = 0; seqid < n_seqs; seqid++) {
        if (seq_len[seqid] != static_cast<int64_t>(h.target_len[seqid])) {
            return false;
        }
        if (seq_name[seqid].compare(h.target_name[seqid]) != 0) {
            return false;
        }
    }

    return true;
}

template <typename T>
int stats_t<T>::import_stats(const std::string &stats_fn) {
    int error_code = 0;

    FILE *fo = fopen(stats_fn.c_str(), "r");
    if (fo == nullptr) {
        error("File with statistics could not be opened ('%s').\n",
              stats_fn.c_str());
        return -1;
    }

    /* number of seqs */
    int32_t n_seqs_loaded;
    fread(&n_seqs_loaded, sizeof(int32_t), 1, fo);

    if (n_seqs_loaded != n_seqs) {
        error(
            "Numbers of sequences in stats and SAM/BAM do not correspond "
            "%" PRId32 "!=%" PRId32 ").\n",
            n_seqs_loaded, n_seqs);
        return -1;
    }

    for (int seqid = 0; seqid < n_seqs; seqid++) {
        /* sequence */

        single_seq_serial_t seq_ser;
        fread(&seq_ser, sizeof(single_seq_serial_t), 1, fo);

        if (seq_ser.seq_active != seq_active[seqid]) {
            error(
                "Active sequences in stats and SAM/BAM do not correspond "
                "(seqid %" PRId32 ").\n",
                seqid);
            return -1;
        }

        if (seq_ser.seq_len != seq_len[seqid]) {
            error(
                "Sequence lengths in stats and SAM/BAM do not correspond "
                "(seqid %" PRId32 ", %" PRId64 "!=%" PRId64 ").\n",
                seqid, seq_ser.seq_len, seq_len[seqid]);
            return -1;
        }

        if (seq_name[seqid].compare(seq_ser.seq_name) != 0) {
            error(
                "Sequence names in stats and SAM/BAM do not correspond "
                "(seqid %" PRId32 ", '%s'!='%s').\n",
                seqid, seq_ser.seq_name, seq_name[seqid].c_str());
            return -1;
        }

        fread(seq_stats[seqid], sizeof(T), seq_len[seqid], fo);
    }

    error_code = fclose(fo);
    if (error_code != 0) {
        error("File with statistics could not be closed ('%s').\n",
              stats_fn.c_str());
        return -1;
    }

    return 0;
}

template <typename T>
int stats_t<T>::export_stats(const std::string &stats_fn) const {
    int error_code = 0;

    FILE *fo = fopen(stats_fn.c_str(), "w+");
    if (fo == nullptr) {
        error("File with statistics could not be opened ('%s').\n",
              stats_fn.c_str());
        return -1;
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
        written += fwrite(seq_stats[seqid], sizeof(T), seq_len[seqid], fo);
        if (written != 1 + static_cast<uint64_t>(seq_len[seqid])) {
            error("Problem with writting to the file with statistics ('%s').\n",
                  stats_fn.c_str());
            return -1;
        }
    }

    error_code = fclose(fo);
    if (error_code != 0) {
        error("File with statistics could not be closed ('%s').\n",
              stats_fn.c_str());
        return -1;
    }
    return 0;
}

template <typename T>
int stats_t<T>::call_consensus(const VcfFile &vcf_file,
                               PileupFile &pileup_file) {
    pos_stats_uncompr_t psu;

    for (int32_t seqid = 0; seqid < n_seqs; seqid++) {
        for (int64_t pos = 0; pos < seq_len[seqid]; pos++) {
            psu.decompress(seq_stats[seqid][pos]);
            call_consensus_position(vcf_file, pileup_file, seqid, pos, psu);
            psu.compress(seq_stats[seqid][pos]);
        }
    }

    return 0;
}

template <typename T>
int stats_t<T>::call_consensus_position(const VcfFile &vcf_file,
                                        PileupFile &pileup_file, int32_t seqid,
                                        int64_t pos, pos_stats_uncompr_t &psu) {
    const char old_base_nt256 = nt16_nt256[psu.nt16];
    const char new_base_nt256 = cons_call_maj(psu, params->min_coverage_upd,
                                              params->majority_threshold);

    if (old_base_nt256 != new_base_nt256) {
        params->n_upd += 1;
        psu.nt16 = nt256_nt16[static_cast<int16_t>(new_base_nt256)];
    }

    if (old_base_nt256 != new_base_nt256 || params->verbose) {
        vcf_file.print_substitution(seq_name[seqid], pos, old_base_nt256,
                                    new_base_nt256, psu);
    }

    pileup_file.print_position(seq_name[seqid], pos, psu);

    return 0;
}

template <typename T>
std::string stats_t<T>::debug_str_counters(int32_t seqid, int64_t pos) const {
    pos_stats_uncompr_t psu;
    psu.decompress(seq_stats[seqid][pos]);
    std::stringstream ss;
    ss << "[" << nt16_nt256[psu.nt16] << "]"
       << "(" << psu.counters[0] << "," << psu.counters[1] << ","
       << psu.counters[2] << "," << psu.counters[3] << ")";
    return ss.str();
}

template <typename T>
void stats_t<T>::debug_print_counters() const {
    for (int seqid = 0; seqid < n_seqs; seqid++) {
        fprintf(stderr, "%s\n", seq_name[seqid]);
        for (int64_t pos = 0; pos < seq_len[seqid]; pos++) {
            fprintf(stderr, "%8" PRId64 " %04x \n", pos, seq_stats[seqid][pos]);
        }
    }
}

}  // namespace ococo
