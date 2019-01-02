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

#include <htslib/faidx.h>
#include <htslib/khash.h>
#include <htslib/kseq.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>

#include <zlib.h>
#include <cassert>
#include <cmath>
#include <ctime>
#include <sstream>

/***********************
 *** Main statistics ***
 ***********************/

namespace ococo {

KSEQ_INIT(gzFile, gzread);

template <typename T, int counter_size, int refbase_size>
struct stats_t {
    static_assert(8 * sizeof(T) >= 4 * counter_size + refbase_size,
                  "Too large counter size (does not fit into the main type).");

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
    int call_consensus(FILE *vcf_file, FILE *out_pileup_file);
    int call_consensus_position(FILE *vcf_file, FILE *out_pileup_file,
                                int32_t seqid, int64_t pos);
    int call_consensus_position_uncompressed(FILE *vcf_file,
                                             FILE *out_pileup_file,
                                             int32_t seqid, int64_t pos,
                                             pos_stats_uncompr_t &psu);

    // Load header and data from a FASTA file and initialize statistics.
    int load_fasta(const std::string &fasta_fn);
    int save_fasta(const std::string &fasta_fn) const;

    int print_vcf_header(FILE *vcf_file, std::string cmd,
                         std::string fasta) const;
    int print_vcf_substitution(FILE *vcf_file, int32_t seqid, int64_t pos,
                               char old_base, char new_base,
                               const pos_stats_uncompr_t &psu) const;

    int print_pileup_line(FILE *out_pileup_file, int32_t seqid, int64_t pos,
                          const pos_stats_uncompr_t &psu) const;

    /*************************
     * Statistics & counters *
     *************************/

    inline int set_nucl_nt256(int32_t seqid, int64_t pos, const char &nt256);
    inline int get_nucl_nt256(int32_t seqid, int64_t pos, char &nt256) const;

    static T compress_position_stats(const pos_stats_uncompr_t &psu);
    static void decompress_position_stats(T psc, pos_stats_uncompr_t &psu);
    static T increment(T psc, nt4_t nt4, int32_t &cov_est);
    static void increment_uncompressed(nt4_t nt4, pos_stats_uncompr_t &psu);

    /***********************
     * Debuging & checking *
     ***********************/

    // Check if everything was initialized.
    bool check_allocation() const;

    // Check if a BAM header corresponds to the stats.
    bool check_headers_bam_hdr(const bam_hdr_t &h) const;

    void debug_print_counters() const;

    std::string debug_str_counters(int32_t seqid, int64_t pos) const;
};

template <typename T, int counter_size, int refbase_size>
stats_t<T, counter_size, refbase_size>::stats_t(ococo::params_t *params,
                                                bam_hdr_t &h)
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

template <typename T, int counter_size, int refbase_size>
stats_t<T, counter_size, refbase_size>::~stats_t() {
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

template <typename T, int counter_size, int refbase_size>
int stats_t<T, counter_size, refbase_size>::load_fasta(
    const std::string &fasta_fn) {
    gzFile fp;
    kseq_t *seq;
    int l;
    fp  = gzopen(fasta_fn.c_str(), "r");
    seq = kseq_init(fp);

    if (fp == nullptr) {
        ococo::error("File '%s' could not be opened.\n", fasta_fn.c_str());
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
            seq_stats[seqid][pos] = compress_position_stats(psu);
        }
    }
    kseq_destroy(seq);
    gzclose(fp);
    return 0;
}

template <typename T, int counter_size, int refbase_size>
int stats_t<T, counter_size, refbase_size>::save_fasta(
    const std::string &fasta_fn) const {
    assert(check_allocation());

    FILE *fasta_file = nullptr;
    fasta_file       = fopen(fasta_fn.c_str(), "w+");

    char fasta_buffer[fasta_line_l];
    for (int s = 0; s < n_seqs; s++) {
        // printf("%s\n",seq_name[s]);
        if (!seq_comment[s].empty()) {
            fprintf(fasta_file, ">%s %s\n", seq_name[s].c_str(),
                    seq_comment[s].c_str());
        } else {
            fprintf(fasta_file, ">%s\n", seq_name[s].c_str());
        }

        for (int64_t i = 0, j = 0; i < seq_len[s]; i++, j++) {
            get_nucl_nt256(s, i, fasta_buffer[j]);

            if (j == fasta_line_l - 1 || i == seq_len[s] - 1) {
                fwrite(fasta_buffer, 1, j + 1, fasta_file);
                fwrite("\n", 1, 1, fasta_file);
                j = -1;
            }
        }
    }

    fclose(fasta_file);

    return 0;
}

template <typename T, int counter_size, int refbase_size>
bool stats_t<T, counter_size, refbase_size>::check_allocation() const {
    if (seq_active == nullptr || seq_len == nullptr || seq_stats == nullptr ||
        seq_name == nullptr || seq_comment == nullptr) {
        return false;
    }

    for (int i = 0; i < n_seqs; i++) {
        if (seq_stats[i] == nullptr) {
            return false;
        }
    }

    return true;
}

template <typename T, int counter_size, int refbase_size>
bool stats_t<T, counter_size, refbase_size>::check_headers_bam_hdr(
    const bam_hdr_t &h) const {
    if (!check_allocation()) {
        return false;
    }

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

template <typename T, int counter_size, int refbase_size>
int stats_t<T, counter_size, refbase_size>::import_stats(
    const std::string &stats_fn) {
    assert(check_allocation());

    int error_code = 0;

    FILE *fo = fopen(stats_fn.c_str(), "r");
    if (fo == nullptr) {
        ococo::error("File with statistics could not be opened ('%s').\n",
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
        ococo::error("File with statistics could not be closed ('%s').\n",
                     stats_fn.c_str());
        return -1;
    }

    return 0;
}

template <typename T, int counter_size, int refbase_size>
int stats_t<T, counter_size, refbase_size>::export_stats(
    const std::string &stats_fn) const {
    assert(check_allocation());

    int error_code = 0;

    FILE *fo = fopen(stats_fn.c_str(), "w+");
    if (fo == nullptr) {
        ococo::error("File with statistics could not be opened ('%s').\n",
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
            ococo::error(
                "Problem with writting to the file with statistics ('%s').\n",
                stats_fn.c_str());
            return -1;
        }
    }

    error_code = fclose(fo);
    if (error_code != 0) {
        ococo::error("File with statistics could not be closed ('%s').\n",
                     stats_fn.c_str());
        return -1;
    }
    return 0;
}

template <typename T, int counter_size, int refbase_size>
int stats_t<T, counter_size, refbase_size>::call_consensus(
    FILE *vcf_file, FILE *out_pileup_file) {
    assert(check_allocation());

    for (int32_t seqid = 0; seqid < n_seqs; seqid++) {
        for (int64_t pos = 0; pos < seq_len[seqid]; pos++) {
            call_consensus_position(vcf_file, out_pileup_file, seqid, pos);
        }
    }

    return 0;
}

template <typename T, int counter_size, int refbase_size>
int stats_t<T, counter_size, refbase_size>::
    call_consensus_position_uncompressed(FILE *vcf_file, FILE *out_pileup_file,
                                         int32_t seqid, int64_t pos,
                                         pos_stats_uncompr_t &psu) {
    const char old_base_nt256 = nt16_nt256[psu.nt16];
    const char new_base_nt256 = cons_call_maj(psu, *params);

    if (old_base_nt256 != new_base_nt256) {
        params->n_upd += 1;
        psu.nt16 = nt256_nt16[static_cast<int16_t>(new_base_nt256)];
    }

    if (vcf_file != nullptr) {
        if (old_base_nt256 != new_base_nt256 || params->verbose) {
            print_vcf_substitution(vcf_file, seqid, pos, old_base_nt256,
                                   new_base_nt256, psu);
        }
    }

    if (out_pileup_file != nullptr) {
        print_pileup_line(out_pileup_file, seqid, pos, psu);
    }

    return 0;
}

template <typename T, int counter_size, int refbase_size>
int stats_t<T, counter_size, refbase_size>::call_consensus_position(
    FILE *vcf_file, FILE *out_pileup_file, int32_t seqid, int64_t pos) {
    pos_stats_uncompr_t psu;
    decompress_position_stats(seq_stats[seqid][pos], psu);
    call_consensus_position_uncompressed(vcf_file, out_pileup_file, seqid, pos,
                                         psu);
    seq_stats[seqid][pos] = compress_position_stats(psu);

    return 0;
}

template <typename T, int counter_size, int refbase_size>
T stats_t<T, counter_size, refbase_size>::compress_position_stats(
    const pos_stats_uncompr_t &psu) {
    T psc = 0;

    // remove if you want to support ambiguous nucleotides
    assert(bitsset_table256[psu.nt16] != 2);

    // 1. incorporate counters
    for (int32_t i = 0; i < 4; i++) {
        psc <<= counter_size;
        psc |= psu.counters[i] & right_full_mask<T, counter_size>();
    }

    // 2. incorporate ref base
    psc <<= refbase_size;
    psc |= psu.nt16 & right_full_mask<T, refbase_size>();

    // 3. if not exact, invert the base bits
    if (psu.bitshifted) {
        psc ^= right_full_mask<T, refbase_size>();
    }

    return psc;
}

template <typename T, int counter_size, int refbase_size>
void stats_t<T, counter_size, refbase_size>::decompress_position_stats(
    T psc, pos_stats_uncompr_t &psu) {
    // 1. reference base(s) (before correction)
    psu.nt16 = psc & right_full_mask<T, refbase_size>();
    psc >>= refbase_size;

    // 2. are the values exact?
    int nones = bitsset_table256[psu.nt16];
    assert(nones != 2);
    if (nones == 1) {
        psu.bitshifted = false;
    } else {
        if (nones == 3) {
            psu.bitshifted = true;
            // if not exact, invert base bits
            psu.nt16 ^= right_full_mask<T, refbase_size>();
        }
    }

    // 3. count of individual nucleotides and the sum
    psu.sum = 0;
    for (int32_t i = 3; i >= 0; i--) {
        psu.counters[i] = psc & right_full_mask<T, counter_size>();
        psu.sum += psu.counters[i];
        psc >>= counter_size;
    }
}

template <typename T, int counter_size, int refbase_size>
int stats_t<T, counter_size, refbase_size>::print_vcf_header(
    FILE *vcf_file, std::string cmd, std::string fasta) const {
    assert(check_allocation());
    assert(vcf_file != nullptr);

    std::time_t tt = std::time(nullptr);
    tm *tm         = localtime(&tt);

    fprintf(vcf_file,
            "##fileformat=VCFv4.3\n"
            "##fileDate=%04d%02d%02d\n"
            "##source=Ococo\n",
            tm->tm_year + 1900, tm->tm_mon + 1, tm->tm_mday);

    if (!cmd.empty()) {
        fprintf(vcf_file, "##ococo_command=%s\n", cmd.c_str());
    }
    fprintf(vcf_file, "##ococo_stats_datatype_size=%zubits\n", 8 * sizeof(T));
    fprintf(vcf_file, "##ococo_counter_size=%dbits\n", counter_size);

    if (!fasta.empty()) {
        fprintf(vcf_file, "##reference=%s\n", fasta.c_str());
    }

    for (int seqid = 0; seqid < n_seqs; seqid++) {
        fprintf(vcf_file, "##contig=<ID=%s,length=%" PRId64 ">\n",
                seq_name[seqid].c_str(), seq_len[seqid]);
    }

    fprintf(vcf_file,
            "##INFO=<ID=AF,Number=A,Type=Float,Description="
            "\"Allele frequency for the ALT allele.\">\n");
    fprintf(vcf_file,
            "##INFO=<ID=CS,Number=4,Type=Integer,Description="
            "\"Values of A,C,G,T counters.\">\n");
    fprintf(vcf_file,
            "##INFO=<ID=COV,Number=1,Type=Integer,Description="
            "\"Coverage\">\n");
    fprintf(
        vcf_file,
        "##INFO=<ID=EX,Number=1,Type=Integer,Description="
        "\"1 if the coverage and counter values are exact (no bitshift made), 0 otherwise\">\n");
    fprintf(vcf_file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

    return 0;
}

template <typename T, int counter_size, int refbase_size>
int stats_t<T, counter_size, refbase_size>::print_vcf_substitution(
    FILE *vcf_file, int32_t seqid, int64_t pos, char old_base, char new_base,
    const pos_stats_uncompr_t &psu) const {
    assert(check_allocation());
    assert(vcf_file != nullptr);

    float alt_freq =
        1.0 * psu.counters[nt256_nt4[static_cast<int16_t>(new_base)]] / psu.sum;

    fprintf(vcf_file,
            "%s\t%" PRId64 "\t.\t%c\t%c\t100\tPASS\tAF=%.2f;CS=%" PRId32
            ",%" PRId32 ",%" PRId32 ",%" PRId32 ";COV=%" PRId32 ";EX=%s\n",
            seq_name[seqid].c_str(), pos + 1, old_base, new_base,
            round(alt_freq * 100.0) / 100, psu.counters[0], psu.counters[1],
            psu.counters[2], psu.counters[3], psu.sum,
            (psu.bitshifted ? "0" : "1"));

    return 0;
}

template <typename T, int counter_size, int refbase_size>
int stats_t<T, counter_size, refbase_size>::print_pileup_line(
    FILE *out_pileup_file, int32_t seqid, int64_t pos,
    const pos_stats_uncompr_t &psu) const {
    assert(check_allocation());
    assert(out_pileup_file != nullptr);

    // todo: fix situation when depth is larger (use the printing buffer more
    // timess)

    const int32_t max_depth = 1000;

    if (psu.sum >= max_depth) {
        ococo::error("Too high coverage at position %" PRId64
                     ". Pileup does not support coverage higher than %" PRId32
                     ".",
                     pos, max_depth);
        return -1;
    }

    char bases[max_depth];
    char qualities[max_depth];

    char ref_nt256 = nt16_nt256[psu.nt16];

    if (psu.sum == 0) {
        return 0;
    }

    int32_t j = 0;

    for (int32_t nt4 = 0; nt4 < 4; nt4++) {
        const char filling_char =
            nt4_nt16[nt4] == psu.nt16 ? '.' : nt4_nt256[nt4];
        for (int32_t i = 0; i < psu.counters[nt4]; i++, j++) {
            bases[j]     = filling_char;
            qualities[j] = '~';
        }
    }

    bases[j]     = '\0';
    qualities[j] = '\0';

    fprintf(out_pileup_file, "%s\t%" PRId64 "\t%c\t%" PRId32 "\t%s\t%s\n",
            seq_name[seqid].c_str(), pos + 1, ref_nt256, psu.sum, bases,
            qualities);

    return 0;
}

template <typename T, int counter_size, int refbase_size>
std::string stats_t<T, counter_size, refbase_size>::debug_str_counters(
    int32_t seqid, int64_t pos) const {
    pos_stats_uncompr_t psu;
    decompress_position_stats(seq_stats[seqid][pos], psu);
    std::stringstream ss;
    ss << "[" << nt16_nt256[psu.nt16] << "]"
       << "(" << psu.counters[0] << "," << psu.counters[1] << ","
       << psu.counters[2] << "," << psu.counters[3] << ")";
    return ss.str();
}

template <typename T, int counter_size, int refbase_size>
void stats_t<T, counter_size, refbase_size>::debug_print_counters() const {
    for (int seqid = 0; seqid < n_seqs; seqid++) {
        fprintf(stderr, "%s\n", seq_name[seqid]);
        for (int64_t pos = 0; pos < seq_len[seqid]; pos++) {
            fprintf(stderr, "%8" PRId64 " %04x \n", pos, seq_stats[seqid][pos]);
        }
    }
}

template <typename T, int counter_size, int refbase_size>
inline int stats_t<T, counter_size, refbase_size>::set_nucl_nt256(
    int32_t seqid, int64_t pos, const char &nt256) {
    pos_stats_uncompr_t psu;
    decompress_position_stats(seq_stats[seqid][pos], psu);
    psu.nt16              = nt256_nt16[static_cast<int16_t>(nt256)];
    seq_stats[seqid][pos] = compress_position_stats(psu);
    return 0;
}

template <typename T, int counter_size, int refbase_size>
inline int stats_t<T, counter_size, refbase_size>::get_nucl_nt256(
    int32_t seqid, int64_t pos, char &nt256) const {
    pos_stats_uncompr_t psu;
    decompress_position_stats(seq_stats[seqid][pos], psu);
    nt256 = nt16_nt256[psu.nt16];
    if (nt256 == '=') {
        nt256 = 'N';
    }
    return 0;
}

template <typename T, int counter_size, int refbase_size>
T stats_t<T, counter_size, refbase_size>::increment(T psc, nt4_t nt4,
                                                    int32_t &cov_est) {
    assert(nt4 < 4);

    pos_stats_uncompr_t psu;
    decompress_position_stats(psc, psu);

    increment_uncompressed(nt4, psu);
    cov_est = psu.sum;

    return compress_position_stats(psu);
}

template <typename T, int counter_size, int refbase_size>
void stats_t<T, counter_size, refbase_size>::increment_uncompressed(
    nt4_t nt4, pos_stats_uncompr_t &psu) {
    assert(nt4 < 4);

    if (psu.counters[nt4] == right_full_mask<uint16_t, counter_size>()) {
        psu.counters[0] >>= 1;
        psu.counters[1] >>= 1;
        psu.counters[2] >>= 1;
        psu.counters[3] >>= 1;
        psu.bitshifted = true;
    }

    psu.counters[nt4]++;

    psu.sum =
        psu.counters[0] + psu.counters[1] + psu.counters[2] + psu.counters[3];
}
}  // namespace ococo
