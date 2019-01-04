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

#include <ctime>
#include "io.h"

namespace ococo {

int print_vcf_header(FILE *file) {
    std::time_t tt = std::time(nullptr);
    tm *tm         = localtime(&tt);

    fprintf(file,
            "##fileformat=VCFv4.3\n"
            "##fileDate=%04d%02d%02d\n"
            "##source=Ococo\n",
            tm->tm_year + 1900, tm->tm_mon + 1, tm->tm_mday);

    // todo: pass as a list of pairs

    /*if (!cmd.empty()) {
        fprintf(file, "##ococo_command=%s\n", cmd.c_str());
    }
    fprintf(file, "##ococo_stats_datatype_size=%zubits\n", 8 * sizeof(T));
    fprintf(file, "##ococo_C=%dbits\n", C);

    if (!fasta.empty()) {
        fprintf(file, "##reference=%s\n", fasta.c_str());
    }

    for (int seqid = 0; seqid < n_seqs; seqid++) {
        fprintf(file, "##contig=<ID=%s,length=%" PRId64 ">\n",
                seq_name[seqid].c_str(), seq_len[seqid]);
    }*/

    fprintf(file,
            "##INFO=<ID=AF,Number=A,Type=Float,Description="
            "\"Allele frequency for the ALT allele.\">\n");
    fprintf(file,
            "##INFO=<ID=CS,Number=4,Type=Integer,Description="
            "\"Values of A,C,G,T counters.\">\n");
    fprintf(file,
            "##INFO=<ID=COV,Number=1,Type=Integer,Description="
            "\"Coverage\">\n");
    fprintf(file,
            "##INFO=<ID=EX,Number=1,Type=Integer,Description="
            "\"1 if the coverage and counter values are exact (no bitshift "
            "made), 0 otherwise\">\n");
    fprintf(file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

    return 0;
}

int print_vcf_substitution(FILE *file, const std::string &seq_name, int64_t pos,
                           char old_base, char new_base,
                           const pos_stats_uncompr_t &psu) {
    const float alt_freq =
        1.0 * psu.counters[nt256_nt4[static_cast<int16_t>(new_base)]] / psu.sum;

    fprintf(file,
            "%s\t%" PRId64 "\t.\t%c\t%c\t100\tPASS\tAF=%.2f;CS=%" PRId32
            ",%" PRId32 ",%" PRId32 ",%" PRId32 ";COV=%" PRId32 ";EX=%s\n",
            seq_name.c_str(), pos + 1, old_base, new_base,
            round(alt_freq * 100.0) / 100, psu.counters[0], psu.counters[1],
            psu.counters[2], psu.counters[3], psu.sum,
            (psu.bitshifted ? "0" : "1"));

    return 0;
}

int print_pileup_line(FILE *file, const std::string &seq_name, int64_t pos,
                      const pos_stats_uncompr_t &psu) {
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

    // todo: check that j does not exceeds buffer size
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

    fprintf(file, "%s\t%" PRId64 "\t%c\t%" PRId32 "\t%s\t%s\n",
            seq_name.c_str(), pos + 1, ref_nt256, psu.sum, bases, qualities);

    return 0;
}

void fatal_error(const char *format, ...) {
    va_list args;
    va_start(args, format);
    fprintf(stderr, "[ococo:fatal-error]: ");
    vfprintf(stderr, format, args);
    va_end(args);
}

void error(const char *format, ...) {
    va_list args;
    va_start(args, format);
    fprintf(stderr, "[ococo:error]: ");
    vfprintf(stderr, format, args);
    va_end(args);
}

void warning(const char *format, ...) {
    va_list args;
    va_start(args, format);
    fprintf(stderr, "[ococo:warning]: ");
    vfprintf(stderr, format, args);
    va_end(args);
}

void info(const char *format, ...) {
    va_list args;
    va_start(args, format);
    fprintf(stderr, "[ococo]: ");
    vfprintf(stderr, format, args);
    va_end(args);
}

bool file_exists(const std::string &fn) {
    FILE *file;

    file = fopen(fn.c_str(), "r");
    if (file) {
        fclose(file);
        return true;
    }
    return false;
}

}  // namespace ococo
