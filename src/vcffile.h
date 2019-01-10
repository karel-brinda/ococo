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

#include "counters.h"
#include "misc.h"
#include "params.h"

namespace ococo {

struct VcfFile {
    std::string fn;
    FILE *file;
    Params params;

    VcfFile(std::string fn, Params params)
        : fn(fn), file(nullptr), params(params) {
        if (!fn.empty()) {
            info("Opening the VCF stream ('%s').\n", fn.c_str());

            if (fn == std::string("-")) {
                file = stdout;
            } else {
                file = fopen(fn.c_str(), "w+");
                if (file == nullptr) {
                    fatal_error("Problem with opening the VCF file '%s'.\n",
                                fn.c_str());
                }
            }

            print_header();
        }
    }

    ~VcfFile() {
        if (file != nullptr && fn != "-") {
            int error_code = fclose(file);
            if (error_code != 0) {
                error("Output VCF file could not be closed.\n");
            }
        }
    }

    void print_header() const {
        if (file == nullptr) {
            return;
        }
        std::time_t tt = std::time(nullptr);
        tm *tm         = localtime(&tt);

        fprintf(file,
                "##fileformat=VCFv4.3\n"
                "##fileDate=%04d%02d%02d\n"
                "##source=Ococo\n",
                tm->tm_year + 1900, tm->tm_mon + 1, tm->tm_mday);

        if (!params.command.empty()) {
            fprintf(file, "##ococo_command=%s\n", params.command.c_str());
        }
        /*fprintf(file, "##ococo_stats_datatype_size=%zubits\n", 8 * sizeof(T));
        fprintf(file, "##ococo_C=%dbits\n", C);

        if (!fasta.empty()) {
            fprintf(file, "##reference=%s\n", fasta.c_str());
        }

        for (int seqid = 0; seqid < n_seqs; seqid++) {
            fprintf(file, "##contig=<ID=%s,length=%" PRId64 ">\n",
                    seq_name[seqid].c_str(), seq_len[seqid]);
        }*/

        fprintf(file,
                "##INFO=<ID=C,Number=4,Type=Integer,Description="
                "\"A,C,G,T counters.\">\n");
        fprintf(file,
                "##INFO=<ID=COV,Number=1,Type=Integer,Description="
                "\"Coverage\">\n");
        fprintf(file,
                "##INFO=<ID=AF,Number=A,Type=Float,Description="
                "\"Allele frequency for the ALT allele.\">\n");
        fprintf(
            file,
            "##INFO=<ID=EX,Number=1,Type=Integer,Description="
            "\"Values are exact (1=yes, 0=no), i.e., no bitshift made\">\n");
        fprintf(file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    }

    void print_substitution(const std::string &seq_name, int64_t pos,
                            char old_base, char new_base,
                            const pos_stats_uncompr_t &psu) const {
        if (file == nullptr) {
            return;
        }

        const float alt_freq =
            1.0 * psu.counters[nt256_nt4[static_cast<int16_t>(new_base)]] /
            psu.sum;

        fprintf(file, "%s\t%" PRId64 "\t.\t%c\t%c\t100\tPASS\t",
                seq_name.c_str(),  //
                pos + 1,           //
                old_base,          //
                new_base           //
        );

        fprintf(file,
                "C="
                "%" PRId32       //
                ",%" PRId32      //
                ",%" PRId32      //
                ",%" PRId32      //
                ";COV=%" PRId32  //
                ";AF=%.2f"       //
                ";EX=%s\n",
                //
                psu.counters[0],                //
                psu.counters[1],                //
                psu.counters[2],                //
                psu.counters[3],                //
                psu.sum,                        //
                round(alt_freq * 100.0) / 100,  //
                (psu.bitshifted ? "0" : "1")    //
        );
    }
};

}  // namespace ococo