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

#include "params.h"

#include <getopt.h>
#include <iostream>

namespace ococo {

/****************************
 *** Consensus parameters ***
 ****************************/

void params_t::init_default_values() {
    verbose               = false;
    counters_str          = "ococo32";
    counter_configuration = OCOCO32;
    mode                  = BATCH;
    mode_str              = "batch";
    min_mapq              = default_q;
    min_baseq             = default_Q;
    min_coverage_upd      = default_c;
    majority_threshold    = default_M;
    coverage_filter       = default_C;

    out_fasta_file = nullptr;

    n_upd = 0;

    correctly_initialized = true;
    return_code           = 0;
}

params_t::params_t() { init_default_values(); }

params_t::params_t(int argc, const char **argv) {
    init_default_values();
    parse_commandline(argc, argv);
}

void params_t::print_help() {
    print_version();

    std::cerr
        <<
        // clang-format off

        // "---------------------------------------------------------------------------------"
           "Usage:   ococo -i <SAM/BAM file> [other options]\n\n"
        // "---------------------------------------------------------------------------------"
           "Input:\n"
           "  -i, --input FILE      input SAM/BAM file (- for standard input)\n"
           "  -f, --fasta-ref FILE  initial FASTA reference (otherwise seq of N's is used)\n"
           "  -s, --stats-in FILE   input statistics\n\n"
        // "---------------------------------------------------------------------------------"
           "Output:\n"
           "  -O, --output FILE     output SAM/BAM file (- for standard output)\n"
           "  -F, --fasta-cons FILE FASTA file with the consensus\n"
           "  -S, --stats-out FILE  output statistics\n"
           "  -V, --vcf-cons FILE   VCF file with updates of consensus (- for standard output)\n"
           "  -P, --pileup FILE     truncated pileup (- for standard output)\n\n"
        // "---------------------------------------------------------------------------------"
           "Variant and consensus calling:\n"
           "  -x, --counters STR    counter configuration: [ococo32]\n"
           "                           - ococo16 (3b/counter, 16b/position)\n"
           "                           - ococo32 (7b/counter, 32b/position)\n"
           "                           - ococo64 (15b/counter, 64b/position)\n"
           "  -m, --mode STR        mode: [batch]\n"
           "                           - real-time (updates reported immediately)\n"
           "                           - batch (updates reported after end of algn stream)\n"
           "  -c, --min-cov INT     minimum coverage required for an update [" << default_c <<"]\n"
           "  -M, --maj-thres FLOAT majority threshold [" << default_M << "]\n\n"
       // "---------------------------------------------------------------------------------"
           "Filtering and coverage normalization:\n"
           "  -q, --min-MQ INT      skip alignments with mapping quality smaller than INT [" << default_q << "]\n"
           "  -Q, --min-BQ INT      skip bases with base quality smaller than INT [" << default_Q <<"]\n"
           "  -C, --cov-filt INT    skip alignments when coverage greater than INT [" << default_C << "]\n\n"
       // "---------------------------------------------------------------------------------"
           "Others:\n"
           "  -v, --version         print version information and exit\n"
           "  -h, --help            print this message and exit\n"
           //"  --log FILE            auxiliary log file\n"
           "  --verbose             verbose mode (report every update of a counter)\n\n"
        // "---------------------------------------------------------------------------------"
           "Examples:\n"
           "   ococo -i test.bam -f test.fa -m real-time -V -\n"
           "   ococo -x ococo64 -i test.bam -f test.fa -P - -V variants.vcf\n\n"
        // "---------------------------------------------------------------------------------"
           "Note:\n"
           "   For more details, see the manual page 'man ./ococo.1'.\n"
        // "---------------------------------------------------------------------------------"
        // clang-format on
        << std::endl;
}

void params_t::parse_commandline(int argc, const char **argv) {
    /* Parse cmd parameters */
    std::stringstream cmd;
    for (int32_t i = 0; i < argc; i++) {
        cmd << argv[i];
        if (i != argc - 1) {
            cmd << " ";
        }
    }
    command = cmd.str();

    if (argc == 1) {
        print_help();
        exit(1);
    }

    const struct option lopts[] = {
        {"version", no_argument, nullptr, 'v'},
        {"help", no_argument, nullptr, 'h'},
        //
        {"input", required_argument, nullptr, 'i'},
        {"fasta-ref", required_argument, nullptr, 'f'},
        {"stats-in", required_argument, nullptr, 's'},
        //
        {"fasta-cons", required_argument, nullptr, 'F'},
        {"stats-out", required_argument, nullptr, 'S'},
        {"vcf-cons", required_argument, nullptr, 'V'},
        {"pileup", required_argument, nullptr, 'P'},
        {"output", required_argument, nullptr, 'O'},
        {"cov-filt", required_argument, nullptr, 'C'},
        {"log", required_argument, nullptr, 'L'},
        {"verbose", no_argument, nullptr, 'W'},
        //
        {"counters", required_argument, nullptr, 'x'},
        {"mode", required_argument, nullptr, 'm'},
        {"min-MQ", required_argument, nullptr, 'q'},
        {"min-BQ", required_argument, nullptr, 'Q'},
        {"ref-weight", required_argument, nullptr, 'w'},
        {"min-cov", required_argument, nullptr, 'c'},
        {"min-coverage", required_argument, nullptr, 'c'},  // deprec
        {"maj-thres", required_argument, nullptr, 'M'},
        {"majority-threshold", required_argument, nullptr, 'M'},  // deprec
        //
        {nullptr, 0, nullptr, 0}};

    int c;
    using std::string;
    while ((c = getopt_long(argc, (char *const *)argv,
                            "vhi:f:s:F:S:V:P:L:W:x:m:t:q:Q:w:c:M:O:C:", lopts,
                            nullptr)) >= 0) {
        switch (c) {
            case 'v': {
                print_version();
                exit(0);
                break;
            }
            case 'h': {
                print_help();
                exit(0);
                break;
            }
            case 'i': {
                in_sam_fn = optarg;
                break;
            }
            case 'f': {
                in_fasta_fn = optarg;
                break;
            }
            case 's': {
                in_stats_fn = optarg;
                break;
            }
            case 'F': {
                out_fasta_fn = optarg;
                break;
            }
            case 'S': {
                out_stats_fn = optarg;
                break;
            }
            case 'V': {
                out_vcf_fn = optarg;
                break;
            }
            case 'P': {
                out_pileup_fn = optarg;
                break;
            }
            case 'O': {
                out_sam_fn = optarg;
                break;
            }
            case 'C': {
                coverage_filter = atoi(optarg);
                break;
            }
            case 'L': {
                out_log_fn = optarg;
                break;
            }
            case 'W': {
                verbose = true;
                break;
            }
            case 'x': {
                counters_str = optarg;

                if (counters_str.compare("ococo8") == 0) {
                    counter_configuration = OCOCO8;
                    counters_str_descr =
                        "ococo8 (8 bits per position, 1bits per nucleotide "
                        "counter)";
                } else if (counters_str.compare("ococo16") == 0) {
                    counter_configuration = OCOCO16;
                    counters_str_descr =
                        "ococo16 (16 bits per position, 3bits per nucleotide "
                        "counter)";
                } else if (counters_str.compare("ococo32") == 0) {
                    counter_configuration = OCOCO32;
                    counters_str_descr =
                        "ococo32 (32 bits per position, 7bits per nucleotide "
                        "counter)";
                } else if (counters_str.compare("ococo64") == 0) {
                    counter_configuration = OCOCO64;
                    counters_str_descr =
                        "ococo64 (64 bits per position, 15bits per nucleotide "
                        "counter)";
                } else {
                    error(
                        "Unknown counter configuration '%s'. Possible modes "
                        "are 'ococo8', 'ococo16', 'ococo32', and 'ococo64'.\n",
                        counters_str.c_str());
                    exit(1);
                }

                break;
            }
            case 'm': {
                mode_str = optarg;

                if (mode_str.compare("batch") == 0) {
                    mode = mode_t::BATCH;
                } else if (mode_str.compare("real-time") == 0) {
                    mode = mode_t::REALTIME;
                } else {
                    error(
                        "Unknown mode '%s'. Possible modes are 'batch' and "
                        "'real-time'.\n",
                        mode_str.c_str());
                    exit(1);
                }

                break;
            }
            case 'q': {
                min_mapq = atoi(optarg);
                break;
            }
            case 'Q': {
                min_baseq = atoi(optarg);
                break;
            }
            case 'c': {
                min_coverage_upd = atoi(optarg);
                break;
            }
            case 'M': {
                majority_threshold = atof(optarg);
                break;
            }
            case '?': {
                error("Unknown error");
                exit(1);
                break;
            }
        }
    }
    if (in_sam_fn.size() == 0) {
        error("SAM/BAM file must be specified (option '-i').\n");
        exit(1);
    }

    if (out_pileup_fn.size() != 0 && out_pileup_fn.compare(out_vcf_fn) == 0) {
        error(
            "Pileup and VCF files cannot be the same (both currently '%s').\n",
            out_pileup_fn.c_str());
        exit(1);
    }

    info("Ococo starting: %s\n", counters_str_descr.c_str());
}
}  // namespace ococo
